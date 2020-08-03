#include <array>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"

#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include "TelescopeTrack.hpp"
#include "TelescopeDetectorElement.hpp"

#include "myrapidjson.h"

using namespace Acts::UnitLiterals;

void Telescope::BuildGeometry(
                              Acts::GeometryContext& nominal_gctx,
                              std::shared_ptr<const Acts::TrackingGeometry>& geo_world,
                              std::vector<std::shared_ptr<Acts::DetectorElementBase>>& element_col,
                              std::vector<std::shared_ptr<const Acts::Surface>>& surface_col,
                              std::vector<Acts::LayerPtr>& layer_col,
                              const rapidjson::GenericValue<rapidjson::UTF8<>, rapidjson::CrtAllocator> &js_opt
                              ){

  // Boundaries of the surfaces (ALPIDE real SIZE: 29.941760 * 13.762560_mm*mm)
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>
   (Acts::RectangleBounds(30_mm/2.0, 14_mm/2.0));

  // Material of the surfaces
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>
    (Acts::MaterialProperties(95.7, 465.2, 28.03, 14., 2.32e-3, 80_um));

  // Set translation vectors

  for(const auto& js_l: js_opt.GetArray()){

    // Acts::Vector3D center3d;
    double cx = js_l["centerX"].GetDouble();
    double cy = js_l["centerY"].GetDouble();
    double cz = js_l["centerZ"].GetDouble();
    Acts::Vector3D translation(cx,cy,cz);

    //Acts::Rotation3D
    double rx = js_l["rotX"].GetDouble();
    double ry = js_l["rotY"].GetDouble();
    double rz = js_l["rotZ"].GetDouble();
    // The rotation around global z axis
    Acts::AngleAxis3D rotZ(rz, Acts::Vector3D::UnitZ());
    // The rotation around global y axis
    Acts::AngleAxis3D rotY(ry, Acts::Vector3D::UnitY());
    // The rotation around global x axis
    Acts::AngleAxis3D rotX(rx, Acts::Vector3D::UnitX());
    Acts::Rotation3D rotation = rotZ * rotY * rotX;

    //Acts::Transform3D
    auto trafo = std::make_shared<Acts::Transform3D>
      (Acts::Translation3D(translation) * rotation);

    // Create the detector element
    auto detElement = std::make_shared<Telescope::TelescopeDetectorElement>
      (trafo, rBounds, 1_um, surfaceMaterial);
    auto detSurface = detElement->surface().getSharedPtr();
    surface_col.push_back(detSurface);
    element_col.push_back(detElement);

    std::unique_ptr<Acts::SurfaceArray> surArray(new Acts::SurfaceArray(detSurface));
    auto detLayer = Acts::PlaneLayer::create(trafo, rBounds, std::move(surArray), 1.5_cm);
    layer_col.push_back(detLayer);

    auto mutableSurface = std::const_pointer_cast<Acts::Surface>(detSurface);
    mutableSurface->associateLayer(*detLayer);
  }



  // Build tracking volume
  float tracker_halfX= 0.1_m;
  float tracker_halfY= 0.1_m;
  float tracker_halfZ= 1.2_m;
  auto tracker_cuboid = std::make_shared<Acts::CuboidVolumeBounds>(tracker_halfX, tracker_halfY, tracker_halfZ);

  auto trafo_no_move = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  trafo_no_move->translation() = Acts::Vector3D(0., 0., 0.);

  auto layer_array_binned = Acts::LayerArrayCreator({}).
    layerArray(nominal_gctx, layer_col,
               - tracker_halfZ, tracker_halfZ,
               Acts::BinningType::arbitrary,
               Acts::BinningValue::binZ);

  auto trackVolume = Acts::TrackingVolume::create
    (trafo_no_move, tracker_cuboid, nullptr, std::move(layer_array_binned), nullptr, {}, "Tracker");



  // Build world volume
  // assumming single tracker in the world
  auto tracker_array_binned = Acts::TrackingVolumeArrayCreator({}).
    trackingVolumeArray( nominal_gctx, { trackVolume }, Acts::BinningValue::binZ );
  auto mtvpWorld = Acts::TrackingVolume::create
    (trafo_no_move, tracker_cuboid, tracker_array_binned , "World");

  geo_world.reset(new Acts::TrackingGeometry(mtvpWorld));
}
