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

#include "TelescopeDetectorElement.hpp"
#include "TelescopeTrack.hpp"

#include "myrapidjson.h"

using namespace Acts::UnitLiterals;

void Telescope::BuildGeometry(
    Acts::GeometryContext &nominal_gctx,
    std::shared_ptr<const Acts::TrackingGeometry> &geo_world,
    std::vector<std::shared_ptr<Telescope::TelescopeDetectorElement>>
        &element_col,
    const std::map<size_t, std::array<double, 6>> &opts, double widthX,
    double heightY, double thickZ) {

  std::vector<std::shared_ptr<const Acts::Surface>> surface_col;
  std::vector<Acts::LayerPtr> layer_col;

  // Set translation vectors
  // for(const auto& js_l: js_opt.GetArray()){
  for (const auto &l : opts) {
    size_t ln = l.first;
    // Acts::Vector3D center3d;
    // double cx = js_l["centerX"].GetDouble();
    // double cy = js_l["centerY"].GetDouble();
    // double cz = js_l["centerZ"].GetDouble();
    double cx = (l.second)[0];
    double cy = (l.second)[1];
    double cz = (l.second)[2];
    Acts::Vector3D translation(cx, cy, cz);

    // Acts::Rotation3D
    // double rx = js_l["rotX"].GetDouble();
    // double ry = js_l["rotY"].GetDouble();
    // double rz = js_l["rotZ"].GetDouble();
    double rx = (l.second)[3];
    double ry = (l.second)[4];
    double rz = (l.second)[5];
    // The rotation around global z axis
    Acts::AngleAxis3D rotZ(rz, Acts::Vector3D::UnitZ());
    // The rotation around global y axis
    Acts::AngleAxis3D rotY(ry, Acts::Vector3D::UnitY());
    // The rotation around global x axis
    Acts::AngleAxis3D rotX(rx, Acts::Vector3D::UnitX());
    Acts::Rotation3D rotation = rotZ * rotY * rotX;

    // Acts::Transform3D
    auto trafo = std::make_shared<Acts::Transform3D>(
        Acts::Translation3D(translation) * rotation);

    // Create the detector element
    // Boundaries of the surfaces (ALPIDE real SIZE: 29.941760
    // * 13.762560_mm*mm)
    auto detElement = std::make_shared<Telescope::TelescopeDetectorElement>(
        ln, trafo, widthX, heightY, thickZ);

    surface_col.push_back(detElement->surface().getSharedPtr());
    layer_col.push_back(detElement->layer());
    element_col.push_back(detElement);
  }

  // Build tracking volume
  float tracker_halfX = 0.1_m;
  float tracker_halfY = 0.1_m;
  float tracker_halfZ = 1.0_m;
  auto tracker_cuboid = std::make_shared<Acts::CuboidVolumeBounds>(
      tracker_halfX, tracker_halfY, tracker_halfZ);

  auto layer_array_binned = Acts::LayerArrayCreator({}).layerArray(
      nominal_gctx, layer_col, -tracker_halfZ, tracker_halfZ,
      Acts::BinningType::arbitrary, Acts::BinningValue::binZ);

  auto trackVolume = Acts::TrackingVolume::create(
      Acts::Transform3D::Identity(), tracker_cuboid, nullptr,
      std::move(layer_array_binned), nullptr, {}, "Tracker");

  // Build world volume
  // assumming single tracker in the world
  auto tracker_array_binned =
      Acts::TrackingVolumeArrayCreator({}).trackingVolumeArray(
          nominal_gctx, {trackVolume}, Acts::BinningValue::binZ);
  auto mtvpWorld = Acts::TrackingVolume::create(Acts::Transform3D::Identity(),
                                                tracker_cuboid,
                                                tracker_array_binned, "World");

  geo_world.reset(new Acts::TrackingGeometry(mtvpWorld));
}
