// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TelescopeDetectorElement.h, Acts project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"

#include "Acts/Surfaces/PlanarBounds.hpp"

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"

#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"

namespace Acts {

class Surface;
class PlanarBounds;
class ISurfaceMaterial;
}  // namespace Acts

namespace Telescope {

/// @class TelescopeDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class TelescopeDetectorElement : public Acts::DetectorElementBase {
 public:

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  TelescopeDetectorElement(
      std::shared_ptr<const Acts::Transform3D> transform,
      double widthX, double heightY, double thickZ)
    : Acts::DetectorElementBase(),
      m_elementTransform(std::move(transform)),
      m_elementThickness(thickZ){

    auto pBounds = std::make_shared<Acts::RectangleBounds>(widthX/2.0, heightY/2.0);

    auto material = std::make_shared<Acts::HomogeneousSurfaceMaterial>
    (Acts::MaterialProperties(95.7, 465.2, 28.03, 14., 2.32e-3, m_elementThickness));

    m_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(pBounds, *this);
    m_surface->assignSurfaceMaterial(material);

    std::unique_ptr<Acts::SurfaceArray> surArray(new Acts::SurfaceArray(m_surface));
    m_layer = Acts::PlaneLayer::create(m_elementTransform, pBounds, std::move(surArray), 5 * Acts::UnitConstants::mm);
    m_surface->associateLayer(*m_layer);
  }

  ///  Destructor
  ~TelescopeDetectorElement() override { /*nop */
  }

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Acts::Transform3D& transform(
      const Acts::GeometryContext& gctx) const override;

  /// Return the nominal local to global transform
  ///
  /// @note the geometry context will hereby be ignored
  const Acts::Transform3D& nominalTransform(
      const Acts::GeometryContext& gctx) const;

  /// Return local to global transform associated with this identifier
  ///
  /// @param alignedTransform is a new transform
  /// @oaram iov is the batch for which it is meant
  void addAlignedTransform(std::unique_ptr<Acts::Transform3D> alignedTransform);


  std::shared_ptr<const Acts::Layer> layer() const {
    return m_layer;
  }

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override;

 private:
  /// the transform for positioning in 3D space
  std::shared_ptr<const Acts::Transform3D> m_elementTransform;
  // the aligned transforms
  std::unique_ptr<Acts::Transform3D> m_alignedTransforms;
  /// the surface represented by it
  std::shared_ptr<Acts::Surface> m_surface;
  std::shared_ptr<Acts::Layer> m_layer;
  /// the element thickness
  double m_elementThickness{0.};
};

inline const Acts::Transform3D& TelescopeDetectorElement::transform(
    const Acts::GeometryContext& gctx) const {
  // Check if a different transform than the nominal exists
  if (m_alignedTransforms) {
    return (*m_alignedTransforms);
  }
  // Return the standard transform if not found
  return nominalTransform(gctx);
}

inline const Acts::Transform3D& TelescopeDetectorElement::nominalTransform(
    const Acts::GeometryContext& /*gctx*/) const {
  return *m_elementTransform;
}

inline void TelescopeDetectorElement::addAlignedTransform(
    std::unique_ptr<Acts::Transform3D> alignedTransform) {
  m_alignedTransforms = std::move(alignedTransform);
}

inline const Acts::Surface& TelescopeDetectorElement::surface() const {
  return *m_surface;
}

inline double TelescopeDetectorElement::thickness() const {
  return m_elementThickness;
}
}  // namespace Telescope
