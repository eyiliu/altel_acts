#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/GeometryObjectSorter.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"

#include "TelEvent.hpp"


#include "myrapidjson.h"
namespace TelActs {

class TelSourceLink;
using PixelSourceLink = TelSourceLink;

class TelSourceLink {
public:

  TelSourceLink(const Acts::Surface &surface, Acts::Vector2D values,
                Acts::BoundMatrix cov)
      : m_values(values), m_cov(cov), m_surface(&surface) {
  }

  TelSourceLink(const Acts::PlaneLayer &planeLayer, std::shared_ptr<altel::TelMeasHit> hitMeas);
  TelSourceLink(std::shared_ptr<altel::TelMeasHit> hitMeas,
                const std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>> &mapDetId2PlaneLayer);

  /// Must be default_constructible to satisfy SourceLinkConcept.
  TelSourceLink() = default;
  TelSourceLink(TelSourceLink &&) = default;
  TelSourceLink(const TelSourceLink &) = default;
  TelSourceLink &operator=(TelSourceLink &&) = default;
  TelSourceLink &operator=(const TelSourceLink &) = default;

  constexpr const Acts::Surface &referenceSurface() const { return *m_surface; }

  Acts::FittableMeasurement<TelSourceLink> operator*() const {
    return Acts::Measurement<TelSourceLink, Acts::BoundIndices,
                             Acts::eBoundLoc0, Acts::eBoundLoc1>{
        m_surface->getSharedPtr(), *this, m_cov.topLeftCorner<2, 2>(),
        m_values[0], m_values[1]};
  }

  const Acts::Vector2D& value() const{
    return m_values;
  }

  Acts::Vector2D& value(){
    return m_values;
  }

  const size_t& detectorId() const{
    return m_detId;
  }

  size_t& detectorId(){
    return m_detId;
  }

  void setCovariance(const Acts::BoundMatrix &cov) { m_cov = cov; }

  Acts::Vector3D globalPosition(const Acts::GeometryContext &gctx) const {
    Acts::Vector3D mom(1, 1, 1);
    Acts::Vector3D global = m_surface->localToGlobal(gctx, m_values, mom);
    return global;
  }

  // static std::vector<TelSourceLink> CreateSourceLinks(const JsonValue &js, const std::vector<std::shared_ptr<TelElement>> eles);

  std::shared_ptr<altel::TelMeasHit> m_hitMeas;
private:
  Acts::Vector2D m_values;
  size_t m_detId;
  Acts::BoundMatrix m_cov;
  // need to store pointers to make the object copyable
  const Acts::Surface *m_surface;
  friend bool operator==(const TelSourceLink &lhs,
                         const TelSourceLink &rhs) {
    return lhs.m_values.isApprox(rhs.m_values) and
           lhs.m_cov.isApprox(rhs.m_cov);
  }
};

} // namespace Telescope
