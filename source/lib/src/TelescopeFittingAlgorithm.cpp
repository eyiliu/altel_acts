// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeFittingAlgorithm.hpp"

#include <stdexcept>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

Telescope::TelescopeFittingAlgorithm::TelescopeFittingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TelescopeFittingAlgorithm", level),
      m_cfg(std::move(cfg)) {

}

FW::ProcessCode Telescope::TelescopeFittingAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  // Read input data
  const auto& trackSourcelinks
    = ctx.eventStore.get<std::vector<Telescope::PixelSourceLink>>(m_cfg.inputSourcelinks);

  // Prepare the output data with MultiTrajectory
  std::vector<PixelMultiTrajectory> trajectories;

  // We can have empty tracks which must give empty fit results
  if (!trackSourcelinks.empty()) {
    // Set initial parameters for the particle track

    // @Below is what is used when the detector aligned along global x in the first beginning. 
    // But this is not working when it's aligned along global z. Needs investigation of the reason.
    //Acts::BoundSymMatrix cov;
    //cov << std::pow(10_mm, 2), 0., 0., 0., 0., 0., 0., std::pow(10_mm, 2), 0.,
    //    0., 0., 0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 0.0001, 0., 0., 0.,
    //    0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 1.;
    //Acts::Vector3D rPos(-120_mm, 0, 0);
    //Acts::Vector3D rMom(4_GeV, 0, 0);
    double resX2 = m_cfg.seedResX * m_cfg.seedResX;
    double resY2 = m_cfg.seedResY * m_cfg.seedResY;
    double resPhi2 = m_cfg.seedResPhi * m_cfg.seedResPhi;
    double resTheta2 = m_cfg.seedResTheta * m_cfg.seedResTheta;

    Acts::BoundSymMatrix cov;
    cov <<
      resX2,0., 0., 0., 0., 0.,
      0., resY2, 0., 0., 0., 0.,
      0., 0., resPhi2, 0., 0., 0.,
      0., 0., 0., resTheta2, 0., 0.,
      0., 0., 0., 0., 0.0001, 0.,
      0., 0., 0., 0.,     0., 1.;

    const Acts::Vector3D global0 = trackSourcelinks.at(0).globalPosition(ctx.geoContext);
    const Acts::Vector3D global1 = trackSourcelinks.at(1).globalPosition(ctx.geoContext);
    Acts::Vector3D distance = global1 - global0;
    const double phi = Acts::VectorHelpers::phi(distance);
    const double theta = Acts::VectorHelpers::theta(distance);
    Acts::Vector3D rPos = global0 - distance / 2;
    Acts::Vector3D rMom(m_cfg.beamEnergy * sin(theta) * cos(phi),
                        m_cfg.beamEnergy * sin(theta) * sin(phi),
                        m_cfg.beamEnergy * cos(theta));

    Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(cov, rPos, rMom, 1., 0);

    // Set the KalmanFitter options
    // TODO, remove target surface, note the vector3d
    auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>
      (Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions
      (ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
       Acts::VoidOutlierFinder(), &(*pSurface));

    ACTS_DEBUG("Invoke fitter");
    auto result = m_cfg.fit(trackSourcelinks, rStart, kfOptions);
    if (result.ok()) {
      // Get the fit output object
      const auto& fitOutput = result.value();
      // The track entry indices container. One element here.
      std::vector<size_t> trackTips;
      trackTips.reserve(1);
      trackTips.emplace_back(fitOutput.trackTip);
      // The fitted parameters container. One element (at most) here.
      IndexedParams indexedParams;
      if (fitOutput.fittedParameters) {
        const auto& params = fitOutput.fittedParameters.value();
        ACTS_VERBOSE("Fitted paramemeters for track " << ctx.eventNumber);
        ACTS_VERBOSE("  position: " << params.position().transpose());
        ACTS_VERBOSE("  momentum: " << params.momentum().transpose());
        // Push the fitted parameters to the container
        indexedParams.emplace(fitOutput.trackTip, std::move(params));
      } else {
        ACTS_DEBUG("No fitted paramemeters for track " << ctx.eventNumber);
      }
      // Create a PixelMultiTrajectory
      trajectories.emplace_back(std::move(fitOutput.fittedStates),
                                std::move(trackTips), std::move(indexedParams));
    } else {
      ACTS_WARNING("Fit failed for track " << ctx.eventNumber << " with error"
                   << result.error());
      trajectories.push_back(PixelMultiTrajectory());
    }
  }
  else{
    ACTS_WARNING("Empty event " << ctx.eventNumber << " found.");
    trajectories.push_back(PixelMultiTrajectory());
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}
