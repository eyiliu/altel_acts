#include "TelescopeJsonTrackReader.hpp"

#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <vector>
#include <ratio>
#include <chrono>

#include <Acts/Utilities/Units.hpp>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"

#include "PixelSourceLink.hpp"

Telescope::TelescopeJsonTrackReader::TelescopeJsonTrackReader(
    const Telescope::TelescopeJsonTrackReader::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(0, SIZE_MAX),
      m_logger(Acts::getDefaultLogger("TelescopeJsonTrackReader", lvl))
{
  m_jsa.reset(new Telescope::JsonAllocator);
  m_jsgen.reset(new Telescope::JsonGenerator(m_cfg.inputDataFile));
  m_jsdoc.reset(new Telescope::JsonDocument(m_jsa.get()));

  m_cov_hit = Acts::BoundMatrix::Zero();
  m_cov_hit(0, 0) =m_cfg.resX*m_cfg.resX;
  m_cov_hit(1, 1) =m_cfg.resY*m_cfg.resY;
}

std::string Telescope::TelescopeJsonTrackReader::name() const {
  return "TelescopeJsonTrackReader";
}

std::pair<size_t, size_t> Telescope::TelescopeJsonTrackReader::availableEvents() const {
  return m_eventsRange;
}

FW::ProcessCode Telescope::TelescopeJsonTrackReader::read(const FW::AlgorithmContext& ctx) {
  Telescope::JsonValue evpack;
  {
    const std::lock_guard<std::mutex> lock(m_mtx_read);
    m_jsdoc->Populate(*m_jsgen);
    if(!m_jsgen->isvalid){
      return FW::ProcessCode::ABORT;
    }
    m_jsdoc->Swap(evpack);
  }

  std::vector<Telescope::PixelSourceLink> sourcelinks;
  const auto &layers = evpack["layers"];
  for(const auto&layer : layers.GetArray()){
    size_t id_ext = layer["ext"].GetUint();
    auto surface_it = m_cfg.surfaces.find(id_ext);
    if(surface_it == m_cfg.surfaces.end()){
      continue;
    }
    for(const auto&hit : layer["hit"].GetArray()){
      double x_hit = hit["pos"][0].GetDouble();
      double y_hit = hit["pos"][1].GetDouble();
      Acts::Vector2D loc_hit;
      loc_hit << x_hit, y_hit;
      sourcelinks.emplace_back(*(surface_it->second), loc_hit, m_cov_hit);
    }
  }

  ctx.eventStore.add(m_cfg.outputTracks, std::move(sourcelinks));
  return FW::ProcessCode::SUCCESS;
}
