// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/ConstantBField.hpp"
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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

#include "TelescopeDetectorElement.hpp"
#include "TelescopeTrackReader.hpp"
#include "TelescopeTrack.hpp"

#include <memory>
#include "getopt.h"
#include "myrapidjson.h"

using JsonValue = rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;
using JsonDocument = rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -file [jsonfile]   name of data json file
  -gev [energy_gev]  beam energy
  -out [jsonfile]    alignment result

)";

int main(int argc, char* argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",       no_argument,       &do_help,      1  },
     { "verbose",    no_argument,       &do_verbose,   1  },
     { "file",      required_argument, NULL,           'f' },
     { "gev",       required_argument, NULL,           'e' },
     { "out",       required_argument, NULL,           'o' },
     { "geo",       required_argument, NULL,           'g' },
     { 0, 0, 0, 0 }};
  
  std::string datafile_name;
  std::string outputfile_name;
  std::string geofile_name;
  uint32_t energy_opt = 4;
  
  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'h':
      do_help = 1;
      std::fprintf(stdout, "%s\n", help_usage.c_str());
      exit(0);
      break;
    case 'f':
      datafile_name = optarg;
      break;
    case 'e':
      energy_opt = std::stoul(optarg);
      break;
    case 'o':
      outputfile_name = optarg;
      break;
    case 'g':
      geofile_name = optarg;
      break;

      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      fprintf(stderr,"case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr,"case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr,"case ?\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }

  if(datafile_name.empty() || outputfile_name.empty()){
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }


  ///////select datapacks/////////////////////
  JsonValue js_selected_datapack_col(rapidjson::kArrayType);
  size_t n_datapack_select_opt = 20000;
  {
    std::FILE* fp = std::fopen(datafile_name.c_str(), "r");
    if(!fp) {
      std::fprintf(stderr, "File opening failed: %s \n", datafile_name.c_str());
      throw std::system_error(EIO, std::generic_category(), "File opening failed");;
    }

    char readBuffer[UINT16_MAX];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>  doc;
    doc.ParseStream(is);
    std::fclose(fp);

    if(!doc.IsArray() || !doc.GetArray().Size()){
      std::fprintf(stderr, "no, it is not data array\n");
      throw std::system_error(EDOM, std::generic_category(), "File is not valid json array");;
    }

    uint64_t processed_datapack_count = 0;
    for(auto ev_it = doc.Begin(); ev_it != doc.End() && js_selected_datapack_col.Size() < n_datapack_select_opt ; ev_it++ ){
      auto &evpack = *ev_it;
      processed_datapack_count ++;
      auto &frames = evpack["layers"];
      bool is_good_datapack = true;
      for(auto& layer : evpack["layers"].GetArray()){
        uint64_t l_hit_n = layer["hit"].GetArray().Size();
        if(l_hit_n != 1){
          is_good_datapack = false;
          continue;
        }
      }
      if(!is_good_datapack){
        continue;
      }
      js_selected_datapack_col.PushBack(evpack, jsa);
    }

    std::fprintf(stdout,
                 "Select %lu datapacks out of %lu processed datapackes. The data file contains %lu datapacks.\n",
                 js_selected_datapack_col.Size(), processed_datapack_count, doc.GetArray().Size() );

  }
  ///////end of select datapacks/////////////////////

  JsonValue js_geometry(rapidjson::kArrayType);
  {
    std::FILE* fp = std::fopen(geofile_name.c_str(), "r");
    if(!fp) {
      std::fprintf(stderr, "File opening failed: %s \n", geofile_name.c_str());
      throw std::system_error(EIO, std::generic_category(), "File opening failed");;
    }
    char readBuffer[UINT16_MAX];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>  doc;
    doc.ParseStream(is);
    std::fclose(fp);
    js_geometry.CopyFrom<rapidjson::CrtAllocator>(doc["alignment_result"],jsa);
  }

  double beamEnergy = energy_opt * Acts::UnitConstants::GeV;
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  // Setup detector geometry
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> element_col;
  std::vector<std::shared_ptr<const Acts::Surface>> surface_col;
  std::vector<Acts::LayerPtr> layer_col;
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
                                           // {-95_mm, 0., 0.},
                                           // {-57_mm, 0., 0.},
                                           // {-19_mm, 0., 0.},
                                           // {19_mm, 0., 0.},
                                           // {57_mm, 0., 0.},
                                           // {95_mm, 0., 0.}};


  Telescope::BuildGeometry(gctx, trackingGeometry, element_col, surface_col, layer_col, js_geometry);
  // Set up the detector elements to be aligned (fix the first one)
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> dets;
  unsigned int idet = 0;
  for (const auto& det : element_col) {
    idet++;
    // Skip the first detector element
    if (idet == 1) {
      continue;
    }
    dets.push_back(det);
  }
  //



  /////////////////////////////////////

  // The criteria to determine if the iteration has converged. @Todo: to use
  // delta chi2 instead
  double chi2ONdfCutOff = 0.01;
  std::pair<size_t, double> deltaChi2ONdfCutOff = {10,0.000001};
  // The maximum number of iterations
  size_t maxNumIterations = 800;
  // set up the alignment dnf for each iteration
  std::map<unsigned int, std::bitset<6>> iterationState;
  for (unsigned int iIter = 0; iIter < maxNumIterations; iIter++) {
    std::bitset<6> mask(std::string("111111"));
    if (iIter % 3 == 0) {
      mask = std::bitset<6>(std::string("000110"));
    } else if (iIter % 3 == 1) {
      mask = std::bitset<6>(std::string("010000"));
    } else {
      mask = std::bitset<6>(std::string("101000"));
    }
    iterationState.emplace(iIter, mask);
  }

  FW::AlignedTransformUpdater alignedTransformUpdaterFun =
    [](Acts::DetectorElementBase& detElement,
       const Acts::GeometryContext& gctx,
       const Acts::Transform3D& aTransform) {
      Telescope::TelescopeDetectorElement* telescopeDetElement =
        dynamic_cast<Telescope::TelescopeDetectorElement*>(&detElement);
      if (telescopeDetElement) {
        telescopeDetElement->addAlignedTransform
          (std::make_unique<Acts::Transform3D>(aTransform));
        return true;
      }
      return false;
    };
  
  auto alignFun = Telescope::makeAlignmentFunction
    (trackingGeometry, magneticField, Acts::Logging::INFO);

  
  while(1){
    // Setup local covariance
    Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
    cov_hit(0, 0) =150_um * 150_um;
    cov_hit(1, 1) =150_um * 150_um;

    std::vector<std::vector<Telescope::PixelSourceLink>> sourcelinkTracks;
    sourcelinkTracks.reserve(js_selected_datapack_col.Size());
    for(auto ev_it = js_selected_datapack_col.Begin(); ev_it != js_selected_datapack_col.End() ; ev_it++ ){
      const auto &evpack = *ev_it;
      std::vector<Telescope::PixelSourceLink> sourcelinks;
      auto &frames = evpack["layers"];
      for(size_t i= 0; i< 6; i++){
        double x = frames[i]["hit"][0]["pos"][0].GetDouble() - 0.02924*1024/2.0;
        double y = frames[i]["hit"][0]["pos"][1].GetDouble() - 0.02688*512/2.0;
        Acts::Vector2D loc;
        loc << x, y;
        sourcelinks.emplace_back(*surface_col.at(i), loc, cov_hit);
      }
      sourcelinkTracks.push_back(sourcelinks);
    }
    
    // Prepare the initial track parameters collection
    Acts::BoundSymMatrix cov_seed;
    cov_seed <<
        50_um * 50_um,     0.,  0.,   0.,   0.,  0.,
        0., 50_um * 50_um,      0.,   0.,   0.,  0.,
        0.,                 0., 0.1, 0.,   0.,  0.,
        0.,                 0., 0.,   0.1, 0.,  0.,
        0.,                 0., 0.,   0.,   0.001,  0.,
        0.,                 0., 0.,   0.,   0.,  1.;
    
    std::vector<Acts::CurvilinearParameters> initialParameters;
    initialParameters.reserve(js_selected_datapack_col.Size());
    for( const auto& sourcelinks : sourcelinkTracks) {
      // Create initial parameters
      // The position is taken from the first measurement
      const Acts::Vector3D global0 = sourcelinks.at(0).globalPosition(gctx);
      const Acts::Vector3D global1 = sourcelinks.at(1).globalPosition(gctx);
      Acts::Vector3D distance = global1 - global0;

      const double phi = Acts::VectorHelpers::phi(distance);
      const double theta = Acts::VectorHelpers::theta(distance);
      // shift along the beam by 100_mm
      Acts::Vector3D rPos = global0 - distance / 2;
      Acts::Vector3D rMom(beamEnergy * sin(theta) * cos(phi),
                          beamEnergy * sin(theta) * sin(phi),
                          beamEnergy * cos(theta));
      
      Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(cov_seed, rPos, rMom, 1., 0);
      initialParameters.push_back(rStart);
    }

    // Set the KalmanFitter options
    auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
                                                                    Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});
    // Acts::PlaneSurface refSurface(Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions
      (gctx, mctx, cctx, Acts::VoidOutlierFinder(), refSurface.get()); //pSurface default nullptr

    // Set the alignment options
    FW::AlignmentOptions<Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>> alignOptions
      (kfOptions,
       alignedTransformUpdaterFun,
       dets,
       chi2ONdfCutOff,
       deltaChi2ONdfCutOff,
       maxNumIterations,
       iterationState
       );
    
    std::printf("Invoke alignment");
    auto result = alignFun(sourcelinkTracks, initialParameters, alignOptions);

    if (result.ok()) {
      std::cout<<"Alignment finished with deltaChi2 = " << result.value().deltaChi2;
      // Print out the alignment parameters of all detector elements (including those not aligned)
      idet=0;
      std::cout<<"iDet, centerX, centerY, centerZ, rotX, rotY, rotZ"<<std::endl;

      JsonValue js_output(rapidjson::kObjectType);
      JsonValue js_align_result(rapidjson::kArrayType);

      for (const auto& det : element_col) {
        JsonValue js_ele(rapidjson::kObjectType);
        const auto& surface = &det->surface();
        const auto& transform =
            det->transform(gctx);
        const auto& translation = transform.translation();
        const auto& rotation = transform.rotation();
        const Acts::Vector3D rotAngles = rotation.eulerAngles(2, 1, 0);
      //std::cout<<"Rotation marix = \n" << rotation<<std::endl;
        js_ele.AddMember("idet",       JsonValue(idet), jsa);
        js_ele.AddMember("centerX",    JsonValue(translation.x()), jsa);
        js_ele.AddMember("centerY",    JsonValue(translation.y()), jsa);
        js_ele.AddMember("centerZ",    JsonValue(translation.z()), jsa);
        js_ele.AddMember("rotX", JsonValue(rotAngles(2)), jsa);
        js_ele.AddMember("rotY", JsonValue(rotAngles(1)), jsa);
        js_ele.AddMember("rotZ", JsonValue(rotAngles(0)), jsa);
        js_align_result.PushBack(std::move(js_ele), jsa);

        std::cout<<idet<<","<<translation.x()<<","<<translation.y()<<","<<translation.z()<<","<<rotAngles(2)<<","<<rotAngles(1)<<","<<rotAngles(0)<<std::endl;
        idet++;
      }

      js_output.AddMember("alignment_result", std::move(js_align_result), jsa);
      std::FILE* fp = std::fopen(outputfile_name.c_str(), "w");
      char writeBuffer[UINT16_MAX];
      rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));
      rapidjson::PrettyWriter< rapidjson::FileWriteStream> writer(os);
      js_output.Accept(writer);
      std::fclose(fp);

    } else {
      std::cout<<"Alignment failed with " << result.error();
    }

    /////////////////////////////
    break;
  }
  return 0;
}
