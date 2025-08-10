#include "TelMille.hh"

#include "getopt.h"

#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"

#include <iostream>
#include <algorithm>
#include <set>

static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -eudaqFiles  <PATH0 [PATH1]...>   paths to input eudaq raw data files (input)
  -maxEventNumber    <int>           max number of events to be processed
  -maxTrackNumber    <int>           max number of tracks to be processed
  -pedeSteeringFile  <PATH>          path to pede steering file (output)
  -milleBinaryFile   <PATH>          path to mille binary file (output)
  -inputGeometryFile <PATH>          geometry input file
  -resolDefault   <  <float_UV>|<float_U float_V> >
                                    default U/V resolution for all detectors
  -resolDetector  <<int_ID>  <<float_UV>|<float_U float_V> >>
                                    U/V resolution(s) for a specific detector by int_ID

example:
./altelMilleBin -pede pede.txt -mille mille.bin  -eudaqFiles  eudaqRaw/altel_Run069017_200824002945.raw eudaqRaw/altel_Run069018_200824003322.raw -input ../init_geo.json -maxE 1000000 -resolDefault 0.04 -resolDet 1 0.1 0.09
)";

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"eudaqFiles", required_argument, NULL, 'f'},
                              {"inputGeometryFile", required_argument, NULL, 'g'},
                              {"pedeSteeringFile", required_argument, NULL, 'u'},
                              {"milleBinaryFile", required_argument, NULL, 'q'},
                              {"resolDefault", required_argument, NULL, 'r'},
                              {"resolDetector", required_argument, NULL, 's'},
                              {"maxEventNumber", required_argument, NULL, 'm'},
                              {"maxTrackNumber", required_argument, NULL, 'n'},
                              {0, 0, 0, 0}};

  std::vector<std::string> rawFilePathCol;
  std::string inputGeometryFile_path;
  std::string pedeSteeringFile_path;
  std::string milleBinaryFile_path;
  size_t maxTrackNumber = -1;
  size_t maxEventNumber = -1;

  std::map<uint16_t, std::pair<double, double>> mapResolDet;
  double resolDefaultU =0.03;
  double resolDefaultV =0.03;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'f':{
      optind--;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        const char* fileStr = argv[optind];
        rawFilePathCol.push_back(std::string(fileStr));
      }
      break;
    }
    case 'g':
      inputGeometryFile_path = optarg;
      break;
    case 'u':
      pedeSteeringFile_path = optarg;
      break;
    case 'q':
      milleBinaryFile_path = optarg;
      break;
    case 'r':{
      optind--;
      std::vector<double> resolVec;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        resolVec.push_back(std::stod(argv[optind]));
      }
      if(resolVec.size()==1){
        resolDefaultU=resolVec[0];
        resolDefaultV=resolVec[0];
      }
      else if(resolVec.size()==2){
        resolDefaultU=resolVec[0];
        resolDefaultV=resolVec[1];
      }
      break;
    }
    case 's':{
      optind--;
      uint16_t detN = -1;
      std::vector<double> resolVec;
      bool isFirstArg = true;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        if(isFirstArg){
          detN = std::stoi(argv[optind]);
          isFirstArg = false;
        }
        else{
          resolVec.push_back(std::stod(argv[optind]));
        }
      }
      if(resolVec.size()==1){
        mapResolDet[detN]=std::pair<double, double>(resolVec[0], resolVec[0]);
      }
      else if(resolVec.size()==2){
        mapResolDet[detN]=std::pair<double, double>(resolVec[0], resolVec[1]);
      }
      break;
    }
    case 'm':
      maxEventNumber = std::stoull(optarg);
      break;
    case 'n':
      maxTrackNumber = std::stoull(optarg);
      break;
      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      fprintf(stderr, "case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr, "case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr, "case ?\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }

  if (rawFilePathCol.empty() ||
      inputGeometryFile_path.empty()||
      milleBinaryFile_path.empty()||
      pedeSteeringFile_path.empty()
    ) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "%d eudaqFiles:\n", rawFilePathCol.size());
  for(auto &rawfilepath: rawFilePathCol){
    std::fprintf(stdout, "  %s\n", rawfilepath.c_str());
  }
  std::fprintf(stdout, "inputGeometryFile:  %s\n", inputGeometryFile_path.c_str());
  std::fprintf(stdout, "milleBinaryFile:    %s\n", milleBinaryFile_path.c_str());
  std::fprintf(stdout, "pedeSteeringFile:   %s\n", pedeSteeringFile_path.c_str());
  std::fprintf(stdout, "resolDefault:       [%f   %f]\n", resolDefaultU, resolDefaultV);
  std::fprintf(stdout, "resolDetector:\n");
  for(auto &[detN, resolUV]: mapResolDet){
    std::fprintf(stdout, "  det #%d:  [%f   %f]\n", detN, resolUV.first, resolUV.second);
  }


  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(inputGeometryFile_path.c_str());
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", inputGeometryFile_path.c_str() );
    throw;
  }
  // JsonUtils::printJsonValue(jsd_geo, true);

  std::set<uint16_t> geoDetNs;
  for(const auto& geoJS: jsd_geo["geometry"]["detectors"].GetArray()){
    geoDetNs.insert(geoJS["id"].GetUint());
  }
  if(geoDetNs.size() != jsd_geo["geometry"]["detectors"].Size()){
    std::fprintf(stderr, "geometory id, something wrong\n");
    throw;
  }
  for(auto &detN: geoDetNs){
    if(mapResolDet.find(detN)==mapResolDet.end()){
      mapResolDet[detN]=std::pair<double, double>(resolDefaultU, resolDefaultV);
    }
  }


  altel::TelMille telmille;
  telmille.setGeometry(jsd_geo);
  for(auto& [detN, resolUV ]: mapResolDet){
    telmille.setResolution(detN, resolUV.first, resolUV.second);
  }
  telmille.startMilleBinary(milleBinaryFile_path);

  JsonAllocator jsa;

  size_t nTracks = 0;
  size_t nEvents = 0;
  uint32_t rawFileN=0;

  eudaq::FileReaderUP reader;

  while(1){
    if (nTracks >= maxTrackNumber || nEvents >= maxEventNumber ) {
      break;
    }

    if(!reader){
      if(rawFileN<rawFilePathCol.size()){
        std::fprintf(stdout, "processing raw file: %s\n", rawFilePathCol[rawFileN].c_str());
        reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), rawFilePathCol[rawFileN]);
        rawFileN++;
      }
      else{
        std::fprintf(stdout, "processed %d raw files, quit\n", rawFileN);
        break;
      }
    }
    auto eudaqEvent = reader->GetNextEvent();
    if(!eudaqEvent){
      reader.reset();
      continue; // goto for next raw file
    }
    nEvents++;

    std::shared_ptr<altel::TelEvent> telEvent = altel::createTelEvent(eudaqEvent);

    // TODO: TelEvent to json
    JsonValue js_track_filtered(rapidjson::kArrayType);
    std::set<uint16_t> measDetNs;
    bool isMoreThanOneHitPerPlane = false;
    std::fprintf(stdout, "\nEvent #%d ", nEvents);
    for(const auto &aMeasHit: telEvent->measHits()){
      uint16_t detN = aMeasHit->detN();
      double measU = aMeasHit->u();
      double measV = aMeasHit->v();

      std::fprintf(stdout, "[%f, %f, %d] ", measU, measV, detN);
      if(geoDetNs.find(detN) == geoDetNs.end()){
        continue;
      }
      auto [it_insert, isSucess] = measDetNs.insert(detN);
      if(!isSucess){
        isMoreThanOneHitPerPlane = true;
        break;
      }
      JsonValue js_hit(rapidjson::kObjectType);
      js_hit.AddMember("id", detN, jsa);
      js_hit.AddMember("x", measU, jsa);
      js_hit.AddMember("y", measV, jsa);
      js_track_filtered.PushBack(std::move(js_hit), jsa);
    }

    if(isMoreThanOneHitPerPlane || measDetNs.size() != geoDetNs.size()){
      std::cout<< "skipping event "<<nEvents <<std::endl;
      continue;
    }

    nTracks++;
    telmille.fillTrackXYRz(js_track_filtered);
  }
  std::fprintf(stdout, "%i tracks are picked from %i events\n", nTracks, nEvents);

  telmille.endMilleBinary();
  telmille.createPedeStreeringModeXYRz(pedeSteeringFile_path);

  return 0;
}
