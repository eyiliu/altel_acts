#include "getopt.h"
#include "mysystem.hh"
#include "myrapidjson.h"

#include <iostream>
#include <algorithm>
#include <regex>


static const std::string help_usage = R"(
Usage:
  -help                    help message
  -verbose                 verbose flag
  -outputGeometry   [PATH]   alignment result file
  -inputGeometry    [PATH]   geometry input file
  -pedeResultFile   [PATH]   path to pede result file (input)

examples:
./altelGeoUpdate -p millepede.res -i oldGeo.json -o newGeo.json
)";

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"inputGeometry", required_argument, NULL, 'g'},
                              {"outputGeometry", required_argument, NULL, 'o'},
                              {"pedeReaultFile", required_argument, NULL, 'u'},
                             {0, 0, 0, 0}};

  std::string outputGeometry_path;
  std::string inputGeometry_path;
  std::string pedeResultFile_path;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'o':
      outputGeometry_path = optarg;
      break;
    case 'g':
      inputGeometry_path = optarg;
      break;
    case 'u':
      pedeResultFile_path = optarg;
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

  if (inputGeometry_path.empty() ||
      outputGeometry_path.empty()||
      pedeResultFile_path.empty()
    ) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "outputGeometry:     %s\n", outputGeometry_path.c_str());
  std::fprintf(stdout, "inputGeometry:      %s\n", inputGeometry_path.c_str());
  std::fprintf(stdout, "pedeResultFile:     %s\n", pedeResultFile_path.c_str());
  std::fprintf(stdout, "\n");

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(inputGeometry_path.c_str());
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", inputGeometry_path.c_str() );
    throw;
  }
  if( !jsd_geo.HasMember("geometry") ||
      !jsd_geo["geometry"].HasMember("detectors") ||
      !jsd_geo["geometry"]["detectors"].IsArray()){
    std::fprintf(stderr, "Geometry file <%s>: unknown json structure.\n", inputGeometry_path.c_str() );
    JsonUtils::printJsonValue(jsd_geo, true);
    throw;
  }
  auto &js_geo = jsd_geo["geometry"];
  auto &js_dets = js_geo["detectors"];
  for(auto& js_det : js_dets.GetArray()){
    JsonUtils::printJsonValue(js_det, false);
  }

  std::string pederes = JsonUtils::readFile(pedeResultFile_path);
  if(pederes.empty()){
    std::fprintf(stderr, "nothing read from file <%s>.\n", pedeResultFile_path.c_str() );
    throw;
  }
  std::istringstream iss(pederes);

  std::regex e("\\s*([0-9]+)\\s+(-?\\d+\\.\\d+(?:[eE]-?\\d+)?)\\s+(-?\\d+\\.\\d+)\\s+(?:.+)");

  std::string aline;
  while (std::getline(iss, aline)){
    std::smatch mt;
    if(std::regex_match(aline, mt, e)){
      std::string label_str = mt[1].str();
      std::string result_str= mt[2].str();
      size_t label = std::stoul(label_str);
      double result = std::stod(result_str);
      uint32_t id = label/10;
      uint32_t index = label%10;
      bool found= false;
      for(auto& js_det : js_dets.GetArray()){
        size_t the_id = js_det["id"].GetUint();
        if(the_id != id){
          continue;
        }
        found = true;
        if(index==1){
          js_det["center"]["x"]=result+js_det["center"]["x"].GetDouble();
        }
        else if(index==2){
          js_det["center"]["y"]=result+js_det["center"]["y"].GetDouble();
        }
        else if(index==3){
          js_det["rotation"]["z"]=result + js_det["rotation"]["z"].GetDouble();
        }
        else if(index==4){
          js_det["rotation"]["x"]=result + js_det["rotation"]["x"].GetDouble();
        }

        else{
          std::cerr<< "something wrong, index is not found\n";
          throw;
        }
      }
      if(!found){
        std::cerr<< "something wrong, id is not found\n";
        throw;
      }
    }
  }

  std::fprintf(stdout, "\n==== updated geo=======\n");
  for(auto& js_det : js_dets.GetArray()){
    JsonUtils::printJsonValue(js_det, false);
  }
  std::fprintf(stdout, "=======================\n");

  // JsonUtils::printJsonValue(jsd_geo, false);
  std::string jsstr = JsonUtils::stringJsonValue(jsd_geo, true);
  std::FILE *fp = std::fopen(outputGeometry_path.c_str(), "w");
  std::fwrite(jsstr.data(), 1, jsstr.size(), fp);
  std::fclose(fp);
  return 0;
}
