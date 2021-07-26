//#define NO_GRAPHICS

#ifndef NO_GRAPHICS
#include "Viewer.hpp"

MAGNUM_APPLICATION_MAIN(Magnum::Viewer);

#else
#include "WaterSurface.hpp"
#include <iostream>

void help_main() {
  std::cout<<"\n     *** WAVE: Help ***\n"<<std::endl;
  std::cout<<"Synopsis: \n     .\\main <options>\n\nOptions:"<<std::endl;
  std::cout<<"     -l, -load <file>: load configuration file"<<std::endl;
  std::cout<<"     -e, -export <name>: export amplitude grid for each frequencies in the file <name>.obj"<<std::endl;
  std::cout<<"     -i, -import <name>: import amplitude grid for each frequencies in the file <name>.obj"<<std::endl;
  std::cout<<"     -stop <t>: stop animation and exit at time t"<<std::endl;
  std::cout<<"     -em <name>: export heightfields and render files in a set of files <name><frame number>.ong and <name><frame number>.xml"<<std::endl;
  std::cout<<"     -es, -export_step <n>: export every n frames"<<std::endl;
  std::cout<<"     -h, -help: print help\n"<<std::endl;
  exit(0);
}

int main(int argc, char **argv) {
  WaterSurface _surface;
  uint stop_time = 1e6;
  
  for (int i = 1;  i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "-l" || s == "-load") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Loading configuration file:"<<" "<<argv[i+1]<<std::endl;
      _surface.setImportConf(argv[i+1]);
      ++i;
    } else if (s == "-i" || s == "-import") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Importing"<<" "<<argv[i+1]<<std::endl;
      _surface.setImport(argv[i+1]);
      ++i;
    } else if (s == "-d" || s == "-data") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Saving amplitude data in"<<" "<<argv[i+1]<<std::endl;
      _surface.setData(argv[i+1]);
      ++i;
    } else if (s == "-e" || s == "-export") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Exporting"<<" "<<argv[i+1]<<std::endl;
      _surface.setExport(argv[i+1]);
      ++i;
    } else if (s == "-em") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Exporting (mitsuba) "<<" "<<argv[i+1]<<std::endl;
      _surface.setExportMitsuba(argv[i+1]);
      ++i;
    } else if (s == "-stop") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Stop at t = "<<argv[i+1]<<std::endl;
      stop_time = atoi(argv[i+1]);
      ++i;
    } else if (s == "-export_step" || s == "-es") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_main();
      }
      std::cout<<"Export every "<<argv[i+1]<<" steps"<<std::endl;
      _surface.setExportStep(atoi(argv[i+1]));
      ++i;
    } else if (s == "-h" || s == "-help") {
      help_main();
    } else {
      std::cerr<<"\nERROR: Unknown option\n"<<std::endl;
      help_main();
    }
  }
  _surface.reset();

  uint time = 0;
  while (time < stop_time) {
    _surface.update();
    ++time;
  }

  return 0;
}

#endif
