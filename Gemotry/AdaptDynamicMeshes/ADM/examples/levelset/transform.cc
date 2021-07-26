#include <tml.h>

int main(int argc, char** argv)
{
  if (argc!=7)
    { std::cerr << "./transform <input> <output> <scale> <tx> <ty> <tz>\n"; std::exit(1); }

  const char* input  = argv[1];
  const char* output = argv[2];
  double scale = atof(argv[3]);
  TML::R3 t( std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6]) );
  
  TML::TriMeshIO m;
  m.Read(input);
  
  std::vector<TML::R3>& vrts = m.get_vertices();
  std::vector<TML::R3>::iterator vi;
  for (vi=vrts.begin(); vi!=vrts.end(); vi++)
    (*vi) = scale*(*vi) + t;
  
  m.Write(output);

  return 1;
}
