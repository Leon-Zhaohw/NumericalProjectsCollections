#include <tml.h>
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char** argv)
{
  if (argc!=3)
    throw TML::Error("USAGE: ./3ds2off <input> <output>");

  int i, pos, nv, nf;
  TML::TriMeshIO mesh;

  std::ifstream input(argv[1]);
  std::string str;
  TML::R3 p;
  
  std::getline( input, str); // Dumping meshes:
  
  // base ...
  std::getline( input, str); 

  pos = str.find_first_of("=");
  sscanf( str.substr(pos+1).c_str(), "%d", &nv );
  
  pos = str.find_first_of("=",pos+1);
  sscanf( str.substr(pos+1).c_str(), "%d", &nf );
  std::cout << nv << " " << nf << std::endl;
  
  // matrix
  std::getline( input, str);
  std::getline( input, str);
  std::getline( input, str);
  std::getline( input, str);
  std::getline( input, str);
  
  // point list
  std::getline( input, str);
  for (i=0; i<nv; i++)
    {
      input >> p;
      mesh.add(p);
      std::getline( input, str);
    }

  // face list
  std::getline( input, str );
  for (i=0; i<nf; i++)
    {
      std::getline( input, str );
      pos = str.find_first_of("0123456789");

      TML::Polygon t(3);
      sscanf( str.substr(pos).c_str(), "%d %d %d",
	      &(t.vertex(0)), &(t.vertex(1)), &(t.vertex(2)) );
      mesh.add(t);
    }
  
    input.close();

  return mesh.Write(argv[2]);
}
