#include "tml.h"

using namespace TML;

bool TriMeshIO::Read_obj(const char* filename)
{
  std::string str;
  std::ifstream input(filename);
  
  while (!input.eof())
    {
      input >> str;
      
      if ( str.compare("v")==0 )
        {
	  R3 v;
	  input >> v;
	  vertices.push_back(v);
	  continue;
	} 
      
      if ( str.compare("f")==0 )
	{
	  Polygon f(3);
	  for (int i=0; i<3; i++)
	    { 
	      input >> f.vertex(i); f.vertex(i)--; // ind[i]
	      input >> str;                        // //ind[i]
	    }
	  facets.push_back(f);
	  continue;
	}
      
      getline( input, str );
    }
  
  input.close();
  return check_sanity();
}

bool TriMeshIO::Write_obj(const char* filename)
{
  std::vector<Polygon>::iterator fi;
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    if (fi->size()>3) return false;

  std::ofstream output(filename);
  
  std::vector<R3>::iterator vi;
  for (vi=vertices.begin(); vi!=vertices.end(); vi++)
    output << "v " << *vi << std::endl;
  
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    {
      output << "f ";
      for (int i=0; i<fi->size(); i++)
        output << ( fi->vertex(i)+1 ) << "\\" << "\\" << ( fi->vertex(i)+1 ) << " ";
      output << std::endl;
    }
  
  output.close();
  return true;
}
