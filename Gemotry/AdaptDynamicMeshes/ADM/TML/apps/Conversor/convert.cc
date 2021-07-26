#ifndef _CONVERT_
#define _CONVERT_ 1

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>

// #define NDEBUG 1
#include <cassert>

#include <tml.h>

#define MESH_INVERT  1
#define MESH_TRIQUAD 2
#define MESH_BBOX    4

void help();
int parse(int argc, char** argv, std::string& input, std::string& output);

int main(int argc, char** argv)
{
  int   options;
  std::string input, output;
  
  options = parse(argc,argv,input,output);

  TML::TriMeshIO m;
  
  assert( m.Read(input.c_str()) );
  
  if (options & MESH_INVERT ) assert( m.invert() );
  if (options & MESH_TRIQUAD) assert( m.quad2tri() );
  if (options & MESH_BBOX   ) assert( m.normalize_bb() );
  
  assert( m.Write(output.c_str()) );
  
  return 1;
}

////////////////////////////////////////////////////////////////////////////
void help()
{
  std::cerr << "Usage:" << std::endl;
  std::cerr << "\t ./meshconvert "
	    << "-i <input> "
	    << "-o <output> "
	    << "[-tq] [-inv] [-bb]" 
	    << std::endl;
  std::cerr << "where:" << std::endl
	    << "\t -i   - input file" << std::endl
	    << "\t -o   - output file" << std::endl
	    << "\t -tq  - convert quad mesh to tri mesh" << std::endl
	    << "\t -inv - change facets orientation" << std::endl
	    << "\t -bb  - normalize mesh size to the bounding box [-1,1]^3" << std::endl;
}

int parse(int argc, char** argv, std::string& input, std::string& output)
{
  int options = 0;
  int iflag=0, oflag=0;
  
  if (argc==1) { help(); assert(false); }
  
  for (int i=1; i<argc; i++)
    {
      std::string str(argv[i]);

      if ( str[0]=='-' )
	{
	  assert(iflag!=1 && oflag!=1);
	  
	  if ( str.compare("-i")==0 )
	    { iflag = 1; continue; }
 
	  if ( str.compare("-o")==0 )
	    { oflag = 1; continue; }
	  
	  if ( str.compare("-tq")==0 )
	    { options |= MESH_TRIQUAD; continue; }
	  
	  if ( str.compare("-inv")==0 )
	    { options |= MESH_INVERT; continue; }
	
	  if ( str.compare("-bb")==0 )
	    { options |= MESH_BBOX; continue; }
	  
	  help(); 
	  assert(false);
	} else {

	  if (iflag==1)
	    { input.assign(argv[i]);  iflag = 2; continue; }
	  
	  if (oflag==1)
	    { output.assign(argv[i]); oflag = 2; continue; }
	  
	  help(); 
	  assert(false);
	} 
    }

  if (iflag!=2 || oflag!=2)
    { help(); assert(false); }

  //std::cout << "FINAL: " << input << " " << output << std::endl;

  return options;
}
////////////////////////////////////////////////////////////////////////////

#endif
