#include "tml.h"

using namespace TML;

bool TriMeshIO::is_triangulated()
{
  std::vector<Polygon>::iterator fi;
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    if (fi->size()!=3) return false;
  return true;
}

bool TriMeshIO::invert()
{
  int i, j;
  std::vector<Polygon>::iterator fi;
  
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    {
      i = 0; j = 1;
      while (i<j)
	{
	  std::swap(fi->vertex(i),fi->vertex(j));
	  i = (i==0)? fi->size()-1 : i-1;
	  j = (j+1)%fi->size();
	}
    }
  return true;
}

bool TriMeshIO::quad2tri()
{
  std::vector<Polygon>::iterator fi;
  
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    if (fi->size()!=4) return false;
  
  int i, nq = facets.size();
  for (i=0; i<nq; i++)
    {
      fi = facets.begin();
      
      Polygon tri1(3);
      tri1.vertex( 0, fi->vertex(2) );
      tri1.vertex( 1, fi->vertex(0) );
      tri1.vertex( 2, fi->vertex(1) );
      
      Polygon tri2(3);
      tri2.vertex( 0, fi->vertex(0) );
      tri2.vertex( 1, fi->vertex(2) );
      tri2.vertex( 2, fi->vertex(3) );
      
      facets.erase(fi);
      facets.push_back(tri1);
      facets.push_back(tri2);
    }
  return true;
}

bool TriMeshIO::normalize_bb()
{
  std::vector<R3>::iterator vi;
  
  double xmax, ymax, zmax, xmin, ymin, zmin;
  vi = vertices.begin();
  xmax = xmin = vi->x;
  ymax = ymin = vi->y;
  zmax = zmin = vi->z;
  for ( ; vi!=vertices.end(); vi++)
    {
      xmax = std::max(xmax,vi->x);
      xmin = std::min(xmin,vi->x);
      ymax = std::max(ymax,vi->y);
      ymin = std::min(ymin,vi->y);
      zmax = std::max(zmax,vi->z);
      zmin = std::min(zmin,vi->z);
    }

  R3 c( 0.5*(xmax+xmin), 0.5*(ymax+ymin), 0.5*(zmax+zmin) ); 
  double scale = 2. / std::max( (xmax-xmin), std::max( (ymax-ymin), (zmax-zmin) ) );
  for (vi=vertices.begin(); vi!=vertices.end(); vi++)
    *vi = scale * ( (*vi) - c );
  
  return true;
}

bool TriMeshIO::check_sanity()
{
  int i, num_verts = vertices.size();
  std::vector<Polygon>::iterator it;
  for (it=facets.begin(); it!=facets.end(); it++)
    for (i=0; i<it->size(); i++)
      if (it->vertex(i)<0 || it->vertex(i)>=num_verts) 
	{ 
	  std::cout << "TriMeshIO::check_sanity: mesh not valid."<< std::endl;
	  return false;
	}
  return true;
}

int TriMeshIO::File_Format(const char* filename)
{
  std::string str(filename);

  if ( (str.rfind(".off") != std::string::npos) ||
       (str.rfind(".OFF") != std::string::npos)   ) return OFF;
  
  if ( (str.rfind(".smf") != std::string::npos) ||
       (str.rfind(".SMF") != std::string::npos)   ) return SMF;
  
  if ( (str.rfind(".ply") != std::string::npos) ||
       (str.rfind(".PLY") != std::string::npos)   ) return PLY;

  if ( (str.rfind(".ifs") != std::string::npos) ||
       (str.rfind(".IFS") != std::string::npos)   ) return IFS;

  if ( (str.rfind(".obj") != std::string::npos) ||
       (str.rfind(".OBJ") != std::string::npos)   ) return OBJ;
  
  return NONE;
}
  
////////////////////////////////////////////////////////////
bool TriMeshIO::Read_smf(const char* filename)
{
  char buf[256], c;
  
  std::ifstream input(filename);
  
  while (!input.eof())
    {
      input >> c;
      switch (c)
	{
	case 'v':
	  {
	    R3 v;
	    input >> v;
	    vertices.push_back(v);
	  } break;
	case 'f': case 't':
	  {
	    Polygon f(3);
	    for (int i=0; i<3; i++)
	      { input >> f.vertex(i); f.vertex(i)--; }
	    facets.push_back(f);
	  } break;
	case 'q':
	  {
	    Polygon f(4);
	    for (int i=0; i<4; i++)
	      { input >> f.vertex(i); f.vertex(i)--; }
	    facets.push_back(f);
	  } break;
	default:
	  input.getline(buf,256);
	  break;
	};
    }  
  
  input.close();
  return check_sanity();
}

bool TriMeshIO::Write_smf(const char* filename)
{
  std::vector<Polygon>::iterator fi;
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    if (fi->size()>4) return false;

  std::ofstream output(filename);
  output << "begin" << std::endl;
  
  std::vector<R3>::iterator vi;
  for (vi=vertices.begin(); vi!=vertices.end(); vi++)
    output << "v " << *vi << std::endl;
  
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    {
      if (fi->size()==3)      output << "f ";
      else if (fi->size()==4) output << "q ";
      
      for (int i=0; i<fi->size(); i++)
	output << ( fi->vertex(i)+1 ) << " ";
      output << std::endl;
    }
  
  output << "end" << std::endl;
  output.close();
  return true;
}
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
bool TriMeshIO::Read_off(const char* filename)
{
  char buf[256];
  int  nv, nf, n;
  
  std::ifstream input(filename);
  input >> buf;
  input >> nv >> nf >> n;
  
  for (int i=0; i<nv; i++)
    { R3 v; input >> v; vertices.push_back( v ); }

  for (int i=0; i<nf; i++)
    {
      input >> n;
      Polygon f(n);
      for (int j=0; j<n; j++) 
	input >> f.vertex(j); 
      facets.push_back(f);
    }
  
  input.close();
  return check_sanity();
}

bool TriMeshIO::Write_off(const char* filename)
{
  std::ofstream output(filename);
  output << "OFF" << std::endl;
  output << vertices.size() << " " 
	 << facets.size()   << " " 
	 << "0"             << std::endl;

  std::vector<R3>::iterator vi;
  for (vi=vertices.begin(); vi!=vertices.end(); vi++)
    output << *vi << std::endl;
  
  std::vector<Polygon>::iterator fi;
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    {
      output << fi->size() << " ";
      for (int i=0; i<fi->size(); i++)
	output << fi->vertex(i) << " ";
      output << std::endl;
    }
  
  output.close();
  return true;
}
////////////////////////////////////////////////////////////
