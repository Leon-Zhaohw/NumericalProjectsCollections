#include "tml.h"

using namespace TML;

bool TriMeshIO::Read_ifs(const char* filename)
{
  char buf[sizeof(float)], str[256];
  int  len, nv, nt;

  std::ifstream input(filename);
  
  // FILE HEADER
  input.read(buf,sizeof(float));
  len = *((unsigned int*)buf);
  input.read(str,len);

  input.read(buf,sizeof(float));

  input.read(buf,sizeof(float));
  len = *((unsigned int*)buf);
  input.read(str,len);

  // VERTEX HEADER
  input.read(buf,sizeof(float));
  len = *((unsigned int*)buf);
  input.read(str,len);

  input.read(buf,sizeof(float));
  nv = *((unsigned int*)buf);

  // VERTEX LIST
  for (int i=0; i<nv; i++)
    { 
      R3 v;
      input.read(buf,sizeof(float));
      v.x = *((float*)buf);
      input.read(buf,sizeof(float));
      v.y = *((float*)buf);
      input.read(buf,sizeof(float));
      v.z = *((float*)buf);
      vertices.push_back(v);
    }

  // TRI HEADER
  input.read(buf,sizeof(float));
  len = *((unsigned int*)buf);
  input.read(str,len);

  input.read(buf,sizeof(float));
  nt = *((unsigned int*)buf);

  // TRI LIST
  for (int i=0; i<nt; i++)
    { 
      Polygon f(3);
      input.read(buf,sizeof(float));
      f.vertex(0) = *((unsigned int*)buf);
      input.read(buf,sizeof(float));
      f.vertex(1) = *((unsigned int*)buf);
      input.read(buf,sizeof(float));
      f.vertex(2) = *((unsigned int*)buf);
      facets.push_back(f);
    }
  
  input.close();
  return check_sanity();
}

bool TriMeshIO::Write_ifs(const char* filename)
{
  std::vector<Polygon>::iterator fi;
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    if (fi->size()!=3) return false;

  char buf[sizeof(float)];
  std::string str(filename);

  std::ofstream output(filename);

  // FILE HEADER
  *((unsigned int*)buf) = 4;
  output.write(buf,sizeof(float));
  output.write("IFS",4);

  *((float*)buf) = 1.;
  output.write(buf,sizeof(float));
  
  *((unsigned int*)buf) = str.size();
  output.write(buf,sizeof(float));
  output.write(str.c_str(),str.size());
  
  // VERTEX HEADER  
  *((unsigned int*)buf) = 9;
  output.write(buf,sizeof(float));
  output.write("VERTICES",9);

  *((unsigned int*)buf) = vertices.size();
  output.write(buf,sizeof(float));

  // VERTEX LIST
  std::vector<R3>::iterator vi;
  for (vi=vertices.begin(); vi!=vertices.end(); vi++)
    {
      *((float*)buf) = float(vi->x);
      output.write(buf,sizeof(float));
      *((float*)buf) = float(vi->y);
      output.write(buf,sizeof(float));
      *((float*)buf) = float(vi->z);
      output.write(buf,sizeof(float));      
    }

  // TRI HEADER
  *((unsigned int*)buf) = 10;
  output.write(buf,sizeof(float));
  output.write("TRIANGLES",10);

  *((unsigned int*)buf) = facets.size();
  output.write(buf,sizeof(float));

  // TRI LIST
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    for (int i=0; i<fi->size(); i++)
      {
	*((unsigned int*)buf) = fi->vertex(i);
	output.write(buf,sizeof(float));
      }
  
  output.close();
  return true;
}
