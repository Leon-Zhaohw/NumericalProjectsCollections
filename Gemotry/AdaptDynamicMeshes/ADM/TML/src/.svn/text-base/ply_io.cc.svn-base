#include "tml.h"

using namespace TML;

bool TriMeshIO::Read_ply(const char* filename)
{
  p_ply ply = ply_open(filename,NULL);
  if (!ply) return false;
  if (!ply_read_header(ply)) return false;

  ply_set_read_cb( ply, "vertex", "x", vertex_cb, (void*)&vertices, 0 );
  ply_set_read_cb( ply, "vertex", "y", vertex_cb, (void*)&vertices, 1 );
  ply_set_read_cb( ply, "vertex", "z", vertex_cb, (void*)&vertices, 2 );
  
  ply_set_read_cb( ply, "face", "vertex_indices", face_cb, (void*)&facets, 0);
  
  if (!ply_read(ply)) return false;

  if (!ply_close(ply)) return false;

  return check_sanity();
}

int TriMeshIO::vertex_cb(p_ply_argument argument)
{
  long eol;
  void *pdata;
  ply_get_argument_user_data( argument, &pdata, &eol );
  
  if (eol==0)
    {
      R3 v;
      v.x = ply_get_argument_value(argument);
      v.y = v.z = 0.;
      ((std::vector<R3>*)pdata)->push_back( v );
      return 1;
    }

  if (eol==1)
    {
      int n = ((std::vector<R3>*)pdata)->size() - 1;
      (*((std::vector<R3>*)pdata))[n].y = ply_get_argument_value(argument);
      return 1;
    } 

  int n = ((std::vector<R3>*)pdata)->size() - 1;
  (*((std::vector<R3>*)pdata))[n].z = ply_get_argument_value(argument);
  return 1;
}

int TriMeshIO::face_cb(p_ply_argument argument)
{
  long eol;
  void *pdata;
  ply_get_argument_user_data( argument, &pdata, &eol );

  long length, value_index;
  ply_get_argument_property(argument, NULL, &length, &value_index);

  if (value_index==-1)
    {
      int n = int(length);
      Polygon f(n);
      ((std::vector<Polygon>*)pdata)->push_back( f );
      return 1;
    }
  
  int n = ((std::vector<Polygon>*)pdata)->size() - 1;
  (*((std::vector<Polygon>*)pdata))[n].vertex(value_index) = int( ply_get_argument_value(argument) );
  
  return 1;
}

/////////////////////////////////////////////////////////////////////////////////

bool TriMeshIO::Write_ply(const char* filename)
{
  p_ply ply = ply_create(filename, PLY_ASCII, NULL); //PLY_LITTLE_ENDIAN, NULL)
  if (!ply) return false;

  // HEADER
  if (!ply_add_element(ply,"vertex",vertices.size())) return false;
  if (!ply_add_property(ply,"x",PLY_FLOAT,PLY_FLOAT,PLY_FLOAT)) return false;
  if (!ply_add_property(ply,"y",PLY_FLOAT,PLY_FLOAT,PLY_FLOAT)) return false;
  if (!ply_add_property(ply,"z",PLY_FLOAT,PLY_FLOAT,PLY_FLOAT)) return false;

  if (!ply_add_element(ply,"face",facets.size())) return false;
  if (!ply_add_property(ply,"vertex_indices",PLY_LIST,PLY_UCHAR,PLY_INT32)) return false;
  
  if (!ply_write_header(ply)) return false;

  // BODY
  std::vector<R3>::iterator vi;
  for (vi=vertices.begin(); vi!=vertices.end(); vi++)
    {
      if (!ply_write( ply, vi->x )) return false;
      if (!ply_write( ply, vi->y )) return false;
      if (!ply_write( ply, vi->z )) return false;
    }

  std::vector<Polygon>::iterator fi;
  for (fi=facets.begin(); fi!=facets.end(); fi++)
    {
      if (!ply_write( ply, fi->size() )) return false;
      for (int i=0; i<fi->size(); i++)
	if (!ply_write( ply, double(fi->vertex(i)) )) return false;
    }

  if (!ply_close(ply)) return false;
  
  return true;
}
