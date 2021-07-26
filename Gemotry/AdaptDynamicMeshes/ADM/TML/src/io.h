#ifndef _TRIMESH_IO_
#define _TRIMESH_IO_ 1

namespace TML
{
  class TriMeshIO
    {
    protected:
      std::vector<R3>      vertices;
      std::vector<Polygon> facets;

      enum { NONE, OFF, SMF, PLY, IFS, OBJ };
      
    public:
      TriMeshIO() : vertices(), facets() { }
      ~TriMeshIO() { vertices.clear(); facets.clear(); }
      
      std::vector<R3>&      get_vertices() { return vertices; }
      std::vector<Polygon>& get_facets()   { return facets; }

      inline void add(R3& v)      { vertices.push_back(v); }
      inline void add(Polygon& f) { facets.push_back(f); }

      bool Read(const char* filename)
	{
	  switch (File_Format(filename))
	    {
	    case OFF: return Read_off(filename);
	    case SMF: return Read_smf(filename);
	    case IFS: return Read_ifs(filename);
	    case PLY: return Read_ply(filename);
	    case OBJ: return Read_obj(filename);
	    };
	  throw Error("TriMeshIO::Read: invalid format");
	  return false;
	}
      
      bool Write(const char* filename)
	{
	  switch (File_Format(filename))
	    {
	    case OFF: return Write_off(filename);
	    case SMF: return Write_smf(filename);
	    case IFS: return Write_ifs(filename);
	    case PLY: return Write_ply(filename);
	    case OBJ: return Write_obj(filename);
	    };
	  throw Error("TriMeshIO::Write: invalid format");
	  return false;
	}
      
      bool is_triangulated();
      bool invert();
      bool quad2tri();
      bool normalize_bb();
      bool check_sanity();

    protected:
      int File_Format(const char* filename);
      
      bool Read_smf(const char* filename);
      bool Write_smf(const char* filename);
      
      bool Read_off(const char* filename);
      bool Write_off(const char* filename);
      
      bool Read_ifs(const char* filename);
      bool Write_ifs(const char* filename);
      
      bool Read_ply(const char* filename);
      bool Write_ply(const char* filename);
      
      bool Read_obj(const char* filename);
      bool Write_obj(const char* filename);
      
      static int vertex_cb(p_ply_argument argument);
      static int face_cb(p_ply_argument argument);
    };
}

#endif
