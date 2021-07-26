#ifndef _TRIMESH_
#define _TRIMESH_ 1

namespace TML
{
  template <class V=Vertex, class E=Edge, class F=Facet> 
    class TriMesh
    {
      public:
      typedef std::set<V*>                       VertexContainer;
      typedef typename VertexContainer::iterator VertexIter;
      
      typedef std::set<E*>                     EdgeContainer;
      typedef typename EdgeContainer::iterator EdgeIter;
      
      typedef std::set<F*>                      FacetContainer;
      typedef typename FacetContainer::iterator FacetIter;
      
      protected:
      FacetContainer  fc_; 
      EdgeContainer   ec_; 
      VertexContainer vc_; 
      double          bbmax_[3], bbmin_[3];
      
      typedef std::map<const Ipair, Halfedge*> HedgeMap;
            
      public:
      TriMesh( const char* filename )
      : fc_(), ec_(), vc_()
      {
	bbmax_[0] = bbmax_[1] = bbmax_[2] = 0.;
	bbmin_[0] = bbmin_[1] = bbmin_[2] = 0.;
	
	TriMeshIO reader;
	if (!reader.Read(filename))    throw Error("TriMesh::Constructor: file corrupted");
	if (!reader.is_triangulated()) throw Error("TriMesh::not triangulated");
	set_mesh( reader.get_vertices(), reader.get_facets() );
	
	set_vertex_attr();
	set_edge_attr();
	set_facet_attr();
      }      

      TriMesh(std::vector<R3>& vrts, std::vector<Polygon>& fcts)
      : fc_(), ec_(), vc_()
      {
	bbmax_[0] = bbmax_[1] = bbmax_[2] = 0.;
	bbmin_[0] = bbmin_[1] = bbmin_[2] = 0.;
	
	set_mesh( vrts, fcts );
	
	set_vertex_attr();
	set_edge_attr();
	set_facet_attr();
      }      

      virtual ~TriMesh()
      {
	for (VertexIter vi=vc_.begin(); vi!=vc_.end(); vi++) delete(*vi);
	vc_.clear();
	for (EdgeIter ei=ec_.begin(); ei!=ec_.end(); ei++)   delete(*ei);
	ec_.clear();
	for (FacetIter fi=fc_.begin(); fi!=fc_.end(); fi++)  delete(*fi);
	fc_.clear();
      }     
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

      void set_mesh(std::vector<R3>& vrts, std::vector<Polygon>& fcts)
      {
	V** verts = new V*[vrts.size()];
	
	int i;
	std::vector<R3>::iterator vi;
	for (i=0, vi=vrts.begin(); vi!=vrts.end(); i++, vi++)
	{
	  verts[i] = new V(*vi);
	  add_vertex(verts[i]);
	  base_cast(verts[i])->set_id(i);
	}
	
	HedgeMap hedges;
	std::vector<Polygon>::iterator fi;
	for (fi=fcts.begin(); fi!=fcts.end(); fi++)
	{ put_face( fi->vertex(2), fi->vertex(0), fi->vertex(1), verts, &hedges); }
	
	link_mesh();
	delete(verts);
      }
      
      void get_mesh(std::vector<R3>& vrts, std::vector<Polygon>& fcts)
      {
	vrts.clear();
	fcts.clear();
	
	vrts.assign( vc_.size(), R3() );
	for (VertexIter vi=vc_.begin(); vi!=vc_.end(); vi++)
	{ vrts[ (*vi)->id() ] = (*vi)->pos(); }
	
	for (FacetIter fi=fc_.begin(); fi!=fc_.end(); fi++)
	{
	  Polygon p(3);
	  p.vertex(0,(*fi)->vertex(1)->id());
	  p.vertex(1,(*fi)->vertex(2)->id());
	  p.vertex(2,(*fi)->vertex(0)->id());
	  fcts.push_back(p);
	}
      }
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      
      Vertex* base_cast(V* v) { return dynamic_cast<Vertex*>(v); }
      Facet*  base_cast(F* f) { return dynamic_cast<Facet*>(f);  }
      Edge*   base_cast(E* e) { return dynamic_cast<Edge*>(e);   }

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

      bool save_mesh( const char* filename )
      {
	TriMeshIO writer;
	get_mesh( writer.get_vertices(), writer.get_facets() );
	return writer.Write(filename);
      }

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      
      // Basic Methods
      VertexIter verts_begin(){ return vc_.begin(); };
      VertexIter verts_end()  { return vc_.end(); };
      
      EdgeIter edges_begin() { return ec_.begin(); };
      EdgeIter edges_end()   { return ec_.end(); };
      
      FacetIter facets_begin() { return fc_.begin(); };
      FacetIter facets_end()   { return fc_.end(); };
      
      inline int num_verts()    { return vc_.size(); };
      inline int num_facets()   { return fc_.size(); };
      inline int num_edges()    { return ec_.size(); };
      inline int num_halfedge() { return 2*num_edges(); }
      inline int num_bdry_edges() 
      {
	int n = 0;
	for (EdgeIter ei=edges_begin(); ei!=edges_end(); ei++)
	{ if ( (*ei)->is_bdry() ) n++; }
	return n;
      }
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

      V* add_vertex()
      {
	V* v = new V();
	add_vertex(v);
	return v;
      }

      bool add_vertex(V* v) 
      { std::pair<VertexIter, bool> r = vc_.insert(v); return r.second; };      

      void del_vertex(V* v) 
      { vc_.erase(v); delete(v); };
      
      Halfedge* add_edge(V* v0, V* v1)
      {
	E* e = new E( base_cast(v0), base_cast(v1) );
	add_edge(e);
	return e->hedge(0);
      }
      
      bool add_edge(E* e)
      { std::pair<EdgeIter, bool> r = ec_.insert(e);  return r.second;};

      void del_edge(E* e) 
      { ec_.erase(e); delete(e); }
      
      F* add_facet(Halfedge* e0, Halfedge* e1, Halfedge* e2)
      {
	F* f = new F(e0, e1, e2);
	add_facet(f);
	return f;
      }

      bool add_facet(F* f)
      { std::pair<FacetIter, bool> r = fc_.insert(f); return r.second;}

      void del_facet(F* f) 
      { fc_.erase(f); delete(f); }
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

      void link_mesh()
      {
	for (FacetIter fi=fc_.begin(); fi!=fc_.end(); fi++) 
	{ (*fi)->link_star_verts(); }

	for (EdgeIter ei=ec_.begin(); ei!=ec_.end(); ei++) 
	{
	  if ( !((*ei)->is_bdry()) ) continue;
	  
	  Halfedge* h = (*ei)->hedge(0);
	  if ( h->is_bdry() && h->mate()->is_bdry() ) 
	    throw Error("TriMesh::link_mesh");
	  else if ( h->is_bdry() )
	    h->org()->set_star(h->mate());
	  else 
	    h->dst()->set_star(h);
	}
      }
      
      Halfedge* get_hedge(int i0, int i1, V* verts[], HedgeMap* hedges)
      {
	bool mate = false;
	if (i0 > i1) 
	{
	  std::swap<int>(i0, i1);
	  mate = true;
	}
	Halfedge* he; 
	HedgeMap::iterator ei  = hedges->find( Ipair(i0,i1) );
	if (ei == hedges->end())
	{ (*hedges)[Ipair(i0,i1)] = he = add_edge(verts[i0], verts[i1]); }
	else
	{ he = (*ei).second; }
	return (mate)? he->mate() : he ;
      }
      
      void put_face(int i0, int i1, int i2, V* verts[], HedgeMap* hedges)
      {
	Halfedge* e0 = get_hedge(i1, i2, verts, hedges);
	Halfedge* e1 = get_hedge(i2, i0, verts, hedges);
	Halfedge* e2 = get_hedge(i0, i1, verts, hedges);
	add_facet(e0, e1, e2);
      }
     
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

      // Geometric Methods
      void set_vertex_attr()
      { 
	/*
	  int i = 0;
	  VertexIter vi;
	  for (i=0, vi=vc_.begin(); vi!=vc_.end(); i++, vi++) (*vi)->set_attr(i); 
	*/
      }
      
      void set_edge_attr()
      { for (EdgeIter ei=ec_.begin(); ei!=ec_.end(); ei++) (*ei)->set_attr(); }
      
      void set_facet_attr()
      { 
	int i;
	FacetIter fi;
	for (i=0, fi=fc_.begin(); fi!=fc_.end(); i++, fi++) (*fi)->set_attr(i); 
      }

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      
      void bounding_box(double& xmax,double& ymax,double& zmax,
			double& xmin,double& ymin,double& zmin)
      {
	bounding_box();
	xmax = bbmax_[0]; ymax = bbmax_[1]; zmax = bbmax_[2];
	xmin = bbmin_[0]; ymin = bbmin_[1]; zmin = bbmin_[2];
      }

      void bounding_box()
      {
	R3 g;
	
	VertexIter vi = vc_.begin();
	g = (*vi)->pos();
	bbmin_[0] = bbmax_[0] = g.x;
	bbmin_[1] = bbmax_[1] = g.y;
	bbmin_[2] = bbmax_[2] = g.z;
	vi++;
	
	for ( ; vi!=vc_.end(); vi++)
	{
	  g = (*vi)->pos();
	  bbmin_[0] = std::min( bbmin_[0],g.x );
	  bbmin_[1] = std::min( bbmin_[1],g.y );
	  bbmin_[2] = std::min( bbmin_[2],g.z );
	  bbmax_[0] = std::max( bbmax_[0],g.x );
	  bbmax_[1] = std::max( bbmax_[1],g.y );
	  bbmax_[2] = std::max( bbmax_[2],g.z );
	}
      }
      
    };
}

#endif
