#ifndef POLYGONMESH_H
#define POLYGONMESH_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "defs_polymesh.h" // For the definition of REAL

////////////////////////////////////////////////////////////////////////////////

// Add the ability to mark a vertex
template < class Refs, class P>
class PolygonMesh_vertex : public CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, P> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef typename CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, P>      Base;
    typedef typename CGAL::Tag_true                       Supports_vertex_halfedge;
    typedef typename CGAL::Tag_true                       Supports_vertex_point;
    typedef P                                    Point;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Halfedge              Halfedge;
    typedef typename Refs::Face                  Face;

    ////////////////////////////////////////////////////////////

protected:
    RealPoint2 _tex;     // Texture coordinates
    bool isEON;
    bool isPolar;
    bool _flag;
    int _index;
    RealPoint _normal;

public:

    // For Tianyun's code
    int t_index;
    bool computed;
    Vector pos;

	// Constructors
    PolygonMesh_vertex()            : isEON(false), isPolar(false)
    {
        // For Tianyun's code
        t_index=-1; computed=false;
    }
    PolygonMesh_vertex(const P &pp) :
        CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, P>(pp), isEON(false), isPolar(false)
    {
        // For Tianyun's code
        t_index=-1; computed=false;
    }
    // Copy constructor to prevent the copy of the flag (important!!)
    PolygonMesh_vertex(const PolygonMesh_vertex &v) :
        CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, P>(v), isEON(false), isPolar(false)
    { }

    // Returns the isEOP flag (settable)
          bool &is_eon()       { return isEON; }
    const bool &is_eon() const { return isEON; }

    // Returns the isPolar flag (settable)
          bool &is_polar()       { return isPolar; }
    const bool &is_polar() const { return isPolar; }

    // Returns flag (settable)
          bool &flag()       { return _flag; }
    const bool &flag() const { return _flag; }

    // Returns texture coordinates (settable)
          RealPoint2 &tex()       { return _tex; }
    const RealPoint2 &tex() const { return _tex; }

    // Returns normal (settable)
          RealPoint &normal()       { return _normal; }
    const RealPoint &normal() const { return _normal; }

    // Returns the index (settable)
          int &index()       { return _index; }
    const int &index() const { return _index; }

    // Returns true if the vertex lies on the border.
    bool
    is_border_vertex() const
    {
        Halfedge_const_handle evc, evc_begin;

        evc = evc_begin = this->halfedge();
        do {
            if (evc->is_border())
               return true;
            evc = evc->next_on_vertex();
        } while (evc != evc_begin);

        return false;
    }

    // Compute the normal associated with a vertex
    RealPoint
    compute_normal() const
    {
        Halfedge_const_handle evc, evc_begin;

        RealPoint norm = RealPoint(0,0,0);

        evc = evc_begin = this->halfedge();
        RealPoint p0 = KToRealPoint(this->point());
        do {
            if (!evc->is_border()) {
                RealPoint p1 = KToRealPoint(evc->opposite()->vertex()->point());
                RealPoint p2 = KToRealPoint(evc->next()->vertex()->point());
                norm += cross(p2 - p0, p1 - p0).normalize();
            }
            evc = evc->next_on_vertex();
        } while (evc != evc_begin);

        norm.normalize();

        return norm;
    }

    // Returns true if the vertex lies on the border, and in this case
    // sets h to the border halfedge.
    bool
    is_border_vertex(Halfedge_handle &border_halfedge)
    {
        Halfedge_handle evc, evc_begin;

        evc = evc_begin = this->halfedge();
        do {
            if (evc->is_border()) {
                border_halfedge = evc;
                return true;
            }
            evc = evc->next_on_vertex();
        } while (evc != evc_begin);

        return false;
    }

    // Returns true iff the vertex is convex
    bool
    is_convex() const {
        CGAL_assertion(this->halfedge()->vertex_degree() > 0);
        Vertex_const_handle vh = this->halfedge()->vertex();

        CGAL_assertion(vh->is_polar() && vh->vertex_degree() >= 3);
        int n = vh->vertex_degree();
        RealPoint eon = KToRealPoint(vh->point());
        Halfedge_const_handle eh_begin, eh = vh->halfedge();
        eh_begin = eh;

        RealPoint p0 = KToRealPoint(eh->opposite()->vertex()->point());
        RealPoint p1 = KToRealPoint(eh->opposite()->next()->vertex()->point());
        RealPoint p2 = KToRealPoint(eh->opposite()->prev()->opposite()->next()->vertex()->point());

        REAL vol = thresh(volume2(eon, p0, p1, p2));
        eh = eh->opposite()->prev();
        for (int i = 1; i < n; ++i) {
            p0 = KToRealPoint(eh->opposite()->vertex()->point());
            p1 = KToRealPoint(eh->opposite()->next()->vertex()->point());
            p2 = KToRealPoint(eh->opposite()->prev()->opposite()->next()->vertex()->point());

            if ( vol * thresh(volume2(eon, p0, p1, p2)) < 0 )
                return false;
            eh = eh->opposite()->prev();
        }
        CGAL_assertion(eh_begin == eh);

        return true;
    }


    // Returns true if the vertex is adjacent to a non-quad vertex.
    bool
    adjacent_to_non_quad() const
    {
        Halfedge_const_handle h = this->halfedge();
        do {
            if (h->facet_degree() != 4)
               return true;
            h = h->next()->opposite();
        } while (h != this->halfedge());

        return false;
    }

    // Returns true iff all the incident faces of the vertex are triangles.
    bool
    neighbors_all_tri() const
    {
        Halfedge_const_handle h = this->halfedge();
        do {
            if (h->facet_degree() != 3)
               return false;
            h = h->next()->opposite();
        } while (h != this->halfedge());

        return true;
    }

    // Returns true iff the vertex and its one-ring form a "flat" polar
    // structure.
    bool
    flat_polar() const
    {
        Halfedge_const_handle h = this->halfedge();
        int n = (int)h->vertex_degree();
        RealPoint cent(0,0,0);
        do {
            if (h->facet_degree() != 3) return false;
            cent += KToRealPoint(h->opposite()->vertex()->point());
            h = h->next()->opposite();
        } while (h != this->halfedge());
        cent /= n;

        REAL d = 0;
        do {
            REAL dist = (cent - KToRealPoint(h->opposite()->vertex()->point())).magsq();
            if (dist > d) d = dist;
            h = h->next()->opposite();
        } while (h != this->halfedge());

        RealPoint v = KToRealPoint(this->point());

        return ( ((cent - v).magsq()/d) <= 0.04 ); // 0.2^2
    }

    // Checks to see if the vertex is a non-border vertex and surrounded entirely by
    // triangles.
    bool
    is_polar_vertex() const {
        if (is_border_vertex()) return false;

        Halfedge_const_handle h = this->halfedge();
        Halfedge_const_handle h_begin = h;

        do {
            if (h->facet_degree() != 3)
                return false;
            h = h->next()->opposite();
        } while (h != h_begin);

        return true;
    }

    // Checks to see if the vertex is a non-border vertex and surrounded entirely by
    // triangles and a certain number of "regular" layers.
    bool
    is_polar_vertex_with_regular_facet_layers(int num_regular_layers = 1) const {
        if (is_border_vertex()) return false;

        Halfedge_const_handle h = this->halfedge();
        Halfedge_const_handle h_begin = h;

        do {
            if (h->facet_degree() != 3)
                return false;

            Halfedge_const_handle eh = h;
            for(int i = 0; i < num_regular_layers; ++i) {
                if (eh->opposite()->vertex_degree() != 4)
                    return false;
                if (eh->prev()->opposite()->facet_degree() != 4)
                    return false;
                eh = eh->prev()->opposite()->prev();
            }

            h = h->next()->opposite();
        } while (h != h_begin);

        return true;
    }

    // Checks to see if the vertex is a non-border vertex and surrounded entirely by
    // triangles and a certain number of "regular" links (checks one less
    // facet layer than above).
    bool
    is_polar_vertex_with_regular_links(int num_regular_links = 1) const {
        if (is_border_vertex()) return false;

        Halfedge_const_handle h = this->halfedge();
        Halfedge_const_handle h_begin = h;

        do {
            if (h->facet_degree() != 3)
                return false;

            Halfedge_const_handle eh = h;
            for(int i = 0; i < num_regular_links; ++i) {
                if (i != 0 && eh->facet_degree() != 4)
                    return false;
                if (eh->opposite()->vertex_degree() != 4)
                    return false;
                eh = eh->prev()->opposite()->prev();
            }

            h = h->next()->opposite();
        } while (h != h_begin);

        return true;
    }
};

////////////////////////////////////////////////////////////

template < class Refs, class Pln>
class PolygonMesh_face {
public:
    typedef Refs                                 HalfedgeDS;
    typedef CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Pln>      Base;
    typedef CGAL::Tag_true                       Supports_face_halfedge;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Halfedge              Halfedge;
    // Additional tags and types required by Polyhedron.
    typedef CGAL::Tag_true                       Supports_face_plane;
    typedef Pln                                  Plane;
    // No longer required.
    //typedef Tag_true                             Supports_face_normal;
    //typedef Trts                                 Traits;
    //typedef typename Traits::Normal              Normal;
    //typedef typename Traits::Plane               Plane;
private:
    Halfedge_handle hdg;
    Plane           pln;
    bool            irreg;
    int             _index;
    int             _type;
public:
    PolygonMesh_face() : irreg(false) {}
    PolygonMesh_face( const Plane& g) : irreg(false), pln(g) {}
    Halfedge_handle       halfedge()                        { return hdg; }
    Halfedge_const_handle halfedge() const                  { return hdg; }
    void                  set_halfedge( Halfedge_handle h)  { hdg = h; }
    Plane&                plane()                           { return pln; }
    const Plane&          plane() const                     { return pln; }
    bool       &is_irreg()       { return irreg; }
    const bool &is_irreg() const { return irreg; }

    // Returns the index (settable)
          int &index()       { return _index; }
    const int &index() const { return _index; }

    // Returns an integer type (used, for example, to figure out the rule to
    // be used) (settable)
          int &type()       { return _type; }
    const int &type() const { return _type; }

    int num_eon() const {
        int numEON = 0;
        Halfedge_handle efc, efc_begin;
        //efc = efc_begin = this->halfedge();
        efc = efc_begin = hdg;
        do {
            if (efc->vertex()->is_eon())
                ++numEON;
            efc = efc->next();
        } while (efc != efc_begin);

        return numEON;
    }

    // No longer required.
    //Normal                normal() const { return pln.orthogonal_vector();}

    // Returns true if the face is isolated; otherwise, it returns false
    // and a non-border edge.
    bool
    is_isolated_face(Halfedge_handle &non_border)
    {
        Halfedge_handle efc, efc_begin;
        efc = efc_begin = this->halfedge();
        do {
            if (!efc->is_border_edge()) {
               non_border = Halfedge_handle(efc);
               return false;
            }
            efc = efc->next();
        } while (efc != efc_begin);

        return true;
    }

    // Returns true if the face has a border edge.
    bool
    is_border_face() const
    {
        Halfedge_const_handle efc, efc_begin;
        efc = efc_begin = this->halfedge();
        do {
            if (efc->is_border_edge()) {
               return true;
            }
            efc = efc->next();
        } while (efc != efc_begin);

        return false;
    }

    // Returns true if the face is adjacent to a non-quad vertex.
    bool
    adjacent_to_non_quad() const
    {
        Halfedge_const_handle efc, efc_begin;
        efc = efc_begin = this->halfedge();
        do {
            if (efc->opposite()->facet_degree() != 4) {
               return true;
            }
            efc = efc->next();
        } while (efc != efc_begin);

        return false;
    }

    // Returns true if the face is adjacet to a vertex that is MARKED as
    // polar, and optionally returns the halfedge that points to the polar
    // vertex.
    bool
    is_polar_face(Halfedge_const_handle *hp = NULL) const
    {
	Halfedge_const_handle hf = this->halfedge();

        if ((int)hf->facet_degree() != 3 ||
            !(
                hf->vertex()->is_polar() ||
                hf->next()->vertex()->is_polar() ||
                hf->next()->next()->vertex()->is_polar()
	     ))
	    return false;

	if (hp) {
	    Halfedge_const_handle &h = *hp;
	    h = hf;
	    if (!h->vertex()->is_polar()) {
		h = h->next();
		if (!h->vertex()->is_polar())
		    h = h->next();
	    }
            assert(h->vertex()->is_polar());
	}

        return true;
    }

    // The non-const version
    bool
    is_polar_face(Halfedge_handle *hp = NULL)
    {
	Halfedge_handle hf = this->halfedge();

        if ((int)hf->facet_degree() != 3 ||
            !(
                hf->vertex()->is_polar() ||
                hf->next()->vertex()->is_polar() ||
                hf->next()->next()->vertex()->is_polar()
	     ))
	    return false;

	if (hp) {
	    Halfedge_handle &h = *hp;
	    h = hf;
	    if (!h->vertex()->is_polar()) {
		h = h->next();
		if (!h->vertex()->is_polar())
		    h = h->next();
	    }
	}

        return true;
    }

    // For Tianyun's code.
    Vector	normal;
};

////////////////////////////////////////////////////////////

template < class Refs >
class PolygonMesh_halfedge : public CGAL::HalfedgeDS_halfedge_base< Refs> {
//template < class Refs, class Pln >
//class HalfedgeDS_face_base< Refs, Tag_true, Pln> {
public:
    typedef Refs                                 HalfedgeDS;
    typedef typename Refs::Vertex                                  Vertex;
//  typedef typename Refs::Facet                                   Facet;
    typedef typename Refs::Vertex_handle                           Vertex_handle;
    typedef typename Refs::Halfedge_handle                         Halfedge_handle;
//  typedef typename Refs::Facet_handle                            Facet_handle;
//  typedef typename Refs::Halfedge_around_vertex_circulator       Halfedge_around_vertex_circulator;
//  typedef typename Refs::Halfedge_around_facet_circulator        Halfedge_around_facet_circulator;
    typedef typename Refs::Vertex_const_handle                     Vertex_const_handle;
    typedef typename Refs::Halfedge_const_handle                   Halfedge_const_handle;
//  typedef typename Refs::Facet_const_handle                      Facet_const_handle;
//  typedef typename Refs::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
//  typedef typename Refs::Halfedge_around_facet_const_circulator  Halfedge_around_facet_const_circulator;
    typedef CGAL::Tag_true Supports_halfedge_prev;
    typedef CGAL::Tag_true Supports_halfedge_vertex;
    typedef CGAL::Tag_true Supports_halfedge_face;
private:
    int             _local_index;
    REAL            _blend_ratio;
public:
    // Returns the index (settable)
          int &local_index()       { return _local_index; }
    const int &local_index() const { return _local_index; }

    // Returns the blend ratio (settable)
          REAL &blend_ratio()       { return _blend_ratio; }
    const REAL &blend_ratio() const { return _blend_ratio; }

    // Returns the blend ratio (settable)
    void set_edge_blend_ratio(REAL blend) { this->blend_ratio() = blend; this->opposite()->blend_ratio() = blend; }

    // Compute the number of times you have to rotate counter-clockwise
    // around a vertex to reach this edge.
    int
    rotation_index_of_edge() const
    {
        Halfedge_const_handle h = this->vertex()->halfedge();
        int i = 0;
        do {
            if (&*h == this) return i;
            h = h->opposite()->prev();
            ++i;
        } while (h != this->vertex()->halfedge());

        assert(false); // Shouldn't happen.
        return -1;
    }

    // For Tianyun's code
    Vector CCedge;  //like in Loop's method
    Vector edge;    // projected
    Vector e310;    // after degree raising
};

////////////////////////////////////////////////////////////

// Polygon mesh items (with new vertex type)
class PolygonMeshItems_3 {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3 Point;
        typedef PolygonMesh_vertex< Refs, Point >                      Vertex;
    };

    template < class Refs, class Traits>
    struct Halfedge_wrapper
    {
        typedef PolygonMesh_halfedge< Refs >                           Halfedge;
    };

    template < class Refs, class Traits>
    struct Face_wrapper
    {
        typedef typename Traits::Plane_3 Plane;
        typedef PolygonMesh_face< Refs, Plane>                         Face;
    };
};

////////////////////////////////////////////////////////////////////////////////

// PolygonMesh type
class PolygonMesh : public CGAL::Polyhedron_3<Kernel, PolygonMeshItems_3> {
public:
    // Add array of quads to mesh. There MUST be at least 1 quad (at least a
    // 2x2 array of points).
    // ORIENTATION:
    //     Creates quads with CCW orientation assuming b(0,0) being the
    //     top-left corner, and b(nr-1,0) begin the bottom-left corner.
    void add_to_mesh(const RealPoint *b, int nr, int nc, int stepr, int stepc)
    {
        assert (nr > 1 && nc > 1);

        PolygonMesh &m = *this;

        int i, j;
        PolygonMesh::Halfedge_handle e, e_next_row;

        // Create the first quad
        e = m.make_triangle(); m.split_edge(e);
        e->vertex()->point() = RealToKPoint(b[0    ]); e = e->next();
        e->vertex()->point() = RealToKPoint(b[stepr]); e = e->next();
        e_next_row = e;
        e->vertex()->point() = RealToKPoint(b[stepr+stepc]); e = e->next();
        e->vertex()->point() = RealToKPoint(b[stepc]);

        // Make the first layer of quads
        for (j = 2; j < nc; ++j) {
            e = m.add_vertex_and_facet_to_border(e->next()->opposite(), e->opposite());
            m.split_edge(e);
            e->vertex()->point()             = RealToKPoint(b[        j*stepc]);
            e->opposite()->vertex()->point() = RealToKPoint(b[stepr + j*stepc]);
        }

        // Make all the subsequent layers of quads
        for (i = 2; i < nr; ++i) {
            // Add the first quad of this row
            e_next_row = m.add_vertex_and_facet_to_border(e_next_row->opposite()->prev(), e_next_row->opposite());
            m.split_edge(e_next_row);
            e_next_row->vertex()->point()             = RealToKPoint(b[i*stepr + stepc]);
            e_next_row->opposite()->vertex()->point() = RealToKPoint(b[i*stepr        ]);
            e = e_next_row->next();

            // Create the rest of the layer of quads
            for (j = 2; j < nc; ++j) {
                e = m.add_vertex_and_facet_to_border(e->opposite()->prev()->prev(), e->opposite());
                e->vertex()->point() = RealToKPoint(b[i*stepr + j*stepc]);
                e = e->next();
            }
        }
    }

    // Creates a quad incident to a border edge.
    // ---+ e        ---+----+
    //    | |           |    |
    //    | V           |    |
    // ---+          ---+----+
    //                   <--- (returned)
    Halfedge_handle add_quad_1(Halfedge_handle e) {
        return this->split_edge( this->add_vertex_and_facet_to_border(e->prev(), e) )->opposite();
    }

    // Creates a quad incident to two edges (e and e->next()).
    //   +----+ e      +----+
    //   |    | |      |    |\.
    //   |    | V      |    | \.
    //   +----+        +----+  +
    //    \    \        \    \ | | (returned)
    //     \____\        \____\| V
    Halfedge_handle add_quad_2(Halfedge_handle e) {
        return this->split_edge( this->add_facet_to_border(e->prev(), e->next()) )->opposite();
    }

    // Creates a quad incident to three border edges
    // (e, e->next(), and e->next()->next()).
    //
    //                               ---> (returned)
    //   +----+    +----+      +----+----+----+
    //   |    ||e  |    |      |    |    |    |
    //   |    |V   |    |      |    |    |    |
    //   +----+----+----+      +----+----+----+
    //   |    |    |    |      |    |    |    |
    //   |    |    |    |      |    |    |    |
    //   +----+----+----+      +----+----+----+
    Halfedge_handle add_quad_3(Halfedge_handle e) {
        return this->add_facet_to_border( e->prev(), e->next()->next() )->opposite();
    }

    // Compute normals at every vertex
    void compute_and_store_normals() {
        for (PolygonMesh::Vertex_iterator vi = this->vertices_begin(); vi != this->vertices_end(); ++vi) {
            vi->normal() = vi->compute_normal();
        }
    }

    // Mark all the vertex indices
    void
    mark_all_vertex_indices()
    {
        Vertex_iterator iter;
        int i;
        for (iter = vertices_begin(), i = 0; iter != vertices_end(); ++iter, ++i) {
            iter->index() = i;
        }
        CGAL_assertion(i == (int)size_of_vertices());
    }

    // Mark all the face indices
    void
    mark_all_facet_indices()
    {
        Facet_iterator iter;
        int i;
        for (iter = facets_begin(), i = 0; iter != facets_end(); ++iter, ++i) {
            iter->index() = i;
        }
        CGAL_assertion(i == (int)size_of_facets());
    }

    // Mark all the local edge indices
    void
    mark_local_edge_indices()
    {
        for (Facet_iterator iter = facets_begin(); iter != facets_end(); ++iter) {
            Halfedge_around_facet_circulator it = iter->facet_begin();
            int i = 0;
            do {
                it->local_index() = i;
                ++it; ++i;
            } while (it != iter->facet_begin());
            CGAL_assertion(i == (int)iter->facet_degree());
        }
    }

};

#endif // POLYGONMESH_H

