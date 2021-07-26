#include "adm.h"

using namespace ADM;

A48Vertex* AdaptiveMesh::split(Halfedge* e)
{
  Halfedge *el, *er;
  if (e->facet() == NULL) e = e->mate();

  A48Facet* f = der_cast(e->facet());
  if (f == NULL) throw TML::Error("subdiv edge");
  int lf0 = f->level();

  A48Facet* fm = der_cast(e->mate()->facet());
  int lf1 = (fm)? fm->level() : 0;

  Halfedge* ef1 = e->next();
  Halfedge* ef2 = e->prev();
  Halfedge* efm1 = (fm)? e->mate()->next() : NULL;
  Halfedge* efm2 = (fm)? e->mate()->prev() : NULL;

  A48Vertex* v = bisect(e, &el, &er);
  bisect(f, ef1, ef2, e, er);
  bisect(fm, efm1, efm2, er->mate(), e->mate());
  v->set_star(e);
  v->set_level( std::max(lf0, lf1) + 1 );

#ifdef CHECKER
  // e -> (el,er)
  er->set_mark( e->mark() );
  er->Events = e->Events;
  er->mate()->set_mark( e->mate()->mark() );
  er->mate()->Events = e->mate()->Events;

  split_color( ef1 , v->level(), e->mark() );
  split_color( ef2 , v->level(), e->mark() );
  split_color( efm1, v->level(), e->mate()->mark() );
  split_color( efm2, v->level(), e->mate()->mark() );
#endif
  
  return v;
}

A48Vertex* AdaptiveMesh::bisect(Halfedge* e, Halfedge** el, Halfedge** er)
{
  A48Vertex* v0 = der_cast(e->org());
  A48Vertex* v1 = der_cast(e->dst());

  A48Vertex* m  = add_vertex();
  m->set_level( der_cast(e->edge())->level() + 1 );
  sample( der_cast(e->edge()), m );

  *el = halfedge_reuse( e, base_cast(v0), base_cast(m) );
  *er = add_edge(m, v1);
  if (v1->star_first() == e)
    v1->set_star((*er));
  return m;
}

Halfedge* AdaptiveMesh::bisect(A48Facet* f, Halfedge* e1, Halfedge* e2, Halfedge* el, Halfedge* er)
{
  if (f == NULL) return 0;
  Halfedge* em = add_edge( der_cast(e2->org()), der_cast(er->org()) );
  f->reuse(e1, em, er);
  add_facet( new A48Facet(e2, el, em->mate()) );
  return em;
}

Halfedge* AdaptiveMesh::weld(A48Vertex* w)
{
  int k, n;
  Halfedge *ee, *e[6]; 
  A48Facet* f[6]; 
  A48Vertex* v[6]; 

  if (w->level() == 0) return NULL;
  
  for (n=0, ee=w->star_first(); ee!=NULL; n++, ee=w->star_next(ee)) 
    {
      if (n > 4) throw TML::Error("weld");
      e[n] = ee;
      v[n] = der_cast(ee->org());
      f[n] = der_cast(ee->facet());
    }
  
  if (n != 4 && n != 3) 
    { std::cerr << "can't weld " << n << "\n"; return NULL; }
  
  Halfedge* n2 = e[2]->mate()->next();
  Halfedge* p0 = e[0]->prev();
  Halfedge* n0 = (f[2])? e[0]->mate()->next() : NULL;
  Halfedge* p2 = (f[2])? e[2]->prev() : NULL;
  
#ifdef CHECKER
  weld_color( n2, w->level(), e[0]->mark() );
  weld_color( p0, w->level(), e[0]->mark() );
  weld_color( n0, w->level(), e[2]->mark() );
  weld_color( p2, w->level(), e[2]->mark() );
#endif

  Halfedge* en = e[0]; 
  halfedge_reuse( e[0], base_cast(v[0]), base_cast(v[2]) );
  if (v[2]->star_first() == e[2]->mate())
    v[2]->set_star(en);
  if (v[1]->star_first() == e[1]->mate())
    v[1]->set_star(n2);
  if (n == 4 && (v[3]->star_first() == e[3]->mate()))
    v[3]->set_star(n0);
  
  f[0]->reuse(en, n2, p0);
  if (f[2] != NULL)
    f[1]->reuse(en->mate(), n0, p2);
  else
    del_facet(f[1]);
  
  for (k = 2; k < n; k++) 
    if (f[k] != NULL)
      del_facet(f[k]);
  for (k = 1; k < n; k++)
    del_edge( der_cast(e[k]->edge()) );
  del_vertex(w);
  
  return en;
}

A48Vertex* AdaptiveMesh::split(A48Facet* f)
{
  Halfedge* e0 = f->hedge(0);
  Halfedge* e1 = f->hedge(1);
  Halfedge* e2 = f->hedge(2);

  A48Vertex* v0 = der_cast(f->vertex(0));
  A48Vertex* v1 = der_cast(f->vertex(1));
  A48Vertex* v2 = der_cast(f->vertex(2));

  A48Vertex* vc = add_vertex();
  vc->set_level( f->level() + 1 );
  sample(f, vc);

  Halfedge* e0c = add_edge(v0, vc);
  Halfedge* e1c = add_edge(v1, vc);
  Halfedge* e2c = add_edge(v2, vc);
  f->reuse(e0, e2c, e1c->mate());
  add_facet( new A48Facet(e1, e0c, e2c->mate()) );
  add_facet( new A48Facet(e2, e1c, e0c->mate()) );
  vc->set_star(e0c);
  return vc;
}

Halfedge* AdaptiveMesh::flip(Halfedge* h)
{
  Halfedge* m = h->mate();
  A48Facet* fl = der_cast(h->facet());
  A48Facet* fr = der_cast(m->facet());
  if (fl == NULL || fr == NULL) return h;

  Halfedge* hp = h->prev(); 
  Halfedge* hn = h->next();
  Halfedge* mp = m->prev();
  Halfedge* mn = m->next();

  A48Vertex* v0 = der_cast(h->org());
  A48Vertex* v1 = der_cast(h->dst());
  A48Vertex* vl = der_cast(hp->org());
  A48Vertex* vr = der_cast(mp->org());

  if (v0->star_first() == m)
    v0->set_star(hp);
  if (v1->star_first() == h)
    v1->set_star(mp);

  Halfedge* o = halfedge_reuse( h, base_cast(vl), base_cast(vr) );
  fr->reuse(o, mp, hn);
  fl->reuse(o->mate(), hp, mn);
  return o;
}
