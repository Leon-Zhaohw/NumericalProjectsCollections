#include <queue>

#include "adm.h"

using namespace ADM;

#ifdef CHECKER

/*********************************
 * Code to add in TML/src/edge.h *
 *********************************
  
  #ifdef CHECKER
  class Conflict
  {
  protected:
  int lv, inv;
  public:
  Conflict(int l, int i) : lv(l), inv(i) { }
  int level()  const { return lv;  }
  int invert() const { return inv; }
  };
  #endif
  
  // inside Halfedge class as public member
  #ifdef CHECKER
  std::stack<Conflict> Events;
  #endif

**********************************/

bool ADM::edge_is_subdiv2(Halfedge* h)
{ return ( halfedge_is_subdiv(h) && halfedge_is_subdiv(h->mate()) ); }
 
void ADM::split_color(Halfedge* e, int lv, int color)
{
  if (e==NULL) return;

  if (e->mark()==-1) 
    {
      e->set_mark(color);
      return;
    }
  
  if (lv%2==1)
    {
      e->Events.push( TML::Conflict(lv,(e->mark()==color)? 0 : 1) ); 
      e->set_mark(color);
      return;
    }

  if ( edge_is_subdiv2(e) )
    {
      e->set_mark( (color==0)? 1 : 0 );
      e->mate()->set_mark( (color==0)? 1 : 0 );
    }
}

void ADM::weld_color(Halfedge* e, int lv, int color)
{
  if (e==NULL) return;
  
  if (lv%2==1)
    {
      if (e->Events.empty())
	{
	  e->set_mark(-1);
	  return;
	}
      
      if ( (e->Events.top()).level()!=lv ) 
	throw TML::Error("Conflict with different levels");
      
      if ( (e->Events.top()).invert()==1 )
	e->set_mark( (e->mark()==0)? 1 : 0 );
      e->Events.pop();
      return;
    }

  if ( edge_is_subdiv2(e) )
    {
      if (e->mark()!=color) 
	{
	  e->set_mark(color);
	  e->mate()->set_mark(color);
	} else e->set_mark(-1);
      return;
    }
  
  if (e->mate()->mark()==-1) e->set_mark(-1);
}

void AdaptiveMesh::greedy_color()
{
  Halfedge* h;
  EdgeIter ei;

  // clear data
  for (ei=edges_begin(); ei!=edges_end(); ei++)
    {
      (*ei)->hedge(0)->set_mark(-1);
      while ( !(*ei)->hedge(0)->Events.empty() )
	(*ei)->hedge(0)->Events.pop();
      (*ei)->hedge(1)->set_mark(-1);
      while ( !(*ei)->hedge(1)->Events.empty() )
	(*ei)->hedge(1)->Events.pop();
    }
  
  // finding first subdiv_edge
  for (ei=edges_begin(); ei!=edges_end(); ei++)
    if ( edge_is_subdiv2((*ei)->hedge(0)) )
      {
	h = (*ei)->hedge(0);
	break;
      }
  h->set_mark(0);
  h->mate()->set_mark(0);
  
  int i, j, color;
  Halfedge* h2;
  std::queue<Halfedge*> Q;
  
  // loop
  Q.push(h);
  while (!Q.empty())
    {
      h = Q.front();
      Q.pop();
      color = (h->mark()==0)? 1 : 0;
      
      for (i=0; i<2; i++)
	{
	  if (!h->is_bdry())
	    {
	      h2 = h->next()->mate();
	      for (j=0; j<2; j++)
		{
		  if (!h2->is_bdry())
		    {
		      h2 = der_cast(h2->facet())->subd_edge();
		      if (h2->mark()==-1) 
			{
			  h2->set_mark(color);
			  h2->mate()->set_mark(color);
			  Q.push(h2);
			}
		    }
		  h2 = h->prev()->mate();
		} // for j
	    }
	  h = h->mate();
	} // for i
    } // while
}

#endif
