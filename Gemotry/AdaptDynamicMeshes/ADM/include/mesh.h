#ifndef _ADAPTIVE_MESH_
#define _ADAPTIVE_MESH_ 1

namespace ADM
{
  class AdaptiveMesh : public TML::TriMesh< A48Vertex, A48Edge, A48Facet >
    {
    protected:
      Oracle*        surf_;
      AdaptCriteria* criteria_;

      MxHeap sf_; // Priority queue for simplification front
      MxHeap rf_; // Priority queue for refinement front

      double threshold_, stiffness_, hysteresis_;
      int    smooth_steps_;
      
    public:
      AdaptiveMesh(const char* base_mesh, Oracle* s);
      
      AdaptiveMesh(const char* base_mesh, AdaptCriteria* cr);

      virtual ~AdaptiveMesh() { }

      /*------------------------------------------------------------*/
      
      double threshold()    const { return threshold_; }
      double stiffness()    const { return stiffness_; }
      double hysteresis()   const { return hysteresis_; }
      int    smooth_steps() const { return smooth_steps_; }

      void set_threshold(double t)  { threshold_ = t; }      
      void set_stiffness(double t)  { stiffness_ = t; }   
      void set_hysteresis(double t) { hysteresis_ = t; }   
      void set_smooth_steps(int t)  { smooth_steps_ = t; }
      
      /*------------------------------------------------------------*/
      
      /***************
       * A48 methods *
       ***************/
      
      void adapt_refine(double t);

      void adapt_refine_step(double t);

      void adapt_simplify(double t);

      void adapt_simplify_step(double t);

      A48Vertex* refine(Halfedge* e);
      
      Halfedge* simplify(A48Vertex* w);
      
      A48Vertex* split(A48Facet* f);

      A48Vertex* split(Halfedge* e);

      Halfedge* weld(A48Vertex* w);
      
      Halfedge*	flip(Halfedge* h);
      
      A48Vertex* bisect(Halfedge* e, Halfedge** el, Halfedge** er);

      Halfedge* bisect(A48Facet* f, Halfedge* e1, Halfedge* e2, Halfedge* el, Halfedge* er);

      void new_front();

      void update_ref_front(A48Edge* e);

      void update_ref_front(A48Vertex* v);

      void update_simpl_front(A48Vertex* v);

      void update_simpl_front(A48Edge* e);

      /************
       * tri-quad *
       ************/

      bool is_triquad();

      void make_triquad();

      /*------------------------------------------------------------*/

      void align_base_mesh();

      void adapt();

      void deform();

      void strutural_phase();

      void simplify_phase();

      void simplify_one_step();

      void refine_phase();

      void refine_one_step();

      void smoothing_phase();

      void stat_curvature(double& max_val, double& min_val);

      void avoid_folding(A48Vertex* v, const R3& u, double& beta);

      void update_attr();

      void sample(A48Edge* e, A48Vertex* v);

      void sample(A48Facet* f, A48Vertex* v);

      /*++++++++++++++++*
       * color routines *
       *++++++++++++++++*/
      void greedy_color();
    };
  
  /*++++++++++++++++*
   * color routines *
   *++++++++++++++++*/
  bool edge_is_subdiv2(Halfedge* h);
  
  void split_color(Halfedge* e, int lv, int color);
  
  void weld_color(Halfedge* e, int lv, int color);
}

#endif
