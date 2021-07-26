#ifndef _PRIMITIVES_H_
#define _PRIMITIVES_H_

#include <list>
#include "cost.h"
#include "sample.h"

//---------------CLASS MY_VERTEX_BASE---------------------
template <class Kernel, class Vbb>
class My_vertex_base : public Vbb
{
public:
    typedef Vbb Base;
    typedef CSample<Kernel> Sample;
    typedef typename Kernel::Point_2 Point;
    typedef typename Base::Face_handle Face_handle;
    
    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Base::template Rebind_TDS<TDS2>::Other Vb2;
        typedef My_vertex_base<Kernel,Vb2> Other;
    };
        
private:
    int   m_id;
    bool  m_pinned;
    Sample* m_sample;
    Point m_relocated;
    
public:
    My_vertex_base()
    : Base()
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }
    
    My_vertex_base(const Point & p)
    : Base(p)
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }
    
    My_vertex_base(Face_handle f)
    : Base(f)
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }
    
    My_vertex_base(const Point & p, Face_handle f)
    : Base(p, f)
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }    
    
    virtual ~My_vertex_base() { }
    
    int  id() const { return m_id; }
    int& id() { return m_id; }
    
    bool  pinned() const { return m_pinned; }
    bool& pinned() { return m_pinned; }
    
    Sample* get_sample() const { return m_sample; }
    void set_sample(Sample* sample) { m_sample = sample; }
    
    const Point& relocated() const { return m_relocated; }
    Point& relocated() { return m_relocated; }
};

//---------------CLASS MY_FACE_BASE---------------------
template <class Kernel, class Fbb>
class My_face_base : public Fbb
{
public:
    typedef Fbb Base;
    typedef typename Base::Vertex_handle Vertex_handle;
    typedef typename Base::Face_handle   Face_handle;
    
    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Base::template Rebind_TDS<TDS2>::Other Fb2;
        typedef My_face_base<Kernel,Fb2> Other;
    };
    
    typedef typename Kernel::FT FT;
    typedef CCost<FT> Cost;
    typedef CSample<Kernel> Sample;
    typedef std::list<Sample*> Sample_list;
    
private:
    Sample_list m_samples[3];
    FT m_mass[3];

    Cost m_cost0[3];
    Cost m_cost1[3];
    int  m_plan[3];

public:
    My_face_base()
    : Base()
    {
        init();
    }
    
    My_face_base(Vertex_handle v1,
                 Vertex_handle v2,
                 Vertex_handle v3)
    : Base(v1,v2,v3)
    {
        init();
    }
    
    My_face_base(Vertex_handle v1,
                 Vertex_handle v2,
                 Vertex_handle v3,
                 Face_handle f1,
                 Face_handle f2,
                 Face_handle f3)
    : Base(v1,v2,v3,f1,f2,f3)
    {
        init();
    }
    
    My_face_base(Face_handle f)
    : Base(f)
    {
        m_samples[0] = f->samples(0);
        m_samples[1] = f->samples(1);
        m_samples[2] = f->samples(2);

        m_mass[0] = f->mass(0);
        m_mass[1] = f->mass(1);
        m_mass[2] = f->mass(2);
        
        m_cost0[0] = f->vertex_cost(0);
        m_cost0[1] = f->vertex_cost(1);
        m_cost0[2] = f->vertex_cost(2);

        m_cost1[0] = f->edge_cost(0);
        m_cost1[1] = f->edge_cost(1);
        m_cost1[2] = f->edge_cost(2);
        
        m_plan[0] = f->plan(0);
        m_plan[1] = f->plan(1);
        m_plan[2] = f->plan(2);        
    }
    
    virtual ~My_face_base()
    {
        clean_all_samples();
    }
    
    void init()
    {
        m_mass[0] = 0.0; 
        m_mass[1] = 0.0;
        m_mass[2] = 0.0;
        
        m_cost0[0] = Cost();
        m_cost0[1] = Cost();
        m_cost0[2] = Cost();

        m_cost1[0] = Cost();
        m_cost1[1] = Cost();
        m_cost1[2] = Cost();
        
        m_plan[0] = 0;
        m_plan[1] = 0;
        m_plan[2] = 0;        
    }
    
    const int plan(int edge) const { return m_plan[edge]; }
    int& plan(int edge) { return m_plan[edge]; }
    
    const FT& mass(int edge) const { return m_mass[edge]; }
    FT& mass(int edge) { return m_mass[edge]; }

    const Cost& vertex_cost(int edge) const { return m_cost0[edge]; }
    Cost& vertex_cost(int edge) { return m_cost0[edge]; }

    const Cost& edge_cost(int edge) const { return m_cost1[edge]; }
    Cost& edge_cost(int edge) { return m_cost1[edge]; }

    const Cost& cost(int edge) const
    {
        if (plan(edge) == 0) return vertex_cost(edge);
        return edge_cost(edge);
    }

    const bool ghost(int edge) const 
    {
        if (mass(edge) == 0.0) return true;
        if (plan(edge) == 0) return true;
        return false;
    }

    const Sample_list& samples(int edge) const { return m_samples[edge]; }
    Sample_list& samples(int edge) { return m_samples[edge]; }

    void add_sample(int edge, Sample* sample) 
    { 
        m_samples[edge].push_back(sample);
    }
    
    void clean_samples(int edge)
    {
        m_samples[edge].clear();
    }
    
    void clean_all_samples() 
    {
        for (int i = 0; i < 3; ++i)
            clean_samples(i);
    }
};

//---------------STRUCT LESS VERTEX_HANDLE---------------------
template <class T>
struct less_Vertex_handle
{
    bool operator() (const T& a, const T& b) const
    {
        return (a->id() < b->id());
    }
};

//---------------STRUCT LESS EDGE---------------------
template <class T>
struct less_Edge
{
    void get_vertices_id(const T& a, int& i, int& j) const
    {
        i = a.first->vertex( (a.second+1)%3 )->id();
        j = a.first->vertex( (a.second+2)%3 )->id();
        if (i > j) std::swap(i, j);
    }
    
    bool operator() (const T& a, const T& b) const
    {
        int a0, a1;
        get_vertices_id(a, a0, a1);
        
        int b0, b1;
        get_vertices_id(b, b0, b1);
        
        if (a0 < b0) return true;
        if (a0 > b0) return false;
        if (a1 < b1) return true;
        return false;
    }    
};

//---------------STRUCT LESS FACE_HANDLE---------------------
template <class T>
struct less_Face_handle
{
    void get_vertices_id(const T& face, int& a, int& b, int& c) const
    {
        a = face->vertex(0)->id();
        b = face->vertex(1)->id();
        c = face->vertex(2)->id();
    }
    
    bool operator() (const T& a, const T& b) const
    {
        int a0, a1, a2;
        get_vertices_id(a, a0, a1, a2);
        
        int b0, b1, b2;
        get_vertices_id(b, b0, b1, b2);
        
        if (a0 < b0) return true;
        if (a0 > b0) return false;
        
        if (a1 < b1) return true;
        if (a1 > b1) return false;
        
        if (a2 < b2) return true;
        return false;
    }
};

#endif // _PRIMITIVES_H_
