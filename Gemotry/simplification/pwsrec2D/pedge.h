#ifndef _PEDGE_H_
#define _PEDGE_H_

//---------------CLASS CPEDGE---------------------
template <class FT, class Edge, class Vertex_handle, class Face_handle>
class CPEdge
{
protected:
    Edge m_edge;
    Vertex_handle m_source;
    Vertex_handle m_target;

    FT  m_before_cost;
    FT  m_after_cost;
    
public:
    CPEdge()
    {
        m_edge = Edge(Face_handle(), 0);
        m_source = Vertex_handle();
        m_target = Vertex_handle();
        
        m_before_cost = 0.0;
        m_after_cost = 0.0;
    }
    
    CPEdge(const CPEdge& pedge) 
    {
        m_edge = pedge.edge();
        m_source = pedge.source();
        m_target = pedge.target();
        
        m_before_cost = pedge.before();
        m_after_cost = pedge.after();
    }
    
    CPEdge(const Edge& edge, 
           const FT before, 
           const FT after)
    {
        m_edge = edge;
        get_vertices();
        
        m_before_cost = before;
        m_after_cost = after;
    }
    
    CPEdge(const Edge& edge,
           const FT priority = 0.0)
    {
        m_edge = edge;
        get_vertices();
        
        m_before_cost = 0.0;
        m_after_cost = priority;
    }
    
    CPEdge(Vertex_handle source, Vertex_handle target)
    {
        m_edge = Edge(Face_handle(), 0);
        m_source = source;
        m_target = target;

        m_before_cost = 0.0;
        m_after_cost = 0.0;
    }
    
    virtual ~CPEdge() { }
    
    CPEdge& operator = (const CPEdge& pedge) 
    {
        m_edge = pedge.edge();
        m_source = pedge.source();
        m_target = pedge.target();

        m_before_cost = pedge.before();
        m_after_cost = pedge.after();
        
        return *this;
    }
    
    bool operator == (const CPEdge& pedge) const
    {
        return (m_source->id() == pedge.source()->id() && 
                m_target->id() == pedge.target()->id());
    }
    
    bool operator < (const CPEdge& pedge) const
    {
        if (m_source->id() < pedge.source()->id()) return true;
        if (m_source->id() > pedge.source()->id()) return false;
        
        if (m_target->id() < pedge.target()->id()) return true;
        return false;
    }
    
    const Edge& edge() const { return m_edge; }
    
    const Vertex_handle& source() const { return m_source; }
    
    const Vertex_handle& target() const { return m_target; }
    
    const FT before() const { return m_before_cost; }
    
    const FT after() const { return m_after_cost; }

    const FT priority() const { return after() - before(); }
    
protected:    
    void get_vertices()
    {
        int index = m_edge.second;
        m_source = m_edge.first->vertex( (index+1)%3 );
        m_target = m_edge.first->vertex( (index+2)%3 );
    }
};

//---------------CLASS CPQUEUE---------------------
template <class T>
class CPQueue : public DynamicPriorityQueue<T> {
public:
    CPQueue() { }

    ~CPQueue() { }

    // a < b
    bool compare(const T& a, const T& b) const 
    {
        return (a.priority() < b.priority());
    }
};

#endif // _PHALFEDGE_H_
