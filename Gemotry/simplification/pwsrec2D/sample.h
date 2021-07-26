#ifndef _SAMPLE_H
#define _SAMPLE_H

//---------------CLASS CSAMPLE---------------------
template <class Kernel>
class CSample
{
public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point;
    
private:
    Point m_point;
    FT m_mass;
    
    FT m_dist2_to_edge;
    FT m_coordinate;
    
    FT m_backup_dist2;
    FT m_backup_coord;
    
public:
    CSample(const Point& point, 
			const FT mass = 1.0)
    {
        m_mass  = mass;
        m_point = point;
        
        m_dist2_to_edge = 0.0;
        m_coordinate    = 0.0;
        
        m_backup_dist2 = 0.0;
        m_backup_coord = 0.0;
    }

    CSample(const CSample& sample)
    {
        m_mass  = sample.mass();
        m_point = sample.point();
        
        m_dist2_to_edge = 0.0;
        m_coordinate    = 0.0;
        
        m_backup_dist2 = 0.0;
        m_backup_coord = 0.0;
    }
    
    ~CSample() { }
    
    const Point& point() const { return m_point; } 
    Point& point() { return m_point; } 
    
    const FT& mass() const { return m_mass; }
    FT& mass() { return m_mass; }
    
    const FT& distance2() const { return m_dist2_to_edge; }
    FT& distance2() { return m_dist2_to_edge; }
    
    const FT& coordinate() const { return m_coordinate; }
    FT& coordinate() { return m_coordinate; }
    
    void backup()
    {
        m_backup_dist2 = m_dist2_to_edge;
        m_backup_coord = m_coordinate;
    }
    
    void restore()
    {
        m_dist2_to_edge = m_backup_dist2;
        m_coordinate = m_backup_coord;
    }
};

//---------------CLASS CPSAMPLE---------------------
template <class	Sample>
class CPSample
{
public:
    typedef typename Sample::FT FT;
    
private:
    Sample* m_sample;
    FT m_priority;
    
public:
    CPSample(Sample* sample, const FT priority = 0.0)
    {
        m_sample   = sample;
        m_priority = priority;
    }
    
    CPSample(const CPSample& psample) 
    {
        m_sample   = psample.sample();
        m_priority = psample.priority();
    }
    
    ~CPSample() { }
    
    CPSample& operator = (const CPSample& psample)
    {
        m_sample   = psample.sample();
        m_priority = psample.priority();
        return *this;
    }
    
    Sample* sample() const { return m_sample; }
    
    const FT priority() const { return m_priority; }
};

//---------------STRUCT GREATER---------------------
template <class T>
struct greater_priority
{
    bool operator() (const T& a, const T& b) const
    {
        return ( a.priority() > b.priority() );
    }
};

#endif // SAMPLE_H
