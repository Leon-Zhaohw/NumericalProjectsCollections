//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Neil Molino, Igor Neverov, Duc Nguyen, Avi Robinson-Mosher,
//     Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX
//#####################################################################
#ifndef __BOX__
#define __BOX__

#include <Geometry/GEOMETRY_FORWARD.h>
#include <Geometry/GEOMETRY_POLICY.h>
#include <Arrays/ARRAYS_FORWARD.h>
#include <Grids/GRID_POLICY.h>
#include <Grids/POLICY_UNIFORM.h>
#include <Math_Tools/clamp.h>
#include <Math_Tools/max.h>
#include <Math_Tools/min.h>
#include <Math_Tools/ZERO.h>
#include <Matrices_And_Vectors/VECTOR_3D.h>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/less.hpp>
namespace PhysBAM{

template<class TV> struct IS_SCALAR_BLOCK<BOX<TV> >:public IS_SCALAR_BLOCK<TV>{};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<BOX<TV>,RW>:public mpl::false_{}; // required since memory format differs from disk format

template<class TV>
class BOX
{
    typedef typename TV::SCALAR T;
    struct UNUSABLE{typedef UNUSABLE UNIFORM_GRID;typedef UNUSABLE ARRAYS_SCALAR;};
    typedef mpl::and_<mpl::less<mpl::int_<0>,mpl::int_<TV::m> >,mpl::less<mpl::int_<TV::m>,mpl::int_<4> > > grid_valid;
    typedef typename IF<grid_valid::value,GRID_POLICY<TV>,UNUSABLE>::TYPE::UNIFORM_GRID T_UNIFORM_GRID;
    typedef typename IF<grid_valid::value,POLICY_UNIFORM<TV,T_UNIFORM_GRID>,UNUSABLE>::TYPE::ARRAYS_SCALAR T_ARRAYS_SCALAR;
public:
    template<class T2> struct REBIND{typedef BOX<T2> TYPE;};
    typedef T SCALAR;
    typedef TV VECTOR_T;
    enum WORKAROUND {d=TV::dimension};

    TV min_corner,max_corner;

    BOX()
        :min_corner(-TV::All_Ones_Vector()),max_corner(TV::All_Ones_Vector())
    {}

    BOX(const T xmin,const T xmax)
        :min_corner(xmin),max_corner(xmax)
    {
        STATIC_ASSERT(d==1);
    }

    BOX(const T xmin,const T xmax,const T ymin,const T ymax)
        :min_corner(xmin,ymin),max_corner(xmax,ymax)
    {
        STATIC_ASSERT(d==2);
    }

    BOX(const T xmin,const T xmax,const T ymin,const T ymax,const T zmin,const T zmax)
        :min_corner(xmin,ymin,zmin),max_corner(xmax,ymax,zmax)
    {
        STATIC_ASSERT(d==3);
    }

    BOX(const TV& minimum_corner,const TV& maximum_corner)
        :min_corner(minimum_corner),max_corner(maximum_corner)
    {}

    template<class T2> explicit BOX(const BOX<T2>& box)
        :min_corner(TV(box.min_corner)),max_corner(TV(box.max_corner))
    {}

    BOX(const TV& point)
        :min_corner(point),max_corner(point)
    {}

    BOX(const BOX<TV>& box,const FRAME<VECTOR<T,1> >& frame) // allow 1d boxes to be used as oriented boxes
        :min_corner(box.min_corner),max_corner(box.max_corner)
    {
        STATIC_ASSERT(d==1);
    }

    BOX<TV> Axis_Aligned_Bounding_Box() const
    {return *this;}

    static BOX<TV> Unit_Box()
    {return BOX<TV>(TV(),TV::All_Ones_Vector());}

    static BOX<TV> Zero_Box()
    {return BOX<TV>(TV(),TV());}

    static BOX<TV> Empty_Box()
    {return BOX<TV>((T)FLT_MAX*TV::All_Ones_Vector(),-(T)FLT_MAX*TV::All_Ones_Vector());}

    static BOX<TV> Full_Box()
    {return BOX<TV>(-(T)FLT_MAX*TV::All_Ones_Vector(),(T)FLT_MAX*TV::All_Ones_Vector());}

    bool Empty() const
    {return !min_corner.All_Less_Equal(max_corner);}

    bool operator==(const BOX<TV>& r) const
    {return min_corner==r.min_corner && max_corner==r.max_corner;}

    bool operator!=(const BOX<TV>& r) const
    {return !(*this==r);}

    BOX<TV> operator-() const
    {return BOX<TV>(-max_corner,-min_corner);}

    BOX<TV>& operator+=(const BOX<TV>& r)
    {min_corner+=r.min_corner;max_corner+=r.max_corner;return *this;}

    BOX<TV>& operator-=(const BOX<TV>& r)
    {min_corner-=r.max_corner;max_corner-=r.min_corner;return *this;}

    BOX<TV> operator+(const BOX<TV>& r) const
    {return BOX<TV>(min_corner+r.min_corner,max_corner+r.max_corner);}

    BOX<TV> operator-(const BOX<TV>& r) const
    {return BOX<TV>(min_corner-r.max_corner,max_corner-r.min_corner);}

    BOX<TV> operator*(const T a) const
    {return a>=0?BOX<TV>(min_corner*a,max_corner*a):BOX<TV>(max_corner*a,min_corner*a);}

    BOX<TV>& operator*=(const T a)
    {return *this=*this*a;}

    BOX<TV> operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    BOX<TV>& operator/=(const T a)
    {return *this=*this/a;}

    TV Edge_Lengths() const
    {return max_corner-min_corner;}

    TV Center() const
    {return (T).5*(min_corner+max_corner);}

    TV Minimum_Corner() const
    {return min_corner;}

    TV Maximum_Corner() const
    {return max_corner;}

    void Corners(ARRAYS_1D<TV>& corners) const
    {STATIC_ASSERT(d==1);corners.Resize(0,1);corners(0)=min_corner;corners(1)=max_corner;}

    void Corners(ARRAYS_2D<TV>& corners) const
    {STATIC_ASSERT(d==2);corners.Resize(0,1,0,1);
    for(int i=0;i<=1;i++) for(int j=0;j<=1;j++)
        corners(i,j)=TV(i?max_corner.x:min_corner.x,j?max_corner.y:min_corner.y);}

    void Corners(ARRAYS_3D<TV>& corners) const
    {STATIC_ASSERT(d==3);corners.Resize(0,1,0,1,0,1);
    for(int i=0;i<=1;i++) for(int j=0;j<=1;j++) for(int k=0;k<=1;k++)
        corners(i,j,k)=TV(i?max_corner.x:min_corner.x,j?max_corner.y:min_corner.y,k?max_corner.z:min_corner.z);}

    T Size() const // assumes nonnegative entries
    {return Edge_Lengths().Product();}

    T Robust_Size() const
    {return Empty()?(T)0:Size();}

    T Surface_Area() const
    {STATIC_ASSERT(d==3);VECTOR<T,3> size(Edge_Lengths());return 2*(size.x*(size.y+size.z)+size.y*size.z);}

    void Reset_Bounds(const TV& point)
    {min_corner=max_corner=point;}

    void Enlarge_To_Include_Point(const TV& point)
    {min_corner=TV::Componentwise_Min(min_corner,point);max_corner=TV::Componentwise_Max(max_corner,point);}

    void Enlarge_Nonempty_Box_To_Include_Point(const TV& point)
    {assert(!Empty());for(int i=1;i<=d;i++) if(point(i)<min_corner(i)) min_corner(i)=point(i);else if(point(i)>max_corner(i)) max_corner(i)=point(i);}

    void Enlarge_Nonempty_Box_To_Include_Points(const TV& p1,const TV& p2)
    {Enlarge_Nonempty_Box_To_Include_Point(p1);Enlarge_Nonempty_Box_To_Include_Point(p2);}

    void Enlarge_Nonempty_Box_To_Include_Points(const TV& p1,const TV& p2,const TV& p3)
    {Enlarge_Nonempty_Box_To_Include_Point(p1);Enlarge_Nonempty_Box_To_Include_Point(p2);Enlarge_Nonempty_Box_To_Include_Point(p3);}

    template<class T_ARRAY>
    void Enlarge_Nonempty_Box_To_Include_Points(const T_ARRAY& points)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    for(int i=1;i<=points.Size();i++) Enlarge_Nonempty_Box_To_Include_Point(points(i));}

    void Enlarge_To_Include_Box(const BOX<TV>& box)
    {min_corner=TV::Componentwise_Min(min_corner,box.min_corner);max_corner=TV::Componentwise_Max(max_corner,box.max_corner);}

    void Change_Size(const T delta)
    {min_corner-=delta;max_corner+=delta;}

    void Change_Size(const TV& delta)
    {min_corner-=delta;max_corner+=delta;}

    BOX<TV> Thickened(const T thickness_over_two) const
    {return BOX<TV>(min_corner-thickness_over_two,max_corner+thickness_over_two);}

    static BOX<TV> Combine(const BOX<TV>& box1,const BOX<TV>& box2)
    {return BOX<TV>(TV::Componentwise_Min(box1.min_corner,box2.min_corner),TV::Componentwise_Max(box1.max_corner,box2.max_corner));}

    static BOX<TV> Intersect(const BOX<TV>& box1,const BOX<TV>& box2) // assumes nonnegative entries
    {return BOX<TV>(TV::Componentwise_Max(box1.min_corner,box2.min_corner),TV::Componentwise_Min(box1.max_corner,box2.max_corner));}

    void Scale_About_Center(const T factor)
    {TV center=(T).5*(min_corner+max_corner),length_over_two=factor*(T).5*(max_corner-min_corner);min_corner=center-length_over_two;max_corner=center+length_over_two;}

    void Scale_About_Center(const TV factor)
    {TV center=(T).5*(min_corner+max_corner),length_over_two=factor*(T).5*(max_corner-min_corner);min_corner=center-length_over_two;max_corner=center+length_over_two;}

    void Scale_About_Center(const T x_factor,const T y_factor)
    {STATIC_ASSERT(d==2);Scale_About_Center(TV(x_factor,y_factor));}

    void Scale_About_Center(const T x_factor,const T y_factor,const T z_factor)
    {STATIC_ASSERT(d==3);Scale_About_Center(TV(x_factor,y_factor,z_factor));}

    bool Lazy_Inside(const TV& location) const
    {return location.All_Greater_Equal(min_corner) && location.All_Less_Equal(max_corner);}

    bool Lazy_Inside_Half_Open(const TV& location) const
    {return location.All_Greater_Equal(min_corner) && location.All_Less(max_corner);}

    bool Inside(const TV& location,const T thickness_over_two) const
    {return Thickened(-thickness_over_two).Lazy_Inside(location);}

    bool Inside(const TV& location,const ZERO thickness_over_two) const
    {return Lazy_Inside(location);}

    bool Lazy_Outside(const TV& location) const
    {return !Lazy_Inside(location);}

    bool Outside(const RAY<TV>& ray,const T thickness_over_two=0) const
    {STATIC_ASSERT(d==3);return Thickened(thickness_over_two).Lazy_Outside(ray);}

    bool Outside(const TV& location,const T thickness_over_two) const
    {return Thickened(thickness_over_two).Lazy_Outside(location);}

    bool Outside(const TV& location,const ZERO thickness_over_two) const
    {return Lazy_Outside(location);}

    bool Boundary(const TV& location,const T thickness_over_two) const
    {bool strict_inside=location.All_Greater(min_corner+thickness_over_two) && location.All_Less(max_corner-thickness_over_two);
    return !strict_inside && !Outside(location,thickness_over_two);}

    TV Clamp(const TV& location) const
    {return clamp(location,min_corner,max_corner);}

    T Clamp(const T& location) const
    {STATIC_ASSERT(d==1);return Clamp(TV(location)).x;}

    void Enlarge_By_Sign(const TV& v)
    {for(int i=1;i<=d;i++) if(v(i)>0) max_corner(i)+=v(i);else min_corner(i)+=v(i);}

    void Project_Points_Onto_Line(const TV& direction,T& line_min,T& line_max) const
    {line_min=line_max=TV::Dot_Product(direction,min_corner);TV e=direction*(max_corner-min_corner);
    for(int i=1;i<=d;i++) if(e(i)>0) line_max+=e(i);else line_min+=e(i);}

    TV Point_From_Normalized_Coordinates(const TV& weights) const
    {return min_corner+weights*(max_corner-min_corner);}

    bool Contains(const BOX<TV>& box) const
    {return min_corner.All_Less_Equal(box.min_corner) && max_corner.All_Greater_Equal(box.max_corner);}

    bool Lazy_Intersection(const BOX<TV>& box) const
    {return min_corner.All_Less_Equal(box.max_corner) && max_corner.All_Greater_Equal(box.min_corner);}

    bool Intersection(const BOX<TV>& box,const T thickness_over_two) const
    {return Thickened(thickness_over_two).Lazy_Intersection(box);}

    bool Intersection(const BOX<TV>& box,const ZERO thickness_over_two) const
    {return Lazy_Intersection(box);}

    bool Intersection(const BOX<TV>& box) const
    {return Lazy_Intersection(box);}

    T Intersection_Area(const BOX<TV>& box) const
    {return Intersect(*this,box).Robust_Size();}

    static BOX<TV> Bounding_Box(const TV& p1,const TV& p2)
    {BOX<TV> box(p1);box.Enlarge_Nonempty_Box_To_Include_Point(p2);return box;}

    static BOX<TV> Bounding_Box(const TV& p1,const TV& p2,const TV& p3)
    {BOX<TV> box(p1);box.Enlarge_Nonempty_Box_To_Include_Points(p2,p3);return box;}

    static BOX<TV> Bounding_Box(const TV& p1,const TV& p2,const TV& p3,const TV& p4)
    {BOX<TV> box(p1);box.Enlarge_Nonempty_Box_To_Include_Points(p2,p3,p4);return box;}

    template<class T_ARRAY>
    static BOX<TV> Bounding_Box(const T_ARRAY& points)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    if(!points.Size()) return Empty_Box();
    BOX<TV> box(points(1));for(int i=2;i<=points.Size();i++) box.Enlarge_Nonempty_Box_To_Include_Point(points(i));return box;}

    BOX<VECTOR<T,d-1> > Get_Horizontal_Box() const
    {return BOX<VECTOR<T,d-1> >(min_corner.Horizontal_Vector(),max_corner.Horizontal_Vector());}

    BOX<VECTOR<T,d-1> > Get_Vertical_Box() const
    {STATIC_ASSERT(d==2);return BOX<VECTOR<T,d-1> >(min_corner.Vertical_Vector(),max_corner.Vertical_Vector());}

    BOX<VECTOR<T,d-1> > Remove_Dimension(int dimension) const
    {return BOX<VECTOR<T,d-1> >(min_corner.Remove_Index(dimension),max_corner.Remove_Index(dimension));}

    const BOX& Bounding_Box() const // for templatization purposes
    {return *this;}

    VECTOR<T,TV::dimension-1> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,TV::dimension-1>();}

    template<class RW>
    void Read(std::istream& input_stream)
    {for(int i=1;i<=d;i++) Read_Binary<RW>(input_stream,min_corner(i),max_corner(i));}

    template<class RW>
    void Write(std::ostream& output_stream) const
    {for(int i=1;i<=d;i++) Write_Binary<RW>(output_stream,min_corner(i),max_corner(i));}

//#####################################################################
    T Halfspace_Intersection_Size(const PLANE<T>& halfspace,VECTOR<T,3>* centroid=0) const;
    T Halfspace_Intersection_Size(const LINE_2D<T>& line,VECTOR<T,2>* centroid=0) const;
    T Halfspace_Intersection_Size(const VECTOR<T,1>& hyperplane,VECTOR<T,1>* center=0) {STATIC_ASSERT(d==1);PHYSBAM_NOT_IMPLEMENTED();} // TODO: unify hyperplane types
    bool Lazy_Intersection(RAY<VECTOR<T,3> >& ray,const T box_enlargement=0) const;
    bool Lazy_Outside(const RAY<VECTOR<T,3> >& ray) const;
    bool Intersection(RAY<VECTOR<T,3> >& ray,const T thickness_over_two=0) const;
    bool Intersection(RAY<VECTOR<T,2> >& ray,const T thickness_over_two=0,const T segment_intersect_epsilon=0) const;
    bool Intersection(RAY<VECTOR<T,1> >& ray,const T thickness_over_two=0) const;
    bool Intersection(const TRIANGLE_3D<T>& triangle,const T thickness_over_two=0) const;
    bool Intersection(const SEGMENT_2D<T>& segment,const T thickness_over_two=0) const;
    bool Intersection(const SPHERE<TV>& sphere,const T thickness_over_two=0) const;
    bool Get_Intersection_Range(const RAY<VECTOR<T,3> >& ray,T& start_t,T& end_t) const;
    TV Normal(const int aggregate) const;
    TV Normal(const TV& X) const;
    TV Surface(const TV& X) const;
    T Signed_Distance(const TV& X) const;
    void Calculate_Signed_Distance_Function(const T_UNIFORM_GRID& grid,T_ARRAYS_SCALAR& phi) const;
    static std::string Name();
    TRIANGULATED_SURFACE<T>* Generate_Triangles() const;
//#####################################################################
};
template<class TV>
inline BOX<TV> operator+(const TV& a,const BOX<TV>& b)
{return BOX<TV>(a+b.min_corner,a+b.max_corner);}

template<class TV>
inline BOX<TV> operator-(const TV& a,const BOX<TV>& b)
{return BOX<TV>(a-b.max_corner,a-b.min_corner);}

template<class TV> inline BOX<TV> operator*(const typename TV::SCALAR a,const BOX<TV>& box)
{return box*a;}

template<class TV>
inline std::ostream& operator<<(std::ostream& output_stream,const BOX<TV>& box)
{for(int i=1;i<=TV::dimension;i++) output_stream<<"("<<box.min_corner(i)<<","<< box.max_corner(i)<<(i<TV::dimension?") x ":")");return output_stream;}

template<class TV>
inline std::istream& operator>>(std::istream& input_stream,BOX<TV>& box)
{for(int i=1;i<=TV::dimension;i++) input_stream>>box.min_corner(i)>>box.max_corner(i);return input_stream;}
}
#endif
