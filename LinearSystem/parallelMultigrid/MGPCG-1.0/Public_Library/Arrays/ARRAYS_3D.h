//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Duc Nguyen, Craig Schroeder, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAYS_3D  
//#####################################################################
#ifndef __ARRAYS_3D__
#define __ARRAYS_3D__

#include <Arrays/ARRAYS_ND_BASE.h>
#include <Grids/GRID_3D.h>
#include <Matrices_And_Vectors/VECTOR.h>
#include <Utilities/EXCEPTIONS.h>
#include <boost/format.hpp>
namespace PhysBAM{

template<class T>
class ARRAYS_3D<T,1>:public ARRAYS_ND_BASE<T,ARRAYS_3D<T> >
{
public:
    enum WORKAROUND1 {dimension=3};
    enum WORKAROUND2 {length=1};
    template<class T2> struct REBIND{typedef ARRAYS_3D<T2> TYPE;};
    template<int length2> struct REBIND_LENGTH{typedef ARRAYS_3D<T,length2> TYPE;};
    typedef T ELEMENT;

    using ARRAYS_ND_BASE<T,ARRAYS_3D>::array; // one-dimensional data storage

    int m_start,m_end,n_start,n_end,mn_start,mn_end; // starting and ending values for x, y and z direction
    int m,n,mn; // used in overloaded () operator
private:
    int n_times_mn,n_times_mn_plus_mn,n_times_mn_minus_mn,n_times_mn_plus_one,n_times_mn_minus_one,mn_plus_one,mn_minus_one,n_times_mn_plus_mn_plus_one,n_times_mn_minus_mn_plus_one,
        n_times_mn_plus_mn_minus_one,n_times_mn_minus_mn_minus_one;
    T* base_pointer;
public:

    ARRAYS_3D()
        :m_start(1),m_end(0),n_start(1),n_end(0),mn_start(1),mn_end(0),m(0),n(0),mn(0)
    {
        Calculate_Acceleration_Constants();
    }

    ARRAYS_3D(const BOX<VECTOR<int,3> >& domain,const bool initialize_using_default_constructor=true)
        :m_start(domain.min_corner.x),m_end(domain.max_corner.x),n_start(domain.min_corner.y),n_end(domain.max_corner.y),mn_start(domain.min_corner.z),mn_end(domain.max_corner.z),
        m(m_end-m_start+1),n(n_end-n_start+1),mn(mn_end-mn_start+1)
    {
        assert(m>=0 && n>=0 && mn>=0);int size=m*n*mn;
        {RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    ARRAYS_3D(const int m_start_input,const int m_end_input,const int n_start_input,const int n_end_input,const int mn_start_input,const int mn_end_input,
        const bool initialize_using_default_constructor=true)
        :m_start(m_start_input),m_end(m_end_input),n_start(n_start_input),n_end(n_end_input),mn_start(mn_start_input),mn_end(mn_end_input),
        m(m_end-m_start+1),n(n_end-n_start+1),mn(mn_end-mn_start+1)
    {
        assert(m>=0 && n>=0 && mn>=0);int size=m*n*mn;
        {RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    template<class T2>
    ARRAYS_3D(const GRID_3D<T2>& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :m_start(1-ghost_cells),m_end(grid.m+ghost_cells),n_start(1-ghost_cells),n_end(grid.n+ghost_cells),mn_start(1-ghost_cells),mn_end(grid.mn+ghost_cells),
        m(m_end-m_start+1),n(n_end-n_start+1),mn(mn_end-mn_start+1)
    {
        {int size=m*n*mn;RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    ARRAYS_3D(const ARRAYS_3D& old_array,const bool initialize_with_old_array=true)
        :ARRAYS_ND_BASE<T,ARRAYS_3D>(old_array.array.Size()),m_start(old_array.m_start),m_end(old_array.m_end),n_start(old_array.n_start),n_end(old_array.n_end),mn_start(old_array.mn_start),
        mn_end(old_array.mn_end),m(old_array.m),n(old_array.n),mn(old_array.mn)
    {
        Calculate_Acceleration_Constants();
        if(initialize_with_old_array) array=old_array.array;
    }

    BOX<VECTOR<int,3> > Domain_Indices() const
    {return BOX<VECTOR<int,3> >(m_start,m_end,n_start,n_end,mn_start,mn_end);}

    void Clean_Memory()
    {Resize(1,0,1,0,1,0,false,false);}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(int i=1;i<=array.Size();i++) delete array(i);Clean_Memory();}

    void Calculate_Acceleration_Constants()
    {T* array_pointer=array.Get_Array_Pointer();
    n_times_mn=n*mn;n_times_mn_plus_mn=n_times_mn+mn;n_times_mn_minus_mn=n_times_mn-mn;n_times_mn_plus_one=n_times_mn+1;n_times_mn_minus_one=n_times_mn-1;mn_plus_one=mn+1;
    mn_minus_one=mn-1;n_times_mn_plus_mn_plus_one=n_times_mn_plus_mn+1;n_times_mn_minus_mn_plus_one=n_times_mn_minus_mn+1;n_times_mn_plus_mn_minus_one=n_times_mn_plus_mn-1;
    n_times_mn_minus_mn_minus_one=n_times_mn_minus_mn-1;
    base_pointer=array_pointer?array_pointer-((m_start*n+n_start)*mn+mn_start):0;}

    ARRAYS_3D& operator=(const ARRAYS_3D& source)
    {if(array.Size()!=source.array.Size()){
        delete[] array.Get_Array_Pointer();
        RAW_ARRAY<T> new_array(source.array.Size(),new T[source.array.Size()]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);}
    else if(this==&source) return *this;
    m=source.m;m_start=source.m_start;m_end=source.m_end;n=source.n;n_start=source.n_start;n_end=source.n_end;mn=source.mn;
    mn_start=source.mn_start;mn_end=source.mn_end;
    Calculate_Acceleration_Constants();
    array=source.array;return *this;}

    T& operator()(const int i,const int j,const int ij)
    {assert(m_start<=i && i<=m_end);assert(n_start<=j && j<=n_end);assert(mn_start<=ij && ij<=mn_end);
    return base_pointer[(i*n+j)*mn+ij];}

    const T& operator()(const int i,const int j,const int ij) const
    {assert(m_start<=i && i<=m_end);assert(n_start<=j && j<=n_end);assert(mn_start<=ij && ij<=mn_end);
    return base_pointer[(i*n+j)*mn+ij];}

    T& operator()(const VECTOR<int,3>& index)
    {assert(m_start<=index.x && index.x<=m_end);assert(n_start<=index.y && index.y<=n_end);assert(mn_start<=index.z && index.z<=mn_end);
    return base_pointer[(index.x*n+index.y)*mn+index.z];}

    const T& operator()(const VECTOR<int,3>& index) const
    {assert(m_start<=index.x && index.x<=m_end);assert(n_start<=index.y && index.y<=n_end);assert(mn_start<=index.z && index.z<=mn_end);
    return base_pointer[(index.x*n+index.y)*mn+index.z];}

    T& operator()(const int k,const int i,const int j,const int ij)
    {assert(k==1);
    return (*this)(i,j,ij);}

    const T& operator()(const int k,const int i,const int j,const int ij) const
    {assert(k==1);
    return (*this)(i,j,ij);}

    T& operator()(const int k,const VECTOR<int,3>& index)
    {assert(k==1);
    return (*this)(index);}

    const T& operator()(const int k,const VECTOR<int,3>& index) const
    {assert(k==1);
    return (*this)(index);}

    bool Valid_Index(const VECTOR<int,3>& index) const
    {return m_start<=index.x && index.x<=m_end && n_start<=index.y && index.y<=n_end && mn_start<=index.z && index.z<=mn_end;}

    bool Valid_Index(const int i,const int j,const int ij) const
    {return m_start<=i && i<=m_end && n_start<=j && j<=n_end && mn_start<=ij && ij<=mn_end;}

    int Standard_Index(const int i,const int j,const int ij) const
    {assert(m_start<=i && i<=m_end);assert(n_start<=j && j<=n_end);assert(mn_start<=ij && ij<=mn_end);
    return ((i-m_start)*n+(j-n_start))*mn+(ij-mn_start)+1;}

    int Standard_Index(const VECTOR<int,3>& index) const
    {assert(m_start<=index.x && index.x<=m_end);assert(n_start<=index.y && index.y<=n_end);assert(mn_start<=index.z && index.z<=mn_end);
    return ((index.x-m_start)*n+(index.y-n_start))*mn+(index.z-mn_start)+1;}

    int I_Plus_One(const int index) const
    {return index+n_times_mn;}

    int I_Minus_One(const int index) const
    {return index-n_times_mn;}

    int J_Plus_One(const int index) const
    {return index+mn;}

    int J_Minus_One(const int index) const
    {return index-mn;}

    int IJ_Plus_One(const int index) const
    {return index+1;}

    int IJ_Minus_One(const int index) const
    {return index-1;}

    int I_Plus_One_J_Plus_One(const int index) const
    {return index+n_times_mn_plus_mn;}

    int I_Plus_One_J_Minus_One(const int index) const
    {return index+n_times_mn_minus_mn;}

    int I_Minus_One_J_Plus_One(const int index) const
    {return index-n_times_mn_minus_mn;}

    int I_Minus_One_J_Minus_One(const int index) const
    {return index-n_times_mn_plus_mn;}

    int I_Plus_One_IJ_Plus_One(const int index) const
    {return index+n_times_mn_plus_one;}

    int I_Plus_One_IJ_Minus_One(const int index) const
    {return index+n_times_mn_minus_one;}

    int I_Minus_One_IJ_Plus_One(const int index) const
    {return index-n_times_mn_minus_one;}

    int I_Minus_One_IJ_Minus_One(const int index) const
    {return index-n_times_mn_plus_one;}

    int J_Plus_One_IJ_Plus_One(const int index) const
    {return index+mn_plus_one;}

    int J_Plus_One_IJ_Minus_One(const int index) const
    {return index+mn_minus_one;}

    int J_Minus_One_IJ_Plus_One(const int index) const
    {return index-mn_minus_one;}

    int J_Minus_One_IJ_Minus_One(const int index) const
    {return index-mn_plus_one;}

    int I_Plus_One_J_Plus_One_IJ_Plus_One(const int index) const
    {return index+n_times_mn_plus_mn_plus_one;}

    int I_Plus_One_J_Plus_One_IJ_Minus_One(const int index) const
    {return index+n_times_mn_plus_mn_minus_one;}

    int I_Plus_One_J_Minus_One_IJ_Plus_One(const int index) const
    {return index+n_times_mn_minus_mn_plus_one;}

    int I_Plus_One_J_Minus_One_IJ_Minus_One(const int index) const
    {return index+n_times_mn_minus_mn_minus_one;}

    int I_Minus_One_J_Plus_One_IJ_Plus_One(const int index) const
    {return index-n_times_mn_minus_mn_minus_one;}

    int I_Minus_One_J_Plus_One_IJ_Minus_One(const int index) const
    {return index-n_times_mn_minus_mn_plus_one;}

    int I_Minus_One_J_Minus_One_IJ_Plus_One(const int index) const
    {return index-n_times_mn_plus_mn_minus_one;}

    int I_Minus_One_J_Minus_One_IJ_Minus_One(const int index) const
    {return index-n_times_mn_plus_mn_plus_one;}

    void Resize(int m_start_new,int m_end_new,int n_start_new,int n_end_new,int mn_start_new,int mn_end_new,const bool initialize_new_elements=true,
        const bool copy_existing_elements=true,const T& initialization_value=T())
    {if(m_start_new==m_start && m_end_new==m_end && n_start_new==n_start && n_end_new==n_end && mn_start_new==mn_start && mn_end_new==mn_end) return;
    int m_new=m_end_new-m_start_new+1,n_new=n_end_new-n_start_new+1,mn_new=mn_end_new-mn_start_new+1;assert(m_new>=0 && n_new>=0 && mn_new>=0);
    int size_new=m_new*n_new*mn_new;RAW_ARRAY<T> array_new(size_new,new T[size_new]);
    if(initialize_new_elements) array_new.Fill(initialization_value);
    if(copy_existing_elements){
        int m1=PhysBAM::max(m_start,m_start_new),m2=PhysBAM::min(m_end,m_end_new),n1=PhysBAM::max(n_start,n_start_new),n2=PhysBAM::min(n_end,n_end_new),
           mn1=PhysBAM::max(mn_start,mn_start_new),mn2=PhysBAM::min(mn_end,mn_end_new);
        for(int i=m1;i<=m2;i++) for(int j=n1;j<=n2;j++) for(int ij=mn1;ij<=mn2;ij++)
            array_new((((i-m_start_new)*n_new+(j-n_start_new))*mn_new+(ij-mn_start_new))+1)=array((((i-m_start)*n+(j-n_start))*mn+(ij-mn_start))+1);}
    m_start=m_start_new;m_end=m_end_new;m=m_end-m_start+1;n_start=n_start_new;n_end=n_end_new;n=n_end-n_start+1;mn_start=mn_start_new;mn_end=mn_end_new;mn=mn_end-mn_start+1;
    delete[] array.Get_Array_Pointer();RAW_ARRAY<T>::Exchange_Arrays(array,array_new);Calculate_Acceleration_Constants();}

    template<class T2>
    void Resize(const GRID_3D<T2>& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(1-ghost_cells,grid.m+ghost_cells,1-ghost_cells,grid.n+ghost_cells,1-ghost_cells,grid.mn+ghost_cells,initialize_new_elements,copy_existing_elements,initialization_value);}

    void Resize(const BOX<VECTOR<int,3> >& box,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z,initialize_new_elements,copy_existing_elements,initialization_value);}

    // TODO: this function is extremely broken, since all the functions in the class depends on size being correct
    void Resize_In_Place(const int m_start_new,const int m_end_new,const int n_start_new,const int n_end_new,const int mn_start_new,const int mn_end_new)
    {if(array.Size()>=(m_end_new-m_start_new+1)*(n_end_new-n_start_new+1)*(mn_end_new-mn_start_new+1)){
        m_start=m_start_new;m_end=m_end_new;m=m_end-m_start+1;n_start=n_start_new;n_end=n_end_new;n=n_end-n_start+1;mn_start=mn_start_new;mn_end=mn_end_new;mn=mn_end-mn_start+1;
        Calculate_Acceleration_Constants();}
    else Resize(m_start_new,m_end_new,n_start_new,n_end_new,mn_start_new,mn_end_new,false,false);}

    void Resize_In_Place(const BOX<VECTOR<int,3> >& box)
    {Resize_In_Place(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z);}

    void Clamp(int& i,int& j,int& ij) const
    {i=clamp(i,m_start,m_end);j=clamp(j,n_start,n_end);ij=clamp(ij,mn_start,mn_end);}

    void Clamp_End_Minus_One(int& i,int& j,int& ij) const
    {i=clamp(i,m_start,m_end-1);j=clamp(j,n_start,n_end-1);ij=clamp(ij,mn_start,mn_end-1);}

    void Clamp_End_Minus_Two(int& i,int& j,int& ij) const
    {i=clamp(i,m_start,m_end-2);j=clamp(j,n_start,n_end-2);ij=clamp(ij,mn_start,mn_end-3);}

    void Clamp_End_Minus_Three(int& i,int& j,int& ij) const
    {i=clamp(i,m_start,m_end-3);j=clamp(j,n_start,n_end-3);ij=clamp(ij,mn_start,mn_end-3);}

    void Clamp_Interior(int& i,int& j,int& ij) const
    {i=clamp(i,m_start+1,m_end-1);j=clamp(j,n_start+1,n_end-1);ij=clamp(ij,mn_start+1,mn_end-1);}

    void Clamp_Interior_End_Minus_One(int& i,int& j,int& ij) const
    {i=clamp(i,m_start+1,m_end-2);j=clamp(j,n_start+1,n_end-2);ij=clamp(ij,mn_start+1,mn_end-2);}

    static void Get(ARRAYS_3D& new_copy,const ARRAYS_3D& old_copy)
    {if(&old_copy!=&new_copy) ARRAYS_3D::Put(old_copy,new_copy,new_copy.m_start,new_copy.m_end,new_copy.n_start,new_copy.n_end,new_copy.mn_start,new_copy.mn_end);}

    static void Shifted_Get(ARRAYS_3D& new_copy,const ARRAYS_3D& old_copy,const VECTOR<int,3>& shift)
    {if(shift==VECTOR<int,3>()) Get(new_copy,old_copy);
    else for(int i=new_copy.m_start;i<=new_copy.m_end;i++)for(int j=new_copy.n_start;j<=new_copy.n_end;j++)for(int ij=new_copy.mn_start;ij<=new_copy.mn_end;ij++)
        new_copy(i,j,ij)=old_copy(i+shift.x,j+shift.y,ij+shift.z);}

    static void Limited_Shifted_Get(ARRAYS_3D& new_copy,const ARRAYS_3D& old_copy,const VECTOR<int,3>& shift)
    {int m_start=PhysBAM::max(new_copy.m_start,old_copy.m_start-shift.x),m_end=PhysBAM::min(new_copy.m_end,old_copy.m_end-shift.x);
    int n_start=PhysBAM::max(new_copy.n_start,old_copy.n_start-shift.y),n_end=PhysBAM::min(new_copy.n_end,old_copy.n_end-shift.y);
    int mn_start=PhysBAM::max(new_copy.mn_start,old_copy.mn_start-shift.z),mn_end=PhysBAM::min(new_copy.mn_end,old_copy.mn_end-shift.z);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++)
        new_copy(i,j,ij)=old_copy(i+shift.x,j+shift.y,ij+shift.z);}

    static void Put(const ARRAYS_3D& old_copy,ARRAYS_3D& new_copy)
    {if(&old_copy!=&new_copy) ARRAYS_3D::Put(old_copy,new_copy,old_copy.m_start,old_copy.m_end,old_copy.n_start,old_copy.n_end,old_copy.mn_start,old_copy.mn_end);}

    static void Shifted_Put(const ARRAYS_3D& old_copy,ARRAYS_3D& new_copy,const VECTOR<int,3>& shift)
    {if(shift==VECTOR<int,3>()) Put(old_copy,new_copy);
    else for(int i=old_copy.m_start;i<=old_copy.m_end;i++)for(int j=old_copy.n_start;j<=old_copy.n_end;j++)for(int ij=old_copy.mn_start;ij<=old_copy.mn_end;ij++)
        new_copy(i+shift.x,j+shift.y,ij+shift.z)=old_copy(i,j,ij);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_3D& old_copy,ARRAYS_3D& new_copy)
    {ARRAYS_3D::Put(constant,old_copy,new_copy,old_copy.m_start,old_copy.m_end,old_copy.n_start,old_copy.n_end,old_copy.mn_start,old_copy.mn_end);}

private:
    static void Put(const ARRAYS_3D& old_copy,ARRAYS_3D& new_copy,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end)
    {assert(old_copy.m_start<=m_start&&m_end<=old_copy.m_end&&old_copy.n_start<=n_start&&n_end<=old_copy.n_end&&old_copy.mn_start<=mn_start&&mn_end<=old_copy.mn_end);
    assert(new_copy.m_start<=m_start&&m_end<=new_copy.m_end&&new_copy.n_start<=n_start&&n_end<=new_copy.n_end&&new_copy.mn_start<=mn_start&&mn_end<=new_copy.mn_end);
    int old_index=old_copy.Standard_Index(m_start,n_start,mn_start),new_index=new_copy.Standard_Index(m_start,n_start,mn_start);
    const int inner_loop_copy_size=(mn_end-mn_start+1);
    for(int i=m_start;i<=m_end;i++){
        int old_save_index=old_index,new_save_index=new_index;
        for(int j=n_start;j<=n_end;j++){
            int old_save_index=old_index,new_save_index=new_index;
            for(int t=1;t<=inner_loop_copy_size;t++) new_copy.array(new_index++)=old_copy.array(old_index++);
            old_index=old_copy.J_Plus_One(old_save_index);new_index=new_copy.J_Plus_One(new_save_index);}
        old_index=old_copy.I_Plus_One(old_save_index);new_index=new_copy.I_Plus_One(new_save_index);}}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_3D& old_copy,ARRAYS_3D& new_copy,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,
        const int mn_end)
    {assert(old_copy.m_start<=m_start&&m_end<=old_copy.m_end&&old_copy.n_start<=n_start&&n_end<=old_copy.n_end&&old_copy.mn_start<=mn_start&&mn_end<=old_copy.mn_end);
    assert(new_copy.m_start<=m_start&&m_end<=new_copy.m_end&&new_copy.n_start<=n_start&&n_end<=new_copy.n_end&&new_copy.mn_start<=mn_start&&mn_end<=new_copy.mn_end);
    int old_index=old_copy.Standard_Index(m_start,n_start,mn_start),new_index=new_copy.Standard_Index(m_start,n_start,mn_start);
    const int inner_loop_copy_size=(mn_end-mn_start+1);
    for(int i=m_start;i<=m_end;i++){
        int old_save_index=old_index,new_save_index=new_index;
        for(int j=n_start;j<=n_end;j++){
            int old_save_index=old_index,new_save_index=new_index;
            for(int t=1;t<=inner_loop_copy_size;t++) new_copy.array(new_index++)=constant*old_copy.array(old_index++);
            old_index=old_copy.J_Plus_One(old_save_index);new_index=new_copy.J_Plus_One(new_save_index);}
        old_index=old_copy.I_Plus_One(old_save_index);new_index=new_copy.I_Plus_One(new_save_index);}}

public:
    template<class TCONSTANT,class TGRID>
    static void Put_Ghost(const TCONSTANT constant,ARRAYS_3D& x,const GRID_3D<TGRID>& grid,const int ghost_cells)
    {for(int j=1-ghost_cells;j<=grid.n+ghost_cells;j++) for(int ij=1-ghost_cells;ij<=grid.mn+ghost_cells;ij++) for(int s=1;s<=ghost_cells;s++) x(1-s,j,ij)=x(grid.m+s,j,ij)=constant;
    for(int i=1;i<=grid.m;i++) for(int ij=1-ghost_cells;ij<=grid.mn+ghost_cells;ij++) for(int s=1;s<=ghost_cells;s++) x(i,1-s,ij)=x(i,grid.n+s,ij)=constant;
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int s=1;s<=ghost_cells;s++) x(i,j,1-s)=x(i,j,grid.mn+s)=constant;}

    static void Exchange_Arrays(ARRAYS_3D& a,ARRAYS_3D& b)
    {RAW_ARRAY<T>::Exchange_Arrays(a.array,b.array);
    exchange(a.m_start,b.m_start);exchange(a.m_end,b.m_end);exchange(a.m,b.m);
    exchange(a.n_start,b.n_start);exchange(a.n_end,b.n_end);exchange(a.n,b.n);
    exchange(a.mn_start,b.mn_start);exchange(a.mn_end,b.mn_end);exchange(a.mn,b.mn);
    a.Calculate_Acceleration_Constants();b.Calculate_Acceleration_Constants();}

    template<class T2>
    static bool Equal_Dimensions(const ARRAYS_3D& a,const ARRAYS_3D<T2>& b)
    {return a.m_start==b.m_start && a.m_end==b.m_end && a.n_start==b.n_start && a.n_end==b.n_end && a.mn_start==b.mn_start && a.mn_end==b.mn_end;}

    static bool Equal_Dimensions(const ARRAYS_3D& a,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end)
    {return a.m_start==m_start && a.m_end==m_end && a.n_start==n_start && a.n_end==n_end && a.mn_start==mn_start && a.mn_end==mn_end;}

    // note that these functions move the *contents* of the grid, not the grid itself
    void Move_Contents_Left(const int increment=1)
    {for(int i=m_start;i<=m_end-increment;i++) for(int j=n_start;j<=n_end;j++) for(int ij=mn_start;ij<=mn_end;ij++) (*this)(i,j,ij)=(*this)(i+increment,j,ij);}

    void Move_Contents_Right(const int increment=1)
    {for(int i=m_end;i>=m_start+increment;i--) for(int j=n_start;j<=n_end;j++) for(int ij=mn_start;ij<=mn_end;ij++) (*this)(i,j,ij)=(*this)(i-increment,j,ij);}

    void Move_Contents_Down(const int increment=1)
    {for(int i=m_start;i<=m_end;i++) for(int j=n_start;j<=n_end-increment;j++) for(int ij=mn_start;ij<=mn_end;ij++) (*this)(i,j,ij)=(*this)(i,j+increment,ij);}

    void Move_Contents_Up(const int increment=1)
    {for(int i=m_start;i<=m_end;i++) for(int j=n_end;j>=n_start+increment;j--) for(int ij=mn_start;ij<=mn_end;ij++) (*this)(i,j,ij)=(*this)(i,j-increment,ij);}

    void Move_Contents_Forward(const int increment=1)
    {for(int i=m_start;i<=m_end;i++) for(int j=n_start;j<=n_end;j++) for(int ij=mn_start;ij<=mn_end-increment;ij++) (*this)(i,j,ij)=(*this)(i,j,ij+increment);}

    void Move_Contents_Backward(const int increment=1)
    {for(int i=m_start;i<=m_end;i++) for(int j=n_start;j<=n_end;j++) for(int ij=mn_end;ij>=mn_start+increment;ij--) (*this)(i,j,ij)=(*this)(i,j,ij-increment);}

    template<class RW>
    void Read(std::istream& input)
    {Read_With_Length<RW>(input,1);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_With_Length<RW>(output,1);}

protected:
    template<class RW>
    void Read_With_Length(std::istream& input,const int length2)
    {int read_length;Read_Binary<RW>(input,read_length,m_start,m_end,n_start,n_end,mn_start,mn_end);
    if(read_length!=length2) throw READ_ERROR(str(boost::format("Read length %d not equal to %d")%read_length%length2));
    m=m_end-m_start+1;
    if(m<0) throw READ_ERROR(str(boost::format("Invalid negative array size m = %d")%m));
    n=n_end-n_start+1;
    if(n<0) throw READ_ERROR(str(boost::format("Invalid negative array size n = %d")%n));
    mn=mn_end-mn_start+1;
    if(mn<0) throw READ_ERROR(str(boost::format("Invalid negative array size mn = %d")%mn));
    delete array.Get_Array_Pointer();int size=m*n*mn;RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(array,new_array);
    Read_Binary_Array<RW>(input,array.Get_Array_Pointer(),array.Size());Calculate_Acceleration_Constants();}

    template<class RW>
    void Write_With_Length(std::ostream& output,const int length2) const
    {Write_Binary<RW>(output,length2,m_start,m_end,n_start,n_end,mn_start,mn_end);Write_Binary_Array<RW>(output,array.Get_Array_Pointer(),array.Size());}

//#####################################################################
};

template<class T,int length_input>
class ARRAYS_3D:public ARRAYS_3D<VECTOR<T,length_input> >
{
    STATIC_ASSERT(length_input>0);
public:
    enum WORKAROUND {length=length_input};
    template<class T2> struct REBIND{typedef ARRAYS_3D<T2,length> TYPE;};
    template<int length2> struct REBIND_LENGTH{typedef ARRAYS_3D<T,length2> TYPE;};

    typedef ARRAYS_3D<VECTOR<T,length> > BASE;
    typedef T ELEMENT_OF_ELEMENT;
    using BASE::m;using BASE::n;using BASE::mn;using BASE::m_start;using BASE::m_end;using BASE::n_start;using BASE::n_end;using BASE::mn_start;using BASE::mn_end;
    using BASE::operator();using BASE::Valid_Index;using BASE::array;

public:
    ARRAYS_3D()
    {}

    ARRAYS_3D(const BOX<VECTOR<int,3> >& domain,const bool initialize_using_default_constructor=true)
        :ARRAYS_3D<VECTOR<T,length> >(domain,initialize_using_default_constructor)
    {}

    ARRAYS_3D(const int m_start_input,const int m_end_input,const int n_start_input,const int n_end_input,const int mn_start_input,const int mn_end_input,
        const bool initialize_using_default_constructor=true)
        :ARRAYS_3D<VECTOR<T,length> >(m_start_input,m_end_input,n_start_input,n_end_input,mn_start_input,mn_end_input,initialize_using_default_constructor)
    {}

    template<class T2>
    ARRAYS_3D(const GRID_3D<T2>& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :ARRAYS_3D<VECTOR<T,length> >(grid,ghost_cells,initialize_using_default_constructor)
    {}

    ARRAYS_3D(const ARRAYS_3D& old_array,const bool initialize_with_old_array=true)
        :ARRAYS_3D<VECTOR<T,length> >(old_array,initialize_with_old_array)
    {}

    T& operator()(const int k,const int i,const int j,const int ij)
    {assert(1<=k && k<=length);
    return (*this)(i,j,ij)[k];}

    const T& operator()(const int k,const int i,const int j,const int ij) const
    {assert(1<=k && k<=length);
    return (*this)(i,j,ij)[k];}

    T& operator()(const int k,const VECTOR<int,3>& index)
    {assert(1<=k && k<=length);
    return (*this)(index)[k];}

    const T& operator()(const int k,const VECTOR<int,3>& index) const
    {assert(1<=k && k<=length);
    return (*this)(index)[k];}

    bool Valid_Index(const int k,const int i,const int j,const int ij) const
    {return 1<=k && k<=length && Valid_Index(i,j,ij);}

    static void Extract_Dimension(const ARRAYS_3D& old_array,ARRAYS_3D<T>& extracted_array,int dimension)
    {extracted_array.Resize(old_array.m_start,old_array.m_end,old_array.n_start,old_array.n_end,old_array.mn_start,old_array.mn_end,false,false);
    for(int i=old_array.m_start;i<=old_array.m_end;i++) for(int j=old_array.n_start;j<=old_array.n_end;j++) for(int ij=old_array.mn_start;ij<=old_array.mn_end;ij++)
        extracted_array(i,j,ij)=old_array(dimension,i,j,ij);}

    template<class T2>
    void Resize(const GRID_3D<T2>& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {ARRAYS_3D<VECTOR<T,length> >::Resize(grid,ghost_cells,initialize_new_elements,copy_existing_elements);}

    template<class T2>
    void Resize(const GRID_3D<T2>& grid,const int ghost_cells,const bool initialize_new_elements,const bool copy_existing_elements,const T& initialization_value)
    {VECTOR<T,length> initialization_vector;initialization_vector.Fill(initialization_value);
    ARRAYS_3D<VECTOR<T,length> >::Resize(grid,ghost_cells,initialize_new_elements,copy_existing_elements,initialization_vector);}

    void Fill(const T& constant)
    {array.Flattened().Fill(constant);}

    template<class TCONSTANT,class TGRID>
    static void Put_Ghost(const TCONSTANT constant,ARRAYS_3D& x,const GRID_3D<TGRID>& grid,const int ghost_cells)
    {for(int j=1-ghost_cells;j<=grid.n+ghost_cells;j++) for(int ij=1-ghost_cells;ij<=grid.mn+ghost_cells;ij++) for(int s=1;s<=ghost_cells;s++) for(int k=1;k<=length;k++)
        x(k,1-s,j,ij)=x(k,grid.m+s,j,ij)=constant;
    for(int i=1;i<=grid.m;i++) for(int ij=1-ghost_cells;ij<=grid.mn+ghost_cells;ij++) for(int s=1;s<=ghost_cells;s++) for(int k=1;k<=x.length;k++)
        x(k,i,1-s,ij)=x(k,i,grid.n+s,ij)=constant;
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int s=1;s<=ghost_cells;s++) for(int k=1;k<=x.length;k++) x(k,i,j,1-s)=x(k,i,j,grid.mn+s)=constant;}

    template<class RW>
    void Read(std::istream& input)
    {ARRAYS_3D<VECTOR<T,length> >::template Read_With_Length<RW>(input,length);}

    template<class RW>
    void Write(std::ostream& output) const
    {ARRAYS_3D<VECTOR<T,length> >::template Write_With_Length<RW>(output,length);}

private:
    // old forms of functions which should not be called!
    ARRAYS_3D(const int length_new,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end,const bool initialize_using_default_constructor=true);
    void Resize(int length_new,int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const bool initialize_new_elements=true,const bool copy_existing_elements=true,
        const T& initialization_value=T());
//#####################################################################
};

template<class T> inline std::ostream& operator<<(std::ostream& output,const ARRAYS_3D<T>& a)
{for(int i=a.m_start;i<=a.m_end;i++){for(int j=a.n_start;j<=a.n_end;j++){for(int ij=a.mn_start;ij<=a.mn_end;ij++)output<<a(i,j,ij)<<" ";output<<std::endl;}
    output<<"------------------------------------------"<<std::endl;}
return output;}
}
#endif
