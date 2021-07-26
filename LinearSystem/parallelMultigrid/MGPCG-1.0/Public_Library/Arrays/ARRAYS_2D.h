//#####################################################################
// Copyright 2002-2007, Robert Bridson, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Duc Nguyen, Craig Schroeder, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAYS_2D
//#####################################################################
#ifndef __ARRAYS_2D__
#define __ARRAYS_2D__

#include <Arrays/ARRAYS_ND_BASE.h>
#include <Grids/GRID_2D.h>
#include <Matrices_And_Vectors/VECTOR.h>
#include <Utilities/EXCEPTIONS.h>
#include <boost/format.hpp>
namespace PhysBAM{

template<class T>
class ARRAYS_2D<T,1>:public ARRAYS_ND_BASE<T,ARRAYS_2D<T> >
{
public:
    enum WORKAROUND1 {dimension=2};
    enum WORKAROUND2 {length=1};
    template<class T2> struct REBIND{typedef ARRAYS_2D<T2,length> TYPE;};
    template<int length2> struct REBIND_LENGTH{typedef ARRAYS_2D<T,length2> TYPE;};
    typedef T ELEMENT;

    using ARRAYS_ND_BASE<T,ARRAYS_2D>::array; // one-dimensional data storage

    int m_start,m_end,n_start,n_end; // starting and ending values for x and y direction
    int m,n; // used in overloaded () operator
private:
    int n_plus_one,n_minus_one;
    T* base_pointer;
public:

    ARRAYS_2D()
        :m_start(1),m_end(0),n_start(1),n_end(0),m(0),n(0)
    {
        Calculate_Acceleration_Constants();
    }

    ARRAYS_2D(const BOX<VECTOR<int,2> >& domain,const bool initialize_using_default_constructor=true)
        :m_start(domain.min_corner.x),m_end(domain.max_corner.x),n_start(domain.min_corner.y),n_end(domain.max_corner.y),m(m_end-m_start+1),n(n_end-n_start+1)
    {
        assert(m>=0 && n>=0);int size=m*n;
        {RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    ARRAYS_2D(const int m_start_input,const int m_end_input,const int n_start_input,const int n_end_input,const bool initialize_using_default_constructor=true)
        :m_start(m_start_input),m_end(m_end_input),n_start(n_start_input),n_end(n_end_input),m(m_end-m_start+1),n(n_end-n_start+1)
    {
        assert(m>=0 && n>=0);int size=m*n;
        {RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    template<class T2>
    ARRAYS_2D(const GRID_2D<T2>& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :m_start(1-ghost_cells),m_end(grid.m+ghost_cells),n_start(1-ghost_cells),n_end(grid.n+ghost_cells),m(m_end-m_start+1),n(n_end-n_start+1)
    {
        {int size=m*n;RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    ARRAYS_2D(const ARRAYS_2D& old_array,const bool initialize_with_old_array=true)
        :ARRAYS_ND_BASE<T,ARRAYS_2D>(old_array.array.Size()),m_start(old_array.m_start),m_end(old_array.m_end),n_start(old_array.n_start),n_end(old_array.n_end),m(old_array.m),n(old_array.n)
    {
        Calculate_Acceleration_Constants();
        if(initialize_with_old_array) array=old_array.array;
    }

    BOX<VECTOR<int,2> > Domain_Indices() const
    {return BOX<VECTOR<int,2> >(m_start,m_end,n_start,n_end);}

    void Clean_Memory()
    {Resize(1,0,1,0,false,false);}

    void Calculate_Acceleration_Constants()
    {T* array_pointer=array.Get_Array_Pointer();
    n_plus_one=n+1;n_minus_one=n-1;
    base_pointer=array_pointer?array_pointer-(m_start*n+n_start):0;}

    ARRAYS_2D& operator=(const ARRAYS_2D& source)
    {if(array.Size()!=source.array.Size()){
        delete[] array.Get_Array_Pointer();
        RAW_ARRAY<T> new_array(source.array.Size(),new T[source.array.Size()]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);}
    else if(this==&source) return *this;
    m=source.m;m_start=source.m_start;m_end=source.m_end;n=source.n;n_start=source.n_start;n_end=source.n_end;
    Calculate_Acceleration_Constants();
    array=source.array;return *this;}

    T& operator()(const int i,const int j)
    {assert(m_start<=i && i<=m_end);assert(n_start<=j && j<=n_end);
    return base_pointer[i*n+j];}

    const T& operator()(const int i,const int j) const
    {assert(m_start<=i && i<=m_end);assert(n_start<=j && j<=n_end);
    return base_pointer[i*n+j];}

    T& operator()(const VECTOR<int,2>& index)
    {assert(m_start<=index.x && index.x<=m_end);assert(n_start<=index.y && index.y<=n_end);
    return base_pointer[index.x*n+index.y];}

    const T& operator()(const VECTOR<int,2>& index) const
    {assert(m_start<=index.x && index.x<=m_end);assert(n_start<=index.y && index.y<=n_end);
    return base_pointer[index.x*n+index.y];}

    T& operator()(const int k,const int i,const int j)
    {assert(k==1);
    return (*this)(i,j);}

    const T& operator()(const int k,const int i,const int j) const
    {assert(k==1);
    return (*this)(i,j);}

    T& operator()(const int k,const VECTOR<int,2>& index)
    {assert(k==1);
    return (*this)(index);}

    const T& operator()(const int k,const VECTOR<int,2>& index) const
    {assert(k==1);
    return (*this)(index);}

    bool Valid_Index(const VECTOR<int,2>& index) const
    {return m_start<=index.x && index.x<=m_end && n_start<=index.y && index.y<=n_end;}

    bool Valid_Index(const int i,const int j) const
    {return m_start<=i && i<=m_end && n_start<=j && j<=n_end;}

    int Standard_Index(const int i,const int j) const
    {assert(m_start<=i && i<=m_end);assert(n_start<=j && j<=n_end);
    return (i-m_start)*n+(j-n_start)+1;}

    int Standard_Index(const VECTOR<int,2>& index) const
    {assert(m_start<=index.x && index.x<=m_end);assert(n_start<=index.y && index.y<=n_end);
    return (index.x-m_start)*n+(index.y-n_start)+1;}

    int I_Plus_One(const int index) const
    {return index+n;}

    int I_Minus_One(const int index) const
    {return index-n;}

    int J_Plus_One(const int index) const
    {return index+1;}

    int J_Minus_One(const int index) const
    {return index-1;}

    int I_Plus_One_J_Plus_One(const int index) const
    {return index+n_plus_one;}

    int I_Plus_One_J_Minus_One(const int index) const
    {return index+n_minus_one;}

    int I_Minus_One_J_Plus_One(const int index) const
    {return index-n_minus_one;}

    int I_Minus_One_J_Minus_One(const int index) const
    {return index-n_plus_one;}

    void Resize(int m_start_new,int m_end_new,int n_start_new,int n_end_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true,
        const T& initialization_value=T())
    {if(m_start_new==m_start && m_end_new==m_end && n_start_new==n_start && n_end_new==n_end) return;
    int m_new=m_end_new-m_start_new+1,n_new=n_end_new-n_start_new+1;assert(m_new>=0 && n_new>=0);
    int size_new=m_new*n_new;RAW_ARRAY<T> array_new(size_new,new T[size_new]);
    if(initialize_new_elements) array_new.Fill(initialization_value);
    if(copy_existing_elements){
        int m1=PhysBAM::max(m_start,m_start_new),m2=PhysBAM::min(m_end,m_end_new),n1=PhysBAM::max(n_start,n_start_new),n2=PhysBAM::min(n_end,n_end_new);
        for(int i=m1;i<=m2;i++) for(int j=n1;j<=n2;j++)
            array_new(((i-m_start_new)*n_new+(j-n_start_new))+1)=array(((i-m_start)*n+(j-n_start))+1);}
    m_start=m_start_new;m_end=m_end_new;m=m_end-m_start+1;n_start=n_start_new;n_end=n_end_new;n=n_end-n_start+1;
    delete[] array.Get_Array_Pointer();RAW_ARRAY<T>::Exchange_Arrays(array,array_new);Calculate_Acceleration_Constants();}

    template<class T2>
    void Resize(const GRID_2D<T2>& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(1-ghost_cells,grid.m+ghost_cells,1-ghost_cells,grid.n+ghost_cells,initialize_new_elements,copy_existing_elements,initialization_value);}

    void Resize(const BOX<VECTOR<int,2> >& box,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,initialize_new_elements,copy_existing_elements,initialization_value);}

    // specialize this; can't just call Resize, it doesn't get the copy right!
    // TODO: merge this into Resize, as the case where copy_existing_elements is false
    // TODO: this function is extremely broken, since all the functions in the class depends on size being correct
    void Resize_In_Place(const int m_start_new,const int m_end_new,const int n_start_new,const int n_end_new)
    {if(array.Size()>=(m_end_new-m_start_new+1)*(n_end_new-n_start_new+1)){
        m_start=m_start_new;m_end=m_end_new;m=m_end-m_start+1;n_start=n_start_new;n_end=n_end_new;n=n_end-n_start+1;
        Calculate_Acceleration_Constants();}
    else Resize(m_start_new,m_end_new,n_start_new,n_end_new,false,false);}

    void Resize_In_Place(const BOX<VECTOR<int,2> >& box)
    {Resize_In_Place(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y);}

    void Clamp(int& i,int& j) const
    {i=clamp(i,m_start,m_end);j=clamp(j,n_start,n_end);}

    void Clamp_End_Minus_One(int& i,int& j) const
    {i=clamp(i,m_start,m_end-1);j=clamp(j,n_start,n_end-1);}

    void Clamp_End_Minus_Two(int& i,int& j) const
    {i=clamp(i,m_start,m_end-2);j=clamp(j,n_start,n_end-2);}

    void Clamp_End_Minus_Three(int& i,int& j) const
    {i=clamp(i,m_start,m_end-3);j=clamp(j,n_start,n_end-3);}

    void Clamp_Interior(int& i,int& j) const
    {i=clamp(i,m_start+1,m_end-1);j=clamp(j,n_start+1,n_end-1);}

    void Clamp_Interior_End_Minus_One(int& i,int& j) const
    {i=clamp(i,m_start+1,m_end-2);j=clamp(j,n_start+1,n_end-2);}

    static void Get(ARRAYS_2D& new_copy,const ARRAYS_2D& old_copy)
    {if(&old_copy!=&new_copy) ARRAYS_2D::Put(old_copy,new_copy,new_copy.m_start,new_copy.m_end,new_copy.n_start,new_copy.n_end);}

    static void Shifted_Get(ARRAYS_2D& new_copy,const ARRAYS_2D& old_copy,const VECTOR<int,2>& shift)
    {if(shift==VECTOR<int,2>()) Get(new_copy,old_copy);
    else for(int i=new_copy.m_start;i<=new_copy.m_end;i++) for(int j=new_copy.n_start;j<=new_copy.n_end;j++)new_copy(i,j)=old_copy(i+shift.x,j+shift.y);}

    static void Limited_Shifted_Get(ARRAYS_2D& new_copy,const ARRAYS_2D& old_copy,const VECTOR<int,2>& shift)
    {int m_start=PhysBAM::max(new_copy.m_start,old_copy.m_start-shift.x),m_end=PhysBAM::min(new_copy.m_end,old_copy.m_end-shift.x);
    int n_start=PhysBAM::max(new_copy.n_start,old_copy.n_start-shift.y),n_end=PhysBAM::min(new_copy.n_end,old_copy.n_end-shift.y);
    for(int i=m_start;i<=m_end;i++) for(int j=n_start;j<=n_end;j++)
        new_copy(i,j)=old_copy(i+shift.x,j+shift.y);}

    static void Put(const ARRAYS_2D& old_copy,ARRAYS_2D& new_copy)
    {if(&old_copy!=&new_copy) ARRAYS_2D::Put(old_copy,new_copy,old_copy.m_start,old_copy.m_end,old_copy.n_start,old_copy.n_end);}

    static void Shifted_Put(const ARRAYS_2D& old_copy,ARRAYS_2D& new_copy,const VECTOR<int,2>& shift)
    {if(shift==VECTOR<int,2>()) Put(old_copy,new_copy);
    else for(int i=old_copy.m_start;i<=old_copy.m_end;i++) for(int j=old_copy.n_start;j<=old_copy.n_end;j++)new_copy(i+shift.x,j+shift.y)=old_copy(i,j);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_2D& old_copy,ARRAYS_2D& new_copy)
    {ARRAYS_2D::Put(constant,old_copy,new_copy,old_copy.m_start,old_copy.m_end,old_copy.n_start,old_copy.n_end);}

private:
    static void Put(const ARRAYS_2D<T,1>& old_copy,ARRAYS_2D<T,1>& new_copy,const int m_start,const int m_end,const int n_start,const int n_end)
    {assert(old_copy.m_start<=m_start&&m_end<=old_copy.m_end&&old_copy.n_start<=n_start&&n_end<=old_copy.n_end);
    assert(new_copy.m_start<=m_start&&m_end<=new_copy.m_end&&new_copy.n_start<=n_start&&n_end<=new_copy.n_end);
    int old_index=old_copy.Standard_Index(m_start,n_start),new_index=new_copy.Standard_Index(m_start,n_start);
    const int inner_loop_copy_size=(n_end-n_start+1);
    for(int i=m_start;i<=m_end;i++){
        int old_save_index=old_index,new_save_index=new_index;
        for(int t=1;t<=inner_loop_copy_size;t++) new_copy.array(new_index++)=old_copy.array(old_index++);
        old_index=old_copy.I_Plus_One(old_save_index);new_index=new_copy.I_Plus_One(new_save_index);}}

    template<class T2>
    static void Put(T2 constant,const ARRAYS_2D& old_copy,ARRAYS_2D& new_copy,const int m_start,const int m_end,const int n_start,const int n_end)
    {assert(old_copy.m_start<=m_start&&m_end<=old_copy.m_end&&old_copy.n_start<=n_start&&n_end<=old_copy.n_end);
    assert(new_copy.m_start<=m_start&&m_end<=new_copy.m_end&&new_copy.n_start<=n_start&&n_end<=new_copy.n_end);
    int old_index=old_copy.Standard_Index(m_start,n_start),new_index=new_copy.Standard_Index(m_start,n_start);
    const int inner_loop_copy_size=(n_end-n_start+1);
    for(int i=m_start;i<=m_end;i++){
        int old_save_index=old_index,new_save_index=new_index;
        for(int t=1;t<=inner_loop_copy_size;t++) new_copy.array(new_index++)=constant*old_copy.array(old_index++);
        old_index=old_copy.I_Plus_One(old_save_index);new_index=new_copy.I_Plus_One(new_save_index);}}

public:
    template<class TCONSTANT,class TGRID>
    static void Put_Ghost(const TCONSTANT constant,ARRAYS_2D& x,const GRID_2D<TGRID>& grid,const int ghost_cells)
    {for(int j=1-ghost_cells;j<=grid.n+ghost_cells;j++) for(int s=1;s<=ghost_cells;s++) x(1-s,j)=x(grid.m+s,j)=constant;
    for(int i=1;i<=grid.m;i++) for(int s=1;s<=ghost_cells;s++) x(i,1-s)=x(i,grid.n+s)=constant;}

    static void Exchange_Arrays(ARRAYS_2D& a,ARRAYS_2D& b)
    {RAW_ARRAY<T>::Exchange_Arrays(a.array,b.array);
    exchange(a.m_start,b.m_start);exchange(a.m_end,b.m_end);exchange(a.m,b.m);
    exchange(a.n_start,b.n_start);exchange(a.n_end,b.n_end);exchange(a.n,b.n);
    a.Calculate_Acceleration_Constants();b.Calculate_Acceleration_Constants();}

    template<class T2>
    static bool Equal_Dimensions(const ARRAYS_2D& a,const ARRAYS_2D<T2>& b)
    {return a.m_start==b.m_start && a.m_end==b.m_end && a.n_start==b.n_start && a.n_end==b.n_end;}

    static bool Equal_Dimensions(const ARRAYS_2D& a,const int m_start,const int m_end,const int n_start,const int n_end)
    {return a.m_start==m_start && a.m_end==m_end && a.n_start==n_start && a.n_end==n_end;}

    void Move_Contents_Left(const int increment=1)
    {for(int i=m_start;i<=m_end-increment;i++) for(int j=n_start;j<=n_end;j++) (*this)(i,j)=(*this)(i+increment,j);}

    void Move_Contents_Right(const int increment=1)
    {for(int i=m_end;i>=m_start+increment;i--) for(int j=n_start;j<=n_end;j++) (*this)(i,j)=(*this)(i-increment,j);}

    void Move_Contents_Down(const int increment=1)
    {for(int i=m_start;i<=m_end;i++) for(int j=n_start;j<=n_end-increment;j++) (*this)(i,j)=(*this)(i,j+increment);}

    void Move_Contents_Up(const int increment=1)
    {for(int i=m_start;i<=m_end;i++) for(int j=n_end;j>=n_start+increment;j--) (*this)(i,j)=(*this)(i,j-increment);}

    template<class RW>
    void Read(std::istream& input)
    {Read_With_Length<RW>(input,1);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_With_Length<RW>(output,1);}

protected:
    template<class RW>
    void Read_With_Length(std::istream& input,const int length2)
    {int read_length;Read_Binary<RW>(input,read_length,m_start,m_end,n_start,n_end);
    if(read_length!=length2) throw READ_ERROR(str(boost::format("Read length %d not equal to %d")%read_length%length2));
    m=m_end-m_start+1;
    if(m<0) throw READ_ERROR(str(boost::format("Invalid negative array size m = %d")%m));
    n=n_end-n_start+1;
    if(n<0) throw READ_ERROR(str(boost::format("Invalid negative array size n = %d")%n));
    delete[] array.Get_Array_Pointer();int size=m*n;RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(array,new_array);
    Read_Binary_Array<RW>(input,array.Get_Array_Pointer(),array.Size());Calculate_Acceleration_Constants();}

    template<class RW>
    void Write_With_Length(std::ostream& output,const int length2) const
    {Write_Binary<RW>(output,length2,m_start,m_end,n_start,n_end);Write_Binary_Array<RW>(output,array.Get_Array_Pointer(),array.Size());}

//#####################################################################
};

template<class T,int length_input>
class ARRAYS_2D:public ARRAYS_2D<VECTOR<T,length_input> >
{
    STATIC_ASSERT(length_input>0);
public:
    enum WORKAROUND {length=length_input};
    template<class T2> struct REBIND{typedef ARRAYS_2D<T2,length> TYPE;};
    template<int length2> struct REBIND_LENGTH{typedef ARRAYS_2D<T,length2> TYPE;};

    typedef ARRAYS_2D<VECTOR<T,length> > BASE;
    typedef T ELEMENT_OF_ELEMENT;
    using BASE::m;using BASE::n;using BASE::m_start;using BASE::m_end;using BASE::n_start;using BASE::n_end;
    using BASE::operator();using BASE::Valid_Index;using BASE::Standard_Index;using BASE::array;

public:
    ARRAYS_2D()
    {}

    ARRAYS_2D(const BOX<VECTOR<int,2> >& domain,const bool initialize_using_default_constructor=true)
        :ARRAYS_2D<VECTOR<T,length> >(domain,initialize_using_default_constructor)
    {}

    ARRAYS_2D(const int m_start_input,const int m_end_input,const int n_start_input,const int n_end_input,const bool initialize_using_default_constructor=true)
        :ARRAYS_2D<VECTOR<T,length> >(m_start_input,m_end_input,n_start_input,n_end_input,initialize_using_default_constructor)
    {}

    template<class T2>
    ARRAYS_2D(const GRID_2D<T2>& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :ARRAYS_2D<VECTOR<T,length> >(grid,ghost_cells,initialize_using_default_constructor)
    {}

    ARRAYS_2D(const ARRAYS_2D& old_array,const bool initialize_with_old_array=true)
        :ARRAYS_2D<VECTOR<T,length> >(old_array,initialize_with_old_array)
    {}

    T& operator()(const int k,const int i,const int j)
    {return (*this)(i,j)[k];}

    const T& operator()(const int k,const int i,const int j) const
    {return (*this)(i,j)[k];}

    T& operator()(const int k,const VECTOR<int,2>& index)
    {return (*this)(index)[k];}

    const T& operator()(const int k,const VECTOR<int,2>& index) const
    {return (*this)(index)[k];}

    bool Valid_Index(const int k,const int i,const int j) const
    {return 1<=k && k<=length && Valid_Index(i,j);}

    static void Extract_Dimension(const ARRAYS_2D& old_array,ARRAYS_2D<T>& extracted_array,int dimension)
    {extracted_array.Resize(old_array.m_start,old_array.m_end,old_array.n_start,old_array.n_end,false,false);
    for(int i=old_array.m_start;i<=old_array.m_end;i++) for(int j=old_array.n_start;j<=old_array.n_end;j++) extracted_array(i,j)=old_array(dimension,i,j);}

    template<class T2>
    void Resize(const GRID_2D<T2>& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {ARRAYS_2D<VECTOR<T,length> >::Resize(grid,ghost_cells,initialize_new_elements,copy_existing_elements);}

    template<class T2>
    void Resize(const GRID_2D<T2>& grid,const int ghost_cells,const bool initialize_new_elements,const bool copy_existing_elements,const T& initialization_value)
    {VECTOR<T,length> initialization_vector;initialization_vector.Fill(initialization_value);
    ARRAYS_2D<VECTOR<T,length> >::Resize(grid,ghost_cells,initialize_new_elements,copy_existing_elements,initialization_vector);}

    void Fill(const T& constant)
    {array.Flattened().Fill(constant);}

    template<class TCONSTANT,class TGRID>
    static void Put_Ghost(const TCONSTANT constant,ARRAYS_2D& x,const GRID_2D<TGRID>& grid,const int ghost_cells)
    {for(int j=1-ghost_cells;j<=grid.n+ghost_cells;j++) for(int s=1;s<=ghost_cells;s++) for(int k=1;k<=length;k++) x(k,1-s,j)=x(k,grid.m+s,j)=constant;
    for(int i=1;i<=grid.m;i++) for(int s=1;s<=ghost_cells;s++) for(int k=1;k<=length;k++) x(k,i,1-s)=x(k,i,grid.n+s)=constant;}

    template<class RW>
    void Read(std::istream& input)
    {ARRAYS_2D<VECTOR<T,length> >::template Read_With_Length<RW>(input,length);}

    template<class RW>
    void Write(std::ostream& output) const
    {ARRAYS_2D<VECTOR<T,length> >::template Write_With_Length<RW>(output,length);}

private:
    // old forms of functions which should not be called!
    ARRAYS_2D(const int length_new,const int m_start,const int m_end,const int n_start,const int n_end,const bool initialize_using_default_constructor=true);
    void Resize(int length_new,int m_start,int m_end,int n_start,int n_end,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T());
//#####################################################################
};

template<class T> inline std::ostream& operator<<(std::ostream& output,const ARRAYS_2D<T>& a)
{for(int j=a.n_end;j>=a.n_start;j--){for(int i=a.m_start;i<=a.m_end;i++) output<<a(i,j)<<" ";output<<std::endl;}return output;}
}
#endif
