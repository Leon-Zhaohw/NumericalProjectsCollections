//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Koltakov, Frank Losasso, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAYS_1D
//#####################################################################
#ifndef __ARRAYS_1D__
#define __ARRAYS_1D__

#include <Arrays/ARRAYS_ND_BASE.h>
#include <Grids/GRID_1D.h>
#include <Matrices_And_Vectors/VECTOR.h>
#include <Utilities/EXCEPTIONS.h>
#include <boost/format.hpp>
namespace PhysBAM{

template<class T>
class ARRAYS_1D<T,1>:public ARRAYS_ND_BASE<T,ARRAYS_1D<T> >
{
public:
    enum WORKAROUND1 {dimension=1};
    enum WORKAROUND2 {length=1};
    template<class T2> struct REBIND{typedef ARRAYS_1D<T2> TYPE;};
    template<int length2> struct REBIND_LENGTH{typedef ARRAYS_1D<T,length2> TYPE;};
    typedef T ELEMENT;

    using ARRAYS_ND_BASE<T,ARRAYS_1D>::array; // one-dimensional data storage

    int m_start,m_end; // starting and ending values for x direction
    int m; // used in overloaded () operator
private:
    T* base_pointer;
public:

    ARRAYS_1D()
        :m_start(1),m_end(0),m(0)
    {
        Calculate_Acceleration_Constants();
    }

    ARRAYS_1D(const BOX<VECTOR<int,1> >& domain,const bool initialize_using_default_constructor=true)
        :m_start(domain.min_corner.x),m_end(domain.max_corner.x),m(m_end-m_start+1)
    {
        assert(m>=0);int size=m;
        {RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    ARRAYS_1D(const int m_start_input,const int m_end_input,const bool initialize_using_default_constructor=true)
        :m_start(m_start_input),m_end(m_end_input),m(m_end-m_start+1)
    {
        assert(m>=0);int size=m;
        {RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    template<class T2>
    ARRAYS_1D(const GRID_1D<T2>& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :m_start(1-ghost_cells),m_end(grid.m+ghost_cells),m(m_end-m_start+1)
    {
        {int size=m;RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);} // allocate a new array
        Calculate_Acceleration_Constants();
        if(initialize_using_default_constructor) array.Fill(T()); // initialize array using default constructor
    }

    ARRAYS_1D(const ARRAYS_1D& old_array,const bool initialize_with_old_array=true)
        :ARRAYS_ND_BASE<T,ARRAYS_1D>(old_array.array.Size()),m_start(old_array.m_start),m_end(old_array.m_end),m(old_array.m)
    {
        Calculate_Acceleration_Constants();
        if(initialize_with_old_array) array=old_array.array;
    }

    BOX<VECTOR<int,1> > Domain_Indices() const
    {return BOX<VECTOR<int,1> >(m_start,m_end);}

    void Clean_Memory()
    {Resize(1,0,false,false);}

    void Calculate_Acceleration_Constants()
    {T* array_pointer=array.Get_Array_Pointer();
    base_pointer=array_pointer?array_pointer-m_start:0;}

    ARRAYS_1D& operator=(const ARRAYS_1D& source)
    {if(array.Size()!=source.array.Size()){
        delete[] array.Get_Array_Pointer();
        RAW_ARRAY<T> new_array(source.array.Size(),new T[source.array.Size()]);RAW_ARRAY<T>::Exchange_Arrays(new_array,array);}
    else if(this==&source) return *this;
    m=source.m;m_start=source.m_start;m_end=source.m_end;
    Calculate_Acceleration_Constants();
    array=source.array;return *this;}

    T& operator()(const int i)
    {assert(m_start<=i && i<=m_end);
    return base_pointer[i];}

    const T& operator()(const int i) const
    {assert(m_start<=i && i<=m_end);
    return base_pointer[i];}

    T& operator()(const VECTOR<int,1>& index)
    {assert(m_start<=index.x && index.x<=m_end);
    return base_pointer[index.x];}

    const T& operator()(const VECTOR<int,1>& index) const
    {assert(m_start<=index.x && index.x<=m_end);
    return base_pointer[index.x];}

    T& operator()(const int k,const int i)
    {assert(k==1);
    return (*this)(i);}

    const T& operator()(const int k,const int i) const
    {assert(k==1);
    return (*this)(i);}

    T& operator()(const int k,const VECTOR<int,1>& index)
    {assert(k==1);
    return (*this)(index.x);}

    const T& operator()(const int k,const VECTOR<int,1>& index) const
    {assert(k==1);
    return (*this)(index.x);}

    bool Valid_Index(const VECTOR<int,1>& index) const
    {return m_start<=index.x && index.x<=m_end;}

    bool Valid_Index(const int i) const
    {return m_start<=i && i<=m_end;}

    int Standard_Index(const int i) const
    {assert(m_start<=i && i<=m_end);
    return i-m_start+1;}

    int Standard_Index(const VECTOR<int,1>& index) const
    {return Standard_Index(index.x);}

    int I_Plus_One(const int index) const
    {return index+1;}

    int I_Minus_One(const int index) const
    {return index-1;}

    void Resize(const int m_start_new,const int m_end_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {if(m_start_new==m_start && m_end_new==m_end) return;
    int m_new=m_end_new-m_start_new+1;assert(m_new>=0);
    int size_new=m_new;RAW_ARRAY<T> array_new(size_new,new T[size_new]);
    if(initialize_new_elements) array_new.Fill(initialization_value);
    if(copy_existing_elements){
        int m1=PhysBAM::max(m_start,m_start_new),m2=PhysBAM::min(m_end,m_end_new);
        int start=(m1-m_start_new)+1,old_start=(m1-m_start)+1,difference=old_start-start,end=(m2-m_start_new)+1;
        for(int k=start;k<=end;k++) array_new(k)=array(k+difference);}
    m_start=m_start_new;m_end=m_end_new;m=m_new;
    delete[] array.Get_Array_Pointer();RAW_ARRAY<T>::Exchange_Arrays(array,array_new);Calculate_Acceleration_Constants();}

    template<class T2>
    void Resize(const GRID_1D<T2>& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(1-ghost_cells,grid.m+ghost_cells,initialize_new_elements,copy_existing_elements,initialization_value);}

    void Resize(const BOX<VECTOR<int,1> >& box,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(box.min_corner.x,box.max_corner.x,initialize_new_elements,copy_existing_elements,initialization_value);}

    void Resize_In_Place(const int m_start_new,const int m_end_new) // TODO: this function is extremely broken, since all the functions in the class depends on size being correct
    {if(array.Size()>=(m_end_new-m_start_new+1)){
        m_start=m_start_new;m_end=m_end_new;m=m_end-m_start+1;Calculate_Acceleration_Constants();}
    else Resize(m_start_new,m_end_new,false,false);}

    void Resize_In_Place(const BOX<VECTOR<int,1> >& box)
    {Resize_In_Place(box.min_corner.x,box.max_corner.x);}

    void Clamp(int& i) const
    {i=clamp(i,m_start,m_end);}

    void Clamp_End_Minus_One(int& i) const
    {i=clamp(i,m_start,m_end-1);}

    void Clamp_End_Minus_Two(int& i) const
    {i=clamp(i,m_start,m_end-2);}

    void Clamp_End_Minus_Three(int& i) const
    {i=clamp(i,m_start,m_end-3);}

    void Clamp_Interior(int& i) const
    {i=clamp(i,m_start+1,m_end-1);}

    void Clamp_Interior_End_Minus_One(int& i) const
    {i=clamp(i,m_start+1,m_end-2);}

    static void Get(ARRAYS_1D& new_copy,const ARRAYS_1D& old_copy)
    {if(&old_copy!=&new_copy) ARRAYS_1D::Put(old_copy,new_copy,new_copy.m_start,new_copy.m_end);}

    static void Shifted_Get(ARRAYS_1D& new_copy,const ARRAYS_1D& old_copy,const VECTOR<int,1>& shift)
    {if(shift==VECTOR<int,1>()) Get(new_copy,old_copy);
    else for(int i=new_copy.m_start;i<=new_copy.m_end;i++) new_copy(i)=old_copy(i+shift.x);}

    static void Limited_Shifted_Get(ARRAYS_1D& new_copy,const ARRAYS_1D& old_copy,const VECTOR<int,1>& shift)
    {int m_start=PhysBAM::max(new_copy.m_start,old_copy.m_start-shift.x),m_end=PhysBAM::min(new_copy.m_end,old_copy.m_end-shift.x);
    for(int i=m_start;i<=m_end;i++)new_copy(i)=old_copy(i+shift.x);}

    static void Put(const ARRAYS_1D& old_copy,ARRAYS_1D& new_copy)
    {if(&old_copy!=&new_copy) ARRAYS_1D::Put(old_copy,new_copy,old_copy.m_start,old_copy.m_end);}

    static void Shifted_Put(const ARRAYS_1D& old_copy,ARRAYS_1D& new_copy,const VECTOR<int,1>& shift)
    {STATIC_ASSERT(length==1);
    if(shift==VECTOR<int,1>()) Put(old_copy,new_copy);
    else for(int i=old_copy.m_start;i<=old_copy.m_end;i++)new_copy(i+shift.x)=old_copy(i);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_1D& old_copy,ARRAYS_1D& new_copy)
    {ARRAYS_1D::Put(constant,old_copy,new_copy,old_copy.m_start,old_copy.m_end);}

    static void Put(const ARRAYS_1D& old_copy,ARRAYS_1D& new_copy,const int m_start,const int m_end)
    {assert(old_copy.m_start<=m_start&&m_end<=old_copy.m_end);assert(new_copy.m_start<=m_start&&m_end<=new_copy.m_end);
    int old_index=old_copy.Standard_Index(m_start),new_index=new_copy.Standard_Index(m_start);
    const int inner_loop_copy_size=(m_end-m_start+1);
    for(int t=1;t<=inner_loop_copy_size;t++) new_copy.array(new_index++)=old_copy.array(old_index++);}

    template<class T2>
    static void Put(T2 constant,const ARRAYS_1D& old_copy,ARRAYS_1D& new_copy,const int m_start,const int m_end)
    {assert(old_copy.m_start<=m_start&&m_end<=old_copy.m_end);assert(new_copy.m_start<=m_start&&m_end<=new_copy.m_end);
    int old_index=old_copy.Standard_Index(1,m_start),new_index=new_copy.Standard_Index(1,m_start);
    const int inner_loop_copy_size=(m_end-m_start+1);
    for(int t=1;t<=inner_loop_copy_size;t++) new_copy.array(new_index++)=constant*old_copy.array(old_index++);}

    template<class TCONSTANT,class TGRID>
    static void Put_Ghost(const TCONSTANT constant,ARRAYS_1D& x,const GRID_1D<TGRID>& grid,const int ghost_cells)
    {for(int s=1;s<=ghost_cells;s++) x(1-s)=x(grid.m+s)=constant;}

    template<class T2>
    void Set_Range(const T2& constant,const int i_start,const int i_end)
    {int start_index=Standard_Index(i_start),end_index=Standard_Index(i_end);
    for(int i=start_index;i<=end_index;i++) array(i)=constant;}

    template<class T2>
    void Set_Range(const T2& constant,const BOX<VECTOR<int,1> >& range)
    {Set_Range(constant,range.min_corner.x,range.max_corner.x);}

    static void Exchange_Arrays(ARRAYS_1D& a,ARRAYS_1D& b)
    {RAW_ARRAY<T>::Exchange_Arrays(a.array,b.array);
    exchange(a.m_start,b.m_start);exchange(a.m_end,b.m_end);exchange(a.m,b.m);
    a.Calculate_Acceleration_Constants();b.Calculate_Acceleration_Constants();}

    template<class T2>
    static bool Equal_Dimensions(const ARRAYS_1D& a,const ARRAYS_1D<T2,length>& b)
    {return a.m_start==b.m_start && a.m_end==b.m_end;}

    static bool Equal_Dimensions(const ARRAYS_1D& a,const int m_start,const int m_end)
    {return a.m_start==m_start && a.m_end==m_end;}

    void Move_Contents_Left(const int increment=1)
    {int jump=increment;for(int i=1;i<=array.Size()-jump;i++) array(i)=array(i+jump);}

    void Move_Contents_Right(const int increment=1)
    {int jump=increment;for(int i=array.Size();i>jump;i--) array(i)=array(i-jump);}

    template<class RW>
    void Read(std::istream& input)
    {Read_With_Length<RW>(input,1);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_With_Length<RW>(output,1);}

protected:
    template<class RW>
    void Read_With_Length(std::istream& input,const int length2)
    {int read_length;Read_Binary<RW>(input,read_length,m_start,m_end);
    if(read_length!=length2) throw READ_ERROR(str(boost::format("Read length %d not equal to %d")%read_length%length2));
    m=m_end-m_start+1;
    if(m<0) throw READ_ERROR(str(boost::format("Invalid negative array size m = %d")%m));
    delete[] array.Get_Array_Pointer();int size=m;RAW_ARRAY<T> new_array(size,new T[size]);RAW_ARRAY<T>::Exchange_Arrays(array,new_array);
    Read_Binary_Array<RW>(input,array.Get_Array_Pointer(),array.Size());Calculate_Acceleration_Constants();}

    template<class RW>
    void Write_With_Length(std::ostream& output,const int length2) const
    {Write_Binary<RW>(output,length2,m_start,m_end);Write_Binary_Array<RW>(output,array.Get_Array_Pointer(),array.Size());}

//#####################################################################
};

template<class T,int length_input>
class ARRAYS_1D:public ARRAYS_1D<VECTOR<T,length_input> >
{
    STATIC_ASSERT(length_input>0);
public:
    enum WORKAROUND {length=length_input};
    template<class T2> struct REBIND{typedef ARRAYS_1D<T2,length> TYPE;};
    template<int length2> struct REBIND_LENGTH{typedef ARRAYS_1D<T,length2> TYPE;};

    typedef ARRAYS_1D<VECTOR<T,length> > BASE;
    typedef T ELEMENT_OF_ELEMENT;
    using BASE::m;using BASE::m_start;using BASE::m_end;
    using BASE::operator();using BASE::Valid_Index;using BASE::Standard_Index;using BASE::array;

public:
    ARRAYS_1D()
    {}

    ARRAYS_1D(const BOX<VECTOR<int,1> >& domain,const bool initialize_using_default_constructor=true)
        :ARRAYS_1D<VECTOR<T,length> >(domain,initialize_using_default_constructor)
    {}

    ARRAYS_1D(const int m_start_input,const int m_end_input,const bool initialize_using_default_constructor=true)
        :ARRAYS_1D<VECTOR<T,length> >(m_start_input,m_end_input,initialize_using_default_constructor)
    {}

    template<class T2>
    ARRAYS_1D(const GRID_1D<T2>& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :ARRAYS_1D<VECTOR<T,length> >(grid,ghost_cells,initialize_using_default_constructor)
    {}

    ARRAYS_1D(const ARRAYS_1D& old_array,const bool initialize_with_old_array=true)
        :ARRAYS_1D<VECTOR<T,length> >(old_array,initialize_with_old_array)
    {}

    T& operator()(const int k,const int i)
    {return (*this)(i)[k];}

    const T& operator()(const int k,const int i) const
    {return (*this)(i)[k];}

    T& operator()(const int k,const VECTOR<int,1>& index)
    {return (*this)(index)[k];}

    const T& operator()(const int k,const VECTOR<int,1>& index) const
    {return (*this)(index)[k];}

    bool Valid_Index(const int k,const int i) const
    {return 1<=k && k<=length && Valid_Index(i);}

    int Standard_Index(const int k,const int i) const
    {assert(1<=k && k<=length);
    return Standard_Index(i);}

    static void Extract_Dimension(const ARRAYS_1D& old_array,ARRAYS_1D<T>& extracted_array,int dimension)
    {extracted_array.Resize(old_array.m_start,old_array.m_end,false,false);
    for(int i=old_array.m_start;i<=old_array.m_end;i++) extracted_array(i)=old_array(dimension,i);}

    template<class T2>
    void Resize(const GRID_1D<T2>& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {ARRAYS_1D<VECTOR<T,length> >::Resize(grid,ghost_cells,initialize_new_elements,copy_existing_elements);}

    template<class T2>
    void Resize(const GRID_1D<T2>& grid,const int ghost_cells,const bool initialize_new_elements,const bool copy_existing_elements,const T& initialization_value)
    {VECTOR<T,length> initialization_vector;initialization_vector.Fill(initialization_value);
    ARRAYS_1D<VECTOR<T,length> >::Resize(grid,ghost_cells,initialize_new_elements,copy_existing_elements,initialization_vector);}

    void Fill(const T& constant)
    {array.Flattened().Fill(constant);}

    template<class T2>
    void Set_Range(const T2& constant,const int i_start,const int i_end)
    {VECTOR<T,length> constant_vector;constant_vector.Fill(constant);
    ARRAYS_1D<VECTOR<T,length> >::template Set_Range<VECTOR<T,length> >(constant_vector,i_start,i_end);}

    template<class T2>
    void Set_Range(const T2& constant,const BOX<VECTOR<int,1> >& range)
    {Set_Range(constant,range.min_corner.x,range.max_corner.x);}

    template<class TCONSTANT,class TGRID>
    static void Put_Ghost(const TCONSTANT constant,ARRAYS_1D& x,const GRID_1D<TGRID>& grid,const int ghost_cells)
    {for(int s=1;s<=ghost_cells;s++) for(int k=1;k<=length;k++) x(k,1-s)=x(k,grid.m+s)=constant;}

    template<class RW>
    void Read(std::istream& input)
    {ARRAYS_1D<VECTOR<T,length> >::template Read_With_Length<RW>(input,length);}

    template<class RW>
    void Write(std::ostream& output) const
    {ARRAYS_1D<VECTOR<T,length> >::template Write_With_Length<RW>(output,length);}

private:
    // old forms of functions which should not be called!
    ARRAYS_1D(const int length_new,const int m_start,const int m_end,const bool initialize_using_default_constructor=true);
    void Resize(int length_new,int m_start,int m_end,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T());
//#####################################################################
};

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAYS_1D<T>& a)
{for(int i=a.m_start;i<=a.m_end;i++) output_stream<<a(i)<<" ";output_stream<<std::endl;return output_stream;}
}
#endif
