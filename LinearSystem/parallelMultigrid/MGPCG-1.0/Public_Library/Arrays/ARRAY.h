//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY
//#####################################################################
#ifndef __ARRAY__
#define __ARRAY__

#include <Arrays/ARRAY_BASE.h>
#include <Math_Tools/exchange.h>
#include <Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Utilities/DEBUG_UTILITIES.h>
#include <Utilities/EXCEPTIONS.h>
#include <boost/format.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<ARRAY<T,ID> >:public mpl::true_{};
template<class T,class ID> struct IS_ARRAY_VIEW<ARRAY<T,ID> >:public mpl::false_{};

template<class T,class ID>
class ARRAY:public RAW_ARRAY<T,ID>
{
public:
    template<class T2> struct REBIND{typedef ARRAY<T2,ID> TYPE;};
    template<int length> struct REBIND_LENGTH:public PhysBAM::REBIND_LENGTH<ARRAY,length>{};
    typedef T ELEMENT;typedef ID INDEX;
private:
    struct UNUSABLE{};
    typedef RAW_ARRAY<T,ID> BASE;
public:
    using BASE::Get_Array_Pointer;using BASE::Valid_Index;using BASE::Size;using BASE::operator();
    using BASE::m; // m is protected in RAW_ARRAY, so make it public

    ARRAY()
        :BASE(ID(),0)
    {}

    explicit ARRAY(const ID m,const bool initialize_using_default_constructor=true)
        :BASE(ID(),0)
    {
        assert(m>=ID());RAW_ARRAY<T,ID> new_memory(m,new T[Value(m)]);RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);
        if(!IS_CLASS<T>::value && initialize_using_default_constructor) Fill(T());
    }

    ARRAY(const ARRAY& source)
        :BASE(ID(),0)
    {
        ID m_new=source.Size();
        RAW_ARRAY<T,ID> new_memory(m_new,new T[Value(m_new)]);
        RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);
        T* array=Get_Array_Pointer();const T* source_array=source.Get_Array_Pointer();
        for(int index=0;index<Value(m_new);index++) array[index]=source_array[index];
    }

    template<class T_ARRAY>
    //explicit ARRAY(const T_ARRAY& source,typename boost::enable_if<IS_SAME<T,typename T_ARRAY::ELEMENT>,UNUSABLE>::type unused=UNUSABLE())
    explicit ARRAY(const T_ARRAY& source,typename boost::enable_if<boost::is_same<T,typename T_ARRAY::ELEMENT> >::type* =0)
        :BASE(ID(),0)
    {
        ID m_new=source.Size();
        RAW_ARRAY<T,ID> new_memory(m_new,new T[Value(m_new)]);
        RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);
        for(ID i(1);i<=m_new;i++) (*this)(i)=source(i);
    }

    ~ARRAY()
    {delete[] Get_Array_Pointer();}

    void Clean_Memory()
    {RAW_ARRAY<T,ID> empty(ID(),0);RAW_ARRAY<T,ID>::Exchange_Arrays(*this,empty);delete[] empty.Get_Array_Pointer();}

    void Remove_All()
    {Resize(0,false,false);}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(ID i(1);i<=Size();i++) delete (*this)(i);Clean_Memory();}

    ARRAY& operator=(const ARRAY& source)
    {ID source_m=source.Size();
    if(!Equal_Dimensions(*this,source)){
        RAW_ARRAY<T,ID> new_memory(source_m,new T[Value(source_m)]);
        RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);
        delete[] new_memory.Get_Array_Pointer();}
    else if(Same_Array(*this,source)) return *this;
    T* array=Get_Array_Pointer();const T* source_array=source.Get_Array_Pointer();
    for(ID index=0;index<Value(source_m);index++) array[index]=source_array[index];return *this;}

    template<class T_ARRAY>
    ARRAY& operator=(const T_ARRAY& source)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    ID source_m=source.Size();
    if(!Equal_Dimensions(*this,source)){
        RAW_ARRAY<T,ID> new_memory(source_m,new T[Value(source_m)]);
        RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);
        delete[] new_memory.Get_Array_Pointer();}
    else if(Same_Array(*this,source)) return *this;
    for(ID i(1);i<=source_m;i++) (*this)(i)=source(i);return *this;}

    void Resize(const ID m_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {ID m_old=Size();if(m_old==m_new) return;
    RAW_ARRAY<T,ID> new_memory(m_new,new T[Value(m_new)]);
    if(copy_existing_elements){ID m_end=PhysBAM::min(m_old,m_new);for(ID index(1);index<=m_end;index++) new_memory(index)=(*this)(index);}
    if(!IS_CLASS<T>::value && initialize_new_elements){ID m_end=PhysBAM::min(m_old,m_new);for(ID index=m_end+1;index<=m_new;index++) new_memory(index)=T();}
    RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);delete[] new_memory.Get_Array_Pointer();}

    // equivalent to Resize, but assumes m_new >= Size()
    void Resize_Grow(const ID m_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {ID m_old=Size();if(m_old==m_new) return;assert(m_old<m_new);
    RAW_ARRAY<T,ID> new_memory(m_new,new T[Value(m_new)]);
    if(copy_existing_elements) for(ID index(1);index<=m_old;index++) new_memory(index)=(*this)(index);
    if(!IS_CLASS<T>::value && initialize_new_elements) for(ID index=m_old+1;index<=m_new;index++) new_memory(index)=T();
    RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);delete[] new_memory.Get_Array_Pointer();}

    void Resize(const ID m_new,const bool initialize_new_elements,const bool copy_existing_elements,const T& initialization_value)
    {ID m_old=Size();if(m_old==m_new) return;
    RAW_ARRAY<T,ID> new_memory(m_new,new T[Value(m_new)]);
    if(copy_existing_elements){ID m_end=PhysBAM::min(m_old,m_new);for(ID index(1);index<=m_end;index++) new_memory(index)=(*this)(index);}
    if(initialize_new_elements){ID m_end=PhysBAM::min(m_old,m_new);for(ID index=m_end+1;index<=m_new;index++) new_memory(index)=initialization_value;}
    RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);delete[] new_memory.Get_Array_Pointer();}

    ID Append(const T& element)
    {assert(m==0 || &element<&(*this)(1) || &element>&(*this)(m)); // make sure element isn't a reference into the current array
    Resize_Grow(m+1);(*this)(m)=element;return m;}

    template<class T_ARRAY>
    void Append_Elements(const T_ARRAY& append_array)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    ID index=m,append_m=append_array.Size();Resize_Grow(m+append_m,false,true);
    for(ID i(1);i<=append_m;i++){index++;(*this)(index)=append_array(i);}}

    void Append_Unique(const T& element)
    {ID index;if(Find(element,index)) return;
    Resize_Grow(m+1);(*this)(m)=element;}

    template<class T_ARRAY>
    void Append_Unique_Elements(const T_ARRAY& append_array)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    ARRAY<int,typename T_ARRAY::INDEX> append(append_array.Size());
    int count=0;
    for(typename T_ARRAY::INDEX j(1);j<=append_array.Size();j++){
        typename T_ARRAY::INDEX index;
        // element occurs duplicated in given array (only copy last instance)
        if(!Find(append_array(j),index) && !append_array.Find(append_array(j),j+1,index)){
            count++;append(j)=1;}}
    ID index=m+1;
    Resize_Grow(m+count);
    for(ID k(1);k<=append.m;k++) if(append(k)) (*this)(index++)=append_array(k);}

    void Remove_Index(const ID index) // maintains ordering of remaining elements
    {assert(Valid_Index(index));
    RAW_ARRAY<T,ID> new_memory(m-1,new T[Value(m)-1]);
    ID i,cut=index;for(i=ID(1);i<cut;i++) new_memory(i)=(*this)(i);
    for(ID j=i+1;i<m;i++,j++) new_memory(i)=(*this)(j);
    RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);delete[] new_memory.Get_Array_Pointer();}

    void Insert(const T& element,const ID index)
    {Resize_Grow(m+1);for(ID i=m;i>index;i--) (*this)(i)=(*this)(i-1);
    (*this)(index)=element;}

    template<class T_ARRAY>
    static ID Argcomp(const T_ARRAY& a,int(*compare)(const void* e1,const void* e2))
    {return argcomp(a,a.Size(),compare);}

    template<class T_ARRAY>
    static ID Argcomp(const T_ARRAY& a,const ID up_to_index,int(*compare)(const void* e1,const void* e2))
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    ID result(1);for(ID index(2);index<=up_to_index;index++) if(compare(&a(index),&a(result))) result=index;return result;}

    template<class T_ARRAY,class T_ARRAY_INT>
    static ID Argcomp(const T_ARRAY& lookup_array,const T_ARRAY_INT& indirection_array,int(*compare)(const void* e1,const void* e2))
    {return argcomp(lookup_array,indirection_array,indirection_array.m,compare);}

    template<class T_ARRAY,class T_ARRAY_INT>
    static ID Argcomp(const T_ARRAY& lookup_array,const T_ARRAY_INT& indirection_array,const ID up_to_index,int(*compare)(const void* e1,const void* e2))
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    STATIC_ASSERT_SAME(int,typename T_ARRAY_INT::ELEMENT);
    ID result(1);for(ID index(2);index<=up_to_index;index++) if(compare(&lookup_array(indirection_array(index)),&lookup_array(indirection_array(result)))) result=index;return result;}

    template<class T_ARRAY>
    static void Heapify(T_ARRAY& a) // largest on top
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    for(ID i=a.Size()/2;i>=ID(1);i--) Heapify(a,i,a.Size());}

    static void Heapify(ARRAY<T,ID>& a,const ID max_index) // largest on top, only does from 1 to max_index
    {for(ID i=max_index/2;i>=ID(1);i--) Heapify(a,i,max_index);}

    template<class T2>
    static void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux) // largest on top
    {for(ID i=a.m/2;i>=ID(1);i--) Heapify(a,aux,i,a.m);}

    template<class T2>
    static void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux,const ID max_index) // largest on top, only does from 1 to max_index
    {for(ID i(Value(max_index/2));i>=ID(1);i--) Heapify(a,aux,i,max_index);}

    template<class T_ARRAY>
    static void Heapify(T_ARRAY& a,ID index,const ID heap_size) // largest on top, only sorts down from index (not up!)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    for(;;){ID left(2*Value(index)),right(2*Value(index)+1),index_of_largest=index;
        if(left<=heap_size && a(left)>a(index_of_largest)) index_of_largest=left;
        if(right<=heap_size && a(right)>a(index_of_largest)) index_of_largest=right;
        if(index_of_largest!=index){exchange(a(index),a(index_of_largest));index=index_of_largest;}else return;}}

    template<class T2>
    static void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux,ID index,const ID heap_size) // largest on top, only sorts down from index (not up!)
    {for(;;){ID left(2*Value(index)),right(2*Value(index)+1),index_of_largest=index;
        if(left<=heap_size && a(left)>a(index_of_largest)) index_of_largest=left;
        if(right<=heap_size && a(right)>a(index_of_largest)) index_of_largest=right;
        if(index_of_largest!=index){exchange(a(index),a(index_of_largest));exchange(aux(index),aux(index_of_largest));index=index_of_largest;}else return;}}

    template<class T2>
    static void Compact_Array_Using_Compaction_Array(ARRAY<T2,ID>& array,const ARRAY<ID,ID>& compaction_array,ARRAY<T2,ID>* temporary_array=0)
    {ID compaction_array_m=compaction_array.Size();
    bool temporary_array_defined=temporary_array!=0;if(!temporary_array_defined) temporary_array=new ARRAY<T2,ID>(compaction_array_m,false);
    ARRAY<T2,ID>::Put(array,*temporary_array);for(ID i(1);i<=compaction_array_m;i++) if(compaction_array(i)>0) array(compaction_array(i))=(*temporary_array)(i);
    if(!temporary_array_defined){delete temporary_array;temporary_array=0;}}

    template<class RW>
    void Read(std::istream& input)
    {Clean_Memory();ID m;Read_Binary<RW>(input,m);
    if(m<ID()) throw READ_ERROR(str(boost::format("Invalid negative array size %d")%m));
    if(!m) return;
    RAW_ARRAY<T,ID> new_memory(m,new T[Value(m)]);
    Read_Binary_Array<RW>(input,new_memory.Get_Array_Pointer(),Value(m));
    RAW_ARRAY<T,ID>::Exchange_Arrays(*this,new_memory);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Prefix<RW>(output,Value(m));}

    template<class RW>
    void Write_Prefix(std::ostream& output,const ID prefix) const
    {PHYSBAM_ASSERT(ID()<=prefix && prefix<=Size());
    Write_Binary<RW>(output,Value(prefix));Write_Binary_Array<RW>(output,Get_Array_Pointer(),Value(prefix));}

private:
    // old forms of functions which should not be called!
    ARRAY(const ID m_start,const ID m_end){}
    void Resize(const ID m_start,const ID m_end){}
    void Resize(const ID m_start,const ID m_end,const bool initialize_new_elements){}
    void Resize(const ID m_start,const ID m_end,const bool initialize_new_elements,const bool copy_existing_elements){}
//#####################################################################
};
}
#endif
