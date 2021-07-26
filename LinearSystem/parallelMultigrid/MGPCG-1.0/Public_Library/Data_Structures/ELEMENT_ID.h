//#####################################################################
// Copyright 2008, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELEMENT_ID
//#####################################################################
#ifndef __ELEMENT_ID__
#define __ELEMENT_ID__

#include <Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Read_Write/READ_WRITE_FUNCTIONS.h>
#include <boost/type_traits/is_fundamental.hpp>
namespace PhysBAM{

namespace ELEMENT_ID_HELPER{
enum {none=0,equality=1,compare=2,increment=4,add_T=8,to_bool=16,negate=32,for_loop=compare|increment,logical=equality|to_bool};
}

template<class ID,class T,int flags>
class ELEMENT_ID
{
    T id_value;
public:
    typedef T VALUE;

    ELEMENT_ID()
        :id_value(0)
    {}

    explicit ELEMENT_ID(T n)
        :id_value(n)
    {}

    T Value() const
    {return id_value;}

    bool operator==(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::equality);return id_value==id.id_value;}

    bool operator!=(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::equality);return id_value!=id.id_value;}

    bool operator<(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value<id.id_value;}

    bool operator>(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value>id.id_value;}

    bool operator<=(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value<=id.id_value;}

    bool operator>=(const ID id) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::compare);return id_value>=id.id_value;}

    ID operator++()
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);id_value++;return static_cast<ID&>(*this);}

    ID operator++(int)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);return ID(id_value++);}

    ID operator--()
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);id_value--;return static_cast<ID&>(*this);}

    ID operator--(int)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::increment);return ID(id_value--);}

    ID operator+(T i) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);return ID(id_value+i);}

    ID& operator+=(T i)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);id_value+=i;return (ID&)*this;}

    ID operator-(T i) const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);return ID(id_value-i);}

    ID& operator-=(T i)
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::add_T);id_value-=i;return (ID&)*this;}

private:
    struct UNUSABLE{void F(){}};
    typedef void (UNUSABLE::*SAFE_BOOL)();
public:

    operator SAFE_BOOL() const // allow conversion to bool without allowing conversion to T
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::to_bool);return id_value?&UNUSABLE::F:0;}

    ID operator-() const
    {STATIC_ASSERT(flags&ELEMENT_ID_HELPER::negate);return ID(-id_value);}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,id_value);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,id_value);}

//#####################################################################
};

template<class ID,class T,int flags> inline std::ostream& operator<<(std::ostream& output,const ELEMENT_ID<ID,T,flags> id)
{return output<<id.Value();}

template<class ID,class T,int flags> inline std::istream& operator>>(std::istream& input,ELEMENT_ID<ID,T,flags>& id)
{T i;input>>i;id=ID(i);return input;}

template<class T> inline typename boost::enable_if<boost::is_fundamental<T>,T>::type
Value(T i)
{return i;}

template<class ID> inline typename ID::VALUE
Value(ID i)
{return i.Value();}

//#####################################################################

#define PHYSBAM_DECLARE_ELEMENT_ID(ID,T,flags)       \
    struct ID:public ELEMENT_ID<ID,T,flags>          \
    {                                                \
        ID(){}                                       \
        explicit ID(T n):ELEMENT_ID<ID,T,flags>(n){} \
    };

PHYSBAM_DECLARE_ELEMENT_ID(INITIAL_SIZE,int,ELEMENT_ID_HELPER::equality);
#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(RIGID_BODY_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T|ELEMENT_ID_HELPER::negate);
PHYSBAM_DECLARE_ELEMENT_ID(JOINT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(COLLISION_BODY_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(FRAGMENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(SUPER_FRAGMENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(PARTITION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(HAIR_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(RIGID_BODY_INDEX,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T|ELEMENT_ID_HELPER::negate);
PHYSBAM_DECLARE_ELEMENT_ID(RIGID_CLUSTER_CONSTITUENT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int RIGID_BODY_ID;
typedef int JOINT_ID;
typedef int COLLISION_BODY_ID;
typedef int FRAGMENT_ID;
typedef int SUPER_FRAGMENT_ID;
typedef int PARTITION_ID;
typedef int HAIR_ID;
typedef int RIGID_BODY_INDEX;
typedef int RIGID_CLUSTER_CONSTITUENT_ID;
#endif

//#####################################################################
}
#endif
