//#####################################################################
// Copyright 2004-2006, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FUNCTIONS
//#####################################################################
// Functions for reading and writing which do the correct thing for objects, pointers, primitive types, etc. In general, use Read/Write_Binary (and Read/Write_Binary_Array) using T for the type
// of the object you're reading/writing and RW the underlying floating point scalar type (float/double).
//#####################################################################
#ifndef __READ_WRITE_FUNCTIONS__
#define __READ_WRITE_FUNCTIONS__

#include <Read_Write/READ_WRITE_FORWARD.h>
#include <Read_Write/TYPED_STREAM.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <boost/utility/enable_if.hpp>
namespace PhysBAM{

//#####################################################################
// Read_Binary
//#####################################################################
template<class RW,class T> inline typename boost::disable_if<HAS_TYPED_READ<T> >::type
Read_Binary(std::istream& input,T& d)
{Read_Write<T,RW>::Read(input,d);}

template<class RW,class T> inline typename boost::enable_if<HAS_TYPED_READ<T> >::type
Read_Binary(std::istream& input,T& d)
{TYPED_ISTREAM typed_input(input,STREAM_TYPE(RW()));d.Read(typed_input);}

template<class T> inline typename boost::disable_if<HAS_TYPED_READ<T> >::type
Read_Binary(TYPED_ISTREAM& input,T& d)
{input.Read(d);}

template<class T> inline typename boost::enable_if<HAS_TYPED_READ<T> >::type
Read_Binary(TYPED_ISTREAM& input,T& d)
{d.Read(input);}

//#####################################################################
// Write_Binary
//#####################################################################
template<class RW,class T> inline typename boost::disable_if<HAS_TYPED_WRITE<T> >::type
Write_Binary(std::ostream& output,const T& d)
{
    Read_Write<T,RW>::Write(output,d);
}

template<class RW,class T> inline typename boost::enable_if<HAS_TYPED_WRITE<T> >::type
Write_Binary(std::ostream& output,const T& d)
{
    TYPED_OSTREAM typed_output(output,STREAM_TYPE(RW()));d.Write(typed_output);
}

template<class T> inline typename boost::disable_if<HAS_TYPED_WRITE<T> >::type
Write_Binary(TYPED_OSTREAM& output,const T& d)
{
    output.Write(d);
}

template<class T> inline typename boost::enable_if<HAS_TYPED_WRITE<T> >::type
Write_Binary(TYPED_OSTREAM& output,const T& d)
{
    d.Write(output);
}

//#####################################################################
// Multiple Argument Read_Binary
//#####################################################################
template<class RW,class T1,class T2>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);}

template<class RW,class T1,class T2,class T3>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);}

template<class RW,class T1,class T2,class T3,class T4>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3,T4& d4)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);Read_Binary<RW>(input,d4);}

template<class RW,class T1,class T2,class T3,class T4,class T5>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);Read_Binary<RW>(input,d4);Read_Binary<RW>(input,d5);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);Read_Binary<RW>(input,d4);Read_Binary<RW>(input,d5);Read_Binary<RW>(input,d6);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);Read_Binary<RW>(input,d4);Read_Binary<RW>(input,d5);Read_Binary<RW>(input,d6);Read_Binary<RW>(input,d7);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7,T8& d8)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);Read_Binary<RW>(input,d4);Read_Binary<RW>(input,d5);Read_Binary<RW>(input,d6);Read_Binary<RW>(input,d7);
Read_Binary<RW>(input,d8);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
inline void Read_Binary(std::istream& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7,T8& d8,T9& d9)
{Read_Binary<RW>(input,d1);Read_Binary<RW>(input,d2);Read_Binary<RW>(input,d3);Read_Binary<RW>(input,d4);Read_Binary<RW>(input,d5);Read_Binary<RW>(input,d6);Read_Binary<RW>(input,d7);
Read_Binary<RW>(input,d8);Read_Binary<RW>(input,d9);}

template<class T1,class T2>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2)
{Read_Binary(input,d1);Read_Binary(input,d2);}

template<class T1,class T2,class T3>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);}

template<class T1,class T2,class T3,class T4>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3,T4& d4)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);Read_Binary(input,d4);}

template<class T1,class T2,class T3,class T4,class T5>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);Read_Binary(input,d4);Read_Binary(input,d5);}

template<class T1,class T2,class T3,class T4,class T5,class T6>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);Read_Binary(input,d4);Read_Binary(input,d5);Read_Binary(input,d6);}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);Read_Binary(input,d4);Read_Binary(input,d5);Read_Binary(input,d6);Read_Binary(input,d7);}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7,T8& d8)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);Read_Binary(input,d4);Read_Binary(input,d5);Read_Binary(input,d6);Read_Binary(input,d7);
Read_Binary(input,d8);}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
inline void Read_Binary(TYPED_ISTREAM& input,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7,T8& d8,T9& d9)
{Read_Binary(input,d1);Read_Binary(input,d2);Read_Binary(input,d3);Read_Binary(input,d4);Read_Binary(input,d5);Read_Binary(input,d6);Read_Binary(input,d7);
Read_Binary(input,d8);Read_Binary(input,d9);}

//#####################################################################
// Multiple Argument Write_Binary
//#####################################################################
template<class RW,class T1,class T2>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);}

template<class RW,class T1,class T2,class T3>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);}

template<class RW,class T1,class T2,class T3,class T4>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);Write_Binary<RW>(output,d4);}

template<class RW,class T1,class T2,class T3,class T4,class T5>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);Write_Binary<RW>(output,d4);Write_Binary<RW>(output,d5);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);Write_Binary<RW>(output,d4);Write_Binary<RW>(output,d5);Write_Binary<RW>(output,d6);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);Write_Binary<RW>(output,d4);Write_Binary<RW>(output,d5);Write_Binary<RW>(output,d6);
Write_Binary<RW>(output,d7);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7,const T8& d8)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);Write_Binary<RW>(output,d4);Write_Binary<RW>(output,d5);Write_Binary<RW>(output,d6);
Write_Binary<RW>(output,d7);Write_Binary<RW>(output,d8);}

template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
inline void Write_Binary(std::ostream& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7,const T8& d8,const T9& d9)
{Write_Binary<RW>(output,d1);Write_Binary<RW>(output,d2);Write_Binary<RW>(output,d3);Write_Binary<RW>(output,d4);Write_Binary<RW>(output,d5);Write_Binary<RW>(output,d6);
Write_Binary<RW>(output,d7);Write_Binary<RW>(output,d8);Write_Binary<RW>(output,d9);}

template<class T1,class T2>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2)
{Write_Binary(output,d1);Write_Binary(output,d2);}

template<class T1,class T2,class T3>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);}

template<class T1,class T2,class T3,class T4>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);Write_Binary(output,d4);}

template<class T1,class T2,class T3,class T4,class T5>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);Write_Binary(output,d4);Write_Binary(output,d5);}

template<class T1,class T2,class T3,class T4,class T5,class T6>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);Write_Binary(output,d4);Write_Binary(output,d5);Write_Binary(output,d6);}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);Write_Binary(output,d4);Write_Binary(output,d5);Write_Binary(output,d6);
Write_Binary(output,d7);}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7,const T8& d8)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);Write_Binary(output,d4);Write_Binary(output,d5);Write_Binary(output,d6);
Write_Binary(output,d7);Write_Binary(output,d8);}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
inline void Write_Binary(TYPED_OSTREAM& output,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7,const T8& d8,const T9& d9)
{Write_Binary(output,d1);Write_Binary(output,d2);Write_Binary(output,d3);Write_Binary(output,d4);Write_Binary(output,d5);Write_Binary(output,d6);
Write_Binary(output,d7);Write_Binary(output,d8);Write_Binary(output,d9);}

//#####################################################################
// Read/Write_Binary_Array
//#####################################################################
// array is C-style (zero-based) array
template<class RW,class T>
inline void Read_Binary_Array(std::istream& input,T* array,const int number_of_elements)
{if(IS_BINARY_IO_SAFE<T,RW>::value) input.read(reinterpret_cast<char*>(array),number_of_elements*sizeof(T));
else for(int i=0;i<number_of_elements;i++) Read_Write<T,RW>::Read(input,array[i]);}

template<class RW,class T>
inline void Write_Binary_Array(std::ostream& output,const T* array,const int number_of_elements)
{if(IS_BINARY_IO_SAFE<T,RW>::value && number_of_elements) output.write(reinterpret_cast<const char*>(array),number_of_elements*sizeof(T));
else for(int i=0;i<number_of_elements;i++) Read_Write<T,RW>::Write(output,array[i]);}
//#####################################################################
}
#include <Read_Write/READ_WRITE.h>
#endif
