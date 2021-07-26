//#####################################################################
// Copyright 2004-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace STRING_UTILITIES
//#####################################################################
#ifndef __STRING_UTILITIES__
#define __STRING_UTILITIES__

#include <Arrays/ARRAYS_FORWARD.h>
#include <string>
#include <sstream>
namespace PhysBAM{

namespace STRING_UTILITIES{

std::string Stripped_Whitespace(const std::string& str);

inline void Strip_Whitespace(std::string& str)
{str=Stripped_Whitespace(str);}

bool Is_Number(const std::string& str); // integer or floating point
int Compare_Strings(const std::string &str1,const std::string &str2,bool case_sensitive=true);
std::string toupper(const std::string& str);
bool Parse_Integer_Range(const std::string& str,LIST_ARRAY<int>& integer_list);

// integer list format: [range],[range],[range],... where each [range] is either a single <number> or <number>-<number>, e.g. "1-3,4,7,10-11"
bool Parse_Integer_List(const std::string& str,LIST_ARRAY<int>& integer_list);

std::string Join(const std::string& separator,const LIST_ARRAY<std::string>& tokens);
void Split(const std::string& str,const std::string& separator,LIST_ARRAY<std::string>& tokens);

template<class T> inline bool String_To_Value(const std::string& str,T& value) // requires operator>>
{std::istringstream string_stream(str);return (string_stream>>value)!=0;}

template<> inline bool String_To_Value<std::string>(const std::string& str,std::string& value)
{value=str;return true;}

template<class T> inline std::string Value_To_String(const T& str)
{std::ostringstream output;output<<str;return output.str();}

template<> inline std::string Value_To_String<std::string>(const std::string& str)
{return str;}

//#####################################################################
}
}
#endif
