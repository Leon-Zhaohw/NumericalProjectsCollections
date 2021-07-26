//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_ARGS  
//##################################################################### 
#ifndef __PARSE_ARGS__
#define __PARSE_ARGS__

#include <Utilities/ARG_DATA.h>
#include <Utilities/NONCOPYABLE.h>
#include <Arrays/LIST_ARRAY.h>
#include <Matrices_And_Vectors/VECTOR.h>
namespace PhysBAM{

class PARSE_ARGS:public NONCOPYABLE
{
private:
    LIST_ARRAY<ARG_DATA> arg_data_list;
    LIST_ARRAY<std::string> extra_arg_list;
    int num_expected_extra_args;
    std::string extra_args_synopsis,extra_args_desc,program_name;
    bool use_help_option;
    void (*extra_usage_callback)();
public:

    PARSE_ARGS() 
        :num_expected_extra_args(0),use_help_option(true),extra_usage_callback(0)
    {}

//#####################################################################
    void Use_Help_Option(bool use_it);
    void Add_Option_Argument(const std::string& arg_str,const std::string& desc="");
    void Add_Integer_Argument(const std::string& arg_str,int default_value,const std::string& val_name="",const std::string& desc="");
    void Add_Double_Argument(const std::string& arg_str,double default_value,const std::string& val_name="",const std::string& desc="");
    void Add_Vector_2D_Argument(const std::string& arg_str,const VECTOR<double,2>& default_value,const std::string& val_name="",const std::string& desc="");
    void Add_Vector_3D_Argument(const std::string& arg_str,const VECTOR<double,3>& default_value,const std::string& val_name="",const std::string& desc="");
    void Add_String_Argument(const std::string& arg_str,const std::string& default_value,const std::string& val_name="",const std::string& desc="");
    void Set_Extra_Arguments(int num,const std::string& synopsis="",const std::string& desc="");
    void Set_Extra_Usage_Callback(void (*extra_usage_callback_input)());
    void Parse(int argc,char* argv[]);
    void Parse(int argc,const char* argv[]); // for backwards compatibility
    bool Get_Option_Value(const std::string& arg_str) const;
    int Get_Integer_Value(const std::string& arg_str) const;
    double Get_Double_Value(const std::string& arg_str) const;
    VECTOR<double,2> Get_Vector_2D_Value(const std::string& arg_str) const;
    VECTOR<double,3> Get_Vector_3D_Value(const std::string& arg_str) const;
    const std::string& Get_String_Value(const std::string& arg_str) const;
    bool Is_Value_Set(const std::string& arg_str) const;
    void Override_String_Value(const std::string& arg_str,const std::string& value);
    int Num_Extra_Args() const;
    const std::string& Extra_Arg(int i) const;    // 1-based indexing
    const std::string& Get_Program_Name() const;
    void Print_Usage(bool do_exit=false) const;
    int Find_Match(const std::string& str) const;
    int Find_Match(const std::string& str,const ARG_DATA::TYPE& type) const;
    static bool Find_And_Remove(const char *str,int& argc,char** argv);
    static void Print_Arguments(int argc,char* argv[]);
//#####################################################################
};
}
#endif
