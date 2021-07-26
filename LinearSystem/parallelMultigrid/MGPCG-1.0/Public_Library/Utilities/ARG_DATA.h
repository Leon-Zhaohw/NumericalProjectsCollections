//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARG_DATA  
//##################################################################### 
#ifndef __ARG_DATA__
#define __ARG_DATA__

#include <string>
#include <Matrices_And_Vectors/VECTOR.h>
namespace PhysBAM{

class ARG_DATA
{
public:
    enum TYPE{OPTION,INTEGER,DOUBLE,VECTOR2,VECTOR3,STRING};
    TYPE type;
    std::string str,val_name,desc;
    bool boolean_value;
    int  integer_default,integer_value;
    double double_default,double_value;
    VECTOR<double,2> vector_2d_default,vector_2d_value;
    VECTOR<double,3> vector_3d_default,vector_3d_value;
    std::string string_default,string_value;
    bool value_set;

    ARG_DATA()
    {}

    ARG_DATA(const std::string& input_str,const std::string& input_desc)
        :type(OPTION),boolean_value(false),integer_default(0),integer_value(0),double_default(0),double_value(0)
    {
        Initialize(input_str,"",input_desc);
    }

    ARG_DATA(const std::string& input_str,const std::string& input_val_name,const std::string& input_desc,int default_value)
        :type(INTEGER),boolean_value(false),integer_default(default_value),integer_value(default_value),double_default(0),double_value(0)
    {
        Initialize(input_str,input_val_name,input_desc,"int");
    }

    ARG_DATA(const std::string& input_str,const std::string& input_val_name,const std::string& input_desc,double default_value)
        :type(DOUBLE),boolean_value(false),integer_default(0),integer_value(0),double_default(default_value),double_value(default_value)
    {
        Initialize(input_str,input_val_name,input_desc,"double");
    }

    ARG_DATA(const std::string& input_str,const std::string& input_val_name,const std::string& input_desc,const VECTOR<double,2>& default_value)
        :type(VECTOR2),boolean_value(false),integer_default(0),integer_value(0),double_default(0),double_value(0),vector_2d_default(default_value),vector_2d_value(default_value)
    {
        Initialize(input_str,input_val_name,input_desc,"vector2");
    }

    ARG_DATA(const std::string& input_str,const std::string& input_val_name,const std::string& input_desc,const VECTOR<double,3>& default_value)
        :type(VECTOR3),boolean_value(false),integer_default(0),integer_value(0),double_default(0),double_value(0),vector_3d_default(default_value),vector_3d_value(default_value)
    {
        Initialize(input_str,input_val_name,input_desc,"vector3");
    }

    ARG_DATA(const std::string& input_str,const std::string& input_val_name,const std::string& input_desc,const std::string& default_value)
        :type(STRING),boolean_value(false),integer_default(0),integer_value(0),double_default(0),double_value(0),string_default(default_value),string_value(default_value)
    {
        Initialize(input_str,input_val_name,input_desc,"string");
    }

    void Initialize(const std::string& input_str,const std::string& input_val_name,const std::string& input_desc,const std::string& default_val_name="")
    {str=input_str;value_set=false;
    if(input_val_name.length())val_name=input_val_name;else if(default_val_name.length())val_name=default_val_name;
    if(!input_desc.length()){if(input_val_name.length())desc=input_val_name;else desc="<no description available>";}
    else desc=input_desc;}

    bool Parse_Value(int argc,char* argv[],int& current_arg)
    {switch(type){
    case OPTION:boolean_value=true;current_arg++;return true;
    case INTEGER:if(current_arg+1<argc){integer_value=atoi(argv[current_arg+1]);current_arg+=2;return true;}else return false;
    case DOUBLE:if(current_arg+1<argc){double_value=atof(argv[current_arg+1]);current_arg += 2;return true;}else return false;
    case VECTOR2:if(current_arg+2<argc){for(int i=1;i<=2;i++)vector_2d_value[i]=atof(argv[current_arg+i]);current_arg+=3;return true;}else return false;
    case VECTOR3:if(current_arg+3<argc){for(int i=1;i<=3;i++)vector_3d_value[i]=atof(argv[current_arg+i]);current_arg+=4;return true;}else return false;
    case STRING:if(current_arg+1<argc){string_value=argv[current_arg+1];current_arg+=2;return true;}else return false;}
    return false;}

    void Print_Synopsis(std::ostream& os) const
    {if(type==OPTION)os<<"["<<str<<"]";else os<<"["<<str<<" <"<<val_name<<">]";}

    void Print_Description(std::ostream& os,int column_width) const
    {os.flags(std::ios::left);os.width(column_width);
    switch(type){
    case OPTION:os<<str<<desc;break;
    case INTEGER:os<<str<<desc<<" (default "<<integer_default<<")";break;
    case DOUBLE:os<<str<<desc<<" (default "<<double_default<<")";break;
    case VECTOR2:os<<str<<desc<<" (default <"<<vector_2d_default<<">)";break;
    case VECTOR3:os<<str<<desc<<" (default <"<<vector_3d_default<<">)";break;
    case STRING:os<<str<<desc;if(string_default.length())os<<" (default "<<string_default<<")";break;}}

//#####################################################################
};
}
#endif
