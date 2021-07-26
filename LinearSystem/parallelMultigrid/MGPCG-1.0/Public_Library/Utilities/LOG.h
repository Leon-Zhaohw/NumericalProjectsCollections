//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG
//##################################################################### 
#ifndef __LOG__
#define __LOG__

#include <Utilities/NONCOPYABLE.h>
#include <Utilities/override.h>
#include <Utilities/TIMER.h>
#include <ostream>
#include <boost/format/format_fwd.hpp>
#include <cassert>
namespace PhysBAM{

class LOG_ENTRY;
class LOG_SCOPE;

namespace LOG{

class LOG_CLASS
{
    friend class ::PhysBAM::LOG_ENTRY;
    friend class ::PhysBAM::LOG_SCOPE;
    friend class LOG_COUT_BUFFER;
    friend class LOG_CERR_BUFFER;
    friend void Reset();
    friend void Dump_Log();

    boost::shared_ptr<TIMER> timer_singleton;
    int timer_id;
    bool suppress_cout;
    bool suppress_cerr;
public:
    bool suppress_timing;
    FILE* log_file;
    int verbosity_level;
    bool log_file_temporary;
    bool xml;

    LOG_ENTRY* root;
    LOG_ENTRY* current_entry;

    LOG_CLASS(const bool suppress_cout,const bool suppress_cerr,const bool suppress_timing,const int verbosity_level,const bool cache_initial_output);
    ~LOG_CLASS();

    static void Push_Scope(const std::string& scope_identifier,const std::string& scope_name);
    static void Pop_Scope();
public:
    static void Time_Helper(const std::string& label);
    void Copy_Log_To_File(const std::string& filename,const bool append);

//##################################################################### 
};

class SCOPE:private NONCOPYABLE
{
    bool active;
public:
    SCOPE()
        :active(false)
    {}

    SCOPE(const std::string& scope_identifier)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,scope_identifier);
    }

    SCOPE(const std::string& scope_identifier,const std::string& scope_name)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,scope_name);
    }

    template<class T1>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,str(boost::format(format)%d1));
    }

    template<class T1,class T2>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,str(boost::format(format)%d1%d2));
    }

    template<class T1,class T2,class T3>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2,const T3& d3)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,str(boost::format(format)%d1%d2%d3));
    }

    template<class T1,class T2,class T3,class T4>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,str(boost::format(format)%d1%d2%d3%d4));
    }

    template<class T1,class T2,class T3,class T4,class T5>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,str(boost::format(format)%d1%d2%d3%d4%d5));
    }

    ~SCOPE()
    {if(active) LOG_CLASS::Pop_Scope();}
    
    void Push(const std::string& scope_identifier,const std::string& scope_name)
    {assert(!active);active=true;LOG_CLASS::Push_Scope(scope_identifier,scope_name);}

    template<class T1>
    void Push(const std::string& scope_identifier,const std::string& format,const T1& d1)
    {Push(scope_identifier,str(boost::format(format)%d1));}

    void Pop()
    {assert(active);active=false;LOG_CLASS::Pop_Scope();}
//##################################################################### 
};

// These next few lines are important to ensure no static data from LOG.cpp is accessed for DLLs
LOG_CLASS* Instance();
std::ostream& cout_Helper();
std::ostream& cerr_Helper();
namespace{
    static std::ostream& cout PHYSBAM_UNUSED =::PhysBAM::LOG::cout_Helper();
    static std::ostream& cerr PHYSBAM_UNUSED =::PhysBAM::LOG::cerr_Helper();
}

void Initialize_Logging(const bool suppress_cout_input=false,const bool suppress_timing_input=false,const int verbosity_level_input=1<<30,const bool cache_initial_output=false);
void Finish_Logging();
void Stop_Time();
template<class T_VALUE> void Stat(const std::string& label,const T_VALUE& value);
void Reset(); 
void Dump_Log();

inline void Time(const std::string& format)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(format);}

template<class T1>
inline void Time(const std::string& format,const T1& d1)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(str(boost::format(format)%d1));}

template<class T1,class T2>
inline void Time(const std::string& format,const T1& d1,const T2& d2)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(str(boost::format(format)%d1%d2));}

template<class T1,class T2,class T3>
inline void Time(const std::string& format,const T1& d1,const T2& d2,const T3& d3)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(str(boost::format(format)%d1%d2%d3));}

template<class T1,class T2,class T3,class T4>
inline void Time(const std::string& format,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(str(boost::format(format)%d1%d2%d3%d4));}

template<class T1,class T2,class T3,class T4,class T5>
inline void Time(const std::string& format,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(str(boost::format(format)%d1%d2%d3%d4%d5));}

//##################################################################### 
}
}
#endif
