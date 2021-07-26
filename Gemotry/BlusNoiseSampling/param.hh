//
// Source code for the paper
// 
// D. Heck and T. Schloemer and O. Deussen, "Blue Noise Sampling with
// Controlled Aliasing", ACM Trans. Graph., 2013, in press
//
// Copyright (C) 2012,2013 Daniel Heck and Thomas Schloemer
//
#ifndef PARAM_HH_INCLUDED
#define PARAM_HH_INCLUDED

#include <exception>
#include <string>
#include <vector>

// Command line parsing. Parameters can be specified on the command 
// line in one of the following forms:
//
//  --param value
//  --param=value
//  param=value

struct Param {
    std::string name;
    std::string value;
    bool set;                   // was specified on command line
    bool used;                  // was queried by program
};

class ParamList {
    std::vector<Param> list;
public:
    ParamList();

    Param *Define(const std::string &key, const std::string &dflt);
    void Parse(int argc, char **argv, std::vector<std::string> &remaining);

    float GetFloat(const std::string &key);
    float GetFloat(const std::string &key, float dflt);
    int GetInt(const std::string &key, int dflt = 0);
    std::string GetString(const std::string &key, std::string dflt = "");
    bool GetBool(const std::string &key, bool dflt = false);

    const Param *UnusedOption() const;

    void Print() const;
private:
    Param *Set(const std::string &key, const std::string &val);
    Param *Find(const std::string &key);
    Param *Insert(const std::string &key);
};

#endif
