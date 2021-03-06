//
// Source code for the paper
// 
// D. Heck and T. Schloemer and O. Deussen, "Blue Noise Sampling with
// Controlled Aliasing", ACM Trans. Graph., 2013, in press
//
// Copyright (C) 2012,2013 Daniel Heck and Thomas Schloemer
//
#include "param.hh"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <string>


static bool IsBool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str == "true" || str == "false";
}

static bool ToBool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str == "true";
}


ParamList::ParamList() {
}

Param *ParamList::Find(const std::string &key) {
    for (int i=0; i<(int)list.size(); ++i) {
        if (list[i].name == key) {
            list[i].used = true;
            return &list[i];
        }
    }
    return NULL;
}

Param *ParamList::Define(const std::string &key, const std::string &dflt) {
    Param *p = Insert(key);
    p->value = dflt;
    p->set = false;
    return p;
}

Param *ParamList::Set(const std::string &key, const std::string &val) {
    Param *p = Insert(key);
    p->value = val;
    p->set = true;
    return p;
}

Param *ParamList::Insert(const std::string &key) {
    if (Param *p = Find(key))
        return p;
    else {
        Param param;
        param.name = key;
        param.set = false;
        param.used = false;
        list.push_back(param);
        return &list[list.size()-1];
    }
}

void ParamList::Parse(int argc, char **argv, std::vector<std::string> &args) {
    for (int index=1; index < argc; index++) {
        std::string arg = argv[index];
        size_t eqpos = arg.find("=");
        size_t beg = arg.find_first_not_of("-");
        std::string name = arg.substr(beg, eqpos-beg);
        if (eqpos != std::string::npos) {
            std::string val = arg.substr(eqpos+1, std::string::npos);
            Set(name, val);
        } else if (!arg.empty() && arg[0] == '-') {
            Param *p = Find(name);
            if (p && IsBool(p->value)) {
                Set(name, "true");
            } else {
                std::string next = ++index < argc ? argv[index] : "";
                Set(name, next);
            }
        } else
            args.push_back(arg);
    }
}

float ParamList::GetFloat(const std::string &key) {
    if (const Param *p = Find(key))
        return atof(p->value.c_str());
    abort();
}

float ParamList::GetFloat(const std::string &key, float dflt) {
    if (const Param *p = Find(key))
        return p->set ? atof(p->value.c_str()) : dflt;
    return dflt;
}

int ParamList::GetInt(const std::string &key, int dflt) {
    if (const Param *p = Find(key))
        return atoi(p->value.c_str());
    return dflt;
}

std::string ParamList::GetString(const std::string &key, std::string dflt) {
    if (const Param *p = Find(key))
        return p->value;
    return dflt;
}

bool ParamList::GetBool(const std::string &key, bool dflt) {
    if (const Param *p = Find(key))
        return ToBool(p->value);
    return dflt;
}

const Param *ParamList::UnusedOption() const {
    for (size_t i=0; i<list.size(); i++) {
        if (!list[i].used && list[i].set)
            return &list[i];
    }
    return NULL;
}

void ParamList::Print() const {
    for (size_t i=0; i<list.size(); i++) {
        const Param *p = &list[i];
        fprintf(stderr, "    %s", p->name.c_str());
        if (!p->value.empty()) 
            fprintf(stderr, " (%s)\n", p->value.c_str());
        else 
            fprintf(stderr, "\n");
    }
}
