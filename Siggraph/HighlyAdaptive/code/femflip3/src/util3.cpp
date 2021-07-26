/*
 *	util3.cpp
 *
 *	Created by Ryoichi Ando on 4/15/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "util3.h"

// Global variables
int nestLevel=0;
int stepNumber=0;
char root_path_string[512];
char *root_path = root_path_string;
uint console_num = 0;

bool ended_with_return=false;
std::stack<FLOAT64> global_timestack;
std::map<std::string,FLOAT64> status_table;