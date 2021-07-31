Thanks for downloading this code.  We hope you find it useful.  We'd really like to see this code get used so if you have any questions don't hesitate to contact us at:

Adam Bargteil (adamb@cs.utah.edu) or
Haimasree Bhattacharya (hb123@cs.utah.edu)

OpenVDB and its denedencies (Boost, OpenEXR, TBB) need to be installed to compile this code. You can get the latest version from http://www.openvdb.org/download/. This code has been tested with version 2.1.0.  You may need to set the path to OpenVDB at the top of the makefile.

You will probably need to edit particleIO.cpp to handle your particle file format(s).

Most of the code has been pretty well tested.  The variable radius code has been less so. We think the knobs we've provided are reasonable, but its hard to say for certain.  Some of the knobs conflict with each other.  In this case the behavior is undefined :)

Cheers.
