These files assume you are using PBRT (http://www.pbrt.org/) version 2.0.0.

Modifications were made to PBRT so that it can read in FIELD_3D files and use them as textures. You will need to copy (in some cases overwrite) the following files into the following directories and rebuild PBRT.

In pbrt/src/core:

  api.cpp
  api.h

In pbrt/src:

  Sconscript

In pbrt/src/textures:

  FIELD_3D.h
  FIELD_3D.cpp
  VEC3.h
  field.cpp
  field.h