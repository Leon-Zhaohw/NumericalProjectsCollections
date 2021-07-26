//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header PARTICLES_FORWARD
//#####################################################################
#ifndef __PARTICLES_FORWARD__
#define __PARTICLES_FORWARD__

namespace PhysBAM{

class PARTICLE_ATTRIBUTE_BASE;
template<class T> class PARTICLE_ATTRIBUTE;
template<class T> class PARTICLE_MASS_ATTRIBUTE;
template<class TV> class PARTICLE_POSITION_ATTRIBUTE;
template<class TV> class PARTICLE_VELOCITY_ATTRIBUTE;

class PARTICLE_BASE;
template<class TV> class PARTICLE;
template<class TV> class SPH_PARTICLE;
template<class TV> class SOLIDS_PARTICLE;
template<class TV> class PARTICLE_LEVELSET_PARTICLE;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLE;
template<class TV> class VORTICITY_PARTICLE;
template<class TV> class RIGID_BODY_PARTICLE;

template<class T_PARTICLE> class PARTICLE_POOL;
template<class T_PARTICLE> class PARTICLE_SUBSET;
template<class T_PARTICLE> class SINGLE_PARTICLE;

template<class TV> class PARTICLE_CONNECTIVITY;
}
#endif
