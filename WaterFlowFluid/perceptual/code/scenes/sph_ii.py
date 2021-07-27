# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL)
# http://www.gnu.org/licenses
#
# Implicit incompressible SPH (IISPH) simulation
#
# ----------------------------------------------------------------------------

import math

guion = False
pause = True

# default solver parameters
params          = {}
params['dim']   = 2                  # dimension
params['sres']  = 2                  # sub-resolution per cell
params['dx']    = 1.0/params['sres'] # particle spacing (= 2 x radius)
params['res']   = 50                 # reference resolution
params['len']   = 1.0                # reference length
params['bnd']   = 4                  # boundary cells
params['gref']  = -9.8               # real-world gravity
params['dens']  = 1000.0             # density
params['avis']  = True               # aritificial viscosity
params['eta']   = 0.1
params['fps']   = 30
params['t_end'] = 5.0
params['sdt']   = None

scaleToManta   = float(params['res'])/params['len']
# NOTE: the original test uses 3.22; but, here it's slightly modified for the sake of convenience in discretization
params['gs']   = [round(float(params['res'])*3.2)+params['bnd']*2, params['res']*3+params['bnd']*2, params['res']+params['bnd']*2 if params['dim']==3 else 1]
params['grav'] = params['gref']*scaleToManta

s             = Solver(name='IISPH', gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = s.frameLength
s.timestep    = s.frameLength

overFld = FlagFluid
overAll = FlagFluid|FlagObstacle

sph  = s.create(SphWorld, delta=params['dx'], density=params['dens'], g=(0,params['grav'],0), eta=params['eta'])
kern = s.create(CubicSpline, h=sph.delta)
print('h = {}, sr = {}'.format(kern.radius(), kern.supportRadius()))

pp = s.create(BasicParticleSystem)

# acceleration data for particle neighbors
gIdxSys  = s.create(ParticleIndexSystem)
gIdx     = s.create(IntGrid)
gCnt     = s.create(IntGrid)
gFlags   = s.create(FlagGrid)
neighbor = s.create(ParticleNeighbors)

pT = pp.create(PdataInt)        # particle type
pV = pp.create(PdataVec3)       # velocity
pF = pp.create(PdataVec3)       # force
pD = pp.create(PdataReal)       # density
pP = pp.create(PdataReal)       # pressure

pDadv  = pp.create(PdataReal)   # density advected
pAii   = pp.create(PdataReal)   # a_ii
pDii   = pp.create(PdataVec3)   # d_ii
pDijPj = pp.create(PdataVec3)   # sum_j(d_ii*pj)

pDtmp = pp.create(PdataReal)
pVtmp = pp.create(PdataVec3)

mesh = {}
if params['dim']==3:
	mesh['mesh']     = s.create(Mesh)
	mesh['levelset'] = s.create(LevelsetGrid)

# boundary setup
gFlags.initDomain(params['bnd']-1)

begin = pp.size()
sampleFlagsWithParticles(flags=gFlags, parts=pp, discretization=params['sres'], randomness=0, ftype=FlagObstacle, notiming=True)
end = pp.size()
pT.setConstRange(s=FlagObstacle, begin=begin, end=end, notiming=True)

# obstacle
a   = vec3(0.744*scaleToManta+params['bnd'], 0.161*0.5*scaleToManta+params['bnd'], 0.5*params['gs'][2] if (params['dim']==3) else 0)
b   = vec3(0.161*0.5*scaleToManta, 0.161*0.5*scaleToManta, 0.403*0.5*scaleToManta if (params['dim']==3) else params['gs'][2])
obs = s.create(Box, center=a, size=b)

begin = pp.size()
sampleShapeWithParticles(shape=obs, flags=gFlags, parts=pp, discretization=params['sres'], randomness=0, notiming=True)
end = pp.size()
pT.setConstRange(s=FlagObstacle, begin=begin, end=end, notiming=True)
obs.applyToGrid(grid=gFlags, value=FlagObstacle, respectFlags=gFlags)

# fluid setup: dam
dam_c = [2.606, 0.275, 0.5]
dam_s = [1.228*0.5, 0.55*0.5, 0.5]
a     = vec3(dam_c[0]*scaleToManta+params['bnd'], dam_c[1]*scaleToManta+params['bnd'], dam_c[2]*scaleToManta+params['bnd'] if (params['dim']==3) else 0)
b     = vec3(dam_s[0]*scaleToManta, dam_s[1]*scaleToManta, dam_s[2]*scaleToManta if (params['dim']==3) else params['gs'][2])
fld   = s.create(Box, center=a, size=b)

begin = pp.size()
sampleShapeWithParticles(shape=fld, flags=gFlags, parts=pp, discretization=params['sres'], randomness=0, notiming=True)
end = pp.size()
pT.setConstRange(s=FlagFluid, begin=begin, end=end, notiming=True)

sph.bindParticleSystem(p_system=pp, p_type=pT, p_neighbor=neighbor, notiming=True)
sph.updateSoundSpeed(math.sqrt(2.0*math.fabs(params['grav'])*0.55*scaleToManta/params['eta']), notiming=True)
pD.setConst(s=sph.density, notiming=True)
gridParticleIndex(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, counter=gCnt, notiming=True)
neighbor.update(pts=pp, indexSys=gIdxSys, index=gIdx, radius=kern.supportRadius(), notiming=True)

if guion:
	gui = Gui()
	gui.show()
	if pause: gui.pause()

while (s.timeTotal<params['t_end']): # main loop
	sphComputeDensity(d=pD, k=kern, sph=sph, itype=overFld, jtype=overAll)
	sphComputeConstantForce(f=pF, v=vec3(0, params['grav']*sph.mass, 0), sph=sph, itype=overFld, accumulate=False)
	sphComputeSurfTension(f=pF, k=kern, sph=sph, kappa=0.8, itype=overFld, jtype=overAll, accumulate=True)
	if(params['avis']):
		sphComputeArtificialViscousForce(f=pF, v=pV, d=pD, k=kern, sph=sph, itype=overFld, jtype=overFld, accumulate=True)

	if params['sdt'] is None:
		adt = min(s.frameLength, kern.supportRadius()/sph.c)
		adt = sph.limitDtByVmax(dt=adt, h=kern.supportRadius(), vmax=pV.getMaxAbsValue(), a=0.4)
		s.adaptTimestepByDt(adt)
	else:
		s.adaptTimestepByDt(params['sdt'])

	sphUpdateVelocity(v=pVtmp, vn=pV, f=pF, sph=sph, dt=s.timestep)
	sphComputeIisphDii(dii=pDii, d=pD, k=kern, sph=sph, dt=s.timestep, itype=overFld, jtype=overAll)

	pDadv.setConst(0)
	sphComputeDivergenceSimple(div=pDadv, v=pVtmp, k=kern, sph=sph, itype=overFld, jtype=overAll) # pDadv = div(v)
	pDadv.multConst(s=-s.timestep)                                                                # pDadv = - dt*div(v)
	pDadv.add(pD)                                                                                 # pDadv = pD - dt*div(v)
	pAii.setConst(0)
	sphComputeIisphAii(aii=pAii, d=pD, dii=pDii, k=kern, sph=sph, dt=s.timestep, itype=overFld, jtype=overAll)

	######################################################################
	# solve pressure
	pP.multConst(s=0.5)         # p = 0.5*p_prev
	d_avg, iters, d_err_th = sph.density, 0, sph.density*sph.eta/100.0
	while ((d_avg - sph.density)>d_err_th) or (iters<2):
		sphComputeIisphDijPj(dijpj=pDijPj, d=pD, p=pP, k=kern, sph=sph, dt=s.timestep, itype=overFld, jtype=overAll)

		pDtmp.setConst(0.0)
		sphComputeIisphP(p_next=pDtmp, p=pP, d_adv=pDadv, d=pD, aii=pAii, dii=pDii, dijpj=pDijPj, k=kern, sph=sph, dt=s.timestep, itype=overFld, jtype=overAll)
		pDtmp.clampMin(0.0)
		pP.copyFrom(pDtmp)

		pDtmp.setConst(0.0)
		sphComputeIisphD(d_next=pDtmp, d_adv=pDadv, d=pD, p=pP, dii=pDii, dijpj=pDijPj, k=kern, sph=sph, dt=s.timestep, itype=overFld, jtype=overAll)
		d_avg = pDtmp.sum(t=pT, itype=overFld)/cntPts(t=pT, itype=overFld)

		iters += 1

		# for the safety
		if iters>999:
			print('\tFail to converge: d_avg = {} (<{}), iters = {}'.format(d_avg, d_err_th+sph.density, iters))
			sys.exit(0)

	print('\td_avg = {} (<{}), iters = {}'.format(d_avg, d_err_th+sph.density, iters))
	######################################################################

	sphComputePressureForce(f=pF, p=pP, d=pD, k=kern, sph=sph, accumulate=False)
	sphUpdateVelocity(v=pV, vn=pVtmp, f=pF, sph=sph, dt=s.timestep)

	sphUpdatePosition(x=pp, v=pV, sph=sph, dt=s.timestep)
	gridParticleIndex(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, counter=gCnt)
	neighbor.update(pts=pp, indexSys=gIdxSys, index=gIdx, radius=kern.supportRadius())

	if params['dim']==3 and guion:
		unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=mesh['levelset'], radiusFactor=1.0, ptype=pT, exclude=FlagObstacle)
		extrapolateLsSimple(phi=mesh['levelset'], distance=4, inside=True)
		mesh['levelset'].createMesh(mesh['mesh'])

	s.step()
