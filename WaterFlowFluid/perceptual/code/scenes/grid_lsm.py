# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL)
# http://www.gnu.org/licenses
#
# Level set simulation
#
# ----------------------------------------------------------------------------

guion = False
pause = True

upls         = 2                # up-scaling factor for surface tracking; don't use a factor greater than 2. (Filtering has not been implemented yet.)
fastMarching = True             # use the fast-marching method
lsBcNeumann  = False            # level-set boundary condition: False means Dirichlet

# default solver parameters
params               = {}
params['dim']        = 2        # dimension
params['res']        = 50       # reference resolution
params['len']        = 1.0      # reference length
params['bnd']        = 4        # boundary cells
params['gref']       = -9.8     # real-world gravity
params['cgaccuracy'] = 1e-3     # cg solver's threshold
params['gfm']        = True     # 2nd order fluid-empty BC
params['fps']        = 30       # frames per second
params['t_end']      = 5.0      # quit simulation
params['sdt']        = None     # fix timestep size

# scale unit in regard to the manta world
scaleToManta   = float(params['res'])/params['len']
# NOTE: the original test uses 3.22; but, here it's slightly modified for the sake of convenience in discretization
params['gs']   = [round(float(params['res'])*3.2)+params['bnd']*2, params['res']*3+params['bnd']*2, params['res']+params['bnd']*2 if params['dim']==3 else 1]
params['grav'] = params['gref']*scaleToManta

s             = Solver(name="LS", gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = s.frameLength
s.timestep    = s.frameLength

# prepare grids and particles
gFlags  = s.create(FlagGrid)
gV      = s.create(MACGrid)
gP      = s.create(RealGrid)
gPhi    = s.create(LevelsetGrid)
gPhiSld = s.create(LevelsetGrid)

if upls>1:
	xl_gs            = vec3(params['gs'][0]*upls, params['gs'][1]*upls, params['gs'][2]*upls if params['dim']==3 else 1)
	xl_s             = Solver(name='Upscale LS',  gridSize=xl_gs, dim=params['dim'])
	xl_s.frameLength = xl_s.frameLength
	xl_gFlags        = xl_s.create(FlagGrid)
	xl_gV            = xl_s.create(MACGrid)
	xl_gPhi          = xl_s.create(LevelsetGrid)
	xl_gPhiSld       = xl_s.create(LevelsetGrid)

phiMesh = xl_s.create(LevelsetGrid) if upls>1 else s.create(LevelsetGrid)
mesh    = xl_s.create(Mesh) if upls>1 else s.create(Mesh)

paramSolvePressure = dict(flags=gFlags, vel=gV, pressure=gP, cgAccuracy=params['cgaccuracy'])
if params['gfm']:               # for the free-surface boundary condition
	paramSolvePressure.update(phi=gPhi)

# boundary setup
gFlags.initDomain(params['bnd']-1)
bndBox = s.create(Box, p0=vec3(0), p1=vec3(params['gs'][0], params['gs'][1], params['gs'][2]))
inBox  = s.create(Box, p0=vec3(params['bnd'], params['bnd'], params['bnd'] if params['dim']==3 else 0), p1=vec3(params['gs'][0]-params['bnd'], params['gs'][1]-params['bnd'], (params['gs'][0]-params['bnd']) if params['dim']==3 else 1))
gPhiSld.join(bndBox.computeLevelset(notiming=True), notiming=True)
gPhiSld.subtract(inBox.computeLevelset(notiming=True), notiming=True)
if upls>1:
	xl_gFlags.initDomain(params['bnd']*upls-1)
	bndBox = xl_s.create(Box, p0=vec3(0), p1=xl_gs)
	inBox  = xl_s.create(Box, p0=vec3(params['bnd'], params['bnd'], params['bnd'] if params['dim']==3 else 0)*upls, p1=vec3((params['gs'][0]-params['bnd'])*upls, (params['gs'][1]-params['bnd'])*upls, (params['gs'][0]-params['bnd'])*upls if params['dim']==3 else 1))
	xl_gPhiSld.join(bndBox.computeLevelset(notiming=True), notiming=True)
	xl_gPhiSld.subtract(inBox.computeLevelset(notiming=True), notiming=True)

# obstacle
a   = vec3(0.744*scaleToManta+params['bnd'], 0.161*0.5*scaleToManta+params['bnd'], 0.5*params['gs'][2] if (params['dim']==3) else 0)
b   = vec3(0.161*0.5*scaleToManta, 0.161*0.5*scaleToManta, 0.403*0.5*scaleToManta if (params['dim']==3) else params['gs'][2])
obs = s.create(Box, center=a, size=b)
obs.applyToGrid(grid=gFlags, value=FlagObstacle, respectFlags=gFlags)
gPhiSld.join(obs.computeLevelset(notiming=True), notiming=True)
if upls>1:
	a   = vec3((0.744*scaleToManta+params['bnd'])*upls, (0.161*0.5*scaleToManta+params['bnd'])*upls, 0.5*params['gs'][2]*upls if (params['dim']==3) else 0)
	b   = vec3(0.161*0.5*scaleToManta*upls, 0.161*0.5*scaleToManta*upls, 0.403*0.5*scaleToManta*upls if (params['dim']==3) else params['gs'][2])
	obs = xl_s.create(Box, center=a, size=b)
	obs.applyToGrid(grid=xl_gFlags, value=FlagObstacle, respectFlags=xl_gFlags)
	xl_gPhiSld.join(obs.computeLevelset(notiming=True), notiming=True)

# fluid setup: dam
dam_c = [2.606, 0.275, 0.5]
dam_s = [1.228*0.5, 0.55*0.5, 0.5]
a     = vec3(dam_c[0]*scaleToManta+params['bnd'], dam_c[1]*scaleToManta+params['bnd'], dam_c[2]*scaleToManta+params['bnd'] if (params['dim']==3) else 0)
b     = vec3(dam_s[0]*scaleToManta, dam_s[1]*scaleToManta, dam_s[2]*scaleToManta if (params['dim']==3) else params['gs'][2])
fld   = s.create(Box, center=a, size=b)
gPhi.join(fld.computeLevelset(), notiming=True)
fld.applyToGrid(grid=gFlags, value=FlagFluid, respectFlags=gFlags)
if upls>1:
	a     = vec3(dam_c[0]*scaleToManta+params['bnd'], dam_c[1]*scaleToManta+params['bnd'], dam_c[2]*scaleToManta+params['bnd'] if (params['dim']==3) else 0)*upls
	b     = vec3(dam_s[0]*scaleToManta*upls, dam_s[1]*scaleToManta*upls, dam_s[2]*scaleToManta*upls if (params['dim']==3) else params['gs'][2])
	fld   = xl_s.create(Box, center=a, size=b)
	xl_gPhi.join(fld.computeLevelset(), notiming=True)
	fld.applyToGrid(grid=xl_gFlags, value=FlagFluid, respectFlags=xl_gFlags)

# lsm special
if (params['dim']==3) and guion:
	phiMesh.copyFrom(xl_gPhi if upls>1 else gPhi)
	phiMesh.subtract(xl_gPhiSld if upls>1 else gPhiSld)
	phiMesh.createMesh(mesh)

if guion:
	gui = Gui()
	gui.show()
	if pause: gui.pause()

while (s.timeTotal<params['t_end']): # main loop
	if params['sdt'] is None: s.adaptTimestep(gV.getMaxValue())
	else: s.adaptTimestepByDt(params['sdt'])

	if upls>1:
		xl_s.timestep = s.timestep*upls

		if fastMarching:
			xl_gPhi.reinitMarching(flags=xl_gFlags, velTransport=xl_gV, maxTime=8.0)
		else :
			extrapolateLsSimple(phi=xl_gPhi, distance=5*upls, inside=False)
			extrapolateLsSimple(phi=xl_gPhi, distance=5*upls, inside=True)

		advectSemiLagrange(flags=xl_gFlags, vel=xl_gV, grid=xl_gPhi, order=1)
		if lsBcNeumann: xl_gPhi.setBoundNeumann(boundaryWidth=params['bnd']*upls-1)
		else: xl_gPhi.subtract(xl_gPhiSld)
		xl_gFlags.updateFromLevelset(levelset=xl_gPhi)
		gFlags.minifyFrom(flags=xl_gFlags, scale=vec3(upls)) # NOTE: or use gFlags.updateFromLevelset(xl_gPhi)

		interpolateGrid(source=xl_gPhi, target=gPhi)
		gPhi.multConst(1.0/upls)

	else:
		if fastMarching:
			gPhi.reinitMarching(flags=gFlags, velTransport=gV)
		else:
			extrapolateLsSimple(phi=gPhi, distance=5, inside=False)
			extrapolateLsSimple(phi=gPhi, distance=5, inside=True)

		advectSemiLagrange(flags=gFlags, vel=gV, grid=gPhi, order=1)
		if lsBcNeumann: gPhi.setBoundNeumann(boundaryWidth=params['bnd']-1)
		else: gPhi.subtract(gPhiSld)
		gFlags.updateFromLevelset(gPhi)

	advectSemiLagrange(flags=gFlags, vel=gV, grid=gV, order=1, boundaryWidth=params['bnd'])
	addGravityNoScale(flags=gFlags, vel=gV, gravity=vec3(0, params['grav'], 0))

	# solve pressure
	setWallBcs(flags=gFlags, vel=gV)
	solvePressure(**paramSolvePressure)
	setWallBcs(flags=gFlags, vel=gV)
	extrapolateMACSimple(flags=gFlags, vel=gV)

	if upls>1: interpolateMACGrid(source=gV, target=xl_gV)

	s.step()

	# lsm special
	if (params['dim']==3) and guion:
		phiMesh.copyFrom(xl_gPhi if upls>1 else gPhi)
		phiMesh.subtract(xl_gPhiSld if upls>1 else gPhiSld)
		phiMesh.createMesh(mesh)
