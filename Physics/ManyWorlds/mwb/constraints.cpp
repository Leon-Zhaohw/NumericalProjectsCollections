#include "stdafx.h"

#include "constraints.h"
#include "mechModel.h"
#include "simulationTree.h"
#include "sphere_tree.h"

#include "twigg/EulerAngles.h"
#include "twigg/linalg.h"
#include "twigg/renderUtils.h"
#include "twigg/stlext.h"

#include <numeric>

namespace planning {

Objective::Objective(const std::string& name)
	: Named(name)
{
}

Objective::~Objective()
{
}

vl::Vec3f Objective::viewDirection(const vl::Vec3f& position, const CameraSource& source, float time) const
{
	vl::Mat4f invModelView = source.inverseModelViewMatrix(time);
	vl::Mat4f modelView = source.modelViewMatrix(time);
	vl::Mat4f projection = source.projectionMatrix(time);
	vl::Mat4f invProjection = vl::inv( projection );

	// first, transform the position into clip coordinates
	vl::Vec3f clipPos = vl::xform( projection, vl::xform(modelView, position) );

	clipPos[2] = 0.0;
	vl::Vec3f nearPos = vl::xform( invModelView, vl::xform(invProjection, clipPos) );

	return vl::norm(position - nearPos);
}

NormalizedIntegratedObjective::NormalizedIntegratedObjective(const std::string& name)
	: Objective(name)
{
}

/*
std::vector<float> NormalizedIntegratedObjective::evaluate( 
	const SimulationNode* tail,
	const SimulationTree& tree,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	std::vector<float> endTimes( objects.size(), tail->state().time() );
	std::vector<bool> moved( objects.size(), false );

	std::vector<float> result( objects.size(), 0.0f );

	// we are allowed to integrate backwards, so we do so
	const SimulationNode* current = tail;
	while( current->parent() != 0 )
	{
		float timeDiff = current->state().time() - current->parent()->state().time();
		std::vector<float> values = this->evaluateState(
			current->state(),
			objects,
			cam );

		assert( values.size() == result.size() );
		for( size_t i = 0; i < result.size(); ++i )
			result[i] += timeDiff*values[i];

		for( size_t iObject = 0; iObject < current->state().stateCount(); ++iObject )
		{
			moved[iObject] = moved[iObject] || !tree.idle( iObject, current->prevState(), current->state() );
			if( !moved[iObject] )
				endTimes[iObject] = current->state().time();
		}

		current = current->parent();
	}

	for( size_t i = 0; i < result.size(); ++i )
		result[i] /= endTimes[i];

	return result;
}
*/

NormalizedIntegratedObjective::~NormalizedIntegratedObjective()
{
}

float NormalizedIntegratedObjective::combine( float oldValue, float newValue ) const
{
	return (oldValue + newValue);
}

float NormalizedIntegratedObjective::startValue() const
{
	return 0.0f;
}

AngularEnergyObjective::AngularEnergyObjective()
	: NormalizedIntegratedObjective("Angular energy")
{
}

float AngularEnergyObjective::evaluate( 
	const PiecewisePath& path,
	size_t frameRate,
	ConstPhysicsObjectPtr object,
	const CameraSource& source ) const
{
	float result = path.angularEnergy() 
		/ boost::numeric_cast<float>( frameRate );
	return result;
}

FinalOrientationObjective::FinalOrientationObjective()
	: Objective( "Ending orientation" )
{
}

float FinalOrientationObjective::evaluate( 
	const PiecewisePath& path,
	size_t frameRate,
	ConstPhysicsObjectPtr object,
	const CameraSource& source ) const
{
	vl::Mat3d startMat = path.rotation( path.startFrame() ).toRotMatd();
	vl::Mat3d endMat = path.rotation( path.endFrame() ).toRotMatd();

	float result = 0.0f;
	for( int i = 0; i < 3; ++i )
	{
		for( int j = 0; j < 3; ++j )
		{
			const float diff = startMat[i][j] - endMat[i][j];
			result += diff*diff;
		}
	}

	return result;
}

float FinalOrientationObjective::combine( float oldValue, float newValue ) const
{
	return oldValue + newValue;
}

float FinalOrientationObjective::startValue() const
{
	return 0.0f;
}

/*
std::vector<float> AngularEnergyObjective::evaluateState( 
	const Simulation::State& state,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	std::vector<float> result;
	result.reserve( objects.size() );
	for( size_t iObject = 0; iObject < objects.size(); ++iObject )
	{
		NxVec3 angularVel = state.state(iObject).angularVelocity();
		result.push_back( angularVel.dot(angularVel) );
	}

	return result;
}

ViewDependentAngularEnergyObjective::ViewDependentAngularEnergyObjective()
	: NormalizedIntegratedObjective("View-dependent angular energy")
{
}

std::vector<float> ViewDependentAngularEnergyObjective::evaluateState( 
	const Simulation::State& state,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	std::vector<float> result;
	result.reserve( objects.size() );
	for( size_t iObject = 0; iObject < objects.size(); ++iObject )
	{
		NxVec3 V = toNxVec( 
			this->viewDirection( toVec3f(state.state(iObject).position()), cam, state.time() ) );

		NxVec3 angularVel = state.state(iObject).angularVelocity();
		result.push_back( angularVel.cross(V).magnitudeSquared() );
	}

	return result;
}

RunTimeObjective::RunTimeObjective()
	: Objective("Running time")
{
}

std::vector<float> RunTimeObjective::evaluate( 
	const SimulationNode* tail,
	const SimulationTree& tree,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	return std::vector<float>( 1, tail->state().time() );
}

CountCollisionsObjective::CountCollisionsObjective()
	: Objective("Collision count")
{
}

std::vector<float> CountCollisionsObjective::evaluate( 
	const SimulationNode* tail,
	const SimulationTree& tree,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	// we'll count to integers first then convert to float
	std::vector<size_t> counts( objects.size(), 0 );

	// we are allowed to integrate backwards, so we do so
	const SimulationNode* current = tail;

	std::vector<Simulation::Interaction> prevInteractions;
	while( current->parent() != 0 )
	{
		std::vector<Simulation::Interaction> newInteractions = 
			current->interactions();

		for( std::vector<Simulation::Interaction>::const_iterator interactionIter = newInteractions.begin();
			interactionIter != newInteractions.end(); ++interactionIter )
		{
			if( std::binary_search(prevInteractions.begin(), prevInteractions.end(), *interactionIter) )
				continue;

			++counts[ interactionIter->first ];
		}

		current = current->parent();

		newInteractions.swap( prevInteractions );
	}

	size_t num = std::accumulate( counts.begin(), counts.end(), 0 );
	std::vector<float> result( 1, num );
	return result;

}

CollisionCoverageObjective::CollisionCoverageObjective()
	: Objective("Collision coverage")
{
}

std::vector<float> CollisionCoverageObjective::evaluate( 
	const SimulationNode* tail,
	const SimulationTree& tree,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	// we are allowed to integrate backwards, so we do so
	const SimulationNode* current = tail;

	std::vector<Simulation::Interaction> allInteractions;
	while( current->parent() != 0 )
	{
		std::vector<Simulation::Interaction> newInteractions = 
			current->interactions();

		typedef std::vector<Simulation::Interaction>::iterator Iter;
		typedef std::pair<Iter, Iter> IterPr;

		for( std::vector<Simulation::Interaction>::const_iterator interactionIter = newInteractions.begin();
			interactionIter != newInteractions.end(); ++interactionIter )
		{
			IterPr pr = std::equal_range(allInteractions.begin(), allInteractions.end(), *interactionIter);
			if( pr.first == pr.second )
				allInteractions.insert( pr.first, *interactionIter );
		}

		current = current->parent();
	}

	std::vector<float> result( 1, allInteractions.size() );
	return result;
}

StackingObjective::StackingObjective()
	: Objective("Stacking")
{
}

std::vector<float> StackingObjective::evaluate(
	const SimulationNode* tail,
	const SimulationTree& tree,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	std::vector<bool> moved( objects.size(), false );

	typedef std::vector<Simulation::Interaction> InteractionList;
	std::vector<InteractionList> interactions( objects.size() );

	// we are allowed to integrate backwards, so we do so
	const SimulationNode* current = tail;
	while( current->parent() != 0 )
	{
		for( size_t iObject = 0; iObject < current->state().stateCount(); ++iObject )
		{
			moved[iObject] = moved[iObject] || 
				!tree.idle( iObject, current->prevState(), current->state() );
		}

		InteractionList currentInteractions = current->interactions();
		for( InteractionList::const_iterator interactionItr = currentInteractions.begin();
			interactionItr != currentInteractions.end(); ++interactionItr )
		{
			if( !moved[ interactionItr->first ] && !moved[ interactionItr->second ] )
			{
				interactions[interactionItr->first].push_back( *interactionItr );
				interactions[interactionItr->second].push_back( *interactionItr );
			}
		}

		current = current->parent();
	}

	std::vector<float> result( objects.size() );
	for( size_t i = 0; i < interactions.size(); ++i )
	{
		std::sort( interactions[i].begin(), interactions[i].end() );
		result[i] = std::distance( interactions[i].begin(),
			std::unique(interactions[i].begin(), interactions[i].end()) );
	}

	return result;
}

YPositionObjective::YPositionObjective()
	: Objective( "Y position" )
{
}

std::vector<float> YPositionObjective::evaluate(
	const SimulationNode* tail,
	const SimulationTree& tree,
	const std::vector<ConstPhysicsObjectPtr>& objects,
	const CameraSource& cam ) const
{
	std::vector<float> result( objects.size() );

	const Simulation::State& state = tail->state();
	for( size_t iObject = 0; iObject < state.stateCount(); ++iObject )
	{
		NxVec3 pos = state.state(iObject).position();
		result[iObject] = pos[1]*pos[1];
	}

	return result;
}
*/

Constraint::Constraint( const std::vector<size_t>& objects )
	: objects_( objects )
{
}

Constraint::~Constraint()
{
}

Constraint::object_iterator Constraint::objects_begin() const
{
	return objects_.begin();
}

Constraint::object_iterator Constraint::objects_end() const
{
	return objects_.end();
}

ProjectiveConstraint::ProjectiveConstraint( 
	const std::vector<size_t>& objects, 
	const vl::Mat4f& matrix, 
	const BoundingBox2f& box )
	: Constraint(objects), matrix_(matrix), box_(box), query_(matrix, box)
{
	/*
	ofstream ofs( "constraints.bin", std::ios::binary | std::ios::app );
	ofs.write( reinterpret_cast<const char*>( &matrix_ ), sizeof( vl::Mat4f ) );
	ofs.write( reinterpret_cast<const char*>( &box_ ), sizeof( BoundingBox2f ) );
	*/
}

ProjectiveConstraint::~ProjectiveConstraint()
{
}

vl::Mat4f ProjectiveConstraint::matrix() const
{
	return this->matrix_;
}

BoundingBox2f ProjectiveConstraint::box() const
{
	return this->box_;
}

bool ProjectiveConstraint::contained( const Path& path ) const
{
	for( Constraint::object_iterator objectItr = this->objects_begin();
		objectItr != this->objects_end(); ++objectItr )
	{
		if( path.obbTree(*objectItr).query( this->query_ ) )
			return true;
	}

	return false;
}

bool ProjectiveConstraint::intersect(const Ray3d& ray, double& t ) const
{
	vl::Mat4f invMatrix = inv(this->matrix_);

	vl::Vec2f min = this->box_.minimum();
	vl::Vec2f max = this->box_.maximum();

	float near_z = -0.9999999;
	float far_z = 1.0;
	boost::array<vl::Vec3f, 4> closePoints = {{
		vl::xform( invMatrix, vl::Vec3f(min[0], min[1], near_z) ),
		vl::xform( invMatrix, vl::Vec3f(min[0], max[1], near_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], max[1], near_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], min[1], near_z) ) }};
	boost::array<vl::Vec3f, 4> farPoints = {{
		vl::xform( invMatrix, vl::Vec3f(min[0], min[1], far_z) ),
		vl::xform( invMatrix, vl::Vec3f(min[0], max[1], far_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], max[1], far_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], min[1], far_z) ) }};

	bool intersected = false;
	t = boost::numeric::bounds<double>::highest();
	for( int i = 0; i < 4; ++i )
	{
		boost::array<vl::Vec3f, 4> quadPoints = 
			{{ closePoints[(i+0)%4], closePoints[(i+1)%4], farPoints[(i+1)%4], farPoints[(i+0)%4] }};
		for( int start = 0; start <= 2; start += 2 )
		{
			std::pair< bool, TriangleIntersection3d > result =
				intersectTriangle( ray, 
					toVec3d(quadPoints[start]), 
					toVec3d(quadPoints[start+1]), 
					toVec3d(quadPoints[(start+2)%4]) );
			if( !result.first )
				continue;

			intersected = true;
			t = std::min( t, result.second.t() );
		}
	}

	return intersected;
	/*
	vl::Mat4d matrix = toMat4d( this->matrix_ );
	vl::Vec3d pos = vl::xform( matrix, ray.getPosition() );
	vl::Vec3d dir = vl::xform( matrix, ray.getPosition() + ray.getDirection() ) - pos;	
	Ray3d clipRay( pos, dir );
	BoundingBox3d box( vl::Vec3d( toVec2d(box_.minimum()), -1.0 ),
		vl::Vec3d( toVec2d(box_.maximum()), 1.0 ) );

	double tmin;
	double tmax;
	bool intersect = box.intersect(clipRay, tmin, tmax);
	if( !intersect )
		return false;

	double tclip = (tmin >= 0) ? tmin : tmax;
	if( tclip < 0.0 )
		return false;

	vl::Vec3d clipPos = clipRay.at( tclip );
	vl::Vec3d globalPos = vl::xform( vl::inv(matrix), clipPos );
	t = vl::dot( globalPos - ray.getPosition(), ray.getDirection() );
	return true;
	*/
}


void ProjectiveConstraint::render(bool selected) const
{
	vl::Mat4f invMatrix = inv(this->matrix_);

	vl::Vec2f min = this->box_.minimum();
	vl::Vec2f max = this->box_.maximum();

	float near_z = -0.99;
	float far_z = 0.99;
	boost::array<vl::Vec3f, 4> closePoints = {{
		vl::xform( invMatrix, vl::Vec3f(min[0], min[1], near_z) ),
		vl::xform( invMatrix, vl::Vec3f(min[0], max[1], near_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], max[1], near_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], min[1], near_z) ) }};
	boost::array<vl::Vec3f, 4> farPoints = {{
		vl::xform( invMatrix, vl::Vec3f(min[0], min[1], far_z) ),
		vl::xform( invMatrix, vl::Vec3f(min[0], max[1], far_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], max[1], far_z) ),
		vl::xform( invMatrix, vl::Vec3f(max[0], min[1], far_z) ) }};

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	GLDisableHandler cullFace(GL_CULL_FACE);

	GLEnableHandler blend( GL_BLEND );

	float alpha = 0.3;
	int lineWidth = 2;

	if( selected )
	{
		alpha = 0.6;
		lineWidth = 4;
	}

	vl::Vec3f color = this->color();
	glColor4f( color[0], color[1], color[2], alpha );
	{
		GLActionHandler quads(GL_QUADS);

		glVertex3fv( closePoints[0].Ref() );
		glVertex3fv( closePoints[1].Ref() );
		glVertex3fv( closePoints[2].Ref() );
		glVertex3fv( closePoints[3].Ref() );

		for( int i = 0; i < 4; ++i )
		{
			glVertex3fv( closePoints[(i+0)%4].Ref() );
			glVertex3fv( closePoints[(i+1)%4].Ref() );
			glVertex3fv( farPoints[(i+1)%4].Ref() );
			glVertex3fv( farPoints[(i+0)%4].Ref() );
		}
	}

	glLineWidth( lineWidth );
	glColor4f( color[0], color[1], color[2], 2.0*alpha );
	{
		GLActionHandler lineLoop(GL_LINE_LOOP);

		for( size_t i = 0; i < 4; ++i )
			glVertex3fv( closePoints[i].Ref() );
	}
	{
		GLActionHandler lines( GL_LINES );
		for( int i = 0; i < 4; ++i )
		{
			glVertex3fv( closePoints[i].Ref() );
			glVertex3fv( farPoints[i].Ref() );
		}
	}
	glLineWidth( 1.0 );
}

PositiveProjectiveConstraint::PositiveProjectiveConstraint( 
	const std::vector<size_t>& objects, 
	const vl::Mat4f& matrix, 
	const BoundingBox2f& box )
	: ProjectiveConstraint( objects, matrix, box )
{
}

PositiveProjectiveConstraint::~PositiveProjectiveConstraint()
{
}

vl::Vec3f PositiveProjectiveConstraint::color() const
{
	return vl::Vec3f( 40.0/255.0, 255.0/255.0, 40.0/255.0 );
}

bool PositiveProjectiveConstraint::satisfies( const Path& path ) const
{
	return this->contained(path);
}

NegativeProjectiveConstraint::NegativeProjectiveConstraint( 
	const std::vector<size_t>& objects, 
	const vl::Mat4f& matrix, 
	const BoundingBox2f& box )
	: ProjectiveConstraint( objects, matrix, box )
{
}

NegativeProjectiveConstraint::~NegativeProjectiveConstraint()
{
}

vl::Vec3f NegativeProjectiveConstraint::color() const
{
	return vl::Vec3f( 255.0/255.0, 40.0/255.0, 40.0/255.0 );
}

bool NegativeProjectiveConstraint::satisfies( const Path& path ) const
{
	return !this->contained(path);
}

} // namespace planning


