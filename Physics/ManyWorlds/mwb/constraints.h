#ifndef __CONSTRAINT_H__

#define __CONSTRAINT_H__

#include "physicsFwd.h"
#include "scene.h"
#include "simulationTree.h"
#include "sphere_tree.h"

#include <boost/utility.hpp>

namespace planning {

class SimulationNode;

class Objective
	: public Named
{
public:
	Objective(const std::string& name);
	virtual ~Objective();

	// An Objective is defined to give a single
	//   value for every object in the scene; this will help the local
	//   optimization significantly
	virtual float evaluate( 
		const PiecewisePath& path,
		size_t frameRate,
		ConstPhysicsObjectPtr object,
		const CameraSource& source ) const = 0;
	
	virtual float combine( float oldValue, float newValue ) const = 0;
	virtual float startValue() const = 0;

protected:
	vl::Vec3f viewDirection(const vl::Vec3f& position, const CameraSource& source, float time) const;
};

// Corresponds to an objective that is integrated across time
class NormalizedIntegratedObjective
	: public Objective
{
public:
	NormalizedIntegratedObjective(const std::string& name);
	virtual ~NormalizedIntegratedObjective();

	float combine( float oldValue, float newValue ) const;
	float startValue() const;

protected:
	/*
	virtual std::vector<float> evaluateState( 
		const Simulation::State& state,
		const std::vector<ConstPhysicsObjectPtr>& objects,
		const CameraSource& cam ) const = 0;
		*/

private:
};

class AngularEnergyObjective
	: public NormalizedIntegratedObjective
{
public:
	AngularEnergyObjective();
	float evaluate( 
		const PiecewisePath& path,
		size_t frameRate,
		ConstPhysicsObjectPtr object,
		const CameraSource& source ) const;
};

class FinalOrientationObjective
	: public Objective
{
public:
	FinalOrientationObjective();
	float evaluate( 
		const PiecewisePath& path,
		size_t frameRate,
		ConstPhysicsObjectPtr object,
		const CameraSource& source ) const;
	
	float combine( float oldValue, float newValue ) const;
	float startValue() const;
};

/*
class ViewDependentAngularEnergyObjective
	: public NormalizedIntegratedObjective
{
public:
	ViewDependentAngularEnergyObjective();

	std::vector<float> evaluate( 
		const std::vector<PiecewisePath> paths,
		const std::vector<ConstPhysicsObjectPtr>& objects,
		const CameraSource& source ) const;
};

class RunTimeObjective
	: public Objective
{
public:
	RunTimeObjective();

	std::vector<float> evaluate( 
		const std::vector<PiecewisePath> paths,
		const std::vector<ConstPhysicsObjectPtr>& objects,
		const CameraSource& source ) const;
};

class CountCollisionsObjective
	: public Objective
{
public:
	CountCollisionsObjective();

	std::vector<float> evaluate( 
		const std::vector<PiecewisePath> paths,
		const std::vector<ConstPhysicsObjectPtr>& objects,
		const CameraSource& source ) const;
};


class StackingObjective
	: public Objective
{
public:
	StackingObjective();

	std::vector<float> evaluate( 
		const std::vector<PiecewisePath> paths,
		const std::vector<ConstPhysicsObjectPtr>& objects,
		const CameraSource& source ) const;
};

class YPositionObjective
	: public Objective
{
public:
	YPositionObjective();

	std::vector<float> evaluate( 
		const std::vector<PiecewisePath> paths,
		const std::vector<ConstPhysicsObjectPtr>& objects,
		const CameraSource& source ) const;
};
*/

class Constraint
{
public:
	Constraint( const std::vector<size_t>& objects );
	virtual ~Constraint();

	virtual bool satisfies( const Path& path ) const = 0;

	virtual void render(bool selected) const = 0;

	typedef std::vector<size_t>::const_iterator object_iterator;
	object_iterator objects_begin() const;
	object_iterator objects_end() const;

	virtual bool intersect(const Ray3d& ray, double& t ) const = 0;

private:
	std::vector<size_t> objects_;
};

// A projective constraint is the standard one drawn in the viewer: it
//   consists of a box on the screen (maybe at some point we will expand this)
//   and the constaint consists of whatever simulations project into the box
class ProjectiveConstraint
	: public Constraint
{
public:
	ProjectiveConstraint( const std::vector<size_t>& objects, const vl::Mat4f& matrix, const BoundingBox2f& box );
	virtual ~ProjectiveConstraint();

	vl::Mat4f matrix() const;
	BoundingBox2f box() const;

	void render(bool selected) const;
	bool intersect(const Ray3d& ray, double& t ) const;

protected:
	bool contained( const Path& path ) const;
	virtual vl::Vec3f color() const = 0;

private:
	vl::Mat4f matrix_;
	BoundingBox2f box_;
	ProjectedInBoxQuery query_;
};

class PositiveProjectiveConstraint
	: public ProjectiveConstraint
{
public:
	PositiveProjectiveConstraint( const std::vector<size_t>& objects, const vl::Mat4f& matrix, const BoundingBox2f& box );
	~PositiveProjectiveConstraint();

	bool satisfies( const Path& path ) const;

protected:
	vl::Vec3f color() const;

private:
};

class NegativeProjectiveConstraint
	: public ProjectiveConstraint
{
public:
	NegativeProjectiveConstraint( const std::vector<size_t>& objects, const vl::Mat4f& matrix, const BoundingBox2f& box );
	~NegativeProjectiveConstraint();

	bool satisfies( const Path& path ) const;

protected:
	vl::Vec3f color() const;

private:
};

} // namespace planning

#endif

