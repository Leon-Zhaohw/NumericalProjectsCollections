#include "stdafx.h"

#include "simulationTree.h"
#include "mechModel.h"

#ifdef CONDOR_SAMPLING

#include "MW.h"
#include "mwSample.h"

#include "parser.h"

const size_t targetMachines = 20;

namespace planning {

Driver_Sample::Driver_Sample( SimulationTreePtr tree )
	: tree_(tree)
{
}

MWReturn Driver_Sample::get_userinfo( int argc, char *argv[] ) 
{
	RMC->set_num_exec_classes ( 1 );

	RMC->set_num_arch_classes( 1 );
	RMC->set_arch_class_attributes ( 0, 
		"((Arch==\"INTEL\") && (Opsys==\"LINUX\") )" );
	RMC->set_num_executables( 1 );
	char exec[] = "/nfs/hn00/cdtwigg/planning/planning.sample";
	RMC->add_executable( 0, 0, exec, "");

	RMC->set_target_num_workers( targetMachines );
	// no checkpointing
	set_checkpoint_frequency ( 0 );

	return OK;
}

Driver_Sample::~Driver_Sample()
{
}

SimulationTreePtr Driver_Sample::tree() const
{
	return this->tree_;
}

MWTask* Driver_Sample::gimme_a_task()
{
	return new Task_Sample;
}

MWReturn Driver_Sample::setup_initial_tasks(int *n_init, MWTask ***init_tasks)
{
	// create one task for each machine
	*n_init = targetMachines;

    *init_tasks = new MWTask *[ *n_init ];
	for ( size_t taskindex = 0 ; taskindex < *n_init; taskindex++ )
		(*init_tasks)[taskindex] = new Task_Sample;

	return OK;
}

MWReturn Driver_Sample::act_on_completed_task( MWTask *task )
{	
	Task_Sample* taskSample = dynamic_cast<Task_Sample*>( task );
	assert( taskSample != 0 );
	tree_->addPath( taskSample->getPath() );
	// should add sample to tree here

	MWTask* replacementTask = new Task_Sample;
	MWDriver::addTask( replacementTask );

	return OK;
}

MWReturn Driver_Sample::pack_worker_init_data( void )
{
	// would be nice to replace this with a binary eventually:
	std::string scene;
	{
		MechModelObjectPtr mechModel = tree_->scene()->toMechModel();
		std::ostringstream oss;
		mechModel->dump( oss, std::string() );
		scene.swap( oss.str() );
	}

	unsigned int frameRate = tree_->samplingRate();
	RMC->pack( &frameRate, 1 );

	unsigned int maxImpulseTime = tree_->maxImpulseFrames();
	RMC->pack( &maxImpulseTime, 1 );

	unsigned int length = scene.size();
	RMC->pack( &length, 1 );
	RMC->pack( &scene[0], length );

	return OK;
}

Worker_Sample::Worker_Sample()
{
	workingTask = new Task_Sample;
}

Worker_Sample::~Worker_Sample()
{
	delete workingTask;
	workingTask = 0;
}

MWReturn Worker_Sample::unpack_init_data( void )
{
	unsigned int frameRate;
	RMC->unpack( &frameRate, 1 );

	unsigned int maxImpulseTime;
	RMC->unpack( &maxImpulseTime, 1 );

	unsigned int length;
	RMC->unpack( &length, 1 );
	
	std::string scene( length, ' ' );
	RMC->unpack( &scene[0], length );

	// total hack
	// @todo fix this so it doesn't need to write to file
#ifdef WIN32
	std::string tmpFilename( "c:\\temp\\tmpScene.txt" );
#else
	std::string tmpFilename( "/tmp/tmpScene.txt" );
#endif
	{
		std::ofstream tempFile( tmpFilename.c_str() );
		assert( tempFile );
		tempFile << scene;
	}
	this->scene_ = parseSimulation( tmpFilename );
	this->tree_.reset( new SimulationTree( 
		"condorSampling", scene_, frameRate, maxImpulseTime ) );

	return OK;
}

void Worker_Sample::execute_task( MWTask* task )
{
	assert( this->tree_ );

	boost::shared_ptr< SimulationTree::SimulationList::Acquire > simPtr( tree_->acquireSimulation() );
	SimulationTree::SimulationList::Acquire& simulation = *simPtr;

	std::deque<Simulation::State> states = 
		this->tree_->expand(
			*simulation(),
			tree_->root(),
			std::vector<NxVec3>(),
			std::deque<ActiveImpulse>(),
			false,
			true,
			boost::numeric::bounds<size_t>::highest(),
			0,
			random_,
			0,
			false );
	std::vector<CompressedPath> compressed = compress( states, random_ );

	Task_Sample* taskSample = dynamic_cast<Task_Sample*>( task );
	assert( taskSample != 0 );
	boost::shared_ptr<WrappedCompressedPath> path(
		new WrappedCompressedPath( compressed, states, tree_->samplingRate(), ScoreArray() ) );
	taskSample->setPath( path );
}

MWTask* Worker_Sample::gimme_a_task()
{
	return new Task_Sample;
}

Task_Sample::Task_Sample(boost::shared_ptr<WrappedCompressedPath> path)
	: path_(path)
{
}

Task_Sample::Task_Sample()
{
}

Task_Sample::~Task_Sample()
{
}

void Task_Sample::setPath( boost::shared_ptr<WrappedCompressedPath> path )
{
	this->path_ = path;
}

boost::shared_ptr< WrappedCompressedPath > Task_Sample::getPath() const
{
	return this->path_;
}

void Task_Sample::pack_work( void )
{
	// all the useful info got passed along in pack_worker_init_data
	// so we don't need to do anything here
}

void Task_Sample::unpack_work( void )
{
	// all the useful info got passed along in pack_worker_init_data
	// so we don't need to do anything here
}

void Task_Sample::pack_results( void )
{
	assert( path_ );
	path_->write( RMC );
}

void Task_Sample::unpack_results( void )
{
	path_.reset( new WrappedCompressedPath(RMC) );
}

} // namespace planning

MWWorker* gimme_a_worker()
{
	return new planning::Worker_Sample;
}

#endif

