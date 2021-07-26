#ifndef __MWSAMPLE_H__
#define __MWSAMPLE_H__

#include "physicsFwd.h"
#include "mwdriver.h"

class WrappedCompressedPath;

namespace planning {

class Task_Sample;
class Driver_Sample;

class Driver_Sample
	: public MWDriver
{
public:
	Driver_Sample( SimulationTreePtr tree );
	~Driver_Sample();

	MWReturn get_userinfo( int argc, char *argv[] );

    /// Set up an array of tasks here
    MWReturn setup_initial_tasks( int *, MWTask *** );

    /// What to do when a task finishes:
    MWReturn act_on_completed_task( MWTask * );

    /// Put things in the send buffer here that go to a worker
    MWReturn pack_worker_init_data( void );

	SimulationTreePtr tree() const;
	MWTask* gimme_a_task();

private:
	SimulationTreePtr tree_;
};

class Worker_Sample
	: public MWWorker
{
public:
	/// Default Constructor
	Worker_Sample();

	/// Default Destructor
	~Worker_Sample();

    MWReturn unpack_init_data( void );
    void execute_task( MWTask * );
	MWTask* gimme_a_task();

private:
	boost::shared_ptr<Scene> scene_;
	boost::shared_ptr<SimulationTree> tree_;
	RandomStream random_;
};

class Task_Sample
	: public MWTask
{
public:
	Task_Sample();
	Task_Sample(boost::shared_ptr<WrappedCompressedPath> path);

	~Task_Sample();

    /// Pack the work for this task into the PVM buffer
    void pack_work( void );
    
    /// Unpack the work for this task from the PVM buffer
    void unpack_work( void );
    
    /// Pack the results from this task into the PVM buffer
    void pack_results( void );
    
    /// Unpack the results from this task into the PVM buffer
    void unpack_results( void );

	void setPath( boost::shared_ptr<WrappedCompressedPath> path );
	boost::shared_ptr< WrappedCompressedPath > getPath() const;

private:
    boost::shared_ptr< WrappedCompressedPath > path_;
};

}

#endif

