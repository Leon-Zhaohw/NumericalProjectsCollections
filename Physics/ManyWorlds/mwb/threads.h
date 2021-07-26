// We will need some threading primitives that aren't in boost,
// so this will hopefully make this cross-platform

#ifndef __THREADS_H__
#define __THREADS_H__

#include <algorithm>

#ifndef WIN32
#include <boost/thread/thread.hpp>
#endif

#ifdef WIN32
class ThreadException
	: public std::exception
{
public:
	ThreadException( DWORD error )
		: error_(error) {}

private:
	DWORD error_;
};

class Event
	: boost::noncopyable
{
public:
	Event(bool manualReset);
	~Event();

	HANDLE handle() const;
	void notify_all() const;
	void reset();

private:
	HANDLE handle_;
};
#endif

// based on this article:
// Compound Win32 Synchronization Objects
// Ruediger R. Asche
// http://msdn.microsoft.com/library/default.asp?url=/library/en-us/dndllpro/html/msdn_locktest.asp
class ReaderWriterLock
{
public:
	ReaderWriterLock();
	~ReaderWriterLock();

	class ScopedWriterLock
	{
	public:
		ScopedWriterLock(ReaderWriterLock& lock);
		~ScopedWriterLock();

	private:
#ifdef WIN32
		ReaderWriterLock& lock_;
#else
		boost::mutex::scoped_lock lock_;
#endif
	};

#ifdef WIN32
	class ScopedReaderLock
	{
	public:
		ScopedReaderLock(ReaderWriterLock& lock);
		~ScopedReaderLock();

	private:
		ReaderWriterLock& lock_;
	};
#else
	typedef ScopedWriterLock ScopedReaderLock;
#endif

private:
#ifdef WIN32
	HANDLE hReaderEvent_;
	HANDLE hMutex_;
	HANDLE hWriterMutex_;
	LONG iCounter_;
#else
	// if we aren't in Win32, just use a garden-variety mutex
	boost::mutex mutex_;
#endif
};

class IncrementDecrement
{
public:
	explicit IncrementDecrement( LONG value = 0 );
	LONG increment();
	LONG decrement();
	LONG value() const;

private:
	LONG value_;
};

#endif

