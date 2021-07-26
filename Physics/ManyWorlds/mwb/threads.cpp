#include "stdafx.h"
#include "threads.h"

Event::Event(bool manualReset)
{
	this->handle_ = CreateEvent(
		NULL,                        // not shared with child processes
		manualReset ? TRUE : FALSE,
		FALSE,                       // not signalled initially
		NULL                         // no name
		);

	if( this->handle_ == NULL )
		throw ThreadException( GetLastError() );
}

Event::~Event()
{
	CloseHandle( this->handle_ );
}

HANDLE Event::handle() const
{
	return this->handle_;
}

void Event::notify_all() const
{
	BOOL result = SetEvent( this->handle_ );
	assert( result == TRUE );
}

void Event::reset()
{
	BOOL result = ResetEvent( this->handle_ );
	assert( result != 0 );
}

ReaderWriterLock::ReaderWriterLock()
	: iCounter_(-1)
{
	this->hReaderEvent_ = CreateEvent(NULL,TRUE,FALSE,NULL);
	this->hMutex_ = CreateEvent(NULL,FALSE,TRUE,NULL);
	this->hWriterMutex_ = CreateMutex(NULL,FALSE,NULL);
}

ReaderWriterLock::~ReaderWriterLock()
{
	CloseHandle( this->hReaderEvent_ );
	CloseHandle( this->hMutex_ );
	CloseHandle( this->hWriterMutex_ );
}

ReaderWriterLock::ScopedWriterLock::ScopedWriterLock(ReaderWriterLock& lock)
	: lock_(lock)
{
	// this lock ensures only one writer at a time
	WaitForSingleObject(lock_.hWriterMutex_, INFINITE);

	// take out the mutex on behalf of the writer (me)
	WaitForSingleObject(lock_.hMutex_, INFINITE);
}

ReaderWriterLock::ScopedWriterLock::~ScopedWriterLock()
{
	SetEvent(lock_.hMutex_);
	ReleaseMutex(lock_.hWriterMutex_);
}

ReaderWriterLock::ScopedReaderLock::ScopedReaderLock(ReaderWriterLock& lock)
	: lock_(lock)
{
	// If I'm the first reader, I'll take out the lock on
	// behalf of all readers:
	if( InterlockedIncrement(&lock_.iCounter_) == 0 )
	{
		WaitForSingleObject(lock_.hMutex_, INFINITE);
		SetEvent(lock_.hReaderEvent_);
	}

	WaitForSingleObject(lock_.hReaderEvent_,INFINITE);
}

ReaderWriterLock::ScopedReaderLock::~ScopedReaderLock()
{
	if( InterlockedDecrement(&lock_.iCounter_) < 0 )
	{
		ResetEvent(lock_.hReaderEvent_);
		SetEvent(lock_.hMutex_);
	}
}

IncrementDecrement::IncrementDecrement( LONG value )
	: value_(value)
{
}

LONG IncrementDecrement::increment()
{
	return InterlockedIncrement( &value_ );
}

LONG IncrementDecrement::decrement()
{
	return InterlockedDecrement( &value_ );
}

LONG IncrementDecrement::value() const
{
	return this->value_;
}

