#ifndef __EXCEPTION_H__

#define __EXCEPTION_H__

#include <string>
#include <iosfwd>
#include <boost/spirit/iterator.hpp>

class Exception
{
public:
	Exception( const std::string& message )
		: message_( message ) {}

	Exception( const char* message )
		: message_( message ) {}

	virtual ~Exception() {}

	std::string message() const { return message_; }

private:
	std::string message_;
};

class CreationException : public Exception
{
public:
	CreationException( const std::string& message )
		: Exception( message ) {}

	CreationException( const char* message )
		: Exception( message ) {}

	virtual ~CreationException() {}
};

class IOException : public Exception
{
public:
	IOException( const std::string& message )
		: Exception( message ) {}

	IOException( const char* message )
		: Exception( message ) {}

	virtual ~IOException() {}
};

class SolverException : public Exception
{
public:
	SolverException( const std::string& message )
		: Exception( message ) {}
};

typedef boost::spirit::position_iterator< boost::spirit::file_iterator<> > iterator_t;

class ParserException
	: public Exception
{
public:
	ParserException(std::string message, const iterator_t& itr);
	ParserException(std::string message, std::string filename, unsigned int line, unsigned int col )
		: Exception(message), filename_(filename), line_(line), col_(col) {}

	friend std::ostream & operator<<( std::ostream &output, const ParserException& a);

	virtual ~ParserException() {}

private:
	std::string filename_;
	unsigned int line_;
	unsigned int col_;
};

#endif
