#include "stdafx.h"
#include "twigg/exception.h"

using namespace std;

ParserException::ParserException( std::string message, const iterator_t& itr )
	: Exception( message )
{
	line_ = itr.get_position().line;
	col_ = itr.get_position().column;
	filename_ = itr.get_position().file;
}

ostream& operator<<( ostream &output, const ParserException& e)
{
	output << "Error in file '" << e.filename_ << "', "
		<< "line " << e.line_ << ", column " << 
		e.col_ << ":" << std::endl << e.message();
	return output;
}


