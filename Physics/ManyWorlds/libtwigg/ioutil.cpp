#include "stdafx.h"
#include "twigg/ioutil.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

std::string percentageAsString( size_t current, size_t total )
{
	float percent = 100.0 * 
		(static_cast<float>(current) / static_cast<float>(total));

	std::ostringstream oss;
	oss << std::setiosflags(std::ios::fixed) 
		<< std::setprecision(1) << std::setw(4) 
		<< std::setfill(' ') << percent << "%";
	return oss.str();
}


PercentageUpdater::PercentageUpdater( size_t max, std::string text,
	std::ostream& out )
	:	maxValue_(max), 
		minUpdate_(max/1000), 
		lastUpdate_(0), 
		stringLength_(0),
		out_(out)
{
	out_ << text << " : ";
	update( 0 );
}

PercentageUpdater::PercentageUpdater( size_t max, std::string text )
	:	maxValue_(max), 
		minUpdate_(max/1000), 
		lastUpdate_(0), 
		stringLength_(0),
		out_(std::cout)
{
	out_ << text << " : ";
	update( 0 );
}

	

PercentageUpdater::~PercentageUpdater()
{
	update( "done.  " );
	out_ << std::endl;
}

void PercentageUpdater::update( size_t value )
{
	if( value - lastUpdate_ < minUpdate_ )
		return;

	lastUpdate_ = value;
	update( percentageAsString( value, maxValue_ ) );
}

void PercentageUpdater::update( const std::string& value )
{
	for( size_t i = 0; i < stringLength_; ++i )
		out_ << "\b";

	out_ << value << std::flush;
	stringLength_ = value.size();
}

