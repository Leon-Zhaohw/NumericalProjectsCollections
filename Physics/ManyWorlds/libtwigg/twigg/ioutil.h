#ifndef __IOUTIL_H__
#define __IOUTIL_H__

#include <iosfwd>
#include <string>

std::string percentageAsString( size_t current, size_t total );

template <typename Streamtype>
typename Streamtype::pos_type filelen( Streamtype& ifs )
{
	typename Streamtype::pos_type curPos = ifs.tellg();
	ifs.seekg(0, std::ios::end);
	typename Streamtype::pos_type l = ifs.tellg();
	ifs.seekg(curPos, std::ios::beg);

	return l;
}

class PercentageUpdater
{
public:
	PercentageUpdater( size_t max, std::string text, std::ostream& out );
	PercentageUpdater( size_t max, std::string text );
	~PercentageUpdater();

	void update( size_t value );
	void update( const std::string& value );

private:
	size_t maxValue_;
	size_t minUpdate_;
	size_t lastUpdate_;
	size_t stringLength_;

	std::ostream& out_;
};

#endif

