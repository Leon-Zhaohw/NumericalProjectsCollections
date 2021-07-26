#include "stdafx.h"

#include "parser.h"
#include "scene.h"
#include "constraints.h"

#include "twigg/objfile.h"


#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/escape_char.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/iterator.hpp>
#include <boost/spirit/error_handling/exceptions.hpp>

#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/exception.hpp>

#include "loki/AssocVector.h"


#include <fstream>
#include <map>
#include <typeinfo>
#include <sstream>

namespace planning {

using namespace std;
using namespace boost::spirit;

using boost::tie;

// Can't be "SetString" due to conflict w/ Apple's libraries.
// I guess I should have used my own namespace from the beginning :(
template <typename IteratorT>
class SetAString
{
public:
	explicit SetAString( std::string& str ) : str_(str) {}

	void operator()(IteratorT first, IteratorT last) const
	{
		str_.append(first, last);
	}

	~SetAString() {}
private:
	std::string& str_;
};

class Token
{
public:
	Token(iterator_t itr) : location_(itr) {}
	virtual ~Token() {}
	iterator_t location() const { return location_; }
	virtual std::string toString() const = 0;

private:
	iterator_t location_;
};

typedef boost::shared_ptr<Token> TokenPtr;

class EOFToken : public Token
{
public:
	EOFToken( const iterator_t& itr ) : Token(itr) {}
	virtual ~EOFToken() {}
	virtual std::string toString() const { return "EOF"; }
};

class PreprocessorDirectiveToken : public Token
{
public:
	enum Type
	{
		INCLUDE,		// "include"
	};

	explicit PreprocessorDirectiveToken( Type type, const iterator_t& itr )
		: Token(itr), type_( type ) {}
	virtual ~PreprocessorDirectiveToken() {}
	Type type() const { return type_; }

	virtual std::string toString() const
	{
		switch( type_ )
		{
		case INCLUDE: return "include";
		default: assert(false); return "<ERROR>";
		}
	}

private:
	Type type_;
};

class PunctuationToken : public Token
{
public:
	enum PunctuationType
	{
		LEFT_BRACE,			// '{'
		RIGHT_BRACE,		// '}'
		SEMICOLON,			// ';'
		LEFT_PAREN,			// '('
		RIGHT_PAREN,		// ')'
		COMMA,				// ','
		PERIOD,				// '.'
		PERCENT,			// '%'
		PLUS,				// '+'
		MINUS,				// '-'
		TIMES,				// '*'
		DIVIDE,				// '/'
	};
	explicit PunctuationToken( PunctuationType type, const iterator_t& itr )
		: Token(itr), punctuationType_(type) {}
	virtual ~PunctuationToken() {}
	PunctuationType type() const { return punctuationType_; }

	virtual std::string toString() const
	{
		return punctuationTypeToString( punctuationType_ );
	}

	static std::string punctuationTypeToString( PunctuationType type )
	{
		switch( type )
		{
		case LEFT_BRACE: return "{";
		case RIGHT_BRACE: return "}";
		case SEMICOLON: return ";";
		case LEFT_PAREN: return "(";
		case RIGHT_PAREN: return ")";
		case COMMA: return ",";
		case PERIOD: return ".";
		case PERCENT: return "%";
		case PLUS: return "+";
		case MINUS: return "-";
		case TIMES: return "*";
		case DIVIDE: return "/";

		default: assert(false); return "<ERROR>";
		}
	}

private:
	PunctuationType punctuationType_;
};

class ValueToken : public Token
{
public:
	ValueToken(const iterator_t& itr) : Token(itr) {}
	virtual ~ValueToken() {}
};

class StringValueToken : public ValueToken
{
public:
	explicit StringValueToken( const std::string& val, const iterator_t& itr )
		: ValueToken(itr), value_(val) {}
	std::string value() const { return value_; }
	virtual ~StringValueToken() {}

	virtual std::string toString() const
	{
		return value_;
	}

private:
	std::string value_;
};

class NumericValueToken : public ValueToken
{
public:
	NumericValueToken( const iterator_t& itr )
		: ValueToken(itr) {}
	virtual double toReal() const = 0;

	virtual ~NumericValueToken() {}
};

class IntegerValueToken : public NumericValueToken
{
public:
	explicit IntegerValueToken( int val, const iterator_t& itr )
		: NumericValueToken(itr), value_(val)			{}
	bool isIntegral() const		{ return true; }
	double toReal() const		{ return static_cast<double>( value_ ); }
	int value() const			{ return value_; }

	virtual ~IntegerValueToken() {}

	virtual std::string toString() const
	{
		std::ostringstream oss;
		oss << value_;
		return oss.str();
	}

private:
	int value_;
};

class RealValueToken : public NumericValueToken
{
public:
	explicit RealValueToken( double val, const iterator_t& itr )
		: NumericValueToken(itr), value_(val)			{}
	bool isIntegral() const		{ return false; }
	double toReal() const		{ return value(); }
	double value() const		{ return value_; }

	virtual ~RealValueToken() {}

	virtual std::string toString() const
	{
		std::ostringstream oss;
		oss << value_;
		return oss.str();
	}

private:
	double value_;
};

// stolen from the Spirit examples.
//  Here's our comment rule
struct skip_grammar : public grammar<skip_grammar>
{
    template <typename ScannerT>
    struct definition
    {
        definition(skip_grammar const& //self
						)
        {
            skip
				=   *(space_p
				|	comment_p("#")				// shell-script-style comment
                |   comment_p("//")				// C++-style comment
                |   comment_p("/*", "*/")		// C-style comment
                );
        }

        rule<ScannerT> skip;

        rule<ScannerT> const&
        start() const { return skip; }
    };
};

class TokenStream
{
public:
	virtual boost::shared_ptr< Token > peek() = 0;
	virtual boost::shared_ptr< Token > eat() = 0;
	virtual bool finished() = 0;

	boost::filesystem::path fileInPath( const std::string& filename )
	{
		try
		{
			boost::filesystem::path p( filename );
			if( p.is_complete() )
				return p;
			else
				return path() / p;
		}
		catch( boost::filesystem::filesystem_error& e )
		{
			throw ParserException(
				"Error in filename '" + filename + "': " + e.what(), current() );
		}
	}

protected:
	boost::filesystem::path path() const
	{
		return filename().branch_path();
	}

	virtual iterator_t current() const = 0;

	virtual boost::filesystem::path filename() const = 0;
};

// This class will wrap our Spirit parser in such a way that it 
// provides a steady stream of tokens
class FileTokenStream : public TokenStream, boost::noncopyable
{
	void check( const parse_info<iterator_t>& info, std::string errorMsg )
	{
		if( !info.hit )
			throw ParserException(errorMsg, info.stop);
		current_ = info.stop;
	}

public:
	FileTokenStream( const boost::filesystem::path& filename )
		: first_(filename.native_file_string()), filename_(filename)
	{
		if (!first_)
		{
			std::string error = "Could not open input file: ";
			error.append( filename.native_file_string() );
			throw IOException( error );
		}

		last_ = first_.make_end();
		current_ = iterator_t( first_, last_, filename.native_file_string() );
	}
	
	virtual ~FileTokenStream() {}

	boost::shared_ptr< Token > peek()
	{
		if( currentToken_.get() == 0 )
		{
			{
				// First, we clear out all the junk with our skip_grammar
				parse_info<iterator_t> info = parse( current_, end_, skip_ );
				current_ = info.stop;
			}

			if( finished() )
			{
				currentToken_.reset( new EOFToken(current_) );
			}
			else if( *current_ == '"' )
			{
				std::string val;
				parse_info<iterator_t> info = parse( current_, end_, 
					lexeme_d[ chlit<>('\"') >>
						*( strlit<>("\\\"") | anychar_p - chlit<>('\"'))[SetAString<iterator_t>(val)]
						>> chlit<>('\"')] );
				currentToken_.reset( new StringValueToken( val, current_ ) );
				check( info, "Unterminated string constant" );
			}
			else if( isdigit(*current_) || *current_ == '-' )
			{
				int val;

				// parse a number
				// first, try parsing as an integer
				parse_info<iterator_t> info = parse( current_, end_, 
					int_p[ assign(val) ]);
				if( *(info.stop) == '.' || *(info.stop) == 'e' || *(info.stop) == 'E' )
				{
					double val;

					// is a real value
					parse_info<iterator_t> info = parse( current_, end_, 
						real_p[ assign(val) ]);

					currentToken_.reset( new RealValueToken(val, current_) );
					check( info, "Invalid value" );
				}
				else
				{
					currentToken_.reset( new IntegerValueToken( val, current_ ) );
					check( info, "Invalid value" );
				}
			}
			else if( ispunct(*current_) )
			{
				char val;
				parse_info<iterator_t> info = parse( current_, end_, 
					punct_p[assign(val)] );
				
				if( val == '{' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::LEFT_BRACE, current_ ) );
				else if( val == '}' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::RIGHT_BRACE, current_ ) );
				else if( val == ';' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::SEMICOLON, current_ ) );
				else if( val == '(' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::LEFT_PAREN, current_ ) );
				else if( val == ')' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::RIGHT_PAREN, current_ ) );
				else if( val == ',' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::COMMA, current_ ) );
				else if( val == '.' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::PERIOD, current_ ) );
				else if( val == '+' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::PLUS, current_ ) );
				else if( val == '-' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::MINUS, current_ ) );
				else if( val == '*' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::TIMES, current_ ) );
				else if( val == '/' )
					currentToken_.reset( new PunctuationToken( PunctuationToken::DIVIDE, current_ ) );
				else if( val == '%' )
					// preprocessor declaration
					currentToken_.reset( new PunctuationToken( PunctuationToken::PERCENT, current_ ) );
				else
				{
					std::string errMsg( "Invalid punctuation: '" );
					errMsg.push_back(val);
					errMsg.append("'");
					throw ParserException( errMsg, current_ );
				}

				std::string errMsg( "Invalid punctuation: '" );
				errMsg.push_back(val);
				errMsg.append("'");
				check( info, errMsg );
			}
			else if( isalpha(*current_) )
			{
				// parse an identifier
				std::string val;
				parse_info<iterator_t> info = parse( current_, end_, 
					(alpha_p >> *(alnum_p || ch_p('_')))[SetAString<iterator_t>(val)] );

				currentToken_.reset( new StringValueToken( val, current_ ) );

				check( info, "Internal parser error" );	// shouldn't ever have an error
			}
			else
			{
				throw ParserException("Unknown token", current_);
			}
		}

		return currentToken_;
	}

	boost::shared_ptr< Token > eat()
	{
		TokenPtr returnVal = currentToken_;
		currentToken_.reset();
		peek();
		return returnVal;
	}

	bool finished()
	{
		return current_ == end_;
	}

protected:
	boost::filesystem::path filename() const
	{
		return this->filename_;
	}

	virtual iterator_t current() const
	{
		return current_;
	}

	friend class PreprocessorTokenStream;
private:
	boost::shared_ptr< Token > currentToken_;

	file_iterator<> first_;
	file_iterator<> last_;
	boost::filesystem::path filename_;

	iterator_t current_;
	iterator_t end_;

	skip_grammar skip_;
};

template <typename TokenType>
bool checkToken( TokenPtr token )
{
	if( dynamic_cast<TokenType *>( token.get() ) != 0 )
		return true;
	else
		return false;
}

bool checkPunctuation( TokenPtr token, 
					  PunctuationToken::PunctuationType type )
{
	PunctuationToken* punctToken;
	if( !(punctToken = dynamic_cast<PunctuationToken*>(token.get()) ) ||
		punctToken->type() != type )
		return false;
	else
		return true;
}

// Parser primitives
void eatPunctuation( TokenStream& stream, 
					PunctuationToken::PunctuationType type )
{
	TokenPtr token = stream.eat();
	if( !checkPunctuation( token, type ) )
	{
		std::string errMsg( "Expected: '" );
		errMsg.append( PunctuationToken::punctuationTypeToString(type) );
		errMsg.append( "'" );
		throw ParserException( errMsg, token->location() );
	}
}

template <typename T>
T parsePrimitive( TokenStream& stream );

template <>
std::string parsePrimitive<std::string>( TokenStream& stream )
{
	TokenPtr nextToken = stream.eat();
	StringValueToken* token = dynamic_cast<StringValueToken *>( nextToken.get() );
	if( token == 0 )
	{
		IntegerValueToken* intTok = dynamic_cast<IntegerValueToken *>( nextToken.get() );
		if( intTok == 0 )
			throw ParserException( "Expected: string", nextToken->location() );
		
		std::ostringstream oss;
		oss << intTok->value();
		return oss.str();
	}

	return token->value();
}

template <>
int parsePrimitive<int>( TokenStream& stream )
{
	TokenPtr nextToken = stream.eat();
	IntegerValueToken* token = dynamic_cast<IntegerValueToken *>( nextToken.get() );
	if( token == 0 )
		throw ParserException( "Expected: integer", nextToken->location() );

	return token->value();
}

template <>
size_t parsePrimitive<size_t>( TokenStream& stream )
{
	TokenPtr nextToken = stream.eat();
	IntegerValueToken* token = dynamic_cast<IntegerValueToken *>( nextToken.get() );
	if( token == 0 || token->value() < 0 )
		throw ParserException( "Expected: unsigned integer", nextToken->location() );

	return token->value();
}

template <>
float parsePrimitive<float>( TokenStream& stream )
{
	TokenPtr nextToken = stream.eat();
	NumericValueToken* token = dynamic_cast<NumericValueToken *>( nextToken.get() );
	if( token == 0 )
		throw ParserException( "Expected: real", nextToken->location() );

	return token->toReal();
}

template <>
double parsePrimitive<double>( TokenStream& stream )
{
	TokenPtr nextToken = stream.eat();
	NumericValueToken* token = dynamic_cast<NumericValueToken *>( nextToken.get() );
	if( token == 0 )
		throw ParserException( "Expected: real", nextToken->location() );

	return token->toReal();
}

template <>
bool parsePrimitive<bool>( TokenStream& stream )
{
	TokenPtr nextToken = stream.eat();
	StringValueToken* token = dynamic_cast<StringValueToken *>( nextToken.get() );
	if( token == 0 )
		throw ParserException( "Expected: boolean", nextToken->location() );

	std::string value = token->value();
	if( value == "true" )
		return true;
	else if( value == "false" )
		return false;
	else
		throw ParserException( "Expected: 'true' or 'false'", nextToken->location() );
}

template <typename T, size_t N>
boost::array<T, N> parseTuple( TokenStream& stream )
{
	boost::array<T, N> result;
	eatPunctuation( stream, PunctuationToken::LEFT_PAREN );
	for( size_t i = 0; i < N && !stream.finished(); ++i )
	{
		result[i] = parsePrimitive<T>( stream );
		TokenPtr nextToken = stream.peek();
		if( !checkPunctuation( nextToken, PunctuationToken::COMMA ) )
			break;
		eatPunctuation( stream, PunctuationToken::COMMA );
	}
	eatPunctuation( stream, PunctuationToken::RIGHT_PAREN );
	return result;
}

template <>
vl::Vec2d parsePrimitive<vl::Vec2d>( TokenStream& stream )
{
	TokenPtr currentToken = stream.peek();
	boost::array<double, 2> values = parseTuple<double, 2>( stream );
	return vl::Vec2d( values[0], values[1] );
}

template <>
vl::Vec3d parsePrimitive<vl::Vec3d>( TokenStream& stream )
{
	TokenPtr currentToken = stream.peek();
	boost::array<double, 3> values = parseTuple<double, 3>( stream );
	return vl::Vec3d( values[0], values[1], values[2] );
}

template <>
vl::Vec4d parsePrimitive<vl::Vec4d>( TokenStream& stream )
{
	TokenPtr currentToken = stream.peek();
	boost::array<double, 4> values = parseTuple<double, 4>( stream );
	return vl::Vec4d( values[0], values[1], values[2], values[3] );
}


template <typename T>
std::vector<T> parseVector( TokenStream& stream )
{
	std::vector<T> result;
	eatPunctuation( stream, PunctuationToken::LEFT_PAREN );
	while( !stream.finished() )
	{
		result.push_back( parsePrimitive<T>( stream ) );
		TokenPtr nextToken = stream.peek();
		if( !checkPunctuation( nextToken, PunctuationToken::COMMA ) )
			break;
		eatPunctuation( stream, PunctuationToken::COMMA );
	}
	eatPunctuation( stream, PunctuationToken::RIGHT_PAREN );

	return result;
}

template <>
std::vector<vl::Vec3d> parsePrimitive< std::vector<vl::Vec3d> >( TokenStream& stream )
{
	return parseVector<vl::Vec3d>( stream );
}

template <>
TinyVec<unsigned int, 3> parsePrimitive< TinyVec<unsigned int, 3> >( TokenStream& stream )
{
	TokenPtr currentToken = stream.peek();
	boost::array<int, 3> tuple = parseTuple< int, 3 >( stream );

	for( unsigned int i = 0; i < 3; ++i )
	{
		if( tuple[i] < 0 )
		{
			std::ostringstream oss;
			oss << "Expected: positive values; found " << tuple[i];
			throw ParserException( oss.str(), currentToken->location() );
		}
	}

	return TinyVec<unsigned int, 3>( tuple[0], tuple[1], tuple[2] );
}

template<>
RandomizedProperty parsePrimitive< RandomizedProperty >( TokenStream& stream )
{
	RandomizedProperty result;
	result.name = parsePrimitive<std::string>( stream );
	result.low = parsePrimitive<float>( stream );
	result.high = parsePrimitive<float>( stream );
	return result;
}

/*
template <typename S, typename T>
std::pair<S, T> parsePrimitive< std::pair<S, T> >( TokenStream& stream )
{
	typedef std::pair<S, T> PrType;
	PrType result;
	result.first = parsePrimitive<S>( stream );
	eatPunctuation( stream, PunctuationToken::PERIOD );
	result.second = parsePrimitive<T>( stream );
	return result;
}
*/

std::string parseString( TokenStream& stream )	{ return parsePrimitive<std::string>(stream); }
double parseReal( TokenStream& stream )			{ return parsePrimitive<double>(stream); }
int parseInt( TokenStream& stream )				{ return parsePrimitive<int>(stream); }
bool parseBoolean( TokenStream& stream )		{ return parsePrimitive<bool>(stream); }
vl::Vec2d parseVec2( TokenStream& stream )		{ return parsePrimitive<vl::Vec2d>(stream); }
vl::Vec3d parseVec3( TokenStream& stream )		{ return parsePrimitive<vl::Vec3d>(stream); }
vl::Vec4d parseVec4( TokenStream& stream )		{ return parsePrimitive<vl::Vec4d>(stream); }



// Wraps the token stream and acts as a preprocessor, handling e.g. "#include"s
class PreprocessorTokenStream : public TokenStream
{
public:
	PreprocessorTokenStream( const boost::filesystem::path& filename )
	{
		includeStack_.push_front( 
			boost::shared_ptr<FileTokenStream>(
				new FileTokenStream( filename ) ) );
	}

	virtual ~PreprocessorTokenStream() {}

	virtual boost::shared_ptr< Token > peek()
	{
		includeStack_.front()->peek();
		if( includeStack_.front()->finished() && includeStack_.size() > 1 )
			includeStack_.pop_front();

		TokenPtr token = includeStack_.front()->peek();
		PunctuationToken* puncToken = dynamic_cast<PunctuationToken*>( token.get() );
		if( puncToken != 0 )
		{
			if( puncToken->type() == PunctuationToken::PERCENT )
			{
				// Preprocessor directive
				includeStack_.front()->eat();

				token = includeStack_.front()->eat();
				PreprocessorDirectiveToken* prePToken = dynamic_cast<PreprocessorDirectiveToken*>(token.get());
				if( prePToken != 0 && prePToken->type() == PreprocessorDirectiveToken::INCLUDE )
				{
					try
					{
						std::string filename = parseString( *includeStack_.front() );
						includeStack_.push_front(
							boost::shared_ptr<FileTokenStream>(
								new FileTokenStream( fileInPath( filename )) ) );
						token = includeStack_.front()->peek();
					}
					catch( IOException& e )
					{
						throw ParserException( "Error parsing #include directive: "
							+ e.message(), token->location() );
					}
				}
				else
					throw ParserException( "Invalid preprocessor directive", token->location() );
			}
		}

		return token;
	}

	virtual boost::shared_ptr< Token > eat()
	{
		peek();
		return includeStack_.front()->eat();
	}

	virtual bool finished()
	{
		peek();
		return includeStack_.front()->finished();
	}

protected:
	boost::filesystem::path filename() const
	{
		return includeStack_.front()->filename();
	}

	iterator_t current() const
	{
		return includeStack_.front()->current();
	}

private:
	std::deque< boost::shared_ptr<FileTokenStream> > includeStack_;
};


template <typename T>
class RequiredAttribute
{
public:
	RequiredAttribute(std::string name) : name_(name), valid_(false) {}
	void set( T val ) { value_ = val; valid_ = true; }

	T get(TokenPtr token) const
	{
		if( !valid_ )
		{
			ostringstream errMsg;
			errMsg << "Expected: " << name_;
			throw ParserException(errMsg.str(), token->location() );
		}

		return value_;
	}

	bool isSet() const		{ return valid_; }

private:
	std::string name_;
	bool valid_;
	T value_;
};

class ParseableAttribute : public Named
{
public:
	ParseableAttribute( const std::string& name )
		: Named(name) {}

	void parse( TokenStream& stream )
	{
		location_ = stream.peek();
		parseLocal( stream );
	}

	TokenPtr token() const
	{
		assert( location_ );
		return location_;
	}

	virtual ~ParseableAttribute() {}

protected:
	virtual void parseLocal( TokenStream& stream ) = 0;

private:
	TokenPtr location_;
};

template <typename T>
class RequiredParseableAttribute
	: public ParseableAttribute
{
	RequiredAttribute<T> attribute_;

public:
	RequiredParseableAttribute( const std::string& name )
		: ParseableAttribute(name), attribute_(name) {}

	void parseLocal( TokenStream& stream )
	{
		attribute_.set( parsePrimitive<T>(stream) );
	}

	T get( TokenPtr t ) const
	{
		return attribute_.get( t );
	}
};

template <typename T>
class OptionalParseableAttribute
	: public ParseableAttribute
{
	T attribute_;

public:
	OptionalParseableAttribute( const std::string& name, const T defaultValue )
		: ParseableAttribute(name), attribute_(defaultValue) {}

	void parseLocal( TokenStream& stream )
	{
		attribute_ = parsePrimitive<T>(stream);
	}

	T get() const
	{
		return attribute_;
	}
};

template <typename T>
class MultivaluedParseableAttribute
	: public ParseableAttribute
{
	std::deque<T> values_;

public:
	MultivaluedParseableAttribute( const std::string& name )
		: ParseableAttribute(name) {}

	void parseLocal( TokenStream& stream )
	{
		values_.push_back( parsePrimitive<T>(stream) );
	}

	std::deque<T> get() const
	{
		return values_;
	}
};

template <typename T>
class ParseMultipleAttributes
	: public ParseableAttribute
{
	std::deque<T> attributes;
	std::deque<TokenPtr> tokens;

public:
	ParseMultipleAttributes( const std::string& name )
		: ParseableAttribute(name) {}

	void parseLocal( TokenStream& stream )
	{
		tokens.push_back( stream.peek() );
		attributes.push_back( parsePrimitive<T>(stream) );
	}

	size_t numAttributes() const 
	{
		assert( attributes.size() == tokens.size() );
		return attributes.size();
	}

	TokenPtr getToken( size_t i )
	{
		return tokens.at(i);
	}

	T getValue( size_t i )
	{
		return attributes.at(i);
	}
};

typedef std::deque< ParseableAttribute* > AttributeList;

std::string quote( const std::string& str )
{
	std::string result;
	result.push_back( '\'' );
	result.append( str );
	result.push_back( '\'' );
	return result;
}

std::string expectedList( const AttributeList& attributes )
{
	assert( !attributes.empty() );
	std::string result;
	if( attributes.size() == 1 )
	{
		result.append( "expected: " );
		result.append( quote(attributes.front()->name()) );
	}
	else if( attributes.size() == 2 )
	{
		std::string result;
		result.append( "expected: " );
		result.append( quote(attributes.front()->name()) );
		result.append( " or " );
		result.append( quote(attributes.back()->name()) );
	}
	else
	{
		std::string result;
		result.append( "expected one of: " );
		for( AttributeList::const_iterator iter = attributes.begin(); 
			iter != attributes.end(); ++iter )
		{
			if( iter != attributes.begin() )
				result.append( ", " );
			result.append( quote((*iter)->name()) );
		}
	}

	return result;
}

class DefaultUnknownAttributeHandler
{
public:
	bool operator()( const std::string& attributeName, TokenStream& stream )
	{
		return false;
	}
};

template <typename UnknownAttributeHandler>
std::pair<TokenPtr, TokenPtr> parseUnnamedObject( 
	TokenStream& stream, 
	AttributeList attributes, 
	UnknownAttributeHandler& handler = UnknownAttributeHandler() )
{
	typedef NamedPtrComparator ComparatorType;

	std::sort( attributes.begin(), attributes.end(), ComparatorType() );

	TokenPtr start = stream.peek();
	eatPunctuation( stream, PunctuationToken::LEFT_BRACE );

	while( !stream.finished() )
	{
		TokenPtr nextToken = stream.peek();
		if( dynamic_cast<StringValueToken*>( nextToken.get() ) )
		{
			std::string attributeName = parseString( stream );

			Named named( attributeName );
			typedef AttributeList::const_iterator AttributeIter;
			typedef std::pair<AttributeIter, AttributeIter> AttributeIterPr;
			AttributeIterPr iterPr = std::equal_range( attributes.begin(), 
				attributes.end(), &named, ComparatorType() );

			// shouldn't be more than one:
			assert( std::distance( iterPr.first, iterPr.second ) <= 1 );

			if( iterPr.first == iterPr.second )		// not found
			{
				if( handler(attributeName, stream) )
					continue;

				std::ostringstream oss;
				oss << "Unknown attribute '" << attributeName << "'; " << expectedList(attributes);
				throw ParserException( oss.str(), nextToken->location() );
			}

			(*iterPr.first)->parse( stream );

			eatPunctuation( stream, PunctuationToken::SEMICOLON );
		}
		else
			break;
	}

	TokenPtr end = stream.peek();
	eatPunctuation( stream, PunctuationToken::RIGHT_BRACE );

	return std::make_pair( start, end );
}

template <typename UnknownAttributeHandler>
std::pair< std::string, std::pair<TokenPtr, TokenPtr> > parseNamedObject( 
	TokenStream& stream, 
	const AttributeList& attributes,
	UnknownAttributeHandler& handler )
{
	std::string name = parseString( stream );
	return std::make_pair( name, parseUnnamedObject( stream, attributes, handler ) );
}

MaterialPtr parseMaterial( TokenStream& stream, Scene& scene )
{
	AttributeList attributes;
	RequiredParseableAttribute<float> density( "density" );                  attributes.push_back( &density );
	RequiredParseableAttribute<float> dynamicFriction( "dynamicFriction" );  attributes.push_back( &dynamicFriction );
	RequiredParseableAttribute<float> staticFriction( "staticFriction" );    attributes.push_back( &staticFriction );
	RequiredParseableAttribute<float> restitution( "restitution" );          attributes.push_back( &restitution );
	OptionalParseableAttribute<bool> coulombFriction( "coulomb", true );     attributes.push_back( &coulombFriction );

	OptionalParseableAttribute<float> raleighAlpha( "raleighAlpha", 0.001 ); attributes.push_back( &raleighAlpha );
	OptionalParseableAttribute<float> raleighBeta( "raleighBeta", 0.001 );   attributes.push_back( &raleighBeta );
	OptionalParseableAttribute<float> youngsModulus( "youngsModulus", 1.9E11 );
	                                                                         attributes.push_back( &youngsModulus );
	OptionalParseableAttribute<float> poissonRatio( "poissonRatio", 0.3 );   attributes.push_back( &poissonRatio );

	OptionalParseableAttribute<vl::Vec3d> color( "color", vl::Vec3d(0.8, 0.8, 0.8) );
	                                                                         attributes.push_back( &color );

	std::pair<TokenPtr, TokenPtr> loc;
	std::string name;
	DefaultUnknownAttributeHandler handler;
	tie( name, loc ) = parseNamedObject( stream, attributes, handler );

	MaterialPtr mat(
		new Material( name, 
			density.get(loc.second), 
			dynamicFriction.get(loc.second), 
			staticFriction.get(loc.second), 
			restitution.get(loc.second) ) );
	mat->setColor( toVec3f(color.get()) );
	mat->setCoulombFriction( coulombFriction.get() );
	mat->setRaleighAlpha( raleighAlpha.get() );
	mat->setRaleighBeta( raleighBeta.get() );
	mat->setYoungsModulus( youngsModulus.get() );
	mat->setPoissonRatio( poissonRatio.get() );
	return mat;
}

struct KeysWithLocation
{
	std::string name;
	TokenPtr location;

	typedef Keyframe<float, float> KeyframeType; 
	std::deque<KeyframeType> keys;
};


class KeysHandler
{
public:
	bool operator()( const std::string& attributeName, TokenStream& stream )
	{
		if( attributeName != "keys" )
			return false;

		KeysWithLocation keys;
		keys.location = stream.peek();
		keys.name = parseString(stream);
		eatPunctuation(stream, PunctuationToken::LEFT_BRACE);

		while( !stream.finished() )
		{
			TokenPtr nextToken = stream.peek();
			if( dynamic_cast<StringValueToken*>( nextToken.get() ) )
			{
				std::string keysString = parseString( stream );
				if( keysString != "key" )
					throw ParserException( "Expected: 'key'", nextToken->location() );

				vl::Vec2d value = parseVec2(stream);
				keys.keys.push_back( KeysWithLocation::KeyframeType(value[0], value[1]) );

				eatPunctuation(stream, PunctuationToken::SEMICOLON);
			}
			else
				break;
		}

		keys_.push_back( keys );
		eatPunctuation(stream, PunctuationToken::RIGHT_BRACE);
		return true;
	}

	const std::deque<KeysWithLocation>& keys() const
	{
		return keys_;
	}

private:
	std::deque<KeysWithLocation> keys_;
};


CameraWrapperPtr parseCamera( TokenStream& stream, Scene& scene )
{
	AttributeList attributes;
	// attributes for all cameras:
	RequiredParseableAttribute<std::string> type( "type" );                  attributes.push_back( &type );
	RequiredParseableAttribute<float> near_z( "near_z" );                    attributes.push_back( &near_z );
	RequiredParseableAttribute<float> far_z( "far_z" );                      attributes.push_back( &far_z );
	OptionalParseableAttribute<vl::Vec3d> up( "up", vl::Vec3d(0.0, 1.0, 0.0) );
	                                                                         attributes.push_back( &up );

	// attributes for perspective camera
	RequiredParseableAttribute<float> elevation( "elevation" );              attributes.push_back( &elevation );
	RequiredParseableAttribute<float> azimuth( "azimuth" );                  attributes.push_back( &azimuth );
	RequiredParseableAttribute<float> dolly( "dolly" );                      attributes.push_back( &dolly );
	RequiredParseableAttribute<float> twist( "twist" );                      attributes.push_back( &twist );
	RequiredParseableAttribute<vl::Vec3d> look( "look" );                    attributes.push_back( &look );
	RequiredParseableAttribute<float> fov( "fov" );                          attributes.push_back( &fov );

	// attributes for orthographic camera
	RequiredParseableAttribute<float> port( "port" );                        attributes.push_back( &port );
	RequiredParseableAttribute<float> rotation( "rotation" );                attributes.push_back( &rotation );
	RequiredParseableAttribute<vl::Vec2d> center( "center" );                attributes.push_back( &center );

	KeysHandler handler;

	std::pair<TokenPtr, TokenPtr> loc;
	std::string name;

	tie( name, loc ) = parseNamedObject( stream, attributes, handler );

	boost::shared_ptr<Camera> cam;
	if( type.get(loc.second) == "perspective" )
	{
		boost::shared_ptr<ProjectiveFixedYCamera> projCam(
			new ProjectiveFixedYCamera );
		cam = projCam;

		projCam->setElevation( elevation.get(loc.second) );
		projCam->setAzimuth( azimuth.get(loc.second) );
		projCam->setDolly( dolly.get(loc.second) );
		projCam->setTwist( twist.get(loc.second) );
		projCam->setLookAt( toVec3f(look.get(loc.second)) );
		projCam->setFOV( fov.get(loc.second)  );
	}
	else if( type.get(loc.second) == "orthographic" )
	{
		boost::shared_ptr<OrthoCamera> orthoCam( 
			new OrthoCamera );
		cam = orthoCam;

		orthoCam->setPort( port.get(loc.second) );
		orthoCam->setAzimRot( rotation.get(loc.second) );
		orthoCam->setCenter( toVec2f(center.get(loc.second)) );
	}
	else
	{
		throw ParserException( "Expected: 'perspective' or 'orthographic'", type.token()->location() );
	}

	cam->setZRange( near_z.get(loc.second), far_z.get(loc.second) );
	cam->setUpVector( toVec3f(up.get()) );

	CameraWrapperPtr camera( new CameraWrapper(name, cam) );

	// now, restore all the keyframes
	{
		const std::deque<KeysWithLocation>& keys = handler.keys();
		for( std::deque<KeysWithLocation>::const_iterator fileKeyableIter = keys.begin();
			fileKeyableIter != keys.end(); ++fileKeyableIter )
		{
			boost::shared_ptr<AttributeWithKeys> attribute;

			for( CameraWrapper::attribute_iterator attIter = camera->keyable_begin();
				attIter != camera->keyable_end(); ++attIter )
			{
				if( (*attIter)->keyable().name() == fileKeyableIter->name )
					attribute = *attIter;
			}

			if( !attribute )
			{
				throw ParserException( "Unknown keyable attribute: '" + fileKeyableIter->name + "'", 
					fileKeyableIter->location->location() );
			}

			// found a match, need to add in the appropriate keys
			for( std::deque<KeysWithLocation::KeyframeType>::const_iterator keyItr = fileKeyableIter->keys.begin();
				keyItr != fileKeyableIter->keys.end(); ++keyItr )
			{
				attribute->addKey( *keyItr );
			}
		}
	}

	return camera;
}

class ObjectParser;

class ObjectParserCreator
{
public:
	virtual ~ObjectParserCreator()
	{
	}

	virtual boost::shared_ptr<ObjectParser> parser() const = 0;
};

typedef Loki::AssocVector<std::string, boost::shared_ptr<ObjectParserCreator> > ParseObjectMap;


class ObjectParser
{
public:
	ObjectParser()
		:	randomizedProperties_( "randomize" ), 
			visible_( "visible", true ), 
			selfCollisions_( "selfCollisions", true )
	{
	}

	virtual ~ObjectParser()
	{
	}

	void parse( 
		TokenStream& stream, 
		Scene& scene, 
		SceneGraphElementPtr parent, 
		const ParseObjectMap& map )
	{
		TokenPtr nameToken = stream.peek();
		std::string name = parseString( stream );
		if( scene.object( name ) )
		{
			std::ostringstream oss;
			oss << "Object with name '" << name << "' already exists.";
			throw ParserException( oss.str(), nameToken->location() );
		}

		try
		{
			SceneGraphElementPtr result( this->createObject( name, scene ) );

			typedef NamedPtrComparator ComparatorType;

			AttributeList attributes = this->attributes(result);
			attributes.push_back( &randomizedProperties_ );
			attributes.push_back( &visible_ );
			attributes.push_back( &selfCollisions_ );
			std::sort( attributes.begin(), attributes.end(), ComparatorType() );

			TokenPtr start = stream.peek();
			eatPunctuation( stream, PunctuationToken::LEFT_BRACE );

			while( !stream.finished() )
			{
				TokenPtr nextToken = stream.peek();
				if( dynamic_cast<StringValueToken*>( nextToken.get() ) )
				{
					std::string attributeName = parseString( stream );

					Named named( attributeName );
					typedef AttributeList::const_iterator AttributeIter;
					typedef std::pair<AttributeIter, AttributeIter> AttributeIterPr;
					AttributeIterPr iterPr = std::equal_range( attributes.begin(), 
						attributes.end(), &named, ComparatorType() );

					// shouldn't be more than one:
					assert( std::distance( iterPr.first, iterPr.second ) <= 1 );

					if( iterPr.first == iterPr.second )		// not found
					{
						ParseObjectMap::const_iterator objectMapItr = 
							map.find(attributeName);
						if( objectMapItr == map.end() )
						{
							std::ostringstream oss;
							oss << "Unknown attribute '" << attributeName << "'; " << expectedList(attributes);
							throw ParserException( oss.str(), nextToken->location() );
						}
						else
						{
							boost::shared_ptr<ObjectParser> parser = 
								objectMapItr->second->parser();
							parser->parse( stream, scene, result, map );
						}
					}
					else
					{
						(*iterPr.first)->parse( stream );

						eatPunctuation( stream, PunctuationToken::SEMICOLON );
					}
				}
				else
					break;
			}

			TokenPtr end = stream.peek();
			eatPunctuation( stream, PunctuationToken::RIGHT_BRACE );

			this->applyAttributes( result, scene, end );
			parent->addChild( result );
			result->setParent( parent );

			result->setVisible( visible_.get() );
			result->setSelfCollisions( selfCollisions_.get() );

			std::vector<std::string> allAttributes = result->floatAttributes();
			for( size_t i = 0; i < randomizedProperties_.numAttributes(); ++i )
			{
				RandomizedProperty att = randomizedProperties_.getValue(i);
				if( att.high < att.low )
					throw ParserException( "Low value and high value inverted.", 
						randomizedProperties_.getToken(i)->location() );

				if( std::find( allAttributes.begin(), allAttributes.end(), 
					att.name ) == allAttributes.end() )
				{
					throw ParserException( "No such attribute: " + att.name, 
						randomizedProperties_.getToken(i)->location() );
				}

				result->setRandomized( att );
			}

			scene.addObject( result );
		}
		catch( ParserException& )
		{
			throw;
		}
		catch( Exception& e )
		{
			throw ParserException( "Unable to create object: " + e.message(), 
				nameToken->location() );
		}
	}

protected:
	virtual AttributeList attributes(SceneGraphElementPtr object) = 0;
	virtual SceneGraphElement* createObject(const std::string& name, Scene& scene) = 0;
	virtual void applyAttributes( SceneGraphElementPtr element, Scene& scene, TokenPtr t ) const = 0;

private:
	ParseMultipleAttributes<RandomizedProperty> randomizedProperties_;
	OptionalParseableAttribute<bool> visible_;
	OptionalParseableAttribute<bool> selfCollisions_;
};

class TransformedSceneGraphElementParser
	: public ObjectParser
{
public:
	virtual ~TransformedSceneGraphElementParser()
	{
	}

	TransformedSceneGraphElementParser()
		:	translate_("translate", vl::vl_0), 
			rotate_("rotate", vl::vl_0), 
			scale_("scale", vl::vl_1),
			pivot_("pivot", vl::vl_0)
	{
	}

protected:
	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result;
		result.push_back( &translate_ );
		result.push_back( &rotate_ );
		result.push_back( &scale_ );
		result.push_back( &pivot_ );
		return result;
	}

	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		boost::shared_ptr<TransformedSceneGraphElement> elt = 
			boost::dynamic_pointer_cast<TransformedSceneGraphElement>(element);
		elt->setTranslate( toVec3f(this->translate_.get()) );
		elt->setRotate( toVec3f(this->rotate_.get()) );
		elt->setPivot( toVec3f(this->pivot_.get()) );

		vl::Vec3f s = toVec3f(this->scale_.get());
		for( vl::Int i = 0; i < 3; ++i )
		{
			if( s[i] <= 0.0f )
			{
				throw ParserException( "Scale must be positive.", 
					this->scale_.token()->location() );
			}
		}

		elt->setScale( toVec3f(this->scale_.get()) );
	}

private:
	OptionalParseableAttribute<vl::Vec3d> translate_;
	OptionalParseableAttribute<vl::Vec3d> rotate_;
	OptionalParseableAttribute<vl::Vec3d> scale_;
	OptionalParseableAttribute<vl::Vec3d> pivot_;
};

template <typename ObjectType>
class PhysicsObjectParser
	: public TransformedSceneGraphElementParser
{
public:
	virtual ~PhysicsObjectParser()
	{
	}

	PhysicsObjectParser()
		:	linearVelocity_("linearVelocity", vl::vl_0), 
			angularVelocity_("angularVelocity", vl::vl_0), 
			static_("static", false),
			hasMass_("hasMass", true),
			hasAudio_("hasAudio", true),
			collides_("collides", true),
			materialName_("material", "")
	{
	}

protected:
	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result = TransformedSceneGraphElementParser::attributes(object);
		result.push_back( &linearVelocity_ );
		result.push_back( &angularVelocity_ );
		result.push_back( &static_ );
		result.push_back( &materialName_ );
		result.push_back( &hasMass_ );
		result.push_back( &hasAudio_ );
		result.push_back( &collides_ );
		return result;
	}

	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		TransformedSceneGraphElementParser::applyAttributes( element, scene, t );

		boost::shared_ptr<PhysicsObject> elt = 
			boost::dynamic_pointer_cast<PhysicsObject>(element);
		elt->setLinearVelocity( toVec3f(linearVelocity_.get()) );
		elt->setAngularVelocity( toVec3f(angularVelocity_.get()) );
		elt->setStatic( static_.get() );
		elt->setHasMass( hasMass_.get() );
		elt->setHasAudio( hasAudio_.get() );
		elt->setCollides( collides_.get() );

		std::string materialName = this->materialName_.get();
		if( !materialName.empty() )
		{
			MaterialPtr material = scene.material( materialName );
			if( !material )
			{
				throw ParserException( "No such material: " + materialName, 
					materialName_.token()->location() );
			}

			elt->setMaterial( material );
		}
	}

	PhysicsObject* createObject(const std::string& name, Scene& scene)
	{
		return new ObjectType(name, scene.defaultMaterial());
	}

private:
	OptionalParseableAttribute<vl::Vec3d> linearVelocity_;
	OptionalParseableAttribute<vl::Vec3d> angularVelocity_;
	OptionalParseableAttribute<bool> static_;
	OptionalParseableAttribute<bool> hasMass_;
	OptionalParseableAttribute<bool> hasAudio_;
	OptionalParseableAttribute<bool> collides_;

	OptionalParseableAttribute<std::string> materialName_;
};

class MeshObjectParser
	: public PhysicsObjectParser<MeshObject>
{
public:
	MeshObjectParser()
		:	objFile_("objFile"), convexHulls_("hull")
	{
	}

	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result = PhysicsObjectParser<MeshObject>::attributes(object);
		result.push_back( &objFile_ );
		result.push_back( &convexHulls_ );
		return result;
	}

	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		PhysicsObjectParser<MeshObject>::applyAttributes( element, scene, t );

		std::string objFile = this->objFile_.get(t);

		boost::shared_ptr<MeshObject> elt = 
			boost::dynamic_pointer_cast<MeshObject>(element);
		elt->setMesh( scene.triMesh( objFile ) );

		ConvexHullList hulls;
		std::deque< std::vector<vl::Vec3d> > values = convexHulls_.get();
		for( std::deque< std::vector<vl::Vec3d> >::const_iterator itr = values.begin();
			itr != values.end(); ++itr )
		{
			hulls.push_back( elt->convexHullForPoints(*itr) );
		}
		elt->setConvexHulls( hulls );
	}

	PhysicsObject* createObject(const std::string& name, Scene& scene)
	{
		return new MeshObject(name, scene.defaultMaterial());
	}

private:
	RequiredParseableAttribute<std::string> objFile_;
	MultivaluedParseableAttribute< std::vector<vl::Vec3d> > convexHulls_;
};

class CylinderObjectParser
	: public PhysicsObjectParser<CylinderObject>
{
public:
	CylinderObjectParser()
		:	lengthRatio_("lengthRatio", 1.0)
	{
	}

	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result = PhysicsObjectParser<CylinderObject>::attributes(object);
		result.push_back( &lengthRatio_ );
		return result;
	}

	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		PhysicsObjectParser<CylinderObject>::applyAttributes( element, scene, t );

		boost::shared_ptr<CylinderObject> elt = 
			boost::dynamic_pointer_cast<CylinderObject>(element);
		elt->setLengthRatio( lengthRatio_.get() );
	}

private:
	OptionalParseableAttribute<double> lengthRatio_;
};

class ConeObjectParser
	: public PhysicsObjectParser<ConeObject>
{
public:
	ConeObjectParser()
		:	height_("height", 1.0)
	{
	}

	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result = PhysicsObjectParser<ConeObject>::attributes(object);
		result.push_back( &height_ );
		return result;
	}

	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		PhysicsObjectParser<ConeObject>::applyAttributes( element, scene, t );

		boost::shared_ptr<ConeObject> elt = 
			boost::dynamic_pointer_cast<ConeObject>(element);
		elt->setHeight( height_.get() );
	}

private:
	OptionalParseableAttribute<double> height_;
};

class CappedCylinderObjectParser
	: public PhysicsObjectParser<CappedCylinderObject>
{
public:
	CappedCylinderObjectParser()
		:	lengthRatio_("lengthRatio", 1.0)
	{
	}

	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result = PhysicsObjectParser<CappedCylinderObject>::attributes(object);
		result.push_back( &lengthRatio_ );
		return result;
	}

	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		PhysicsObjectParser<CappedCylinderObject>::applyAttributes( element, scene, t );

		boost::shared_ptr<CappedCylinderObject> elt = 
			boost::dynamic_pointer_cast<CappedCylinderObject>(element);
		elt->setLengthRatio( lengthRatio_.get() );
	}

private:
	OptionalParseableAttribute<double> lengthRatio_;
};

class GroupParser
	: public TransformedSceneGraphElementParser
{
public:
	virtual ~GroupParser()
	{
	}

protected:
	TransformGroup* createObject(const std::string& name, Scene& scene)
	{
		return new TransformGroup(name);
	}
};

template <typename JointType>
class JointParser
	: public TransformedSceneGraphElementParser
{
public:
	JointParser()
		: firstObject_("first"), secondObject_("second"), cfm_("cfm", 1e-5), 
			jointFriction_("jointFriction", 0.0f)
	{
	}

	virtual ~JointParser()
	{
	}

protected:
	void applyAttributes( boost::shared_ptr<SceneGraphElement> element, Scene& scene, TokenPtr t ) const
	{
		boost::shared_ptr<Joint> elt = 
			boost::dynamic_pointer_cast<Joint>(element);
		TransformedSceneGraphElementParser::applyAttributes( element, scene, t );

		std::string firstName = firstObject_.get(t);
		PhysicsObjectPtr first = boost::dynamic_pointer_cast<PhysicsObject>(scene.object( firstName ));
		if( !first )
			throw ParserException( "No such object: '" + firstName + "'", 
				firstObject_.token()->location() );

        if( elt->numObjects() > 1 )
        {
		    std::string secondName = secondObject_.get(t);
		    PhysicsObjectPtr second = boost::dynamic_pointer_cast<PhysicsObject>(scene.object( secondName ));
		    if( !second )
			    throw ParserException( "No such object: '" + secondName + "'", 
				    secondObject_.token()->location() );

            elt->setObjects( PhysicsObjectPair(first, second) );
        }
        else
        {
            elt->setObjects( first );
        }

		elt->setCFM( cfm_.get() );
		for( std::deque< AttributePtr >::const_iterator attItr = otherAttributes_.begin();
			attItr != otherAttributes_.end(); ++attItr )
		{
			elt->setAttribute( (*attItr)->name(), (*attItr)->get() );
		}
	}

	Joint* createObject(const std::string& name, Scene& scene)
	{
		return new JointType(name);
	}

	AttributeList attributes(SceneGraphElementPtr object)
	{
		AttributeList result = TransformedSceneGraphElementParser::attributes(object);
		result.push_back( &firstObject_ );
		result.push_back( &secondObject_ );
		result.push_back( &cfm_ );
		result.push_back( &jointFriction_ );

		boost::shared_ptr<const Joint> joint = 
			boost::dynamic_pointer_cast<const Joint>( object );
		assert( joint );
		std::vector<std::string> jointAttributes = joint->jointAttributes();
		for( std::vector<std::string>::const_iterator jointAttItr = jointAttributes.begin();
			jointAttItr != jointAttributes.end(); ++jointAttItr )
		{
			float defaultValue = joint->getAttribute( *jointAttItr );
			otherAttributes_.push_back( 
				AttributePtr(new OptionalParseableAttribute<double>(
					*jointAttItr, defaultValue) ) );
			result.push_back( otherAttributes_.back().get() );
		}

		return result;
	}

private:
	RequiredParseableAttribute<std::string> firstObject_;
	RequiredParseableAttribute<std::string> secondObject_;
	OptionalParseableAttribute<float> cfm_;
	OptionalParseableAttribute<float> jointFriction_;

	typedef boost::shared_ptr<OptionalParseableAttribute<double> > AttributePtr;
	std::deque< AttributePtr > otherAttributes_;
};

template <typename T>
class GenericObjectParserCreator
	: public ObjectParserCreator
{
public:
	boost::shared_ptr<ObjectParser> parser() const
	{
		boost::shared_ptr<ObjectParser> result(new T);
		return result;
	}
};

template <typename T>
boost::shared_ptr<ObjectParserCreator> createParserCreator()
{
	boost::shared_ptr<ObjectParserCreator> result( 
		new GenericObjectParserCreator<T> );
	return result;
}

ParseObjectMap initParseObjects()
{
	ParseObjectMap result;
	result.insert( std::make_pair( "box", createParserCreator< PhysicsObjectParser<BoxObject> >() ) );
	result.insert( std::make_pair( "sphere", createParserCreator< PhysicsObjectParser<SphereObject> >() ) );
	result.insert( std::make_pair( "plane", createParserCreator< PhysicsObjectParser<PlaneObject> >() ) );
	result.insert( std::make_pair( "capsule", createParserCreator< CappedCylinderObjectParser >() ) );
	result.insert( std::make_pair( "cylinder", createParserCreator< CylinderObjectParser >() ) );
	result.insert( std::make_pair( "cone", createParserCreator< ConeObjectParser >() ) );
	result.insert( std::make_pair( "mesh", createParserCreator< MeshObjectParser >() ) );

	result.insert( std::make_pair( "ballJoint", createParserCreator< JointParser<BallAndSocketJoint> >() ) );
	result.insert( std::make_pair( "sliderJoint", createParserCreator< JointParser<SliderJoint> >() ) );
	result.insert( std::make_pair( "hingeJoint", createParserCreator< JointParser<HingeJoint> >() ) );
	result.insert( std::make_pair( "universalJoint", createParserCreator< JointParser<UniversalJoint> >() ) );
	result.insert( std::make_pair( "hinge2Joint", createParserCreator< JointParser<Hinge2Joint> >() ) );
	result.insert( std::make_pair( "planeJoint", createParserCreator< JointParser<PlaneJoint> >() ) );
	result.insert( std::make_pair( "group", createParserCreator< GroupParser >() ) );

	result.insert( std::make_pair( "fused", createParserCreator< PhysicsObjectParser<CombinedObject> >() ) );


	return result;
}

boost::shared_ptr<Scene> parseSimulation( const std::string& filename )
{
	boost::filesystem::path fn( filename );

	try
	{
		static ParseObjectMap parseObjectMap = initParseObjects();

		boost::shared_ptr<Scene> scene( new Scene() );
		PreprocessorTokenStream stream(fn);

		size_t cameraNo = 0;
		size_t materialNo = 0;

		while( !stream.finished() )
		{
			std::string object = parseString( stream );

			if( object == "material" )
			{
				MaterialPtr mat = parseMaterial( stream, *scene );
				if( materialNo == 0 )
				{
					MaterialPtr defaultMat = scene->defaultMaterial();
					*defaultMat = *mat;
					// only safe b/c this is the very first material:
					defaultMat->setName( mat->name() );
				}
				else
				{
					scene->addMaterial( mat );
				}

				++materialNo;
				continue;
			}

			if( object == "camera" )
			{
				CameraWrapperPtr cam = parseCamera( stream, *scene );
				if( cameraNo == 0 )
				{
					scene->removeCamera( scene->camera() );
					scene->addCamera( cam );
					scene->setCamera( cam );
				}
				else
				{
					scene->addCamera( cam );
				}

				++cameraNo;
				continue;
			}

			ParseObjectMap::const_iterator objectMapItr = 
				parseObjectMap.find(object);
			if( objectMapItr != parseObjectMap.end() )
			{
				boost::shared_ptr<ObjectParser> parser = 
					objectMapItr->second->parser();
				parser->parse( stream, *scene, scene->root(), parseObjectMap );
				continue;
			}

			std::ostringstream oss;
			oss << "Expected one of: ";
			for( objectMapItr = parseObjectMap.begin();
				objectMapItr != parseObjectMap.end(); ++objectMapItr )
			{
				oss << quote(objectMapItr->first) << ", ";
			}

			oss << "'camera', or " << quote("material");

			throw ParserException( oss.str(), stream.peek()->location() ); 
		}

		if( !stream.finished() )
			throw ParserException( "Expected: object declaration.", stream.peek()->location() );

		return scene;
	}
	catch( const boost::filesystem::filesystem_error & ex )
	{
		throw IOException( std::string("Error opening file '") + filename + 
			std::string(": ") + ex.what() );
	}
}

} // namespace planning
