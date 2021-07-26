// This server based on the forking socket server described
// in http://www.scit.wlv.ac.uk/~jphb/comms/sockets.metask.html

// described in http://www.lowtek.com/sockets/select.html

#include "stdafx.h"

#include "parser.h"
#include "scene.h"
#include "simulationTree.h"
#include "mechModel.h"

#include "twigg/exception.h"

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/static_assert.hpp>
#include <boost/crc.hpp>

#include <signal.h>
#include <iostream>
#include <algorithm>

#include "pvm3.h"
#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/poll.h>
#include <sys/select.h>
#include <netinet/in.h>
#include <netdb.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>

#include <exception>
#include <string>

#define HOSTDELETE 12
#define HOSTSUSPEND 13
#define HOSTRESUME 14
#define TASKEXIT 15
#define HOSTADD 16
#define INIT_DATA 17
#define SCENE_ERROR 18
#define RESULT_TAG 30

#define DEBUG_PVM
//#define VERBOSE
//#define DEBUG_NETWORK

// no more than 1MB of data ever:
const size_t maxSize = 1024*1024;

#define ntohll(x) (((_int64)(ntohl((int)((x << 32) >> 32))) << 32) | \
                     (unsigned int)ntohl(((int)(x >> 32)))) //By Runner
#define htonll(x) ntohll(x)

namespace planning
{
void handle( int wsd );
void pvmClient(int parent);
}

bool die = false;
void handleKill( int value )
{
	std::cout << "Caught interrupt." << std::endl;
	die = true;
}

void my_pvm_perror( const char* function, const char* file, int line )
{
	std::ostringstream oss;
	oss << "Error in " << function << "; line " << line << " in file " << file;
	std::string value = oss.str();
	std::vector<char> vecValue;  vecValue.reserve(value.size()+1);
	std::copy(value.begin(), value.end(), std::back_inserter(vecValue) );
	vecValue.push_back(0);
	pvm_perror( &vecValue[0] );
}

#define PVM_PERROR(function) my_pvm_perror(#function, __FILE__, __LINE__);

std::vector<char> fullPath( 1024 );

class SocketException
	: std::exception
{
public:
	SocketException( const char* call )
		: call_( call ) {}

	const char* call() const
	{
		return this->call_;
	}
	
private:
	const char* call_;
};

class PVMException
	: std::exception
{
public:
	PVMException( const char* call, int info )
		: call_(call), info_(info) {}

private:
	const char* call_;
	int info_;
};

class SocketWrapper
{
public:
	SocketWrapper( int socket )
		: socket_(socket)
	{
		if( socket_ < 0 )
		{
			throw SocketException("socket");
		}
	}

	~SocketWrapper()
	{
		close( socket_ );
	}

	int socket() const
	{
		return socket_;
	}

private:
	int socket_;
};

namespace po = boost::program_options;

class PVMAutoExit
{
public:
	PVMAutoExit()  {}
	~PVMAutoExit() { pvm_exit(); }
};

class PVMHost
{
public:
	PVMHost( const char* name, int hostId )
		: hostname_( name ), hostId_( hostId )
	{
	}

	bool operator==( const PVMHost& other ) const
	{
		return (this->hostId_ == other.hostId_);
	}

	bool operator<( const PVMHost& other ) const
	{
		return (this->hostId_ < other.hostId_);
	}

	int hostId() const
	{
		return this->hostId_;
	}

	std::string name() const
	{
		return this->hostname_;
	}

private:
	std::string hostname_;
	int hostId_;
};

class PVMTask
{
public:
	PVMTask( int tid, PVMHost host )
		: tid_(tid), host_(host) {}

	~PVMTask()
	{
		int info = pvm_kill(tid_);
#ifdef VERBOSE
		std::cout << "Ended task: " << this->tid_ << "\n";
#endif
		// can't throw an exception here, so if it doesn't
		// kill it properly we'll just have to ignore it.
		if( info < 0 )
			PVM_PERROR(pvm_kill);
	}

	PVMHost host() const
	{
		return this->host_;
	}

	int tid() const
	{
		return this->tid_;
	}

private:
	int tid_;
	PVMHost host_;
};


std::deque<PVMHost> pvmHosts;
int mytid;

void handleHostAdd();
void handleHostDelete();
void handleSceneError();

int main(int argc, char** argv)
{
	chdir( "/nfs/hn00/cdtwigg/planning" );

/*
	signal( SIGINT, handleKill );
	signal( SIGHUP, handleKill );
	signal( SIGTERM, handleKill );
*/

	realpath( argv[0], &fullPath[0] );

	PVMAutoExit autoExit;
	//pvm_catchout(stdout);

	pvm_mytid();
	if( mytid < 0 )
	{
		PVM_PERROR(pvm_mytid);
		return EXIT_FAILURE;
	}

	int myparent = pvm_parent();
	if ((myparent < 0) && (myparent != PvmNoParent))
	{
		PVM_PERROR(pvm_parent);
		return EXIT_FAILURE;
	}

	if( myparent != PvmNoParent )
	{
		try
		{
			planning::pvmClient(myparent);
			return 0;
		}
		catch( PVMException& e )
		{
			PVM_PERROR("init")
		}
	}

	int defaultPort = 10037;
	std::vector< std::string > defaultClients;
	defaultClients.push_back( "slide.graphics.cs.cmu.edu" );
	defaultClients.push_back( "bend.graphics.cs.cmu.edu" );

	std::ostringstream portNumberStr;
	portNumberStr << "port number (default: " << defaultPort << ")";

	std::ostringstream clientStr;
	clientStr << "client to listen to (default: ";
	for( std::vector< std::string >::const_iterator itr = defaultClients.begin();
		itr != defaultClients.end(); ++itr )
	{
		if( itr != defaultClients.begin() )
			clientStr << ", ";
		clientStr << *itr;
	}
	clientStr << ")";

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("port", po::value<int>(), portNumberStr.str().c_str())
		("client", po::value<std::vector<std::string> >(), clientStr.str().c_str());

	// first, we will add an initial host.  Once this
	//   gets added, we will add more
	{
		int info = pvm_notify( PvmHostAdd, HOSTADD, -1, NULL );
	}
#ifdef DEBUG_PVM
	{
		pvmhostinfo *hostp;
		int nhost, narch;
		int info = pvm_config( &nhost, &narch, &hostp );
		assert( info == 0 );
		for( int i = 0; i < nhost; ++i )
			pvmHosts.push_back( PVMHost(hostp[i].hi_name, hostp[i].hi_tid) );
	}
#else
	{
		char* zero = "0";
		const size_t numInitHosts = 5;

		int infos[numInitHosts];
		char* hosts[numInitHosts];
		for( size_t i = 0; i < numInitHosts; ++i )
			hosts[i] = zero;

		int info = pvm_addhosts( hosts, numInitHosts, infos );
		if( info < 1 )
		{
			std::cout << "Unable to start PVM." << std::endl;
			PVM_PERROR(pvm_addhosts);
			return 1;
		}
	}
#endif

	po::variables_map vm;
	try
	{
		po::store(po::command_line_parser(argc, argv).
			options(desc).run(), vm);
		po::notify(vm);
	}
	catch( po::error& err )
	{
		std::cout << err.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	if( vm.count("help"))
	{
		std::cout << desc << "\n";
		return 0;
	}

	int port = defaultPort;
	if( vm.count("port") )
	{
		if( port < 1024 || port > 65535 )
		{
			std::cout << "Port must be between 1024 and 65535\n";
			return EXIT_FAILURE;
		}

		port = vm["port"].as<int>();
	}

	std::vector< std::string > clients = defaultClients;
	if( vm.count("client") )
	{
		std::vector< std::string > otherClients = 
			vm["client"].as< std::vector< std::string > >();
		clients = otherClients;
	}

	std::vector<unsigned long> addresses;

	for( std::vector< std::string >::const_iterator itr = clients.begin();
		itr != clients.end(); ++itr )
	{
		hostent* hp = 
			gethostbyname2( itr->c_str(), AF_INET );
		if( hp == 0 )
		{
			std::cout << "No such host: " << *itr << "\n";
			return EXIT_FAILURE;
		}

		assert( hp->h_addrtype == AF_INET );
		for( int i = 0; hp->h_addr_list[i] != NULL; ++i )
		{
			struct in_addr* address = ( struct in_addr*)( hp -> h_addr_list[i]);
			addresses.push_back( address->s_addr);
		}
	}

	try
	{
		// write my host name to the appropriate file
		{
			std::vector<char> hostname( 256 );
			int result = gethostname( &hostname[0], hostname.size() );
			assert( result == 0 );
			std::ofstream ofs( "/nfs/hn00/cdtwigg/planning/serverAddress" );
			ofs << &hostname[0];
		}

		SocketWrapper sockWrapper( 
			socket(AF_INET, SOCK_STREAM, 0) );

		/*
		{
			int reuse_addr = 1;
			setsockopt(sockWrapper.socket(), SOL_SOCKET, SO_REUSEADDR, &reuse_addr,
				sizeof(reuse_addr));
		}
		*/


		// set socket to nonblocking
		/*
		{
			int opts = fcntl(sockWrapper.socket(), F_GETFL);
			if( opts < 0 )
				throw SocketException( "fcntl(F_GETFL)" );

			opts = (opts | O_NONBLOCK);
			if( fcntl(sockWrapper.socket(), F_GETFL) < 0 )
				throw SocketException( "fcntl(F_SETFL)" );
		}
		*/

		// Get the address information, and bind it to the socket
		{
			sockaddr_in server_address;
			bzero((char *) &server_address, sizeof(server_address));
			server_address.sin_family = AF_INET;
			server_address.sin_addr.s_addr = INADDR_ANY;
			server_address.sin_port = htons(port);

			if( bind(sockWrapper.socket(), 
				reinterpret_cast<sockaddr*>(&server_address), 
				sizeof(server_address)) < 0 )
			{
				throw SocketException( "bind" );
			}

		}

		// set up queue for incoming connections.
		if( listen(sockWrapper.socket(), 1) == -1 )
		{
			throw SocketException( "listen" );
		}

		std::cout << "Listening on port " << port << ".\n";
		std::cout << "IP addresses accepted: ";
		for( std::vector<unsigned long>::const_iterator itr = addresses.begin();
			itr != addresses.end(); ++itr )
		{
			in_addr address;
			address.s_addr = *itr;
			// note that this is not thread-safe
			std::cout << inet_ntoa(address);
			if( (itr+1) != addresses.end() )
				std::cout << ", ";
		}
		std::cout << std::endl;
	
		// now, start the main server loop
		while( !die )
		{
			int *d;
			int nfds = pvm_getfds(&d);
#ifdef VERBOSE
			std::cout << "Num condor fds: " << nfds << "\n";
#endif

			int fd_max = sockWrapper.socket();
			fd_set socks;
			FD_ZERO(&socks);
			FD_SET(sockWrapper.socket(), &socks);
			for( int i = 0; i < nfds; ++i )
			{
				FD_SET(d[i], &socks);
				fd_max = std::max( fd_max, d[i] );
			}

#ifdef VERBOSE
			std::cout << "Selecting... " << std::flush;
#endif
			if(select(fd_max+1, &socks, NULL, NULL, NULL) == -1)
			{
				perror( "select" );
				// ignore select() errors
				continue;
			}
#ifdef VERBOSE
			std::cout << "saw something." << std::endl;
#endif

			if( FD_ISSET(sockWrapper.socket(), &socks) )
			{
#ifdef VERBOSE
				std::cout << "Accepting... ";
#endif
				sockaddr_in addr;
				socklen_t addlen = sizeof(addr);
				int wsd = accept(sockWrapper.socket(), reinterpret_cast<sockaddr*>(&addr), &addlen);

				if( wsd < 0 )
				{
					std::cout << "Error on accept.\n";
					continue;
				}
#ifdef VERBOSE
				std::cout << "done. \n";
#endif

				std::vector<unsigned long>::iterator itr = 
					std::find(addresses.begin(), addresses.end(), addr.sin_addr.s_addr);
	
				if( addr.sin_family != AF_INET
					|| itr == addresses.end() )
				{
					std::cout << "Dropped illegal connection from host: ";
					std::cout << inet_ntoa(addr.sin_addr) << std::endl;
					close( wsd );
					continue;
				}

				std::cout << "Received connection from host: ";
				std::cout << inet_ntoa(addr.sin_addr) << std::endl;
	
				planning::handle( wsd );
			}
			else
			{
				std::cout << "Receiving message... " << std::flush;
				int buf_id = pvm_nrecv(-1, -1);
				std::cout << "Done." << std::endl;

				if( buf_id == 0 )
				{
					continue;
				}
				else if( buf_id < 0 )
				{
					std::cout << "Error caught in pvm_nrecv; ignoring." << std::endl;
					continue;
				}
				if( buf_id == PvmSysErr )
				{
					PVM_PERROR(pvm_recv);
					throw SocketException( "PVM server failure." );
				}

				int msg_len;
				int msg_tag;
				int msg_src;
				int info = pvm_bufinfo(buf_id, &msg_len, &msg_tag, &msg_src);
				assert( info != PvmNoSuchBuf && info != PvmBadParam );
#ifdef VERBOSE
			std::cout << "Message tag: " << msg_tag << std::endl;
#endif
				if( msg_tag == HOSTADD )
				{
					handleHostAdd();
				}
				else if( msg_tag == HOSTDELETE )
				{
					handleHostDelete();
				}
				else if( msg_tag == TASKEXIT )
				{
					// ignore
				}
				else if( msg_tag == SCENE_ERROR )
				{
					handleSceneError();
				}
				if( msg_tag >= RESULT_TAG )
				{
					// need to kill old process
					// don't want to kill myself though
					if( msg_src != mytid )
						pvm_kill( msg_src );
				}
				else
				{
					// ignore
				}
			}
		}
	}
	catch( SocketException& e )
	{
		perror( e.call() );
		return EXIT_FAILURE;
	}
}

class FDWrapper
{
public:
	FDWrapper( int fd )
		: fd_(fd) {}

	~FDWrapper()
	{
		close( fd_ );
	}

	int fd() const
	{
		return this->fd_;
	}

private:
	int fd_;
};

void readBuf( int connection, char* buf, size_t num )
{
	size_t current = 0;
	while( current < num )
	{
		ssize_t numRead = read( connection, buf + current, num - current );
		if( numRead == 0 )
			throw SocketException( "read: EOF" );
		else if( numRead < 0 )
			throw SocketException( "read: error" );
		current += numRead;
	}
}

uint32_t readUnsigned( int connection )
{
	uint32_t result;
	readBuf( connection, reinterpret_cast<char*>( &result ), sizeof(uint32_t) );
	return ntohl( result );
}

float readFloat( int connection )
{
	uint32_t uintResult = readUnsigned( connection );
	const float* result = reinterpret_cast<const float*>( &uintResult );
	return *result;
}

template <typename T>
T toHostOrder( T value );

template <>
float toHostOrder( float value )
{
	const uint32_t* ptr = reinterpret_cast<const uint32_t*>( &value );
	uint32_t result = ntohl( *ptr );
	const float* res = reinterpret_cast<const float*>( &result );
	return *res;
}

template <>
uint32_t toHostOrder( uint32_t value )
{
	return ntohl( value );
}


template <typename T>
void readFromSocket( int connection, std::vector<T>& values )
{
	uint32_t size = readUnsigned( connection );
#ifdef DEBUG_NETWORK
	std::cout << "Read size: " << size << std::endl;
#endif
	if( size > maxSize )
		throw SocketException( "Array size too big." );

	values.resize( size );

	readBuf( connection, reinterpret_cast<char*>( &values[0] ), size*sizeof(T) );

	std::transform( values.begin(), values.end(), values.begin(), 
		&toHostOrder<T> );
/*
#ifdef DEBUG_NETWORK
	std::cout << "Read values: ";
	for( size_t i = 0; i < values.size(); ++i )
		std::cout << values[i] << " ";
	std::cout << std::endl;
#endif
*/
}

// Actually have to be careful here; anything other than char will
//   need the byte order switched so will have to be treated differently
template <typename T>
void writeToSocket( int connection, const T& value );

void send_uint32( int connection, uint32_t value )
{
	uint32_t networkOrder = htonl( value );
	int result = send( connection, &networkOrder, 4, MSG_NOSIGNAL );
	if( result == 0 )
		throw SocketException( "send: nothing written" );
	else if( result < 0 )
		throw SocketException( "send: error" );
}

template <>
void writeToSocket<uint32_t>( int connection, const uint32_t& value )
{
	send_uint32( connection, value );
	send_uint32( connection, planning::checksum(value) );
}

template <>
void writeToSocket<float>( int connection, const float& value )
{
	const uint32_t* uintVal = reinterpret_cast<const uint32_t*>( &value );
	writeToSocket( connection, *uintVal );
	send_uint32( connection, planning::checksum(value) );
}

template <>
void writeToSocket<std::string>( int connection, const std::string& value )
{
	writeToSocket<uint32_t>( connection, value.size() );
	send_uint32( connection, planning::checksum(value) );
	if( value.empty() )
		return;
	ssize_t result = send( connection, 
		&value[0], value.size(), MSG_NOSIGNAL );
	if( result == 0 )
		throw SocketException( "send: nothing written" );
	else if( result < 0 )
		throw SocketException( "send: error" );
}


template <>
void writeToSocket< std::vector<char> >( int connection, const std::vector<char>& vec )
{
	writeToSocket<uint32_t>( connection, vec.size() );
	send_uint32( connection, planning::checksum(vec) );
	if( vec.empty() )
		return;
	ssize_t result = send( connection, 
		&vec[0], vec.size(), MSG_NOSIGNAL );
	if( result == 0 )
		throw SocketException( "send: nothing written" );
	else if( result < 0 )
		throw SocketException( "send: error" );
}

template <>
void writeToSocket< std::vector<float> >( int connection, const std::vector<float>& vec )
{
	writeToSocket<uint32_t>( connection, vec.size() );
	send_uint32( connection, planning::checksum(vec) );
	if( vec.empty() )
		return;
	std::vector<uint32_t> networkOrder( vec.size() );
	for( size_t i = 0; i < vec.size(); ++i )
		networkOrder[i] = htonl( *reinterpret_cast<const uint32_t*>( &vec[0] + i ) );

	ssize_t result = send( connection, 
		&networkOrder[0], sizeof(uint32_t)*networkOrder.size(), MSG_NOSIGNAL );
	if( result == 0 )
		throw SocketException( "send: nothing written" );
	else if( result < 0 )
		throw SocketException( "send: error" );
}

template <>
void writeToSocket< std::vector<vl::Vec3f> >( int connection, const std::vector<vl::Vec3f>& vec )
{
	writeToSocket<uint32_t>( connection, vec.size() );
	send_uint32( connection, planning::checksum(vec) );
	if( vec.empty() )
		return;
	std::vector<uint32_t> networkOrder( 3*vec.size() );
	for( size_t i = 0; i < vec.size(); ++i )
		for( size_t j = 0; j < 3; ++j )
			networkOrder[3*i+j] = htonl( *reinterpret_cast<const uint32_t*>( vec[i].Ref() + j ) );

	ssize_t result = send( connection, 
		&networkOrder[0], sizeof(uint32_t)*networkOrder.size(), MSG_NOSIGNAL );
	if( result == 0 )
		throw SocketException( "send: nothing written" );
	else if( result < 0 )
		throw SocketException( "send: error" );
}

template <typename T>
void pvm_unpack( T& value );

template <>
void pvm_unpack<unsigned int>( unsigned int& value )
{
	int info = pvm_upkuint( &value, 1, 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkuint", info );
}

template <>
void pvm_unpack<float>( float& value )
{
	int info = pvm_upkfloat( &value, 1, 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkfloat", info );
}

template <>
void pvm_unpack< std::vector<char> >( std::vector<char>& value )
{
	unsigned int size;
	pvm_unpack(size);
	value.resize( size );
	int info = pvm_upkbyte( &value[0], size, 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkbyte", info );
}

template <>
void pvm_unpack< std::string >( std::string& value )
{
	unsigned int size;
	pvm_unpack(size);
#ifdef DEBUG_NETWORK
	std::cout << "unpacked size: " << size << std::endl;
#endif
	value.resize( size, ' ' );
	int info = pvm_upkbyte( &value[0], size, 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkbyte", info );
}

template <>
void pvm_unpack< std::vector<vl::Vec3f> >( std::vector<vl::Vec3f>& value )
{
	unsigned int size;
	pvm_unpack(size);
	value.resize(size);
	int info = pvm_upkfloat( reinterpret_cast<float*>( &value[0] ), 3*value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkfloat", info );
}

template <>
void pvm_unpack< std::vector<float> >( std::vector<float>& value )
{
	unsigned int size;
	pvm_unpack(size);
	value.resize(size);
	int info = pvm_upkfloat( &value[0], value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkfloat", info );
}

template <>
void pvm_unpack< std::vector<unsigned int> >( std::vector<unsigned int>& value )
{
	unsigned int size;
	pvm_unpack(size);
#ifdef DEBUG_NETWORK
	std::cout << "unpacked size: " << size << std::endl;
#endif
	value.resize(size);
	int info = pvm_upkuint( &value[0], value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_upkuint", info );
}


template <typename T>
void pvm_pack( T& value );

template <>
void pvm_pack<unsigned int>( unsigned int& value )
{
	int info = pvm_pkuint( &value, 1, 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkuint", info );
}

template <>
void pvm_pack<float>( float& value )
{
	int info = pvm_pkfloat( &value, 1, 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkfloat", info );
}

template <>
void pvm_pack< std::string >( std::string& value )
{
	unsigned int size = value.size();
	pvm_pack( size );
#ifdef DEBUG_NETWORK
	std::cout << "packed size: " << size << std::endl;
#endif
	int info = pvm_pkbyte( &value[0], value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkbyte", info );
}

template <>
void pvm_pack< std::vector<char> >( std::vector<char>& value )
{
	unsigned int size = value.size();
	pvm_pack<unsigned int>( size );
	int info = pvm_pkbyte( &value[0], value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkbyte", info );
}

template <>
void pvm_pack< std::vector<unsigned int> >( std::vector<unsigned int>& value )
{
	unsigned int size = value.size();
	pvm_pack<unsigned int>( size );
#ifdef DEBUG_NETWORK
	std::cout << "packed size: " << size << std::endl;
#endif
	int info = pvm_pkuint( &value[0], value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkuint", info );
}

template <>
void pvm_pack< std::vector<float> >( std::vector<float>& value )
{
	unsigned int size = value.size();
	pvm_pack<unsigned int>( size );
	int info = pvm_pkfloat( &value[0], value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkfloat", info );
}

template <>
void pvm_pack< std::vector<vl::Vec3f> >( std::vector<vl::Vec3f>& value )
{
	unsigned int size = value.size();
	pvm_pack<unsigned int>( size );
	int info = pvm_pkfloat( reinterpret_cast<float*>( &value[0] ), 
		3*value.size(), 1 );
	if( info < 0 )
		throw PVMException( "pvm_pkfloat", info );
}

void handleSceneError()
{
	std::string error;
	pvm_unpack( error );
	std::cout << "Caught error in scene: " << error << std::endl;
}

void handleHostAdd()
{
	// we got a new host.
	int tid;
	int result = pvm_upkint(&tid, 1, 1);
	// if we get an error here we have major problems, and should probably die
	assert( result >= 0 );

	// this is scary because PVM doesn't provide a way to
	// specify the buffer length
	char hostname[ 1000 ];
	pvm_upkstr( &hostname[0] );
#ifdef VERBOSE
	std::cout << "HOSTADD: " << hostname << std::endl;
#endif

	PVMHost newHost( hostname, tid );
	if( std::find( pvmHosts.begin(), pvmHosts.end(), newHost ) == pvmHosts.end() )
	{
		pvmHosts.push_back( newHost );
		std::cout << "Host " << tid << " added.\n";
	}
	else
	{
		std::cout << "Warning: duplicate add detected of host " << tid << ".\n";
	}
	
	// want to know when this host disappears
	pvm_notify(PvmHostDelete, HOSTDELETE, 1, &tid);

#ifndef DEBUG_PVM
	// also, add another host
	{
		int infos[1];
		char* hosts[] = { "0" };
		int info = pvm_addhosts( hosts, 1, infos );
		if( info < 1 )
		{
			std::cout << "Unable to add another host." << std::endl;
			PVM_PERROR(pvm_addhosts);
		}
	}
#endif
}

void handleHostDelete()
{
	int tid;
	int result = pvm_upkint(&tid, 1, 1);
	assert( result >= 0 );

	// don't need to worry about multiple copies because we handled that above
	std::deque<PVMHost>::iterator itr = 
		std::remove( pvmHosts.begin(), pvmHosts.end(), 
			PVMHost( "", tid ) );
	if( itr == pvmHosts.end() )
	{
		std::cout << "Warning: received HostDelete notification for nonexistent host " << tid << ".\n";
	}
	else
	{
		pvmHosts.erase( itr );
	}
}

namespace planning
{

void handle( int wsd_in )
{
	// we will need to make sure that any messages we get 
	//   are not kicking around from the last connection
	// We will therefore use this unique id to identify them
	// WARNING: not currently thread safe
	static unsigned int idCounter = 0;
	unsigned int myId = idCounter++;

	typedef boost::shared_ptr<PVMTask> PVMTaskPtr;
	typedef std::deque<PVMTaskPtr> PVMTaskList;
	PVMTaskList tasks;

	try
	{
		FDWrapper wsd( wsd_in );

		uint32_t frameRate = readUnsigned( wsd.fd() );

		uint32_t sceneSize = readUnsigned( wsd.fd() );
		std::string sceneDesc( sceneSize, ' ' );
		readBuf( wsd.fd(), &sceneDesc[0], sceneSize );

		uint32_t parentSceneSize = readUnsigned( wsd.fd() );
		std::string parentSceneDesc( parentSceneSize, ' ' );
		readBuf( wsd.fd(), &parentSceneDesc[0], parentSceneSize );

#ifdef DEBUG_NETWORK
		std::cout << "Read scene description." << std::endl;
#endif

		std::vector<unsigned int> fixedObjects;
		readFromSocket( wsd.fd(), fixedObjects );
#ifdef DEBUG_NETWORK
		std::cout << "Read fixed objects: ";
		for( size_t i = 0; i < fixedObjects.size(); ++i )
			std::cout << " " << fixedObjects[i];
		std::cout << std::endl;
#endif

		float startTime = readFloat( wsd.fd() );
#ifdef DEBUG_NETWORK
		std::cout << "Read start time:" << startTime << std::endl;
#endif

		float endTime = readFloat( wsd.fd() );
#ifdef DEBUG_NETWORK
		std::cout << "Read end time:" << endTime << std::endl;
#endif

		uint32_t nPaths = readUnsigned( wsd.fd() );
#ifdef DEBUG_NETWORK
		std::cout << "Read nPaths:" << nPaths << std::endl;
#endif

		std::vector< std::vector<unsigned int> > positionTimes( nPaths );
		std::vector< std::vector<float> > positionCoeffs( nPaths );
		std::vector< std::vector<unsigned int> > rotationTimes( nPaths );
		std::vector< std::vector<float> > rotations( nPaths );
		std::vector< std::vector<float> > rotationCoeffs( nPaths );

		for( size_t i = 0; i < nPaths; ++i )
		{
			readFromSocket( wsd.fd(), positionTimes[i] );
#ifdef DEBUG_NETWORK
			std::cout << "Read position times." << std::endl;
#endif
			readFromSocket( wsd.fd(), positionCoeffs[i] );
#ifdef DEBUG_NETWORK
			std::cout << "from socket: positionTimes: " << positionTimes.size() << "; positionCoeffs: " << positionCoeffs.size() << std::endl;
#endif
			readFromSocket( wsd.fd(), rotationTimes[i] );
			readFromSocket( wsd.fd(), rotations[i] );
			readFromSocket( wsd.fd(), rotationCoeffs[i] );
		}

		while( !die )
		{
			// detect if user requests disconnect
			{
				pollfd fds[1];
				fds[0].fd = wsd.fd();
				fds[0].events = POLLIN;
				int rv = poll( fds, 1, 0 );
				if( fds[0].revents & POLLIN )
				{
					// for now we will assume always that the user
					//   wants to close
					std::cout << "Remote user requested connection closed." << std::endl;
					break;
				}
			}

			// if we have more hosts than jobs, start up new ones
			while( tasks.size() < (pvmHosts.size()-1) )
			{
				std::set<int> usedMachines;
				for( PVMTaskList::const_iterator taskItr = tasks.begin();
					taskItr != tasks.end(); ++taskItr )
				{
					usedMachines.insert( (*taskItr)->host().hostId() );
				}

				// now find an empty machine:
				std::deque<PVMHost>::const_iterator hostItr = pvmHosts.begin();
				while( hostItr != pvmHosts.end() )
				{
					if( usedMachines.find( hostItr->hostId() ) == 
						usedMachines.end() )
					{
						break;
					}

					++hostItr;
				}

				if( hostItr == pvmHosts.end() )
					break;

				int tid;

				char* arch = "0";
				std::cout << "Spawning PVM process: '" << &fullPath[0] << "'..." << std::flush;
#ifdef DEBUG_PVM
				int numt = pvm_spawn( &fullPath[0], NULL, PvmTaskDefault, 0, 1, &tid );
#else
				/*
				std::string hn = hostItr->name();
				std::vector<char> hostname( hn.begin(), hn.end() );
				hostname.push_back( 0 );
				*/
				int numt = pvm_spawn( "/nfs/hn00/cdtwigg/planning/runSampleServer", NULL, PvmTaskArch, "0", 1, &tid );
#endif
				std::cout << "Done." << std::endl;
				PVMTaskPtr task( new PVMTask(tid, *hostItr) );
				if( numt == 1 )
				{
#ifdef VERBOSE
					std::cout << "Sending over scene... \n" << std::flush;
#endif
					int status = pvm_notify(PvmTaskExit, TASKEXIT, 1, &tid);
					if (pvm_pstat(tid) != PvmOk)
					{
						PVM_PERROR(pvm_pstat);
						break;
					}

					// need to send the scene on over
					int bufid = pvm_initsend( PvmDataRaw );
					if( bufid < 0 )
					{
						PVM_PERROR(pvm_initsend);
						break;
					}

					pvm_pack( myId );
#ifdef VERBOSE
					std::cout << "Sent id: " << myId << std::endl;
#endif
					pvm_pack( frameRate );
#ifdef VERBOSE
					std::cout << "Sent frame rate: " << frameRate << std::endl;
#endif
					pvm_pack( sceneDesc );
#ifdef VERBOSE
					std::cout << "Sent scene desc" << std::endl;
#endif
					pvm_pack( parentSceneDesc );
#ifdef VERBOSE
					std::cout << "Sent parent scene desc" << std::endl;
#endif

					pvm_pack( fixedObjects );
#ifdef VERBOSE
					std::cout << "Sent fixed objects" << std::endl;
#endif
					pvm_pack( startTime );
#ifdef VERBOSE
					std::cout << "Sent start time" << std::endl;
#endif
					pvm_pack( endTime );
#ifdef VERBOSE
					std::cout << "Sent end time" << std::endl;
#endif

					pvm_pack( nPaths );
					for( size_t i = 0; i < nPaths; ++i )
					{
						pvm_pack( positionTimes[i] );
						pvm_pack( positionCoeffs[i] );
						pvm_pack( rotationTimes[i] );
						pvm_pack( rotations[i] );
						pvm_pack( rotationCoeffs[i] );
					}
#ifdef VERBOSE
					std::cout << "Sent paths" << std::endl;
#endif

					{
						int info = pvm_send( tid, INIT_DATA );
						if( info < 0 )
						{
							PVM_PERROR(pvm_send);
							break;
						}
					}
#ifdef VERBOSE
					std::cout << "Done. " << std::endl;
#endif
					tasks.push_back( task );
				}
				else
				{
					std::cout << "Failed to spawn task.\n";
					PVM_PERROR(pvm_spawn);
					break;
				}
			}

			// now, check for any messages from PVM
			timeval tv;
			tv.tv_sec = 1;
			tv.tv_usec = 0;
#ifdef VERBOSE
			std::cout << "Receive messages... " << std::flush;
#endif
			int buf_id = pvm_trecv(-1, -1, &tv);
#ifdef VERBOSE
			std::cout << "Done." << std::endl;
#endif

			if( buf_id == PvmSysErr )
			{
				std::cout << "PVM server failure." << std::endl;
				PVM_PERROR(pvm_recv);
				break;
			}
			else if( buf_id <= 0 )
			{
#ifdef VERBOSE
				std::cout << "Timed out." << std::endl;
#endif
				// timed out
				continue;
			}

			int msg_len;
			int msg_tag;
			int msg_src;
			int info = pvm_bufinfo(buf_id, &msg_len, &msg_tag, &msg_src);
			assert( info != PvmNoSuchBuf && info != PvmBadParam );

#ifdef VERBOSE
			std::cout << "Message tag: " << msg_tag << std::endl;
#endif

			if( msg_tag == (RESULT_TAG+myId) )
			{
				// result returned by worker
				// first read it off the PVM wire
				std::vector<CompressedPath> paths;
				std::vector< std::vector<vl::Vec3f> > cmPaths;

				typedef std::pair< float, std::vector<RigidStaticState> > Snap;
				std::vector<Snap> snapshots;
				std::vector<float> metricScores;

				try
				{
					unsigned int nPaths;
					pvm_unpack( nPaths );
					paths.resize( nPaths );
					cmPaths.resize( nPaths );
					for( size_t i = 0; i < nPaths; ++i )
					{
						pvm_unpack( paths[i].positionTimes );
						for( size_t j = 0; j < 3; ++j )
							pvm_unpack( paths[i].positionLinearCoeffs[j] );
						pvm_unpack( paths[i].positionQuadraticCoeffs[1] );

						pvm_unpack( paths[i].rotationTimes );
						for( size_t j = 0; j < 3; ++j )
							pvm_unpack( paths[i].rotationCoeffs[j] );

						pvm_unpack( cmPaths[i] );
					}

					unsigned int nSnaps;
					pvm_unpack( nSnaps );

					snapshots.resize( nSnaps );
					for( size_t iSnapshot = 0; iSnapshot < nSnaps; ++iSnapshot )
					{
						pvm_unpack( snapshots[iSnapshot].first );

						unsigned int numStates;
						pvm_unpack( numStates );
						std::vector<RigidStaticState>& state = snapshots[iSnapshot].second;
						state.resize( numStates );
						int info = pvm_upkfloat( 
							reinterpret_cast<float*>(&state[0]), 
							7*state.size(), 1 );
						if( info < 0 )
							throw PVMException( "pvm_unpkfloat", info );
					}

#ifdef VERBOSE
					std::cout << "recv " << nSnaps << " snapshots\n";
#endif

					pvm_unpack( metricScores );
#ifdef VERBOSE
					std::cout << "recv " << nSnaps << " scores\n";
#endif
				}
				catch( PVMException& e )
				{
					PVM_PERROR(pvm_unpack);
					continue;
				}

				// now write back to the client
				// note that if we had a way to error out a path in the 
				//   middle, we could interleave this with the PVM receive
				//   step.

				// report on the number of used machines
				writeToSocket<uint32_t>( wsd.fd(), tasks.size() );

				writeToSocket<uint32_t>( wsd.fd(), paths.size() );
				for( size_t iPath = 0; iPath < paths.size(); ++iPath )
				{
					writeToSocket( wsd.fd(), paths[iPath].positionTimes );
					for( size_t i = 0; i < 3; ++i )
						writeToSocket( wsd.fd(), paths[iPath].positionLinearCoeffs[i] );
					writeToSocket( wsd.fd(), paths[iPath].positionQuadraticCoeffs[1] );
	
					writeToSocket( wsd.fd(), paths[iPath].rotationTimes );
					for( size_t i = 0; i < 3; ++i )
						writeToSocket( wsd.fd(), paths[iPath].rotationCoeffs[i] );

					writeToSocket( wsd.fd(), cmPaths[iPath] );
				}

				writeToSocket<uint32_t>( wsd.fd(), snapshots.size() );
				for( size_t iSnapshot = 0; iSnapshot < snapshots.size(); ++iSnapshot )
				{
					writeToSocket( wsd.fd(), snapshots[iSnapshot].first );

					const size_t nStates = snapshots[iSnapshot].second.size();
					writeToSocket( wsd.fd(), nStates );

					const uint32_t* start = 
						reinterpret_cast<const uint32_t*>( &(snapshots[iSnapshot].second)[0] );
					std::vector<uint32_t> networkOrder( 7*nStates );
					std::transform( start, start + 7*nStates,
						networkOrder.begin(),
						&htonl );

					ssize_t result = send( wsd.fd(), 
						&networkOrder[0], sizeof(uint32_t)*networkOrder.size(), MSG_NOSIGNAL );
					if( result == 0 )
						throw SocketException( "send: nothing written" );
					else if( result < 0 )
						throw SocketException( "send: error" );
				}

				writeToSocket( wsd.fd(), metricScores );
			}
			else if( msg_tag == HOSTADD )
			{
				handleHostAdd();
			}
			else if( msg_tag == HOSTDELETE )
			{
				handleHostDelete();
			}
			else if( msg_tag == SCENE_ERROR )
			{
				handleSceneError();
			}
			else if( msg_tag == TASKEXIT )
			{
#ifdef VERBOSE
				std::cout << "Received task exit message.\n";
#endif
				// remove the task from the list of tasks
				int tid;
				int result = pvm_upkint(&tid, 1, 1);

				for( PVMTaskList::iterator itr = tasks.begin();
					itr != tasks.end(); )
				{
					if( (*itr)->tid() == tid )
					{
						std::cout << "Erasing task " << tid 
							<< " running on " << (*itr)->host().name() << std::endl;
						itr = tasks.erase( itr );
					}
					else
						++itr;
				}
			}
			else
			{
				if( msg_tag >= RESULT_TAG )
				{
					// need to kill old process
					pvm_kill( msg_src );
					std::cout << "Killed process " << msg_src 
						<< " due to unexpected mesage." << std::endl;
				}
				else
				{
					std::cout << "Warning: ignored tag " << msg_tag << std::endl;
				}
			}


			// @todo fix this?
			//fflush( wsd.fd() );
		}
	}
	catch( SocketException& e )
	{
		perror( e.call() );
	}
	catch( ParserException& e )
	{
		std::ostringstream oss;
		oss << e;

		std::cout << "Caught exception parsing remote scene description: ";
		std::cout << e;
		std::cout << std::endl;
	}
	catch( std::bad_alloc& )
	{
		std::cout << "Error: invalid size received." << std::endl;
	}

	std::cout << "Dropping connection." << std::endl;
}

planning::ScenePtr getScene( const std::string& sceneDesc )
{
	// Okay, we have the entire scene description
	// total hack
	// @todo fix this so it doesn't need to write to file
	// also this is not thread-safe
	std::string tmpFilename( "/tmp/tmpScene.txt" );
	{
		std::ofstream tempFile( tmpFilename.c_str() );
		assert( tempFile );
		tempFile << sceneDesc;
	}
	planning::ScenePtr scene = planning::parseSimulation( tmpFilename );
	return scene;
}

void sendSceneError( int parent, const std::string& error )
{
	int bufid = pvm_initsend( PvmDataRaw );
	if( bufid < 0 )
	{
		PVM_PERROR(pvm_initsend);
		return;
	}

	std::string localError = error;
	pvm_pack( localError );

	int info = pvm_send( parent, SCENE_ERROR );
	if( info < 0 )
	{
		PVM_PERROR(pvm_send);
		return;
	}
}

void pvmClient(int parent)
{
	RandomStream random;

	// first need to grab the scene off the wire
	int buf_id = pvm_recv(parent, INIT_DATA);
	assert( buf_id != PvmDataInPlace );
	if( buf_id == PvmSysErr )
	{
		std::cout << "PVM server failure." << std::endl;
		PVM_PERROR(pvm_recv);
		return;
	}

	unsigned int myId;
	pvm_unpack( myId );
#ifdef VERBOSE
	std::cout << "Recv id: " << myId << std::endl;
#endif

	unsigned int frameRate;
	pvm_unpack( frameRate );
#ifdef VERBOSE
	std::cout << "Recv frame rate: " << frameRate << std::endl;
#endif

	std::string sceneDesc;
	pvm_unpack( sceneDesc );
#ifdef VERBOSE
	std::cout << "Recv scene desc: size " << sceneDesc.size() << std::endl;
#endif

	std::string parentSceneDesc;
	pvm_unpack( parentSceneDesc );
#ifdef VERBOSE
	std::cout << "Recv parent scene desc: size " << parentSceneDesc.size() << std::endl;
#endif

	std::vector<unsigned int> fixedObjects;
	pvm_unpack( fixedObjects );
#ifdef VERBOSE
	std::cout << "Recv fixed objects" << std::endl;
#endif

	float startTime;
	pvm_unpack( startTime );
#ifdef VERBOSE
	std::cout << "Recv start time" << std::endl;
#endif

	float endTime;
	pvm_unpack( endTime );
#ifdef VERBOSE
	std::cout << "Recv end time" << std::endl;
#endif

	unsigned int nPaths;
	pvm_unpack( nPaths );
#ifdef VERBOSE
	std::cout << "Recv npaths" << std::endl;
#endif

	std::vector<PiecewisePath> fixedPaths;
	for( size_t i = 0; i < nPaths; ++i )
	{
		std::vector<unsigned int> positionTimes;
		std::vector<float> positionCoeffs;
		std::vector<unsigned int> rotationTimes;
		std::vector<float> rotations;
		std::vector<float> rotationCoeffs;

		std::cout << "in pvm: positionTimes: " << positionTimes.size() << "; positionCoeffs: " << positionCoeffs.size() << std::endl;

		pvm_unpack( positionTimes );
		pvm_unpack( positionCoeffs );
		pvm_unpack( rotationTimes );
		pvm_unpack( rotations );
		pvm_unpack( rotationCoeffs );

		fixedPaths.push_back( 
			PiecewisePath( 
				positionTimes, 
				positionCoeffs, 
				rotationTimes, 
				rotations, 
				rotationCoeffs ) );
	}
#ifdef VERBOSE
	std::cout << "Recv fixedPaths" << std::endl;
#endif


	boost::shared_ptr<planning::SimulationTree> tree;
	try
	{
		planning::ScenePtr scene = getScene( sceneDesc );
#ifdef VERBOSE
		std::cout << "created scene." << std::endl;
#endif

		if( !parentSceneDesc.empty() )
		{
#ifdef VERBOSE
		std::cout << "Creating parent tree... " << std::flush;
#endif
			planning::ScenePtr parentScene = getScene( parentSceneDesc );
			boost::shared_ptr<planning::SimulationTree> parentTree(
				new SimulationTree("parentTree", parentScene, frameRate) );
#ifdef VERBOSE
		std::cout << "done. " << std::flush;
#endif

			std::sort( fixedObjects.begin(), fixedObjects.end() );
			std::vector<ConstPhysicsObjectPtr> dynamicObjects = parentTree->dynamicObjects();
			std::deque<ConstPhysicsObjectPtr> activeObjectPtrs;
			for( size_t i = 0; i < dynamicObjects.size(); ++i )
			{
				if( !std::binary_search(fixedObjects.begin(), fixedObjects.end(), i) )
					activeObjectPtrs.push_back( dynamicObjects.at(i) );
			}
	
			tree.reset(
				new SimulationTree("childCondorSampling",
					scene,
					parentTree,
					activeObjectPtrs,
					fixedPaths,
					parentTree->samplingRate(),
					startTime ) );
		}
		else
		{
			tree.reset( 
				new SimulationTree("condorSampling", scene, frameRate) );
		}
#ifdef VERBOSE
		std::cout << "create SimulationTree" << std::endl;
#endif
	}
	catch( ParserException& e )
	{
		std::ostringstream oss;
		oss << e;
		sendSceneError( parent, oss.str() );
		return;
	}
	catch( Exception& e )
	{
		sendSceneError( parent, e.message() );
		return;
	}

	tree->setEndTime( endTime );

	for( size_t i = 0; i < 200 && !die; ++i )
	{
#ifdef VERBOSE
		std::cout << "Generating sample... " << std::flush;
#endif

		Path* path = 
			tree->sample(
				boost::numeric::bounds<size_t>::highest(),
				0,
				random );

#ifdef VERBOSE
		std::cout << "Done." << std::endl;
#endif

		// now send it back over the wire
		try
		{
#ifdef VERBOSE
			std::cout << "Calling pvm_initsend... " << std::flush;
#endif
			int bufid = pvm_initsend( PvmDataRaw );
			if( bufid < 0 )
			{
				PVM_PERROR(pvm_initsend);
				break;
			}
#ifdef VERBOSE
			std::cout << "Done" << std::endl;
#endif

#ifdef VERBOSE
			std::cout << "Packing result... " << std::flush;
#endif
			unsigned int size = path->treeCount();
			pvm_pack( size );
			for( size_t i = 0; i < path->treeCount(); ++i )
			{
				CompressedPath compressed = path->path( i );
				pvm_pack( compressed.positionTimes );
				for( size_t j = 0; j < 3; ++j )
					pvm_pack( compressed.positionLinearCoeffs[j] );
				pvm_pack( compressed.positionQuadraticCoeffs[1] );

				pvm_pack( compressed.rotationTimes );
				for( size_t j = 0; j < 3; ++j )
					pvm_pack( compressed.rotationCoeffs[j] );

				std::vector< vl::Vec3f > points = path->obbTree(i).points();
				pvm_pack( points );
			}

			// now need to pack the snapshots
			unsigned int numSnapshots = path->snapshotCount();
			pvm_pack( numSnapshots );
			for( size_t iSnapshot = 0; iSnapshot < numSnapshots; ++iSnapshot )
			{
				float snapshotTime = path->snapshotTime(iSnapshot);
				pvm_pack( snapshotTime );

				BOOST_STATIC_ASSERT( sizeof(RigidStaticState) == (4*7) );
				std::vector<RigidStaticState> state = path->snapshot( iSnapshot );
				unsigned int numStates = state.size();
				pvm_pack( numStates );

				int info = pvm_pkfloat( reinterpret_cast<float*>(&state[0]), 7*state.size(), 1 );
				if( info < 0 )
					throw PVMException( "pvm_pkfloat", info );
			}
#ifdef VERBOSE
			std::cout << "Sent " << numSnapshots << " snapshots." << std::endl;
#endif
			std::vector<float> metricScores = path->allScores();
			pvm_pack( metricScores );
#ifdef VERBOSE
			std::cout << "Sent " << metricScores.size() << " scores. " << std::endl;
#endif

#ifdef VERBOSE
			std::cout << "Done" << std::endl;
#endif

#ifdef VERBOSE
			std::cout << "Sending result... " << std::flush;
#endif
			{
				int info = pvm_send( parent, RESULT_TAG+myId );
				if( info < 0 )
				{
					PVM_PERROR(pvm_send);
					break;
				}
			}
#ifdef VERBOSE
			std::cout << "Done" << std::endl;
#endif
		}
		catch( PVMException& e )
		{
			PVM_PERROR(pvm_unpack);
			continue;
		}
	}
}

} // namespace planning


