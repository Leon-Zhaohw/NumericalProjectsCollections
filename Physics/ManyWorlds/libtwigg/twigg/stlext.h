#ifndef __STLEXT_H__

#define __STLEXT_H__

#include <algorithm>

#ifdef WIN32
#include <hash_set>
#include <hash_map>
#else
#include <ext/algorithm>
#include <ext/functional>
#include <ext/hash_set>
#include <ext/hash_map>
#endif

#ifdef WIN32
#define STLEXT stdext::
#else
#define STLEXT __gnu_cxx::
#endif


namespace std
{

#ifdef WIN32



template <typename ForwardIterator>
bool is_sorted( const ForwardIterator first, const ForwardIterator last )
{
	if( first == last )
		return true;

	ForwardIterator prev = first;
	ForwardIterator next = first;
	for( ++next; next != last; prev = next, ++next )
	{
		if( std::less<ForwardIterator::value_type>()( *next, *prev ) )
			return false;
	}

	return true;
}

template <class _Pair>
struct select1st : public std::unary_function<_Pair, typename _Pair::first_type> {
  typename _Pair::first_type& operator()(_Pair& __x) const {
    return __x.first;
  }
  const typename _Pair::first_type& operator()(const _Pair& __x) const {
    return __x.first;
  }
};

template <class _Pair>
struct select2nd : public std::unary_function<_Pair, typename _Pair::second_type>
{
  typename _Pair::second_type& operator()(_Pair& __x) const {
    return __x.second;
  }
  const typename _Pair::second_type& operator()(const _Pair& __x) const {
    return __x.second;
  }
};
#else

using __gnu_cxx::select1st;
using __gnu_cxx::select2nd;

using __gnu_cxx::is_sorted;



/*
template <typename ForwardIterator>
bool is_sorted( const ForwardIterator first, const ForwardIterator last )
{
	return __gnu_cxx::is_sorted(first, last);
}
*/

#endif

}

#ifndef WIN32
namespace __gnu_cxx
{
template<> struct hash<std::string>
{
	size_t operator()( const std::string& __s) const
	{
		return __stl_hash_string(__s.c_str());
	}
};
}
#endif

#endif

