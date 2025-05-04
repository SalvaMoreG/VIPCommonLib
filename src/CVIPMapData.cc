
#include "CVIPMapData.h"
#include <string>

#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdlib>

using namespace std;


template <>
void
CVIPMapData<int>::TestSpecialized() const
{
	std::cout << "YES! This is specialized for integer: " << std::endl;
}

template <>
void
CVIPMapData<unsigned int>::TestSpecialized() const
{
	std::cout << "YES! This is specialized for unsigned integer: " << std::endl;
}

template <>
void
CVIPMapData<unsigned int>::Serialize()
{
	if (m_maxValue < USHRT_MAX)
	{
		DoSerialize<unsigned short>();
	}
	else
	{
		DoSerialize<unsigned int>();
	}
}

template <>
void
CVIPMapData<int>::Serialize()
{
	if (m_maxValue < USHRT_MAX)
	{
		DoSerialize<unsigned short>();
	}
	else
	{
		DoSerialize<unsigned int>();
	}
}

template <>
void
CVIPMapData<unsigned int>::Deserialize( int in_setNumber )
{
	if (m_maxValue < USHRT_MAX)
	{
// cout << "unsigned short: " << endl;
		DoDeserialize<unsigned short>( in_setNumber );
	}
	else
	{
// cout << "unsigned int: " << endl;
		DoDeserialize<unsigned int>( in_setNumber );
	}
}

template <>
void
CVIPMapData<int>::Deserialize( int in_setNumber )
{
	if (m_maxValue < USHRT_MAX)
		DoDeserialize<unsigned short>( in_setNumber );
	else
		DoDeserialize<unsigned int>( in_setNumber );
}


