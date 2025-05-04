
#ifndef __CC_VIP_MAPDATA_H___
#define __CC_VIP_MAPDATA_H___

#include "CVIPFieldOfView.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cassert>
#include <string>

#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdlib>

#include <climits>

/* =========================================================================================================
This is a complicated code. 
The basic idea is simple: store a map, using a key (type integer) and a value (type whatever). 
So for instance: map<int room-number, float room-size> or map<int room-number, int number-of-windows)
For the "whatever" type (int, float, ...), we use a template class "CVIPMapData". 
So, an object of type CVIPMapData<double> internally will have a map of type <int, std::vector<double> >. 

So far so good.

Additionally (!!!), this class also uses 2 additional template functions "DoSerialize()" and "DoDeserialize()".
We only do this because in case the entire object is of type "CVIPMapData<int>", (so all values are of type integer), 
	we want to write the output either as "unsigned int" or "unsigned short" (in order to have smaller storage files).
The output can only be of type "unsigned short" if the maximum value for the values in the map is smaller than USHRT_MAX.
USHRT_MAX = 65535, defined in the "climits" include file.
So: for this optimization, you have to use the Set value and make sure the 3th parameter is given.
For instance, for LM-OSEM, the FOV-voxel index is always smaller than the total number of FOV voxels. 
So, if the FOV->GetNumberOfVoxels() is smaller than USHRT_MAX, the output can be written in "unsigned short", 
	resulting in filesizes of half the size compared with "unsigned int".

Final warning: 	"Templates are beautiful, but also rather difficult"  (unknown, 23 B.C.)
============================================================================================================= */

template <class TType>
class CVIPMapData 
{
public:
					CVIPMapData();
					CVIPMapData( int in_setSize, const std::string& in_basefilename );
	virtual			~CVIPMapData();

	// Set function
	void			Set( int in_setSize, const std::string& in_basefilename, int in_maxValue = USHRT_MAX+1 );

	// Add function
	// The name "Add" might be misleading. It adds a new (!) entry with key "in_key" and value(s) "in_value(s)".
	// It does NOT add new values to an existing entry. 
	bool			Add( int in_key, const std::vector<TType>& in_values ); // in_key gets a vector of values
	bool			Add( int in_key, const TType& in_value );	// in_key gets a single value
	
	inline void		Finalize() { Serialize(); }
	
	// Get functions cannot be const because they do I/O functions, clearing and filling internal map
	bool			Get( int in_key, std::vector<TType>& io_values );
	bool			Get( int in_key, TType& io_value );	// get first element in the vector	
	
	inline int		GetCurrentSize() const { return m_map.size(); } 
	inline int		GetCurrentSetNumber() const { return m_currentSetNumber; } 	

	void			TestSpecialized() const;		// TESTING ONLY

private:
	// copy and assign constructors, not defined, not used
					CVIPMapData(const CVIPMapData& in_obj);
	CVIPMapData& 	operator= (const CVIPMapData& in_obj);

	// aux.
	void			Serialize();
	void			Deserialize( int in_setNumber );

	template<typename TF>
	void 			DoSerialize();
	template<typename TF>
	void 			DoDeserialize( int in_setNumber );


	int				GetSetNumber(int in_key) const;
	void			GetSetFileName(std::string& io_fname, int in_nr) const;

	// data
	typedef std::map<int, std::vector<TType> >  CVIPMapDataMap;
	CVIPMapDataMap					m_map;
	int								m_setSize;
	int 							m_currentSetNumber;
	std::string 					m_basefilename;
	int								m_maxValue;
};

// #include "CVIPMapData.inl"


template <class TType>
CVIPMapData<TType>::CVIPMapData()
	: m_setSize(0)		
	, m_currentSetNumber(-1)
	, m_basefilename("")
	, m_maxValue(USHRT_MAX+1)
{
	m_map.clear();
}

template <class TType>
CVIPMapData<TType>::CVIPMapData( int in_setSize, const std::string& in_basefilename )
	: m_setSize(in_setSize)		
	, m_currentSetNumber(-1)
	, m_basefilename( in_basefilename )
	, m_maxValue(USHRT_MAX+1)
{
	m_map.clear();
}

template <class TType>
CVIPMapData<TType>::~CVIPMapData()
{
	m_map.clear();
}

template <class TType>
void
CVIPMapData<TType>::Set( int in_setSize, const std::string& in_basefilename, int in_maxValue )
{
	m_setSize = in_setSize;		
	m_basefilename = in_basefilename;
	m_maxValue = in_maxValue;
}

template <class TType>
bool
CVIPMapData<TType>::Add( int in_key, const std::vector<TType>& in_values )
{
	// make sure this is a new entry....
	if ( m_map.find(in_key) != m_map.end() )
	{
		std::cout << "ERROR! Entry already stored. in_key: " << in_key << std::endl;

		typename CVIPMapDataMap::const_iterator iter = m_map.find(in_key);

		std::cout << "a.k.a.: " << iter->first << std::endl;
		return false;
	}
	assert( m_map.find(in_key) == m_map.end() );

	// load correct set into memory map
	int	setNumber = GetSetNumber( in_key );
	if (setNumber != m_currentSetNumber)
	{
		// Serialize previous set
		Serialize();

		// Clean current map
 		m_map.clear();

		// Change "current set number"
		m_currentSetNumber = setNumber;
	}

	// Add vector to new map entry
	m_map[in_key] = in_values;
	return true;
}

template <class TType>
bool
CVIPMapData<TType>::Add( int in_key, const TType& in_value )
{
	std::vector<TType> dumVector;
	dumVector.push_back( in_value );
	return Add( in_key, dumVector );
}

template <class TType>
bool
CVIPMapData<TType>::Get( int in_key, std::vector<TType>& io_values )
{
	// get the correct set
	int	setNumber = GetSetNumber( in_key );
	if (setNumber != m_currentSetNumber)
	{
		// Clean current map
 		m_map.clear();

		// Deserialize new vector map
		Deserialize( setNumber );

		// Change "current set number"
		m_currentSetNumber = setNumber;
	}

	// Get the key entry from the current map
	typename CVIPMapDataMap::const_iterator iter = m_map.find(in_key);
	if ( iter == m_map.end() )
	{
		std::cout << "ERROR, failed getting key: " << in_key << std::endl;
		std::cout << "setNumber: " << setNumber << " m_currentSetNumber: " << m_currentSetNumber << std::endl;
	}
	assert ( iter != m_map.end() );		// in list...

	io_values = (*iter).second;
	int key = (*iter).first;

	assert( in_key == key );

	return true;
}

template <class TType>
bool
CVIPMapData<TType>::Get( int in_key, TType& io_value )
{
	std::vector<TType> values;
	if ( Get( in_key, values ) && values.size() > 0)
	{
		io_value = values[0];
		return true;
	}
	return false;
}

template <class TType>
int
CVIPMapData<TType>::GetSetNumber( int in_key ) const
{
	int setNumber = (m_setSize == 0) ? 0 : in_key / m_setSize;
	return setNumber;
}

template <class TType>
void
CVIPMapData<TType>::GetSetFileName(std::string& io_fname, int in_nr) const
{
	io_fname = m_basefilename; 			// "eventToFovMap_";
	std::stringstream setNrStr;
	setNrStr << in_nr;
	io_fname += setNrStr.str();
	io_fname += ".dat";
}

// <<<<<<<<<<<<<<<<<<<<<< TESTING ONLY
// Explicit instantiation:
// https://msdn.microsoft.com/en-us/library/by56e477.aspx
template <>
void
CVIPMapData<unsigned int>::TestSpecialized() const;

template <>
void
CVIPMapData<int>::TestSpecialized() const;

template <class TType>
void
CVIPMapData<TType>::TestSpecialized() const
{
	std::cout << "This is not specialized (TType): " << std::endl;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// Explicit instantiation:
// https://msdn.microsoft.com/en-us/library/by56e477.aspx
template <>
void CVIPMapData<unsigned int>::Serialize();

template <>
void CVIPMapData<int>::Serialize();

template <class TType>
void CVIPMapData<TType>::Serialize()
{
	// std::cout << "CVIPMapData<TType>, calling TType (float)" << std::endl;
	DoSerialize<TType>();
}

// Explicit instantiation:
// https://msdn.microsoft.com/en-us/library/by56e477.aspx
template <>
void CVIPMapData<unsigned int>::Deserialize( int in_setNumber );

template <>
void CVIPMapData<int>::Deserialize( int in_setNumber );

template <class TType>
void CVIPMapData<TType>::Deserialize( int in_setNumber )
{
	DoDeserialize<TType>( in_setNumber );
}


template <class TType>		// type of class (e.g. unsigned int)
template<typename TF>		// type of method (e.g. short or int)
void CVIPMapData<TType>::DoSerialize()
{
	if (m_map.size() == 0) return;

	std::string fname;
	GetSetFileName(fname, m_currentSetNumber);

	std::ofstream outfile(fname.c_str(), std::ios::binary | std::ios::out );
	assert (outfile.is_open());

	unsigned int size = m_map.size();
	unsigned int key; 
	unsigned int vecsize;
	TF value;
	// WRITE
	outfile.write((char*) &size, sizeof(size));

	typename CVIPMapDataMap::const_iterator iter;
	for (iter = m_map.begin(); iter != m_map.end(); iter++)
	{
		key = (*iter).first;
		// WRITE
		outfile.write((char*) &key, sizeof(key));

		std::vector<TType> vec = (*iter).second;
		vecsize = vec.size();
		// WRITE
		outfile.write((char*) &vecsize, sizeof(vecsize));	

		for (int idx = 0; idx < vecsize; idx++)
		{
			value = vec[idx];
			// WRITE
			outfile.write((char*) &value, sizeof(value));	
		}
	}

	outfile.flush();
	outfile.close();
}

template <class TType>
template<typename TF>
void CVIPMapData<TType>::DoDeserialize( int in_setNumber )
{
	std::string fname;
	GetSetFileName(fname, in_setNumber);

	std::ifstream infile(fname.c_str(), std::ios::binary | std::ios::in );
	if ( !infile.is_open() )
	{
		// No previous data available. Start with empty map...
		std::cout << "ERROR, file is not open! " << fname << std::endl;
		exit(1);
	}

	unsigned int size; 
	unsigned int key; 
	unsigned vecsize;	
	TF value;	

	// READ
	infile.read((char*) &size, sizeof(size));
	std::vector<TType> vec;
	for (int iloop = 0; iloop < size; iloop++)
	{
		// READ
		infile.read((char*) &key, sizeof(key));
		infile.read((char*) &vecsize, sizeof(vecsize));

		vec.clear();
		for (int idx = 0; idx < vecsize; idx++)
		{
			// READ
			infile.read((char*) &value, sizeof(value));
			vec.push_back( value );
		}
		m_map[key] = vec;
	}
	infile.close();
}

#endif
