#pragma once

#ifndef _VIP_SERIAL_READER___
#define _VIP_SERIAL_READER___

#include <fstream>

// NOTE: this class uses template functions. These cannot be defined in a cc file
// because the compiler wouldn't know where to look for a particular implementation
// (when it's called from somewhere). 
// In other words: during compilation all possible implementations (T = int, T = double, 
// T = string, ...) should be defined so they can be called from other places in the code.
// For this reason, at the end of this header there is this line: 
//   #include "CVIPSerialIO.tpp"

class CVIPSerialIO
{
public:
	// ctors and dtor
				CVIPSerialIO(){};
                ~CVIPSerialIO(){};

	// Read and Write
    // Variadic functions Read and Write
    template <typename T, typename... Types>
    void        SerialRead(std::ifstream& infile, bool is_binary, T& var1, Types&... var2);

    template <typename T, typename... Types>
    void        SerialWrite(std::ofstream& outfile, bool is_binary, const T& var1, const Types&... var2);
    
// protected:    
    // at the end of the line 
    // executed when there are no variables left...
    inline void        SerialRead(std::ifstream& infile, bool is_binary)
                        { if (!is_binary) infile.ignore(1024, '\n'); };  
    inline void        SerialWrite(std::ofstream& outfile, bool is_binary)
                        { if (!is_binary) outfile << std::endl; }; 
    
private:
	// copy and assign constructors, not defined, not used
                    CVIPSerialIO(const CVIPSerialIO& in_obj);
	CVIPSerialIO& 	operator= (const CVIPSerialIO& in_obj);
};

#include "CVIPSerialIO.tpp"

#endif
