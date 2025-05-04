
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>	// <---- modern C++ smart pointers

using namespace std;

int main()
{
	cout << "give the list of files: ";
	string listfname; cin >> listfname;
	ifstream inlist( listfname.c_str() );
	if ( !inlist.is_open() )
	{
		cout << "ERROR! List file not open: " << listfname << endl;
		return -1;	
	}

	// Collect all streams
	std::vector< shared_ptr<ifstream> > streams;
	string fname;
	while ( !inlist.eof() )
	{
		inlist >> fname;
		if ( !inlist.eof() )
		{
			auto infile_ptr = std::make_shared<ifstream>( fname.c_str() );
			if ( !infile_ptr->is_open() )
			{
				cout << "ERROR! file not open: " << fname << endl;
				return -1;	
			}
			streams.push_back(infile_ptr);
		}
	}

	// create outfile
	ofstream outfile("MIXED.dat_NEW");

	// Loop over streams
	double z1, y1, x1, E1;
	double z2, y2, x2, E2;
	bool all_eof( false );
	while ( !all_eof )
	{
		all_eof = true;
		for (auto stream_ptr: streams)
		{
			ifstream& stream = *stream_ptr;
			stream  >> z1 >> y1 >> x1 >> E1 
		        	>> z2 >> y2 >> x2 >> E2;
			if ( !stream.eof() )
			{
				all_eof = false;
				outfile << z1 << " " << y1 << " " << x1 << " " << E1 
						<< "    "
			        	<< z2 << " " << y2 << " " << x2 << " " << E2 << endl;
			}
		}	
	}	
	streams.clear();	// here all shared_ptrs will be deleted automatically (cause nobody is using them anymore)

	/*	// No longer necessary to delete anything
	if (all_eof)
	{
		cout << "ready for destruction" << endl;
		for (auto* stream_ptr: streams)
			delete stream_ptr;
		streams.clear();
	}
	*/

	return 0;
}


