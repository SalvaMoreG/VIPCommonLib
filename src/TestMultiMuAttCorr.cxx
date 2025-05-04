
const int NCOLS(10);
const int NROWS(10);
const double BINSIZE(1.0); // in X and Y

const double halfSize_in = 2;   // mm
const double halfSize_out = 4;  // mm

const double mu_in  = 0.00968; // mm^-1
const double mu_out = 0.0173184; // mm^-1

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


double getMuValue( int in_index );
bool IsInSquare(const double& in_x, const double& in_y, const double in_halfSize);
int getIndex(int icol, int irow);
void getXandY( int in_index, double& out_x, double& out_y );
void getColRow(int in_index, int& out_col, int& out_row);
double getDistance(int firstIdx, int lastIdx);
void FindIndicesEachRegion( const std::vector<int> in_idxVec);

// ==========================================

int main()
{
	std::vector<int> idxVec_1{ 40, 41, 42, 43, 44, 45, 46, 47, 48, 49 };
    FindIndicesEachRegion( idxVec_1 );
    
    cout << "  -----------" << endl;
    
	std::vector<int> idxVec_2{ 4, 14, 24, 34, 44, 54, 64, 74, 84, 94 };
    FindIndicesEachRegion( idxVec_2 );

    cout << "  -----------" << endl;
    
	std::vector<int> idxVec_3{ 0, 11, 22, 33, 44, 55, 66, 77, 88, 99 };
    FindIndicesEachRegion( idxVec_3 );
    
    cout << "  -----------" << endl;
    
	std::vector<int> idxVec_4{ 99, 88, 77, 66, 55, 44, 33, 22, 11, 0 };
    FindIndicesEachRegion( idxVec_4 );
   
    cout << "  -----------" << endl;
    
	std::vector<int> idxVec_5{ 49, 58, 67, 76, 85, 94 };
    FindIndicesEachRegion( idxVec_5 );
    
    cout << "  -----------" << endl;
    
	std::vector<int> idxVec_6{ 39, 48, 57, 66, 75, 84, 93 };
    FindIndicesEachRegion( idxVec_6 );
}

// ==========================================

double getMuValue( int in_index )
{
    double x, y;
    getXandY( in_index, x, y );
    double mu(0.0);
    if ( IsInSquare(x, y, halfSize_in) )
        mu = mu_in;
    else if ( IsInSquare(x, y, halfSize_out) )
        mu = mu_out;
    
    return mu;
}

// ==========================================

bool IsInSquare(const double& in_x, const double& in_y, const double in_halfSize)
{
	return ( fabs(in_x) < in_halfSize && fabs(in_y) < in_halfSize );
}

// ==========================================

int getIndex(int icol, int irow)
{
	int index = irow * NCOLS + icol;
	return index;
}

// ==========================================

void getXandY( int in_index, double& out_x, double& out_y )
{
	int icol, irow;
	getColRow(in_index, icol, irow);
    // center of voxel coordinates
	out_x = (-0.5 * (double) (NCOLS) + (double) (icol) + 0.5) * BINSIZE;
	out_y = (-0.5 * (double) (NROWS) + (double) (irow) + 0.5) * BINSIZE;
    // cout << "idx: " << in_index << " x, y: " << out_x << " " << out_y << endl;
}

// ==========================================

void getColRow(int in_index, int& out_col, int& out_row)
{
	out_row = (int) (in_index/NCOLS);
	out_col = in_index % NCOLS;
}

// ==========================================

double getDistance(int in_firstIdx, int in_lastIdx)
{
    double x0, y0;
    getXandY( in_firstIdx, x0, y0 );
    double xend, yend;
    getXandY( in_lastIdx, xend, yend );
    
    double distance = std::sqrt((xend-x0)*(xend-x0) + (yend-y0)*(yend-y0)) + BINSIZE;
    return distance;
}

// ==========================================

void FindIndicesEachRegion( const std::vector<int> in_idxVec)
{
	int firstindex(-1), lastindex(-1);
    int numvoxels(0);
    
    vector<double> muValues;
    double muValue(0.0), prevMu(0.0);
    vector<double> distances;
    double distance(0.0);
    
    for ( unsigned int kiter = 0; kiter < in_idxVec.size(); kiter++ )
    {
        unsigned int index = in_idxVec[kiter]; 
        muValue = getMuValue( index );
        if (muValue > 0.0001 && muValue != prevMu)
        {
            if (muValue != prevMu)
            {
                // cout << "previous mu: " << prevMu << " last index: " << lastindex << endl;
                // cout << "   new mu: " << muValue << " first index: " << index << endl;
                if (prevMu > 0.0001)
                {
                    muValues.push_back(prevMu);
                    distance = getDistance(firstindex, lastindex);
                    distances.push_back(distance);
                    cout << "attenuated indices : [" 
                        << firstindex << ";" << lastindex 
                        << "] mu: " << prevMu
                        << " D: " << distance << endl;
                }
                
                firstindex = index;
            }
        }
        if (muValue > 0.0001)
        {
            lastindex = index;
            prevMu = muValue;
            numvoxels++;
        }
    }
    
    if (lastindex != -1)
    {
        // cout << "previous mu: " << prevMu << " last index: " << lastindex << endl;
        muValues.push_back(prevMu);
        distance = getDistance(firstindex, lastindex);
        distances.push_back(distance);
        cout << "attenuated indices : [" 
             << firstindex << ";" << lastindex 
             << "] mu: " << prevMu
             << " D: " << distance << endl;
    }
    
    cout << endl;
    cout << "indices: ";
    for(int i=0; i < in_idxVec.size(); i++)
        cout << in_idxVec[i] << " ";
    cout << endl;
    
    for (int id = 0; id < distances.size(); id++)
    {
        distance = distances[id];
        cout << "mu: " << muValues[id]
             << " D: " << distance << endl;
    }
    cout << "total #attenuated_voxels: " << numvoxels << endl;
}

// ==========================================


