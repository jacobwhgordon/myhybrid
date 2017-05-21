//include file for hybrid functions

#include <iostream>
#include <vector>
#include <complex> 
#include <fstream>
#include <sstream>
#include <string>

#include "TH1.h"

//AnitaTools libs
#include "../../AnitaTools/include/FFTtools.h"

using namespace std;
using namespace FFTtools;

class myHybrid
{

//Functions from myHybridFunctions.cc
//I guess we need them to be public? and say that.

public:

template <typename T> int sgn(T val);
vector<double> unwrapFunction(vector<double> );
vector<vector<double> > unwrapFunction(vector<vector<double> > );
vector<vector<double> > matchFunction(vector<vector<double> >, vector<vector<double> >) ;
vector<vector<complex<double> > > matchFunction(vector<vector<complex<double> > >, vector<vector<complex<double> > > ); 
vector< vector< complex<double> > > extendS21(vector< vector< complex< double > > >, vector< vector< complex< double > > >, int);
//vector< vector< complex<double> > > extendS21double(vector< vector< complex< double > > >, vector< vector< complex< double > > >, int);
vector<vector<double> > pulsePadder ( vector< vector<double> > , int, int);
vector<vector<double > > myHybridSim(vector<vector<complex<double> > >, vector<vector<complex<double> > >, vector<vector<complex<double> > >);
vector<vector<double > > myHybridSim(vector<vector<complex<double> > >, vector<vector<complex<double> > >, vector<vector<complex<double> > >, vector<vector<complex<double> > >);
vector<vector<double> > timePulse ( vector< vector<double> >);
vector<vector<complex<double> > > freqPulse( vector< vector<double> >);
vector<vector<double> > getdBFunc(vector< vector< complex<double> > > );
vector<vector<double> > getPhaseFunc( vector< vector< complex<double > > >);
vector<vector<double> > getNAData( string );
vector<vector<double> > getPulseData( string );
vector<vector<complex<double> > > getCorrectionData( string );
vector< vector< complex< double > > >  rawToCalibReIm (vector< vector<double> >,vector< vector<double> >,vector< vector<double> >,vector< vector<double> >);
vector<vector<complex<double> > > FFTWtovector ( FFTWComplex * , int, double );
FFTWComplex * vectortoFFTW (vector<vector<complex<double> > >);
FFTWComplex * vectortoFFTW (vector<complex<double> >);
TH1D * histFunction ( vector<vector<double> >, const char*, const char* );
TH1D * histFunction ( vector<vector<complex<double> > >, const char*, const char*, int );
TH1D * histFunction ( TGraph * , const char*, const char* );
void histShift ( TH1D &, int );
void hybridPathCorrectionFunction ( vector<vector<double> > &, vector<vector<double> > & );

};
















































