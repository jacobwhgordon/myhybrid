//include file for hybrid functions

#include <iostream>
#include <vector>
#include <complex> 
#include <fstream>
#include <sstream>
#include <string>

#include "myHybridFunctions.h"

#include "TH1.h"

//AnitaTools libs
#include "FFTtools.h"

using namespace std;
using namespace FFTtools;


const double PI  =3.141592653589793238463;



/////////////////////////////////////////////////////////////////////////////////
/////sgn function ///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//pulled from internet.
template <typename T> int myHybrid::sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}


/////////////////////////////////////////////////////////////////////////////////
/////Unwrap function ////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

vector<double> myHybrid::unwrapFunction(vector<double> data)
{
    int loops = 0;
    vector<double> out;
    out.push_back(data[0]);                           //first data point will have no loops.
    for (int i = 1; i < data.size(); i++)
    {
        if ( abs(data[i] - data[i-1]) > 90 )
        {
            if ( sgn(data[i] - data[i-1]) == 1)
            {loops = loops-1;}
            else if ( sgn(data[i] - data[i-1]) == -1)
            {loops = loops+1;}
        }
        out.push_back(data[i]+loops*360);
    }
    return out;
}

vector<vector<double> > myHybrid::unwrapFunction(vector<vector<double> > data)
{
    int loops = 0;
    vector<vector<double> > out;
    out.push_back(data[0]);                           //first data point will have no loops.
    for (int i = 1; i < data.size(); i++)
    {
        if ( abs(data[i][1] - data[i-1][1]) > 90 )
        {
            if ( sgn(data[i][1] - data[i-1][1]) == 1)
            {loops = loops-1;}
            else if ( sgn(data[i][1] - data[i-1][1]) == -1)
            {loops = loops+1;}
        }
        vector<double> temp;
        temp.push_back(data[i][0]);
        temp.push_back(data[i][1]+loops*360);
        out.push_back(temp);
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////Match function /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//matches the x-axis of data 2 to data 1
//NOTE!!! WE NEED TO DO SOMETHING TO MAKE SURE out IS THE SAME SIZE AS data1
vector<vector<double> > myHybrid::matchFunction(vector<vector<double> > data1, vector<vector<double> > data2) 
{
    double data1Spacing = data1[1][0]-data1[0][0];
    int iEnd = data1.size();                                                     //we want a data 2 point at each data1 point
    //int iEnd = data2[data2.size()-1][0]/data1Spacing;                          //only needed if... something... idk
    vector<vector<double> > out;                                                 //I think it was incase data2 had more points than data1?
    for( int i=0;i<iEnd; i++ )
    {
        int index = 0; 
        double data1here = data1[0][0] + i*data1Spacing;
        //find points above and below current data2 position
        for (int j=0;j < data2.size()-1; j++)
        {
            if( data2[j][0] <= data1here &&  data1here <= data2[j+1][0])
            {
                index = j;
                break;
            }
        }
        double m = (data2[index+1][1]-data2[index][1]) / (data2[index+1][0]-data2[index][0]);
        double b = (data2[index][1] - m*data2[index][0]);
        vector < double > temp;
        temp.push_back(data1here);
        temp.push_back(m*data1here + b);
        out.push_back(temp);
    }
    return out;
}

vector<vector<complex<double> > > myHybrid::matchFunction(vector<vector<complex<double> > >data1, vector<vector<complex<double> > > data2)
{
    double data1Spacing = data1[1][0].real()-data1[0][0].real();
    int iEnd = data1.size();                                                     //we want a data 2 point at each data1 point
    //int iEnd = data2[data2.size()-1][0].real()/data1Spacing;                   //only needed if... something... idk
    vector<vector<complex<double> > > out;
    for( int i=0;i<iEnd; i++ )                                                   //iEnd could just be data1.size()
    {
        int index = 0; 
        double data1here = data1[0][0].real() + i*data1Spacing;                  //data1here could just be data1[i][0]
        //find points above and below current data2 position                     //why did i do it this way...
        for (int j=0;j < data2.size()-1; j++)
        {
            if( data2[j][0].real() <= data1here && data1here <= data2[j+1][0].real() )
            {
                index = j;
                break;
            }
        }
        double mReal = (data2[index+1][1].real()-data2[index][1].real()) / (data2[index+1][0].real()-data2[index][0].real());
        double mImag = (data2[index+1][1].imag()-data2[index][1].imag()) / (data2[index+1][0].real()-data2[index][0].real());
        double bReal = (data2[index][1].real() - mReal*data2[index][0].real());
        double bImag = (data2[index][1].imag() - mImag*data2[index][0].real());
        vector<complex<double> > temp;
        temp.push_back(complex<double>(data1here,0));
        temp.push_back(complex<double>(mReal*data1here + bReal,mImag*data1here + bImag));
        out.push_back(temp);
    }
    //cout << "size of data1: " << data1.size() << endl;
    //cout << "size of data2: " << data2.size() << endl;
    //cout << "size of iEnd: " << iEnd << endl;
    //cout << "size of out: " << out.size() << endl;
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////Extend function ////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//function written to extend out the frequency domain data taken to match the 
//frequency's seen in the pulse FFT
//if fill = 0, fill with zeros, 
//if fill = 1, Fill with freq and phase matched sine function in logistic function envolope. for smoothness.
vector< vector< complex<double> > > myHybrid::extendS21(vector< vector< complex< double > > > data1, vector< vector< complex< double > > > data2, int fill)
{
    double endPoint = data2[data2.size()-1][0].real();
    double data1Spacing = data1[1][0].real() - data1[0][0].real();
    double originalData1End = data1[data1.size()-1][0].real();                     //where data1 ends before being extended
    int data1EndIndex = floor( (endPoint-data1[0][0].real())/data1Spacing )+1;     //the index where data1 should end after being extended
    
    //cout << data1.size() << " " << data2.size() << " " << data1EndIndex << endl;
        
    //we need to find the wavelengths of the re and imag parts, do this by scanning for zeros.
    vector<double> reZeros;
    vector<double> imZeros;
    //we only want to look for zeros in the valid range of our filter, ~100MHz-1.4GHz
    int lowSearchBound = floor( (pow(10,8)-data1[0][0].real())/data1Spacing );
    int highSearchBound = floor( (1.4*pow(10,9)-data1[0][0].real())/data1Spacing );
    for ( int i=lowSearchBound; i<highSearchBound; i++)
    {
        if ( sgn(data1[i][1].real()) != sgn(data1[i+1][1].real()) )
        { reZeros.push_back(data1[i][0].real()); }
        if ( sgn(data1[i][1].imag()) != sgn(data1[i+1][1].imag()) )
        { imZeros.push_back(data1[i][0].real()); }
    }
    //define the wavelengths
    double reLambda = (reZeros[reZeros.size()-1]-reZeros[0])/(0.5*(reZeros.size()-1));
    double imLambda = (imZeros[imZeros.size()-1]-imZeros[0])/(0.5*(imZeros.size()-1));   
     
    //for phase matching we also need to find the LAST zero
    double lastReZero, lastImZero, reShift, imShift;
    for (int i=0; i< data1.size()-1; i++)
    {
        if ( sgn(data1[i][1].real()) != sgn(data1[i+1][1].real()) )
        {
            lastReZero = data1[i][0].real();
            if ( sgn(data1[i+1][1].real()) > 0 )
            { reShift = 0.0*reLambda; } 
            else 
            { reShift = 0.5*reLambda; }
        }
        if ( sgn(data1[i][1].imag()) != sgn(data1[i+1][1].imag()) )
        {
            lastImZero = data1[i][0].real();
            if ( sgn( data1[i+1][1].imag()) > 0 ) 
            { imShift = 0.0*imLambda; } 
            else 
            { imShift = 0.5*imLambda; }
        }
    }
    //now we calculate that the re and imag amplitudes and phase 
    //are at the end of the data 1 data before populating the output.
    double rePhi = (2*PI)/reLambda * (originalData1End-lastReZero+reShift);
    double imPhi = (2*PI)/imLambda * (originalData1End-lastImZero+imShift);
    double reAmp = data1[data1.size()-1][1].real()/sin(rePhi);
    double imAmp = data1[data1.size()-1][1].imag()/sin(imPhi);
    
    //cout << imLambda << " " << lastImZero << " " << imShift << " " << imPhi << " " << imAmp << endl;
    
    //now populate the output
    vector<vector<complex<double> > > out;
    for (int i=0; i<data1EndIndex; i++)
    {
        double data1Here = data1[0][0].real() + i*data1Spacing;
        if (i < data1.size())                                                  //if we have data, use that
        {
            out.push_back(data1[i]);
        }
        else
        {
            if (fill == 0)
            {
                vector<complex<double> > temp;
                temp.push_back(complex<double> (data1Here,0));                //the x coordinate
                temp.push_back(complex<double> (0,0));                        //the y coordinate
                out.push_back(temp);
            } 
            else if (fill == 1)
            {
                double reFunction = reAmp*( 1 - 1/( 1+exp(-10/(endPoint-originalData1End)*(data1Here-originalData1End) + 5) ) )
                                                                     *sin((2*PI/reLambda)*(data1Here-originalData1End) + rePhi);
                double imFunction = imAmp*( 1 - 1/( 1+exp(-10/(endPoint-originalData1End)*(data1Here-originalData1End) + 5) ) )
                                                                     *sin((2*PI/imLambda)*(data1Here-originalData1End) + imPhi);
                //note, this is what i had for mathematica
                //if I want it to decay faster (Like i had in C++ before) then change endPoint->endPoint/2
                //doesnt seem to have a significant effect though.
                vector<complex<double> > temp;
                temp.push_back(complex<double>(data1Here,0));
                temp.push_back(complex<double>(reFunction,imFunction));
                out.push_back(temp);
            }
        }
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////pulsePadder function ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


//this pads data with 0s out to a size of 'size'
//centerType = 0 will center the start of the data at the center of the new window size
//centerType = 1 will center the peak of the data at the center of the new window size
//centerTpye = 2 will put the start of the data at the start of the new window. 
vector<vector<double> > myHybrid::pulsePadder ( vector< vector<double> > data , int size, int centerType)
{
    vector<vector<double> > out;
    int iShift;
    if (centerType == 0)
    { 
        iShift = floor(size/2); 
    } else if (centerType == 1) {        
        //need to find the peak
        double peakI = 0;
        for (int i=0; i<data.size(); i++)
        {
            if ( peakI < data[i][1] )
            { peakI = data[i][1]; }
        }
        iShift = peakI;
    } else if (centerType == 2) { 
        iShift = 0; 
    }
    for (int i=0; i<size ; i++)
    {
        vector<double> temp;
        if ( iShift <= i && i < iShift+data.size() )
        {
            temp.push_back(i*(data[1][0]-data[0][0]));
            temp.push_back(data[i-iShift][1]);
            out.push_back(temp);
        } else {
            temp.push_back(i*(data[1][0]-data[0][0]));
            temp.push_back(0);
            out.push_back(temp);
        }
    }
    return out;
        
}


/////////////////////////////////////////////////////////////////////////////////
/////myHybridSim function ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//for use when the pulse going into A and B (the h and v pulses) are the same
vector<vector<double > > myHybrid::myHybridSim(vector< vector< complex< double > > > pulseData, vector< vector <complex <double> > > path1, vector< vector< complex<double> > > path2)
{
    //cout << "size of pulseData: " << pulseData.size() << endl;
    //cout << "size of path1: " << path1.size() << endl;
    //cout << "size of path2: " << path2.size() << endl;
    double pulseSpacing = pulseData[1][0].real()-pulseData[0][0].real();    //in [Hz]
    vector<complex<double> > outPath1;
    vector<complex<double> > outPath2;
    for ( int i=0; i<pulseData.size(); i++ )
    {
       //cout << pulseData[i][1] << " " << path2[i][1] << endl;
       outPath1.push_back( pulseData[i][1]*path1[i][1] );
       outPath2.push_back( pulseData[i][1]*path2[i][1] );
    }
    //do the inv FFTW
    int outSize = (pulseData.size()-1)*2;
    
    FFTWComplex * outPath1FFT = vectortoFFTW (outPath1);
    double * outPath1FFTTime = doInvFFT(outSize, outPath1FFT); 
    FFTWComplex * outPath2FFT = vectortoFFTW (outPath2);
    double * outPath2FFTTime = doInvFFT(outSize, outPath2FFT);
    
    double outSpacing = 1/(pulseSpacing*outSize)*pow(10,9);                 //should that be outSize,or pulseData.size()
    vector<vector<double> > out;
    for (int i=0; i<outSize ; i++)
    {
        //cout << outPath1FFTTime[i] << " " << outPath2FFTTime[i] << endl;
        vector<double> temp;
        temp.push_back(i*outSpacing);
        temp.push_back(outPath1FFTTime[i]+outPath2FFTTime[i]);              //Fuck... i dont remember, is that a - or a +? 
                                                                            //logic says +, but I seem to remember -?
        out.push_back(temp);
    }
    return out;
}

//for use when the signals going into A/B (the h/v pulses) are different!
vector<vector<double > > myHybrid::myHybridSim(
    vector<vector<complex< double > > > pulse1Data, vector<vector<complex< double > > > pulse2Data, 
    vector<vector<complex <double> > > path1, vector<vector<complex<double> > > path2)
{
    double pulseSpacing = pulse1Data[1][0].real()-pulse1Data[0][0].real();
    vector<complex<double> > outPath1;
    vector<complex<double> > outPath2;
    for ( int i=0; i<pulse1Data.size(); i++ )
    {
       outPath1.push_back( pulse1Data[i][1]*path1[i][1] );
       outPath2.push_back( pulse2Data[i][1]*path2[i][1] );
    }

    //do the inv FFTW
    int outSize = (pulse1Data.size()-1)*2;
    
    FFTWComplex * outPath1FFT = vectortoFFTW (outPath1);
    double * outPath1FFTTime = doInvFFT(outSize, outPath1FFT); 
    FFTWComplex * outPath2FFT = vectortoFFTW (outPath2);
    double * outPath2FFTTime = doInvFFT(outSize, outPath2FFT);
    
    double outSpacing = 1/(pulseSpacing*outSize)*pow(10,9);                 //should that be outSize,or pulseData.size()
    vector<vector<double> > out;
    for (int i=0; i<outSize ; i++)
    {
        //cout << outPath1FFTTime[i] << " " << outPath2FFTTime[i] << endl;
        vector<double> temp;
        temp.push_back(i*outSpacing);
        temp.push_back(outPath1FFTTime[i]+outPath2FFTTime[i]);              //Fuck... i dont remember, is that a - or a +? 
                                                                            //logic says +, but I seem to remember -?
        out.push_back(temp);
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////pulse getting functions ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//puts into units of [ns]
//subtracts off DC voltage (/first bin, which is ~ the DC voltage)
//pads the pulse.
vector<vector<double> > myHybrid::timePulse ( vector< vector<double> > data)
{
    double tPerSample = 4*pow(10,-11);
    vector<vector<double> > prePad;
    for (int i=0; i<data.size(); i++)
    {
        vector<double> temp;
        temp.push_back(data[i][0]*pow(10,9)+30);                            //convert to ns, and make it start at t=0
        temp.push_back(data[i][1]-data[0][1]);                              //subtract off DC
        prePad.push_back(temp);
    }
    vector<vector<double> > out = pulsePadder(prePad, pow(2,14), 0);        //pads to ~ 8 times the size with 0s
    //cout << prePad.size() << " " << out.size() << endl;
    return out;
}

//given raw data, returns freq domain pulse
vector<vector<complex<double> > > myHybrid::freqPulse( vector< vector<double> > data)
{
    vector<vector<double> > timeData = timePulse(data);
    
    double doFFTinput [timeData.size()];
    for (int i=0; i < timeData.size(); i++) 
    { doFFTinput[i] = timeData[i][1]; } 
    int lengthFFT = timeData.size()/2+1;                                                          //this is the size of the FFTWComplex * object we get back
    FFTWComplex * pulseFFT = doFFT( timeData.size() , doFFTinput);
    
    double freqSpacing = 1/((timeData[1][0]-timeData[0][0])*pow(10,-9)*timeData.size());          //because pulseData is in units of [ns]
    vector<vector<complex<double> > > out = FFTWtovector(pulseFFT,lengthFFT,freqSpacing); 
    
    return out;
}

//given timeData-ed data, returns freq domain pulse
vector<vector<complex<double> > > myHybrid::freqPulseNotRaw( vector< vector<double> > timeData)
{
    
    double doFFTinput [timeData.size()];
    for (int i=0; i < timeData.size(); i++) 
    { doFFTinput[i] = timeData[i][1]; } 
    int lengthFFT = timeData.size()/2+1;                                                          //this is the size of the FFTWComplex * object we get back
    FFTWComplex * pulseFFT = doFFT( timeData.size() , doFFTinput);
    
    double freqSpacing = 1/((timeData[1][0]-timeData[0][0])*pow(10,-9)*timeData.size());          //because pulseData is in units of [ns]
    vector<vector<complex<double> > > out = FFTWtovector(pulseFFT,lengthFFT,freqSpacing); 
    
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////real+i*imag -> mag, phase functions ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

vector<vector<double> > myHybrid::getdBFunc(vector<vector<complex<double> > > data)
{
    vector<vector<double> > out;
    for (int i=0; i<data.size(); i++)
    {
        vector<double> temp;
        temp.push_back(data[i][0].real());
        temp.push_back(20*log(abs(data[i][1])));
        out.push_back(temp);
    }
    return out;
}

vector<vector<double> > myHybrid::getPhaseFunc( vector<vector<complex<double> > > data)
{
    vector<double> wrappedPulse;
    for (int i=0; i<data.size(); i++)
    {
        wrappedPulse.push_back( (360/PI)*atan(data[i][1].imag()/data[i][1].real()) );   
        //note we use the wrong correction here, we fix it down below
    }
    vector<vector<double> > out;
    vector<double> unwrappedPulse = unwrapFunction(wrappedPulse);
    for (int i=0; i<data.size(); i++)
    {
        vector<double> temp;
        temp.push_back(data[i][0].real());
        temp.push_back(unwrappedPulse[i]/2);
        out.push_back(temp);
    }
    return out; 
}


/////////////////////////////////////////////////////////////////////////////////
/////get Data functions /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


vector<vector<double> > myHybrid::getNAData( string fileLocation )
{
    vector<vector<double> > data;
    string line;
    ifstream myfile (fileLocation.c_str());
    if (myfile.is_open())
    {
        cout << "file '" << fileLocation << "' opened" << endl;
        int lineNum = 0;
        while( getline(myfile,line) )
        {
            if (lineNum < 3) 
            { 
                lineNum++;
                continue;                                   //skip the first 3 lines
            }
            vector<double> lineVector;
            stringstream s(line);                           //istringstream or stringstream?
            string field;
            int columnNum = 0;                          
            while (getline(s, field, ','))
            {
                if (columnNum > 1) 
                { 
                    columnNum++;
                    continue;                               //skip the 3rd (last) column
                }
                stringstream fs( field );
                double temp = 0.0;
                fs >> temp;
                lineVector.push_back(temp);
            }
            data.push_back(lineVector);
        }
    }
    myfile.close();
    cout << "file '" << fileLocation << "' closed" << endl;
    return data;
}

vector<vector<double> > myHybrid::getPulseData( string fileLocation )
{
    vector<vector<double> > data;
    string line;
    ifstream myfile (fileLocation.c_str());
    if (myfile.is_open())
    {
        cout << "file '" << fileLocation << "' opened" << endl;
        while( getline(myfile,line) )
        {
            vector<double> lineVector;
            stringstream s(line);                           //istringstream or stringstream?
            string field;
            int columnNum = 0;                          
            while (getline(s, field, ','))
            {
                if (columnNum < 3) 
                { 
                    columnNum++;                
                    continue;                               //skip the first 3 column
                }
                stringstream fs( field );
                double temp = 0.0;
                fs >> temp;
                lineVector.push_back(temp);
            }
            data.push_back(lineVector);
        }
    }
    myfile.close();
    cout << "file '" << fileLocation << "' closed" << endl;
    return data;
}

vector<vector<complex<double> > > myHybrid::getCorrectionData( string fileLocation )
{
    vector<vector<complex<double> > > data;
    string line;
    ifstream myfile (fileLocation.c_str());
    if (myfile.is_open())
    {
        cout << "file '" << fileLocation << "' opened" << endl;
        int lineNum = 0;
        while( getline(myfile,line) )
        {
            if (lineNum < 1) 
            { 
                lineNum++;
                continue;                                   //skip the first line
            }
            vector<complex<double> > lineVector;
            stringstream s(line);              
            string field;                       
            while (getline(s, field, ';'))
            {
                stringstream fs( field );
                complex<double> temp = (0.0,0.0);
                fs >> temp;
                lineVector.push_back(temp);
            }
            data.push_back(lineVector);
        }
    }
    myfile.close();
    cout << "file '" << fileLocation << "' closed" << endl;
    return data;
}



/////////////////////////////////////////////////////////////////////////////////
/////Calibration function ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


//unwrap the phase data
//calibrate S21 hybrid data by subtracting off cable corrections
//combines logmag and phase into one complex value.
vector< vector< complex< double > > >  myHybrid::rawToCalibReIm 
    (vector< vector<double> > hybridLogMag, vector< vector<double> > cableLogMag,
     vector< vector<double> > hybridPhase, vector< vector<double> > cablePhase)
{
    vector<vector<double> > unwrappedHybridPhase = unwrapFunction(hybridPhase);
    vector<vector<double> > unwrappedCablePhase = unwrapFunction(cablePhase);
    vector<vector<complex<double> > > out;
    for (int i=0; i<hybridLogMag.size(); i++)
    {
        vector<complex<double> > temp;
        double realTemp = pow(10, (hybridLogMag[i][1]-cableLogMag[i][1])/20) * cos( (PI/180)*(hybridPhase[i][1]-cablePhase[i][1]) ); 
        double imagTemp = pow(10, (hybridLogMag[i][1]-cableLogMag[i][1])/20) * sin( (PI/180)*(hybridPhase[i][1]-cablePhase[i][1]) );
        temp.push_back( complex<double> (hybridLogMag[i][0],0) );
        temp.push_back( complex<double> (realTemp,imagTemp) );
        out.push_back(temp);
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////data format convertion functions ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

vector<vector<complex<double> > > myHybrid::FFTWtovector ( FFTWComplex * data, int size, double spacing)
{
    vector<vector<complex<double> > > out;
    for (int i=0; i<size ; i++)
    {
        vector<complex<double> > temp;
        temp.push_back( complex<double>(i*spacing,0) );
        temp.push_back( complex<double>(data[i].re,data[i].im) );
        out.push_back(temp);
    }
    return out; 
}

FFTWComplex * myHybrid::vectortoFFTW (vector<vector<complex<double> > > data)
{
    FFTWComplex * out = new FFTWComplex[data.size()];
    for (int i=0; i<data.size(); i++)
    {
        out[i].re = data[i][1].real();
        out[i].im = data[i][1].imag();
    }
    return out;
}

FFTWComplex * myHybrid::vectortoFFTW (vector<complex<double> > data)
{
    FFTWComplex * out = new FFTWComplex[data.size()];
    for (int i=0; i<data.size(); i++)
    {
        out[i].re = data[i].real();
        out[i].im = data[i].imag(); 
    }
    return out;
}


/////////////////////////////////////////////////////////////////////////////////
/////Histograming functions function ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

TH1D * myHybrid::histFunction ( vector<vector<double> > data, const char* title, const char* name )
{
    int nBins = data.size();
    double binSize = data[1][0]-data[0][0];
    double lowBound =  data[0][0] - binSize/2;
    double hiBound = data[data.size()-1][0] + binSize/2;
    TH1D * out = new TH1D( title, name, nBins, lowBound, hiBound );
    for (int i=0; i< data.size(); i++)
    {
        out->Fill(data[i][0],data[i][1]);
    }
    return out;
}

//reim == 0 for real part, reim == 1 for imag part.
TH1D * myHybrid::histFunction ( vector<vector<complex<double> > > data, const char* title, const char* name, int reim )
{
    int nBins = data.size();
    double binSize = data[1][0].real()-data[0][0].real();
    double lowBound =  data[0][0].real() - binSize/2;
    double hiBound = data[data.size()-1][0].real() + binSize/2;
    TH1D * out = new TH1D( title, name, nBins, lowBound, hiBound );
    for (int i=0; i< data.size(); i++)
    {
        if (reim == 0) 
        {
            out->Fill(data[i][0].real(),data[i][1].real());
        } 
        else if (reim == 1)
        {
            out->Fill(data[i][0].real(),data[i][1].imag());
        }
        else
        {
            cout << "ERROR: Incorrect usage: reim = 0 for real part, reim = 1 for imag part" << endl;
            return out;
        }
    }
    return out;
}

TH1D * myHybrid::histFunction ( TGraph * data, const char* title, const char* name )
{
    int nBins = data->GetN();
    double binSize = data->GetX()[1]-data->GetX()[0];
    double lowBound =  data->GetX()[0] - binSize/2;
    double hiBound = data->GetX()[nBins-1] + binSize/2;
    //cout << data->GetX()[0] << " " << data->GetX()[nBins-1] << endl;
    //cout << "nBins: " << nBins << endl;
    //cout << "binSize: " << binSize << endl;
    //cout << "lowBound: " << lowBound << endl;
    //cout << "hiBound: " << hiBound << endl;
    TH1D * out = new TH1D( title, name, nBins, lowBound, hiBound );
    for (int i=0; i< nBins; i++)
    {
        out->Fill(data->GetX()[i],data->GetY()[i]);
    }
    return out;
}



/////////////////////////////////////////////////////////////////////////////////
/////Histogram shift function ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//function to apply a shift to the incoming data by a number of 'shift' bins
void myHybrid::histShift ( TH1D & data, int shift )
{
    int numBins = data.GetNbinsX();
    double dataSpacing = data.GetBinContent(1) - data.GetBinContent(0);
    cout << numBins << " " << dataSpacing << endl;
    cout << data.GetBinContent(10) << endl;
    if (shift <= 0)
    {
        for (int i=0;i<numBins;i++)
        {
            if ( (i+shift) < 0 || numBins <= (i+shift) )                                          //my logic might not be right here.
            { 
                data.SetBinContent( i, 0.0 );
            }
            else
            {
                data.SetBinContent( i, data.GetBinContent(i-shift) );
            }
        }
    }
    else
    {
        for (int i=numBins-1;i>=0;i--)
        {
            if ( (i+shift) < 0 || numBins <= (i+shift) )                                          //my logic might not be right here.
            { 
                data.SetBinContent( i, 0.0 );
            }
            else
            {
                data.SetBinContent( i, data.GetBinContent(i-shift) );
            }
        }
    }
    cout << data.GetBinContent(10) << endl;
}

/////////////////////////////////////////////////////////////////////////////////
///// hybrid path correction function ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//we pass this a 2d vector of pulse data, 
//it passes back the SAME 2d vector with our hybid sim having acted on it!
void myHybrid::hybridPathCorrectionFunction ( vector<vector<double> > & APulse, vector<vector<double> > & BPulse)
{
    // 1) FFT pulse into freq domain

    //cout << "t= " << APulse[1][0] << " " << BPulse[1][0] << endl;
    
    //pulse's time coordinates need to be in [ns]!!!
    vector<vector<complex<double> > > fAPulse = freqPulseNotRaw(APulse);                   //in A
    vector<vector<complex<double> > > fBPulse = freqPulseNotRaw(BPulse);                   //in B
    //sould come out in units of [Hz]
    
    //cout << "f= " << fAPulse[1][0] << " " << fBPulse[1][0] << endl;
    
    // 2) load in the path corrections
    
    string loc = "pathCorrections.dat";
    vector<vector<complex<double> > > pathCorrections = getCorrectionData(loc);
    
    
    // 3) do myHybrid sim stuff   
    
    vector<vector<complex<double> > > AC,AD,BC,BD;
    for (int i=0; i<pathCorrections.size(); i++)
    {
        vector<complex<double> > temp;
        temp.push_back(pathCorrections[i][0]);
        temp.push_back(pathCorrections[i][1]);                                      //need to make sure this works.
        AC.push_back(temp);
        temp[1] = pathCorrections[i][2];
        AD.push_back(temp);
        temp[1] = pathCorrections[i][3];
        BC.push_back(temp);
        temp[1] = pathCorrections[i][4];
        BD.push_back(temp);
    }   
    
    //for (int i=0; i<30; i++)
    //{cout << AC[i][1] << " " << AD[i][1] << " " << BC[i][1] << " " << BD[i][1] << endl;}

    
    vector<vector<double> > outC = myHybridSim ( fAPulse, fBPulse, AC, BC );  
    vector<vector<double> > outD = myHybridSim ( fBPulse, fAPulse, BD, AD );  
    
    //NOTE: this outputs data in time units of [ns] with the first point set to t=0 
    
    //cout << "t= " << outD[1][0] << " " << outC[1][0] << endl;    
    
    // 4) Alter APulse and BPulse
    cout << pathCorrections.size() << " " << outC.size() << endl;
    
    for (int i=0; i<APulse.size(); i++)
    {
        APulse[i] = outD[i];
        BPulse[i] = outC[i];
    }
}

/* 
//Not working yet, future how to integrate into anitatools

void myHybrid::hybridFilter ( FilteredAnitaEvent * event)
{
    //Get the waveforms from he root tree
    AnalysisWaveform horz = getwf(event, 1, AnitaPol::kHorizontal);
    AnalysisWaveform vert = getwf(event, 1, AnitaPol::kVertical);


    //whats this do?
    int N = TMath::Min(x->Neven(), y->Neven()); 

    // Put the horz and vert into vectors
    vector<vector<double> > horzVect, vertVect;
    for (int i = 0; i < N; i++)
    {
        vector<double> temp;
        temp.push_back(horz->GetX()[i]);
        temp.push_back(horz->GetY()[i]);
        horzVect.push_back(temp);
        temp[0] = vert->GetX()[i];
        temp[0] = vert->GetY()[i];
        vertVect.push_back(temp);
    }
    
    //For testing, if we made it this far, lets see what our data looks like!
    
    //Next do the hybrid conversion stuff, (ref the above function I guess?)    
    
    hybridPathCorrectionFunction ( * horzVect, * vertVect);

    // Get the transofmed data into the root tree event... or whatever
    
    double new_horz[N] __attribute__((aligned)); 
    double new_vert[N] __attribute__((aligned)); 

    for (int i = 0; i < N; i++) 
    {
      new_horz[i] = horzVect[i][1];
      new_vert[i] = vertVect[i][1];
    }
    
    memcpy(horz->updateEven()->GetY(), new_horz, N * sizeof(double)); 
    memcpy(vert->updateEven()->GetY(), new_vert, N * sizeof(double)); 

    //whats this do?
    x->forceEvenSize(N); 
    y->forceEvenSize(N); 
    
}

*/






























