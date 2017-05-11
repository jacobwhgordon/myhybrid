//Main function to do hybrid stuff

#include <iostream>
#include <vector>
#include <complex> 
#include <fstream>
#include <sstream>
#include <string>

#include <stdio.h>
#include <math.h>

#include <TCanvas.h>

//AnitaTools libs
#include "../AnitaTools/include/FFTtools.h"

#include "myHybridFunctions.h"

using namespace std;
using namespace FFTtools;

const double PI  =3.141592653589793238463;

int main()
{
    myHybrid mh;
    string loc;
    
    TCanvas* myCanv = new TCanvas("myCanv", "title", 600, 600);                    // instantiate, name must be unique
    myCanv->cd();                                                                  // select the canvas for display
    
        
    //import pulse
    cout << endl << "Importing Data" << endl;    
    loc = "data/PulseData/MiniFilter.csv";
    vector<vector<double> > rawPulseData = mh.getPulseData(loc);    
    //this will pad the pulse, and also subtract off the biase in voltage and convert the time to units of [ns]
    vector<vector<double> > pulseData = mh.timePulse(rawPulseData);
    //outcoming 2d vector should be 2^14 entries long.

    //Also need outCData
    loc = "data/PulseData/MiniFilterOutC.csv";
    vector<vector<double> > rawPulseDataOutC = mh.getPulseData(loc);    
    vector<vector<double> > pulseDataOutC = mh.timePulse(rawPulseDataOutC);

    
    //cout << "Pulse Data: (" << pulseData[0].size() << " by " << pulseData.size() << " entry array)" << endl;
    //for (int i; i < pulseData.size(); i++)
    //{
    //    cout << "t=" << pulseData[i][0] << "  V=" << pulseData[i][1] << endl;
    //}
    
    //Import log mag data
    loc = "data/NetworkAnalizer/HYBRID_AC_LOGMAG_S21.CSV";
    vector<vector<double> > logMagACData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/HYBRID_AD_LOGMAG_S21.CSV";
    vector<vector<double> > logMagADData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/HYBRID_BC_LOGMAG_S21.CSV";
    vector<vector<double> > logMagBCData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/HYBRID_BD_LOGMAG_S21.CSV";
    vector<vector<double> > logMagBDData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/CABLE_LOGMAG_S21.CSV";
    vector<vector<double> > logMagCableData = mh.getNAData(loc);
    //import phase data
    loc = "data/NetworkAnalizer/HYBRID_AC_PHASE_S21.CSV";
    vector<vector<double> > phaseACData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/HYBRID_AD_PHASE_S21.CSV";
    vector<vector<double> > phaseADData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/HYBRID_BC_PHASE_S21.CSV";
    vector<vector<double> > phaseBCData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/HYBRID_BD_PHASE_S21.CSV";
    vector<vector<double> > phaseBDData = mh.getNAData(loc);
    loc = "data/NetworkAnalizer/CABLE_PHASE_S21.CSV";
    vector<vector<double> > phaseCableData = mh.getNAData(loc);
    //cout << "LogMag AC Data: (" << logMagACData[0].size() << " by " << logMagACData.size() << " entry array)" << endl;

    cout << "Data Imported" << endl;
    cout << endl << "==================================================" << endl << endl;
    

    

    cout << "Calibrating S21 Data" << endl;
    
    vector<vector<complex<double> > > ACS21 = mh.rawToCalibReIm ( logMagACData, logMagCableData, phaseACData, phaseCableData );
    vector<vector<complex<double> > > ADS21 = mh.rawToCalibReIm ( logMagADData, logMagCableData, phaseADData, phaseCableData );
    vector<vector<complex<double> > > BCS21 = mh.rawToCalibReIm ( logMagBCData, logMagCableData, phaseBCData, phaseCableData );
    vector<vector<complex<double> > > BDS21 = mh.rawToCalibReIm ( logMagBDData, logMagCableData, phaseBDData, phaseCableData );
    
    cout << "S21 Data Calibrated" << endl;
    cout << endl << "==================================================" << endl << endl;
    
    
    
    
    cout << "FFTing the Pulse data" << endl;
    
    //doFFT needs a 1d array, not a 2d vector, so lets make a 1d array
    //this may need to be padded out to be a size of 2^(something).  if you dont do that it might give incorrect answers.
    double doFFTinput [pulseData.size()];
    for (int i=0; i < pulseData.size(); i++) { doFFTinput[i] = pulseData[i][1]; } 
    int lengthFFT = pulseData.size()/2+1; //this is the size of the FFTWComplex * object we get back
    FFTWComplex * pulseFFT = doFFT( pulseData.size() , doFFTinput);
    
    //Lets convert the doFFT output into a 2d complex vector
    double freqSpacing = 1/((pulseData[1][0]-pulseData[0][0])*pow(10,-9)*pulseData.size()); //because pulseData is in units of [ns]
    //cout << "timeSpacing= " << pulseData[1][0]-pulseData[0][0] << " in [ns], freqSpacing= " << freqSpacing << " in [hz]" << endl;
    vector<vector<complex<double> > > pulseFreq = mh.FFTWtovector(pulseFFT,lengthFFT,freqSpacing); 
    
    cout << "Pulse data Tranformed using FFTW" << endl;
    cout << endl << "==================================================" << endl << endl;
    
    cout << "'Extending' S21 Data" << endl;
    
    vector<vector<complex<double> > > extendAC = mh.extendS21 ( ACS21, pulseFreq, 1);
    vector<vector<complex<double> > > extendAD = mh.extendS21 ( ADS21, pulseFreq, 1);
    vector<vector<complex<double> > > extendBC = mh.extendS21 ( BCS21, pulseFreq, 1);
    vector<vector<complex<double> > > extendBD = mh.extendS21 ( BDS21, pulseFreq, 1);
        
    cout << "S21 Data range 'Extended' to match FFTed Pulse Data" << endl;
    cout << endl << "==================================================" << endl << endl;
    
    cout << "Matching S21 and freqPulse data " << endl;
    
    vector<vector<complex<double> > > AC = mh.matchFunction (pulseFreq, extendAC);
    vector<vector<complex<double> > > AD = mh.matchFunction (pulseFreq, extendAD);
    vector<vector<complex<double> > > BC = mh.matchFunction (pulseFreq, extendBC);
    vector<vector<complex<double> > > BD = mh.matchFunction (pulseFreq, extendBD);
    
    cout << "S21 and freqPulse data matched" << endl;
    
    /*
    // outputs for testing that extendS21 and match function is working.  It seems to be!
    cout << " ACS21.size() =       " << ACS21.size();
    cout << " extendAC.size() =    " << extendAC.size();
    cout << " AC.size() =          " << AC.size();
    cout << " pulseFreq.size() =   " << pulseFreq.size() << endl;
    cout << " ACS21 endpoint =     " << ACS21[ACS21.size()-1][0].real() ;
    cout << " extendAC endpoint =  " << extendAC[extendAC.size()-1][0].real();
    cout << " AC endpoint =        " << AC[AC.size()-1][0].real();
    cout << " pulseFreq endpoint = " << pulseFreq[pulseFreq.size()-1][0].real() << endl;
    cout << " ACS21 spacing =      " << ACS21[1][0].real()-ACS21[0][0].real();
    cout << " extendAC spacing =   " << extendAC[1][0].real() - extendAC[0][0].real();
    cout << " AC spacing =         " << AC[1][0].real() - AC[0][0].real();
    cout << " pulseFreq spacing =  " << pulseFreq[1][0].real()-pulseFreq[0][0].real() << endl;
    */
    cout << endl << "==================================================" << endl << endl;
    
    cout << "Applying transfer function to pulse, inverse FFT back into time domain and summing hybrid port outputs" << endl;
    
    vector<vector<double> > outC = mh.myHybridSim ( pulseFreq, AC, BC );  
    vector<vector<double> > outD = mh.myHybridSim ( pulseFreq, BD, AD );  
    //NOTE: this outputs data in time units of [ns] with the first point set to t=0 
    
    cout << "Back in time domain.  Pulse transformed." << endl;
    
    /*
    cout << "size of outC " << outC.size() << endl;
    cout << "size of time pulse " << pulseData.size() << endl;
    cout << "last point of outC " << outC[outC.size()-1][0] << endl;
    cout << "last point of time pulse " << pulseData[pulseData.size()-1][0] << endl;
    cout << "range of outC " << outC[outC.size()-1][0]-outC[0][0] << endl;
    cout << "range of time Pulse " << pulseData[pulseData.size()-1][0]-pulseData[0][0] << endl;
    */
    cout << endl << "==================================================" << endl << endl;
    

    //we want to test the hilbert methoid too.
    
    cout << "Hilbert Transform Method" << endl;
    double pulseX [pulseData.size()];
    double pulseY [pulseData.size()];
    for (int i=0; i < pulseData.size() ; i++)
    {
        pulseX[i] = pulseData[i][0];   
        pulseY[i] = pulseData[i][1];
    }
    
    TGraph * pulseGraph = new TGraph(pulseData.size(), pulseX, pulseY);
    TGraph * hilbert = getHilbertTransform(pulseGraph); 
        
    //convert into hists
    TH1D * hilbertHist = mh.histFunction(hilbert, "name", "hilbertHist");
    TH1D * pulseDataHist = mh.histFunction(pulseData, "name", "pulseData");    

    //note outC and outD are exactly the same with this methoid     
    hilbertHist->Add(pulseDataHist,-1);                                               //subtract the hilbertHist fromthe pulseDataHist
    hilbertHist->Scale(-1*pow(2,-0.5));                                                  //nomalize by 1/root(2), extra neg because the above subtraction was backwords
     
    hilbertHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                  //set y axis
    hilbertHist->GetXaxis()->SetRangeUser(350,450);                                   //set x axis
    hilbertHist->Draw();                                                              // Draw the histogram
    myCanv->SaveAs("plots/hilbertData.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/hilbertData.root");
    


    cout << endl << "==================================================" << endl << endl;





    
    
    //simple Loop for testing stuff
    
    /*
    for (int i=0; i < outC.size() ; i++)
    {
        cout << "i=" << i ;
        //cout << ", input=" << doFFTinput[i] ;
        //cout << ", PulseFFT=" << pulseFFT[i]; 
        //cout << ", PulseFFT.re=" << pulseFFT[i].re; 
        //cout << " Freq: " << ACS21[i][0].real();
        //cout << " Real S21: " << ACS21[i][1].real();
        //cout << " Imag S21: " << ACS21[i][1].imag();
        cout << " Time: " << outC[i][0];
        cout << " Volt: " << outC[i][1];
        cout << endl;
    } 
    */
    
    
    cout << "Plotting!" << endl;

    /*
    TCanvas* myCanv = new TCanvas("myCanv", "title", 600, 600);                    // instantiate, name must be unique
    myCanv->cd();                                                                  // select the canvas for display
    */

    //Plot the pulseData    
    //TH1D * pulseDataHist = mh.histFunction(pulseData, "name", "pulseData");       //this is done up higher now
    pulseDataHist->GetYaxis()->SetRangeUser(-0.2,0.2);                              //set y axis
    pulseDataHist->GetXaxis()->SetRangeUser(350,450);                               //set x axis
    pulseDataHist->Draw();                                                          // Draw the histogram
    myCanv->SaveAs("plots/pulseData.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/pulseData.root");

    //plot the pulseFreq.real()
    TH1D * pulseFreqHist = mh.histFunction(pulseFreq, "name", "pulseFreq", 0);
    //pulseFreqHist->GetYaxis()->SetRangeUser(-0.00015, 0.00015);                   //this was EXACTLY the window in mathematica?
    pulseFreqHist->GetYaxis()->SetRangeUser(-8, 8);                                 //differenece in size of about a factor of 4000? thats... weird
    pulseFreqHist->Draw();                                                          // Draw the histogram
    myCanv->SaveAs("plots/pulseFreq.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/pulseFreq.root");    
    
    //plot the S21 AC            'ACS21'
    TH1D * ACS21Hist = mh.histFunction(ACS21, "name", "raw AC", 0);
    ACS21Hist->GetYaxis()->SetRangeUser(-1,1);                                      //set y axis
    ACS21Hist->Draw();                                                              // Draw the histogram
    myCanv->SaveAs("plots/ACS21.jpg");                                              // saves the plot
    myCanv->SaveAs("plots/ACS21.root");
    
    //plot the S21 AC extended   'extendAC'
    TH1D * extendACHist = mh.histFunction(extendAC, "name", "extended AC", 0);
    extendACHist->GetYaxis()->SetRangeUser(-1,1);                                   //set y axis
    extendACHist->Draw();                                                           // Draw the histogram
    myCanv->SaveAs("plots/extendAC.jpg");                                           // saves the plot
    myCanv->SaveAs("plots/extendAC.root");
        
    //plot the S21 AC matched    'AC'
    TH1D * ACHist = mh.histFunction(AC, "name", "AC", 0);
    ACHist->GetYaxis()->SetRangeUser(-1,1);                                         //set y axis
    ACHist->Draw();                                                                 // Draw the histogram
    myCanv->SaveAs("plots/AC.jpg");                                                 // saves the plot
    myCanv->SaveAs("plots/AC.root");

    //plot the outC   'outC'
    TH1D * outCHist = mh.histFunction(outC, "name", "outC");
    outCHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                   //set y axis
    outCHist->GetXaxis()->SetRangeUser(350,450);                                    //set x axis
    outCHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/outC.jpg");                                               // saves the plot
    myCanv->SaveAs("plots/outC.root");
    
    
    //make a plot with multiple lines on it
    //we want hilbert outC and outCdata all on one plot.
    //hilbert and outC are already made, so just need to create outCData
    TH1D * outCDataHist = mh.histFunction(pulseDataOutC, "name", "pulseDataOutC");  
    outCDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                              //set y axis
    outCDataHist->GetXaxis()->SetRangeUser(350,400);                               //set x axis
    outCDataHist->Draw();                                                          // Draw the histogram
    
    outCHist->SetLineColor(kRed);
    outCHist->Draw("same");
    pulseDataHist->SetLineColor(kBlue);
    pulseDataHist->Draw("same");
    myCanv->SaveAs("plots/outCAll.jpg");
    myCanv->SaveAs("plots/outCAll.root");
    
    
    
    cout << "Done generating plots" << endl;
    cout << endl << "==================================================" << endl << endl;


    //Testing Area!!!//////////////////////////////////////////////

    cout << "Entering Testing area." << endl;
    
    
    
    
    //plot doInvFFT(doFFT(pulseData)) to make sure we get back pulseData.  fftTestPulse
    
    double * fftTest = doInvFFT(pulseData.size(), pulseFFT);
    //make the output
    vector<vector<double > > fftTestPulse;
    double timeSpacing = pulseData[1][0]-pulseData[0][0];
    for (int i=0; i < pulseData.size(); i++)
    {
        vector<double>  temp;
        temp.push_back( i*timeSpacing );
        temp.push_back ( fftTest[i] );
        fftTestPulse.push_back(temp);
    }
    
    TH1D * fftTestHist = mh.histFunction(fftTestPulse, "name", "fftTest");
    outCHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                   //set y axis
    outCHist->GetXaxis()->SetRangeUser(200,450);                                    //set x axis
    outCHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/fftTest.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/fftTest.root");
        
    
    
    
    
    cout << "Leaving Testing area." << endl;
    cout << endl << "==================================================" << endl << endl;
    
    return(0);
}
    
     










































































