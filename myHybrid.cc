/////////////////////////////////////////////////////////////////
//
// Main function to do hybrid stuff
//
//  -takes in pulse data
//  -FFTs it
//  -Applies hybrid transfer function to pulse data from:
//  --Network analyser data
//  --Correct and incorrect hilbert transfer methods
//  --FFT phase shift method (used in iceMC)
//  -Converts back to time domain
//  -creates plots
//
//  functions refernced found in myHybridFunction.h
// 
// by Jacob Gordon
//
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <complex> 
#include <fstream>
#include <sstream>
#include <string>

#include <stdio.h>
#include <math.h>

#include <TCanvas.h>
#include <TLegend.h>

//AnitaTools libs
#include "FFTtools.h"

#include "myHybridFunctions.h"

using namespace std;
using namespace FFTtools;

const double PI = 3.141592653589793238463;

int main()
{
    myHybrid mh;
    string loc;
    
    TCanvas* myCanv = new TCanvas("myCanv", "title", 800, 600);                    // instantiate, name must be unique
    myCanv->cd();                                                                  // select the canvas for display
    
        
    //import pulse
    cout << endl << "Importing Data" << endl;    
    loc = "data/PulseData/MiniFilter.csv";
    vector<vector<double> > rawPulseData = mh.getPulseData(loc);      
      
    //this will pad the pulse, and also subtract off the bias in voltage and convert the time to units of [ns]
    vector<vector<double> > pulseData = mh.timePulse(rawPulseData);
    //outcoming 2d vector should be 2^14 entries long.

    //Also need outCData
    loc = "data/PulseData/MiniFilterOutC.csv";
    vector<vector<double> > rawPulseDataOutC = mh.getPulseData(loc);    
    vector<vector<double> > pulseDataOutC = mh.timePulse(rawPulseDataOutC);

    //cout << rawPulseData.size() << endl;
    cout << "Raw Pulse Data: (" << rawPulseData[0].size() << " by " << rawPulseData.size() << " entry array)" << endl;    
    cout << "Pulse Data: (" << pulseData[0].size() << " by " << pulseData.size() << " entry array)" << endl;

    for (int i=pulseData.size()/2; i < (pulseData.size()/2+100); i++)
    {
        cout << "t=" << pulseData[i][0] << "  V=" << pulseData[i][1] << endl;
    }
    
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

    //cout << " ACS21.size() =       " << ACS21.size();
    //cout << " ACS21 endpoint =     " << ACS21[ACS21.size()-1][0].real() ;
    //cout << " ACS21 spacing =      " << ACS21[1][0].real()-ACS21[0][0].real() << endl;

    
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
    
    //cout << "pulseFreq size: " << pulseFreq.size() << endl;
    //cout << "ACS21     size: " << ACS21.size() << endl;
    
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
    
    /*
    for (int i=0; i<AC.size(); i++)
    {
        cout << "extend:  (" << extendAC[i][0].real() << "," << extendAC[i][1].real() << " + i*" << extendAC[i][1].imag() << ")" << endl; 
        cout << "matched: (" << AC[i][0].real() << "," << AC[i][1].real() << " + i*" << AC[i][1].imag() << ")" << endl; 
    }
    */
    
    cout << "S21 and freqPulse data matched" << endl;
    cout << endl << "==================================================" << endl << endl;
    
    
    
    cout << "Outputting path corrections to file" << endl;
    // we want to output AC,AD,ect to a .dat file for future reference later 
    // (so we dont need to load all 8 of the SA tables.
    // format:  freq, AC, AD, BC, BD
    
    ofstream outputFile;
    outputFile.open ("pathCorrections.dat");
    outputFile << "Frequency [Hz] (re,im); AC correction (re,im); AD (re,im); BC (re,im); BD (re,im)" << endl; 
    for (int i=0; i< AC.size(); i++)
    {
        outputFile << AC[i][0] << "; " << AC[i][1] << "; " << AD[i][1] << "; " << BC[i][1] << "; " << BD[i][1] << endl;
    }
    outputFile.close();
    
    outputFile.open ("na.dat");
    outputFile << "ACS21 (re,im); extendAC (re,im); AC (re,im)" << endl; 
    for (int i=0; i< ACS21.size(); i++)
    {
        outputFile <<  ACS21[i][1] << "; " << extendAC[i][1] << "; " << AC[i][1] << endl;
    }
    outputFile.close();
    
    cout << "AC, AD, BC, BD exported to file" << endl;
    cout << endl << "==================================================" << endl << endl;
    
    
    
    
    
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
    

    //we want to test the hilbert method too.
    
    cout << "Other Transform Methods" << endl;
    //hilbert test
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
    TH1D * hilbertHistWrong = mh.histFunction(hilbert, "name", "hilbertHist");
    TH1D * pulseDataHist = mh.histFunction(pulseData, "name", "pulseData");    
    //note outC and outD are exactly the same with this method     
    hilbertHist->Add(pulseDataHist,1);                                               //subtract the hilbertHist fromthe pulseDataHist
                                                                                     //note, it seems now we are adding them.... ummmmmmmmm
    hilbertHist->Scale(1*pow(2,-0.5));                                               //nomalize by 1/root(2),
                                                                                     
    hilbertHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                  //set y axis
    hilbertHist->GetXaxis()->SetRangeUser(350,400);                                   //set x axis
    hilbertHist->SetOption("HIST C");                                               //set to draw smooth "C"urve and without error bars
    hilbertHist->Draw();                                                              // Draw the histogram
    myCanv->SaveAs("plots/hilbertData.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/hilbertData.root");
    
    //Make plot to show what happens when we do it the way shown in AnitaTools (the wrong way)...
    hilbertHistWrong->Add(pulseDataHist,-1);
    hilbertHistWrong->Scale(1*pow(2,-0.5));
    
    hilbertHistWrong->GetYaxis()->SetRangeUser(-0.2,0.2);                                  //set y axis
    hilbertHistWrong->GetXaxis()->SetRangeUser(350,400);                                   //set x axis
    hilbertHistWrong->SetOption("HIST C");                                                 //set to draw smooth "C"urve and without error bars
    hilbertHistWrong->SetLineColor(kRed);
    hilbertHistWrong->Draw();
    hilbertHist->Draw("same HIST C");
    myCanv->SaveAs("plots/hilbertDataWrong.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/hilbertDataWrong.root");

    //now we also wnat to test the 'standard' method which is used in iceMC
    
    //we need two new variables, r_ice and l_ice, and we need to use the inputs (h and v) in freq domain.  pulseFreq
    vector<vector<complex<double> > > r_ice_freq;
    vector<vector<complex<double> > > l_ice_freq;
    for( int i=0; i<pulseFreq.size(); i++) 
    {
        vector<complex<double> > temp;
        temp.push_back(complex<double>(pulseFreq[i][0].real(),pulseFreq[i][0].imag()));
        temp.push_back(complex<double>(pow(2,-0.5)*(pulseFreq[i][1].real()+pulseFreq[i][1].imag()), pow(2,-0.5)*(pulseFreq[i][1].imag()-pulseFreq[i][1].real())));
        r_ice_freq.push_back(temp);
        temp[1]=complex<double>(pow(2,-0.5)*(pulseFreq[i][1].real()-pulseFreq[i][1].imag()), pow(2,-0.5)*(pulseFreq[i][1].imag()+pulseFreq[i][1].real()));
        l_ice_freq.push_back(temp);
    }
    
    //now fft back;
    vector<vector<double> > r_ice = mh.doInverseFFT( r_ice_freq );
    vector<vector<double> > l_ice = mh.doInverseFFT( l_ice_freq );
        
        
    TH1D * outCDataHist = mh.histFunction(pulseDataOutC, "name", "pulseDataOutC");  
    outCDataHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                  //set y axis
    outCDataHist->GetXaxis()->SetRangeUser(350,400);                                   //set x axis
    outCDataHist->Draw("Hist C");                                                      // Draw the histogram        
    TH1D * iceRHist = mh.histFunction(r_ice, "name", "r_ice");
    iceRHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                      //set y axis
    iceRHist->GetXaxis()->SetRangeUser(350,400);                                       //set x axis
    iceRHist->SetLineColor(kRed);
    mh.histShift(*iceRHist,29);
    iceRHist->Draw("Hist C Same");                                                     // Draw the histogram
    TH1D * iceLHist = mh.histFunction(l_ice, "name", "l_ice");
    iceLHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                      //set y axis
    iceLHist->GetXaxis()->SetRangeUser(350,400);                                       //set x axis
    iceLHist->SetLineColor(kBlue);
    mh.histShift(*iceLHist,36);
    iceLHist->Draw("Hist C Same");
    myCanv->SaveAs("plots/IceMCOut.jpg");                                              // saves the plot
    myCanv->SaveAs("plots/IceMCOut.root");

    
    
    /*
    cout << "pulseData spacing: " << pulseData[1][0]-pulseData[0][0] << endl;
    cout << "pulseX spacing: " << pulseX[1]-pulseX[0] << endl;    
    cout << "pulseGraph spacing: " << pulseGraph->GetX()[1]-pulseGraph->GetX()[0] << endl;
    cout << "hilbert spacing: " << hilbert->GetX()[1]-hilbert->GetX()[0] << endl;
    cout << "hibertHist spacing: " << hilbertHist->GetXaxis()->GetBinCenter(2) - hilbertHist->GetXaxis()->GetBinCenter(1) << endl;  
    */


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
    //pulseDataHist->GetXaxis()->SetRangeUser(360,370);
    pulseDataHist->GetXaxis()->SetRangeUser(350,450);                               //set x axis
    pulseDataHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    pulseDataHist->Draw();                                                          // Draw the histogram
    myCanv->SaveAs("plots/pulseData.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/pulseData.root");

    //plot the pulseFreq.real()
    TH1D * pulseFreqHist = mh.histFunction(pulseFreq, "name", "pulseFreq", 0);
    //pulseFreqHist->GetYaxis()->SetRangeUser(-0.00015, 0.00015);                   //this was EXACTLY the window in mathematica?
    pulseFreqHist->GetYaxis()->SetRangeUser(-8, 8);                                 //differenece in size of about a factor of 4000? thats... weird
    pulseFreqHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    pulseFreqHist->Draw();                                                          // Draw the histogram
    myCanv->SaveAs("plots/pulseFreq.jpg");                                          // saves the plot
    myCanv->SaveAs("plots/pulseFreq.root");    
    
    //plot the S21 AC            'ACS21'
    TH1D * ACS21Hist = mh.histFunction(ACS21, "name", "raw AC", 0);
    ACS21Hist->GetYaxis()->SetRangeUser(-1,1);                                      //set y axis
    ACS21Hist->Draw();                                                              // Draw the histogram
    ACS21Hist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    myCanv->SaveAs("plots/ACS21.jpg");                                              // saves the plot
    myCanv->SaveAs("plots/ACS21.root");
    
    //plot the S21 AC extended   'extendAC'
    TH1D * extendACHist = mh.histFunction(extendAC, "name", "extended AC", 0);
    extendACHist->GetYaxis()->SetRangeUser(-1,1);                                   //set y axis
    extendACHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    extendACHist->Draw();                                                           // Draw the histogram
    myCanv->SaveAs("plots/extendAC.jpg");                                           // saves the plot
    myCanv->SaveAs("plots/extendAC.root");
        
    //plot the S21 AC matched    'AC real'
    TH1D * ACHist = mh.histFunction(AC, "name", "AC real", 0);
    ACHist->GetYaxis()->SetRangeUser(-1,1);                                         //set y axis
    ACHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    ACHist->Draw();                                                                 // Draw the histogram
    myCanv->SaveAs("plots/AC_real.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/AC_real.root");
        
    //plot the S21 AD matched    'AD real'
    TH1D * ADHist = mh.histFunction(AD, "name", "AD real", 0);
    ADHist->GetYaxis()->SetRangeUser(-1,1);                                         //set y axis
    ADHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    ADHist->Draw();                                                                 // Draw the histogram
    myCanv->SaveAs("plots/AD_real.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/AD_real.root");
        
    //plot the S21 BC matched    'BC real'
    TH1D * BCHist = mh.histFunction(BC, "name", "BC real", 0);
    BCHist->GetYaxis()->SetRangeUser(-1,1);                                         //set y axis
    BCHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    BCHist->Draw();                                                                 // Draw the histogram
    myCanv->SaveAs("plots/BC_real.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/BC_real.root");
        
    //plot the S21 BD matched    'BD real'
    TH1D * BDHist = mh.histFunction(BD, "name", "BD real", 0);
    BDHist->GetYaxis()->SetRangeUser(-1,1);                                         //set y axis
    BDHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    BDHist->Draw();                                                                 // Draw the histogram
    myCanv->SaveAs("plots/BD_real.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/BD_real.root");
        
    //plot the S21 AC matched    'AC imag'
    TH1D * ACimHist = mh.histFunction(AC, "name", "AC imag", 1);
    ACimHist->GetYaxis()->SetRangeUser(-1,1);                                       //set y axis
    ACimHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    ACimHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/AC_imag.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/AC_imag.root");
        
    //plot the S21 AD matched    'AD imag'
    TH1D * ADimHist = mh.histFunction(AD, "name", "AD imag", 1);
    ADimHist->GetYaxis()->SetRangeUser(-1,1);                                       //set y axis
    ADimHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    ADimHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/AD_imag.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/AD_imag.root");
        
    //plot the S21 BC matched    'BC imag'
    TH1D * BCimHist = mh.histFunction(BC, "name", "BC imag", 1);
    BCimHist->GetYaxis()->SetRangeUser(-1,1);                                       //set y axis
    BCimHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    BCimHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/BC_imag.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/BC_imag.root");
        
    //plot the S21 BD matched    'BD imag'
    TH1D * BDimHist = mh.histFunction(BD, "name", "BD imag", 1);
    BDimHist->GetYaxis()->SetRangeUser(-1,1);                                       //set y axis
    BDimHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    BDimHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/BD_imag.jpg");                                            // saves the plot
    myCanv->SaveAs("plots/BD_imag.root");

    //plot the outC   'outC'
    TH1D * outCHist = mh.histFunction(outC, "name", "outC");
    outCHist->GetYaxis()->SetRangeUser(-0.2,0.2);                                   //set y axis
    outCHist->GetXaxis()->SetRangeUser(350,450);                                    //set x axis
    outCHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    outCHist->Draw();                                                               // Draw the histogram
    myCanv->SaveAs("plots/outC.jpg");                                               // saves the plot
    myCanv->SaveAs("plots/outC.root");
    
    
    //make a plot with multiple lines on it
    //we want hilbert outC and outCdata all on one plot.
    //hilbert and outC are already made, so just need to create outCData
    //TH1D * outCDataHist = mh.histFunction(pulseDataOutC, "name", "pulseDataOutC");  
    outCDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                              //set y axis
    outCDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    outCDataHist->SetOption("HIST C");                                             //set to draw smooth "C"urve and without error bars
    outCDataHist->Draw();                                                          // Draw the histogram
    
    mh.histShift(*outCHist,-1);
    outCHist->SetLineColor(kRed);
    outCHist->Draw("same HIST C");
    mh.histShift(*hilbertHist, 34);
    hilbertHist->SetLineColor(kGreen);
    hilbertHist->Draw("same HIST C");
    
    
    TLegend * legend = new TLegend(0.6,0.7,0.9,0.9);                      //(x1,y1,x2,y2)
    //legend->SetHeader("The Legend Title","C");                          // option "C" allows to center the header
    legend->AddEntry(outCHist,"Hybrid outC sim (my method)","l");
    legend->AddEntry("outCDataHist","Hybrid outC data","l");
    legend->AddEntry("hilbertHist","Hybrid outC sim (hilbert)","l");
    legend->Draw();
    
    
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
        
    
    //Lets have it print out sum squared... thing, or whatever.  
    //basically V[0]^2+V[1]^2+ect...
    double pulseDataOutCSumSqr=0;                                    //pulseDataOutC
    double outCSumSqr=0;                                             //outC
    double outCHilbertSumSqr=0;                                      //base this on hilbertHist (which is a hist, so its trickyer)
    double l_iceSumSqr=0;
    double r_iceSumSqr=0;

    for (int i=0; i < pulseDataOutC.size(); i++) 
    {
        pulseDataOutCSumSqr = pulseDataOutCSumSqr +pow(pulseDataOutC[i][1],2);
        outCSumSqr = outCSumSqr + pow(outC[i][1],2);
        outCHilbertSumSqr = outCHilbertSumSqr + pow(hilbertHist->GetBinContent(i),2);
        l_iceSumSqr = l_iceSumSqr + pow(l_ice[i][1],2);
        r_iceSumSqr = r_iceSumSqr + pow(r_ice[i][1],2);
        
    }
    cout << "Sum Squared Power:" << endl;
    cout << "For Data: " << pulseDataOutCSumSqr << endl;
    cout << "For Data Driven Method: " << outCSumSqr << endl;
    cout << "For Hilbert Method: " << outCHilbertSumSqr << endl;
    cout << "For IceMC L Method: " << l_iceSumSqr << endl;
    cout << "For IceMC R Method: " << r_iceSumSqr << endl;
    cout << "Data Driven Method/Data: " << outCSumSqr/pulseDataOutCSumSqr << endl;
    cout << "Hilbert Method/Data: " << outCHilbertSumSqr/pulseDataOutCSumSqr << endl; 
    cout << "IceMC L Method/Data: " << l_iceSumSqr/pulseDataOutCSumSqr << endl; 
    cout << "IceMC R Method/Data: " << r_iceSumSqr/pulseDataOutCSumSqr << endl; 
    

    
    cout << "Leaving Testing area." << endl;
    cout << endl << "==================================================" << endl << endl;
    
    return(0);
}









