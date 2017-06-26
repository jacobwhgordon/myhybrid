//Test file to make sure fuctions are working
//
// By Jacob Gordon
// 5/26/17
//
/////////////////////////////////////////////////////////

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
#include "../../AnitaTools/include/FFTtools.h"

#include "myHybridFunctions.h"

using namespace std;
using namespace FFTtools;

const double PI  =3.141592653589793238463;

int main()
{
    myHybrid mh;
    string loc;
    
    //We need outC and outD data to compare too
    loc = "data/PulseData/MiniFilterOutC.csv";
    vector<vector<double> > rawPulseDataOutC = mh.getPulseData(loc);    
    vector<vector<double> > pulseDataOutC = mh.timePulse(rawPulseDataOutC);
    loc = "data/PulseData/MiniFilterOutD.csv";
    vector<vector<double> > rawPulseDataOutD = mh.getPulseData(loc);    
    vector<vector<double> > pulseDataOutD = mh.timePulse(rawPulseDataOutD);
    
    

    // First Test hybridPathCorrectionFunction
    // It takes h+v
    // For now let both h+v be the stadard input pulse
    // I think it needs them to be in [ns], address that later!

    loc = "data/PulseData/MiniFilter.csv";
    vector<vector<double> > rawPulseData = mh.getPulseData(loc);  
    vector<vector<double> > x = mh.timePulse(rawPulseData);                       //Right now x is h (?)
    vector<vector<double> > y = mh.timePulse(rawPulseData);                       //Right now y is v (?)
    vector<vector<double> > y_backup = mh.timePulse(rawPulseData);                       //Right now y is v (?)

    /*
    for (int i=0;i<30;i++)
    {
        cout << x[i][0] << " " << x[i][1] << " " << y[i][1] << endl;
    }
    */
    
    mh.hybridPathCorrectionFunction( x, y );

    /*
    for (int i=0;i<30;i++)
    {
        cout << x[i][0] << " " << x[i][1] << " " << y[i][1] << endl;
    }
    */ 

    //now we have our l and r (or outC and outD)
    //just need to make some hists and plot them
    
    TCanvas* myCanv = new TCanvas("myCanv", "title", 800, 600);                    // instantiate, name must be unique
    myCanv->cd();                                                                  // select the canvas for display


    //make outC plot    
    
    TH1D * outCDataHist = mh.histFunction(pulseDataOutC, "name", "pulseDataOutC");  
    outCDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                            //set y axis
    outCDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    outCDataHist->Draw();                                                          // Draw the histogram

    TH1D * yHist = mh.histFunction(y, "name", "y");  
    yHist->GetYaxis()->SetRangeUser(-0.16,0.16);                                   //set y axis
    yHist->GetXaxis()->SetRangeUser(350,390);                                      //set x axis
    mh.histShift(*yHist,0);
    yHist->SetLineColor(kRed);
    yHist->Draw("same");
    
    myCanv->SaveAs("plots/outCShort.jpg");                                         // saves the plot
    myCanv->SaveAs("plots/outCShort.root");
    
    //make outD plot
    
    TH1D * outDDataHist = mh.histFunction(pulseDataOutD, "name", "pulseDataOutD");  
    outDDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                            //set y axis
    outDDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    outDDataHist->Draw();                                                          // Draw the histogram
    
    TH1D * xHist = mh.histFunction(x, "name", "x");  
    xHist->GetYaxis()->SetRangeUser(-0.16,0.16);                                   //set y axis
    xHist->GetXaxis()->SetRangeUser(350,390);                                      //set x axis
    mh.histShift(*xHist,0);
    xHist->SetLineColor(kRed);
    xHist->Draw("same");
    
    myCanv->SaveAs("plots/outDShort.jpg");                                         // saves the plot
    myCanv->SaveAs("plots/outDShort.root");

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Test the anita data type function thing///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////











    return 0;
}


















