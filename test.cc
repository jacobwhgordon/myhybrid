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
#include "FFTtools.h"

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
//Lets make sure we are getting the right answer for when we only put signal into one of the inputs/////
////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<vector<double> > inA = mh.timePulse(rawPulseData); 
    vector<vector<double> > inB = mh.timePulse(rawPulseData);
    vector<vector<double> > inA_0;
    vector<vector<double> > inB_0;
    for (int i=0; i < inA.size(); i++)
    {
        vector<double> temp;
        temp.push_back(inA[i][0]);
        temp.push_back(0);
        inA_0.push_back(temp);
        inB_0.push_back(temp);
    }
    
    //now we have out inputs, but them in
    
    mh.hybridPathCorrectionFunction( inA, inB_0 );
    mh.hybridPathCorrectionFunction( inA_0, inB );
    
    //so... this is confusing but... 
    //inA is now noBOutD
    //inB_0 is now noBOutC
    //inA_0 is now noAOutD
    //inB is now noAOutC
    
    //lets make some plots and load the data for comparison
    
    loc = "data/PulseData/noAMiniFilterOutC.csv";
    vector<vector<double> > rawNoAOutC = mh.getPulseData(loc);    
    vector<vector<double> > noAOutC = mh.timePulse(rawNoAOutC);
    loc = "data/PulseData/noAMiniFilterOutD.csv";
    vector<vector<double> > rawNoAOutD = mh.getPulseData(loc);    
    vector<vector<double> > noAOutD = mh.timePulse(rawNoAOutD);
    loc = "data/PulseData/noBMiniFilterOutC.csv";
    vector<vector<double> > rawNoBOutC = mh.getPulseData(loc);    
    vector<vector<double> > noBOutC = mh.timePulse(rawNoBOutC);
    loc = "data/PulseData/noBMiniFilterOutD.csv";
    vector<vector<double> > rawNoBOutD = mh.getPulseData(loc);    
    vector<vector<double> > noBOutD = mh.timePulse(rawNoBOutD);
    
    
    //make noAoutC plot    
    
    TH1D * noAOutCDataHist = mh.histFunction(noAOutC, "name", "noAOutC Data");  
    noAOutCDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                            //set y axis
    noAOutCDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    noAOutCDataHist->SetOption("HIST C");     
    noAOutCDataHist->Draw();                                                          // Draw the histogram

    TH1D * inBHist = mh.histFunction(inB, "name", "noAOutC Sim");  
    inBHist->GetYaxis()->SetRangeUser(-0.16,0.16);                                   //set y axis
    inBHist->GetXaxis()->SetRangeUser(350,390);                                      //set x axis
    mh.histShift(*inBHist,0);
    inBHist->SetLineColor(kRed);
    inBHist->SetOption("HIST C");     
    inBHist->Draw("HIST C same");
    
    myCanv->SaveAs("plots/noAoutC.jpg");                                         // saves the plot
    myCanv->SaveAs("plots/noAoutC.root");
    
    //make noAoutD plot    
    
    TH1D * noAOutDDataHist = mh.histFunction(noAOutD, "name", "noAOutD Data");  
    noAOutDDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                            //set y axis
    noAOutDDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    noAOutDDataHist->SetOption("HIST C");     
    noAOutDDataHist->Draw();                                                          // Draw the histogram

    TH1D * inA_0Hist = mh.histFunction(inA_0, "name", "noAOutD Sim");  
    inA_0Hist->GetYaxis()->SetRangeUser(-0.16,0.16);                                   //set y axis
    inA_0Hist->GetXaxis()->SetRangeUser(350,390);                                      //set x axis
    mh.histShift(*inA_0Hist,0);
    inA_0Hist->SetLineColor(kRed);
    inA_0Hist->SetOption("HIST C");     
    inA_0Hist->Draw("HIST C same");
    
    myCanv->SaveAs("plots/noAoutD.jpg");                                         // saves the plot
    myCanv->SaveAs("plots/noAoutD.root");
        
    //make noBoutC plot    
    
    TH1D * noBOutCDataHist = mh.histFunction(noBOutC, "name", "noBOutC Data");  
    noBOutCDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                            //set y axis
    noBOutCDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    noBOutCDataHist->SetOption("HIST C");     
    noBOutCDataHist->Draw();                                                          // Draw the histogram

    TH1D * inB_0Hist = mh.histFunction(inB_0, "name", "noBOutC Sim");  
    inB_0Hist->GetYaxis()->SetRangeUser(-0.16,0.16);                                   //set y axis
    inB_0Hist->GetXaxis()->SetRangeUser(350,390);                                      //set x axis
    mh.histShift(*inB_0Hist,0);
    inB_0Hist->SetLineColor(kRed);
    inB_0Hist->SetOption("HIST C");     
    inB_0Hist->Draw("HIST C same");
    
    myCanv->SaveAs("plots/noBoutC.jpg");                                         // saves the plot
    myCanv->SaveAs("plots/noBoutC.root");
        
    //make noBoutD plot    
    
    TH1D * noBOutDDataHist = mh.histFunction(noBOutD, "name", "noBOutD Data");  
    noBOutDDataHist->GetYaxis()->SetRangeUser(-0.16,0.16);                            //set y axis
    noBOutDDataHist->GetXaxis()->SetRangeUser(350,390);                               //set x axis
    noBOutDDataHist->SetOption("HIST C");     
    noBOutDDataHist->Draw();                                                          // Draw the histogram

    TH1D * inAHist = mh.histFunction(inA, "name", "noBOutD Sim");  
    inAHist->GetYaxis()->SetRangeUser(-0.16,0.16);                                   //set y axis
    inAHist->GetXaxis()->SetRangeUser(350,390);                                      //set x axis
    mh.histShift(*inAHist,0);
    inAHist->SetLineColor(kRed);
    inAHist->SetOption("HIST C");     
    inAHist->Draw("HIST C same");
    
    myCanv->SaveAs("plots/noBoutD.jpg");                                         // saves the plot
    myCanv->SaveAs("plots/noBoutD.root");        


////////////////////////////////////////////////////////////////////////////////////////////////////////
//Test the anita data type function thing///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////











    return 0;
}




















