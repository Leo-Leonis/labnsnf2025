#include "TCanvas.h"
#include "RooRealVar.h"
#include "TFormula.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TRatioPlot.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "TFile.h"
#include "TList.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "TMultiGraph.h"

using namespace RooFit;

double scd(const int a) { return static_cast<double>(a); }

void Lab2ff(){
  gStyle->SetOptStat(111111);
  const int n_bins = 50;
  // (only in debug mode) decides the number of non triple-fff events being
  // recorded and in all these cases the p_n values are printed in the terminal
  const int debug_detections = 1000;
  // maximum value of each pn count
  const int max_pn = 4095;
  // max difference in absolute value between coincident counts (events in which
  // there are more than 2 positive stop signals)
  const int max_diff = 3;
  const int min_x = 0;
  const int max_x = 16500;
  double const R=1.21;
  const double Q=0.975;



  // histogram that contains all true positive events
  TH1D *tp_h = new TH1D();

TFile *result=new TFile("result.root") ;
tp_h=(TH1D*)result->Get("tp_h");
int tp_count=tp_h->GetEntries();
//Fitting
//Warning: cannot call variable t in TFormula (TF1) as it apparently treats it as a fourth
//variable even if no other variables are present and contrary to the reference guide
//Full equation, fitting for tau0, tauc/tau_mu^minus, and b

TF1 *fullForm=new TF1("fullForm","[N]*exp(-x/[tau0])*(1+exp(-x/[tauC])/[R])+[b]",0,17000);

std::cout<<fullForm->GetParNumber("tauC")<<std::endl;
std::cout<<fullForm->GetParNumber("t0")<<std::endl;
std::cout<<fullForm->GetParNumber("R")<<std::endl;
std::cout<<fullForm->GetParNumber("N")<<std::endl;
std::cout<<fullForm->GetParNumber("tau0")<<std::endl;
std::cout<<fullForm->GetParNumber("b")<<"\n"<<std::endl;

//N
fullForm->FixParameter(0,tp_count);
//t0
fullForm->FixParameter(3,4);
//R
fullForm->FixParameter(1,R);
//b
fullForm->SetParLimits(2,0,16000);
//tauC
fullForm->SetParLimits(5,0,16000);
//tau0
fullForm->SetParLimits(4,0,16000);
//Reduced equation, fitting only for tau0 and b.
TF1 *redForm=new TF1("redForm","[N]*exp(-[t0]/[tau0])*exp(-x/[tau0])+[b]",4000,17000);
std::cout<<redForm->GetParNumber("t0")<<std::endl;
std::cout<<redForm->GetParNumber("N")<<std::endl;
std::cout<<redForm->GetParNumber("tau0")<<std::endl;
std::cout<<redForm->GetParNumber("b")<<std::endl;
//N
redForm->FixParameter(0,tp_count);
//tau0
redForm->SetParLimits(3,0,1000);
//b
fullForm->SetParLimits(1,0,1000);
//t0
fullForm->FixParameter(2,4);

  gROOT->SetStyle("Modern");

  fullForm->SetLineColor(kOrange);
  redForm->SetLineColor(kCyan);

    int f=4;
  int b=23;

  gStyle->SetCanvasColor(b);
    gStyle->SetTitleFillColor(b);
  gStyle->SetStatColor(b);
    gStyle->SetFrameLineColor(f);
  gStyle->SetGridColor(f);
  gStyle->SetStatTextColor(f);
  gStyle->SetTitleTextColor(f);
  gStyle->SetLabelColor(f,"xy");
  gStyle->SetTitleColor(f,"xy");
  gStyle->SetAxisColor(f,"xy");
TCanvas *canvas1 = new TCanvas("canvas1", "Final Fitting", 1280, 720);
//TMultiGraph *collection=new TMultiGraph();
tp_h->SetFillColor(kPink);
tp_h->SetMarkerStyle(kFullCircle);
TFitResultPtr fullFit=tp_h->Fit(fullForm,"RS","");
TFitResultPtr redFit=tp_h->Fit(redForm,"RS+","");


TF1 *testForm=new TF1("testForm","[e]*exp(-x/[s])",0,17000);
  testForm->SetLineColor(kViolet);
  testForm->FixParameter(0,tp_count);
TFitResultPtr testFit=tp_h->Fit(testForm,"RS+","");

  
  gStyle->SetPalette(kRainBow); // "kRainBow" is not colourblind friendly!
  canvas1->SetGrid();

  //canvas1->SetLogy(1);

  canvas1->BuildLegend();

gPad->Update();


}