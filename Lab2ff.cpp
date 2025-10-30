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
#include "TLatex.h"


#include <fstream>
#include <iostream>


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
  const double lambdaC=4.4*pow(10,-6);



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
std::cout<<fullForm->GetParNumber("R")<<std::endl;
std::cout<<fullForm->GetParNumber("N")<<std::endl;
std::cout<<fullForm->GetParNumber("tau0")<<std::endl;
std::cout<<fullForm->GetParNumber("b")<<"\n"<<std::endl;

//N
fullForm->FixParameter(0,tp_count);
//R
fullForm->FixParameter(1,R);
//b
//fullForm->SetParLimits(2,0,16000);
fullForm->SetParameter(2,50);
//tauC
//fullForm->SetParLimits(4,0,16000);
fullForm->SetParameter(4,100);
//tau0
//fullForm->SetParLimits(3,0,16000);
fullForm->SetParameter(3,800);
//Reduced equation, fitting only for tau0 and b.
TF1 *redForm=new TF1("redForm","[N]*exp(-x/[tau])+[b]",4000,17000);
std::cout<<redForm->GetParNumber("t0")<<std::endl;
std::cout<<redForm->GetParNumber("N")<<std::endl;
std::cout<<redForm->GetParNumber("tau")<<std::endl;
std::cout<<redForm->GetParNumber("b")<<std::endl;
//N
redForm->FixParameter(0,tp_count);
//tau0
//redForm->SetParLimits(3,0,1000);
redForm->SetParameter(3,200);
//b
//redForm->SetParLimits(1,0,1000);
redForm->SetParameter(1,50);
//t0
redForm->FixParameter(2,4);

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
tp_h->SetFillColor(kPink);
tp_h->SetMarkerStyle(kFullCircle);
TFitResultPtr fullFit=tp_h->Fit(fullForm,"RS","");
TFitResultPtr redFit=tp_h->Fit(redForm,"RS+","");


TF1 *testForm=new TF1("testForm","[0]*exp(-x/[1])",1000,17000);
  testForm->SetLineColor(kViolet);
  testForm->FixParameter(0,tp_count);
  testForm->SetParameter(1,10000);
TFitResultPtr testFit=tp_h->Fit(testForm,"RS+","");

  
  gStyle->SetPalette(kRainBow); // "kRainBow" is not colourblind friendly!
  canvas1->SetGrid();

  //canvas1->SetLogy(1);

  canvas1->BuildLegend();

gPad->Update();

//TODO
//add TLatex symbols for parameters
//Fix cutoff for redform, the equation that ignores muon capture
//Check redForm equation
std::cout<<"Chi^2:"<<std::endl;
std::cout<<"Full fit: "<<fullForm->GetChisquare()<<std::endl;
std::cout<<"Reduced fit: "<<redForm->GetChisquare()<<std::endl;
double tau0=fullForm->GetParameter("tau0");
double tauC=fullForm->GetParameter("tauC");
double Q=(1/(tau0+tauC)-lambdaC)*tau0;
std::cout<<"Huff Factor (Q): "<<std::endl;
std::cout<<Q<<std::endl;

}