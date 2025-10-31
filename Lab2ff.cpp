#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRatioPlot.h"
#include "TStyle.h"

#include <iostream>

double scd(const int a) { return static_cast<double>(a); }

// (LEO) sets makar's default style of the graphs produced
void set_MAKAR_style() {
  gROOT->SetStyle("Modern");

  int f = 4;
  int b = 23;

  gStyle->SetCanvasColor(b);
  gStyle->SetTitleFillColor(b);
  gStyle->SetStatColor(b);
  gStyle->SetFrameLineColor(f);
  gStyle->SetGridColor(f);
  gStyle->SetStatTextColor(f);
  gStyle->SetTitleTextColor(f);
  gStyle->SetLabelColor(f, "xy");
  gStyle->SetTitleColor(f, "xy");
  gStyle->SetAxisColor(f, "xy");
}

// (LEO) Sets Leo's (superior) default style of the graphs produced
void set_LEO_style() {
  // gStyle->SetOptStat(111110);
  gStyle->SetOptStat(0);

  // pad
  gStyle->SetPadLeftMargin(0.12);
  // gStyle->SetPadRightMargin(1);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.1);

  // title
  gStyle->SetTitleFont(62, "");
  gStyle->SetTitleSize(.05, "");

  // axis
  gStyle->SetAxisMaxDigits(4);
  gStyle->SetStripDecimals(kFALSE); // all numbers with same number of decimals
  gStyle->SetLabelSize(.05, "X");
  gStyle->SetLabelSize(.05, "Y");
  gStyle->SetTitleSize(.05, "X");
  gStyle->SetTitleSize(.05, "Y");
  // gStyle->SetTitleOffset(1.3, "X"); // ROOT bug: "x" does y-axis and
  // viceversa
  gStyle->SetTitleOffset(1.3, "Y"); // ROOT bug: "x" does y-axis and viceversa
}

Double_t full_eq(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t N = par[0];  // N
  Double_t t0 = par[0]; // tau_0
  Double_t R = par[1];  // R
  Double_t tc = par[2]; // tau_C
  Double_t b = par[3];  // b

  return N * (TMath::Exp(xx / b) - 1);
}

void Lab2ff() {
  // gStyle->SetOptStat(111111);

  set_LEO_style();

  double const R = 1.21;
  double const lambdaC = 4.4 * pow(10, -6);

  // histogram that contains all true positive events
  TH1D *tp_h = new TH1D();

  TFile *res_file = new TFile("result.root");
  if (res_file->IsZombie()) {
    std::cout << "\033[31;1mLEO_ERROR: could not open file \""
              << res_file->GetName() << "\".\033[0m" << '\n';
    return;
  }
  tp_h = (TH1D *)res_file->Get("tp_h");
  if (!tp_h) {
    std::cout
        << "\033[31;1mLEO_ERROR: could not find histogram \"tp_h\".\033[0m"
        << '\n';
    return;
  }
  int tp_count = tp_h->GetEntries();
  // Fitting
  // Warning: cannot call variable t in TFormula (TF1) as it apparently treats
  // it as a fourth variable even if no other variables are present and contrary
  // to the reference guide Full equation, fitting for tau0, tauc/tau_mu^minus,
  // and b

  TF1 *fullForm = new TF1(
      "fullForm", "[N]*exp(-x/[tau0])*(1+exp(-x/[tauC])/[R])+[b]", 0, 17000);

  std::cout << fullForm->GetParNumber("tauC") << '\n';
  std::cout << fullForm->GetParNumber("R") << '\n';
  std::cout << fullForm->GetParNumber("N") << '\n';
  std::cout << fullForm->GetParNumber("tau0") << '\n';
  std::cout << fullForm->GetParNumber("b") << "\n" << '\n';

  // N
  // fullForm->FixParameter(0, tp_count);
  // R
  fullForm->FixParameter(1, R);
  // b
  // fullForm->SetParLimits(2,0,16000);
  fullForm->SetParameter(2, 50);
  // tauC
  // fullForm->SetParLimits(4,0,16000);
  fullForm->SetParameter(4, 100);
  // tau0
  // fullForm->SetParLimits(3,0,16000);
  fullForm->SetParameter(3, 800);

  // Reduced equation, fitting only for tau0 and b.
  TF1 *redForm = new TF1("redForm", "[N]*exp(-x/[tau])+[b]", 4000, 8000);
  std::cout << redForm->GetParNumber("t0") << '\n';
  std::cout << redForm->GetParNumber("N") << '\n';
  std::cout << redForm->GetParNumber("tau") << '\n';
  std::cout << redForm->GetParNumber("b") << '\n';
  // N
  // redForm->FixParameter(0, tp_count);
  // tau0
  // redForm->SetParLimits(3,0,1000);
  redForm->SetParameter(3, 200);
  // b
  // redForm->SetParLimits(1,0,1000);
  redForm->SetParameter(1, 50);
  // t0
  redForm->FixParameter(2, 4);

  fullForm->SetLineColor(kRed);
  redForm->SetLineColor(kCyan);

  TCanvas *canvas1 = new TCanvas("canvas1", "Final Fitting", 720, 720);
  // tp_h->SetFillColor(kPink);
  tp_h->SetMarkerStyle(kFullCircle);
  tp_h->SetMarkerColor(kBlue);
  TFitResultPtr fullFit = tp_h->Fit(fullForm, "RS", "");
  TFitResultPtr redFit = tp_h->Fit(redForm, "RS+", "");

  TF1 *testForm = new TF1("testForm", "[0]*exp(-x/[1])", 1000, 17000);
  testForm->SetLineColor(kViolet);
  testForm->FixParameter(0, tp_count);
  testForm->SetParameter(1, 10000);
  TFitResultPtr testFit = tp_h->Fit(testForm, "RS+", "");

  canvas1->SetLogy();

  canvas1->BuildLegend(.5, .75, .9, .93);

  gPad->Update();

  // TODO
  // add TLatex symbols for parameters
  // Fix cutoff for redform, the equation that ignores muon capture
  // Check redForm equation
  std::cout << "Chi^2:" << '\n';
  std::cout << '\t' << "Full fit: " << fullForm->GetChisquare() << '\n';
  std::cout << '\t' << "Reduced fit: " << redForm->GetChisquare() << '\n';
  double tau0 = fullForm->GetParameter("tau0");
  double tauC = fullForm->GetParameter("tauC");
  // double Q = (1 / (tau0 + tauC) - lambdaC) * tau0;
  std::cout << "Huff Factor (Q): " << '\n';
  // std::cout << Q << '\n';
  std::cout << (1 / (tau0 + tauC) - lambdaC) * tau0 << '\n';

  canvas1->Print("graphs/Lab2ff/Lab2ff.pdf");
}