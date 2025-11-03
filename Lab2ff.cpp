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
  // gStyle->SetTitleOffset(1, "Y"); // ROOT bug: "x" does y-axis and viceversa
}

Double_t fullEq(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N = par[0];  // N
  Double_t const R = par[1];  // R
  Double_t const b = par[2];  // offset
  Double_t const t0 = par[3]; // tau_0
  Double_t const Lc = par[4]; // Lambda_C
  Double_t const Q = par[5];  // Huff factor

  return N * (Exp(-t / t0) + (1. / R) * Exp((-t) * (Lc + (Q / t0)))) + b;
}

Double_t altFullEq(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N1 = par[0]; // N1
  Double_t const N2 = par[1]; // N2
  Double_t const b = par[2];  // offset
  Double_t const t0 = par[3]; // tau_0
  Double_t const Lc = par[4]; // Lambda_C
  Double_t const Q = par[5];  // Huff factor

  return N1 * (Exp(-t / t0)) + N2 * Exp((-t) * (Lc + (Q / t0))) + b;
}

Double_t redEq(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N = par[0];  // N
  Double_t const b = par[1];  // offset
  Double_t const t0 = par[2]; // tau_0

  return N * Exp(-t / t0) + b;
}

Double_t test_exp(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N = par[0];   // N
  Double_t const tau = par[1]; // tau

  return N * Exp(-t / tau);
}

void Lab2ff() {
  // gStyle->SetOptStat(111111);

  set_LEO_style();

  double const R = 1.21;
  double const lambdaC = 4.4e-3;
  double const Q = 0.975;

  // histogram that contains all true positive events
  TH1D *tp_h = new TH1D();

  double const range_min = 0;
  double const range_max = 16500;

  TFile *res_file = new TFile("result.root");
  if (res_file->IsZombie()) {
    std::cout << "\033[31;1mLEO_ERROR: could not open file \""
              << res_file->GetName() << "\".\033[0m" << '\n';
    return;
  }
  tp_h = (TH1D *)res_file->Get("tp_h_new");
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

  TF1 *fullForm = new TF1("fullForm", fullEq, range_min, range_max, 6);

  ///////                0    1    2     3        4       5
  fullForm->SetParNames("N", "R", "b", "tau0", "labdaC", "Q");
  ///////                0    1    2     3        4       5

  // fullForm->FixParameter(0, tp_count);

  // fullForm->FixParameter(1, R);
  fullForm->SetParameter(1, R);
  fullForm->SetParLimits(1, 0, 2);

  // fullForm->SetParLimits(2,0,16000);
  fullForm->SetParameter(2, 50);

  // fullForm->SetParLimits(3,0,16000);
  fullForm->SetParameter(3, 800);

  // fullForm->SetParLimits(4,0,16000);
  fullForm->SetParameter(4, lambdaC);
  fullForm->SetParLimits(4, 1e-3, 1e-2);

  fullForm->FixParameter(5, Q);
  // fullForm->SetParameter(5, 1);
  // fullForm->SetParLimits(5, 0.7, 1);

  TF1 *altFullForm = new TF1("altFullForm", altFullEq, range_min, range_max, 6);

  ///////                    0     1     2        3         4         5
  altFullForm->SetParNames("N+", "N-", "b_a", "tau0_a", "labdaC_a", "Q_a");
  ///////                    0     1     2        3         4         5

  // altFullForm->FixParameter(0, tp_count);
  altFullForm->SetParameter(0, 2000);

  altFullForm->SetParameter(1, 2000);

  // altFullForm->SetParLimits(2,0,16000);
  altFullForm->SetParameter(2, 50);

  // altFullForm->SetParLimits(3,0,16000);
  altFullForm->SetParameter(3, 800);

  // altFullForm->SetParLimits(4,0,16000);
  altFullForm->SetParameter(4, lambdaC);
  altFullForm->SetParLimits(4, 1e-3, 1e-2);

  // altFullForm->FixParameter(5, Q);
  altFullForm->FixParameter(5, 1);
  // altFullForm->SetParLimits(5, 0.7, 1);

  // Reduced equation, fitting only for tau0 and b.
  TF1 *redForm = new TF1("redForm", redEq, range_min, range_max, 3);

  ///////                0     1     2
  redForm->SetParNames("N_r", "b", "t0_r");
  ///////                0     1     2

  redForm->SetParameter(0, 4000);
  // redForm->SetParLimits(1,0,1000);
  redForm->SetParameter(1, 50);
  // redForm->FixParameter(2, 4);
  redForm->SetParameter(2, 2000);

  fullForm->SetLineColor(kRed);
  altFullForm->SetLineColor(kGreen);
  altFullForm->SetLineStyle(kDashed);
  redForm->SetLineColor(kCyan);

  // TCanvas *canvas1 = new TCanvas("canvas1", "Final Fitting", 720, 720);
  // tp_h->SetFillColor(kPink);
  tp_h->SetMarkerStyle(kFullCircle);
  tp_h->SetFillColor(kBlue);
  tp_h->SetMarkerColor(kBlue);
  TFitResultPtr fullFit = tp_h->Fit(fullForm, "NSR", "");
  std::cout << '\n'
            << '\t' << "tau_C is "
            << 1. /
                   (fullForm->GetParameter(4) + (Q / fullForm->GetParameter(3)))
            << '\n';
  TFitResultPtr altFullFit = tp_h->Fit(altFullForm, "NSR+", "");
  std::cout << '\n'
            << '\t' << "alt_tau_C is "
            << 1. / (altFullForm->GetParameter(4) +
                     (Q / altFullForm->GetParameter(3)))
            << '\n';
  TFitResultPtr redFit = tp_h->Fit(redForm, "NSR+", "", 2000, range_max);

  TF1 *testForm = new TF1("testForm", test_exp, range_min, range_max, 2);
  testForm->SetLineColor(kViolet);
  testForm->SetParNames("N_t", "tau");
  testForm->SetParameter(0, 4000);
  testForm->SetParameter(1, 2000);
  TFitResultPtr testFit = tp_h->Fit(testForm, "NSR+", "", 1000, 8000);

  TCanvas *canvas1 = new TCanvas("canvas1", "Final Fitting", 720, 720);
  canvas1->SetLogy();
  tp_h->Draw();
  fullForm->Draw("same");
  altFullForm->Draw("same");
  canvas1->BuildLegend(.5, .75, .9, .93);

  gPad->Update();

  // TODO
  // add TLatex symbols for parameters (leo: what you mean by that?)
  std::cout << '\n' << "Chi^2:" << '\n';
  std::cout << '\t' << "Full fit: " << fullForm->GetChisquare() << '\n';
  std::cout << '\t' << "Reduced fit: " << redForm->GetChisquare() << '\n';
  // double tau0 = fullForm->GetParameter("tau0");
  // double tauC = fullForm->GetParameter("tauC");
  // std::cout << "Huff Factor (Q): " << '\n';
  // std::cout << (1 / (tau0 + tauC) - lambdaC) * tau0 << '\n';

  canvas1->Print("graphs/Lab2ff/Lab2ff.pdf");
}