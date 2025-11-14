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

#include <iomanip>
#include <iostream>
#include <string>

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
  Double_t const b = par[2];  // background
  Double_t const t0 = par[3]; // tau_0
  Double_t const LC = par[4]; // Lambda_C
  Double_t const Q = par[5];  // Huff factor

  return N * (Exp(-t / t0) + (1. / R) * Exp((-t) * (LC + (Q / t0)))) + b;
}

Double_t fullEq2(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N = par[0];  // N
  Double_t const R = par[1];  // R
  Double_t const b = par[2];  // background
  Double_t const t0 = par[3]; // tau_0
  Double_t const tC = par[4]; // tau_C

  return N * Exp(-t / t0) * (1. + (1. / R) * Exp(-t / tC)) + b;
}

// Double_t fullEq3(Double_t const *x, Double_t const *par) {
//   using TMath::Exp;
//   Float_t const t = x[0];
//   Double_t const N = par[0];  // N
//   Double_t const R = par[1];  // R
//   Double_t const b = par[2];  // background
//   Double_t const t0 = par[3]; // tau_0
//   Double_t const tC = par[4]; // tau_C

//   return N * Exp(-t / t0) * (1. + (1. / R) * Exp(-t / tC)) + b;
// }

Double_t altFullEq(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N1 = par[0]; // N1
  Double_t const N2 = par[1]; // N2
  Double_t const b = par[2];  // background
  Double_t const t0 = par[3]; // tau_0
  Double_t const LC = par[4]; // Lambda_C
  Double_t const Q = par[5];  // Huff factor

  return N1 * (Exp(-t / t0)) + N2 * Exp((-t) * (LC + (Q / t0))) + b;
}

Double_t redEq(Double_t const *x, Double_t const *par) {
  using TMath::Exp;
  Float_t const t = x[0];
  Double_t const N = par[0];  // N
  Double_t const b = par[1];  // background
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

/// @brief Sets the fit mode of the param. R, LambdaC and Q (fix or fit) and
/// prints the result of the fit
/// @param tp_h histogram
/// @param form the function to fit
/// @param R_opt R fit option (1 is fit, 0 is fix)
/// @param LC_opt LambdaC fit option (1 is fit, 0 is fix)
/// @param Q_opt Q fit option (1 is fit, 0 is fix)
void fitNPrint(TH1D *tp_h, TF1 *ff, const bool R_opt, const bool LC_opt,
               const bool Q_opt) {

  double const R = 1.21; // R
  double const R_min = 1.10;
  double const R_max = 1.35;
  double const LC = 4.46e-3; // LambdaC
  double const LC_min = 3e-3;
  double const LC_max = 6e-3;
  double const Q = 0.975; // Q
  double const Q_min = 0.8;
  double const Q_max = 1;

  double const range_min = 0;
  double const range_max = 16500;

  if (!R_opt) {
    ff->FixParameter(1, R);
  } else {
    ff->SetParameter(1, R);
    ff->SetParLimits(1, R_min, R_max);
  }
  if (!LC_opt) {
    ff->FixParameter(4, LC);
  } else {
    ff->SetParameter(4, LC);
    ff->SetParLimits(4, LC_min, LC_max);
  }
  if (!Q_opt) {
    ff->FixParameter(5, Q);
  } else {
    ff->SetParameter(5, Q);
    ff->SetParLimits(5, Q_min, Q_max);
  }

  tp_h->Fit(ff, "QNR", "");
  // calculating tau- and error of tau-
  double tauM;    // tau-
  double tauMErr; // tau- error
  {
    using TMath::Power, TMath::Sqrt;
    double_t const t0ff = ff->GetParameter(3);            // tau0
    double_t const lcff = ff->GetParameter(4);            // lambdaC
    double_t const qff = ff->GetParameter(5);             // Q
    double_t const D2t0ff = Power(ff->GetParError(3), 2); // tau0_err^2
    double_t const D2lcff = Power(ff->GetParError(4), 2); // lambdaC_err^2
    double_t const D2qff = Power(ff->GetParError(5), 2);  // Q_err^2
    double_t const fff = Power(lcff + (qff / t0ff), -4);  // den^(-4)
    tauM = 1. / (lcff + (qff / t0ff));
    tauMErr = Sqrt((fff * D2lcff) + (fff * Power(t0ff, -2) * D2qff) +
                   (fff * Power(qff, 2) * Power(t0ff, -4) * D2t0ff));
  }

  using std::setprecision, std::fixed, std::scientific, std::setw,
      std::to_string;
  std::cout << "FIT" << R_opt << LC_opt << Q_opt; //
  std::cout << '\t' << fixed << setprecision(0) << ff->GetParameter(3) << " ± "
            << ff->GetParError(3); // tau0
  std::cout << '\t' << fixed << setprecision(1) << tauM << " ± " << setw(4)
            << tauMErr; // tau-
  std::cout << '\t' << fixed << setprecision(2) << ff->GetParameter(1);
  if (R_opt)
    std::cout << " ± " << ff->GetParError(1); // R
  else
    std::cout << setw(7) << "";
  std::cout << '\t' << setprecision(2) << scientific << ff->GetParameter(4);
  if (LC_opt)
    std::cout << " ± " << ff->GetParError(4); // LambdaC
  else
    std::cout << setw(11) << "";

  std::cout << '\t' << fixed << setprecision(2) << ff->GetParameter(5);
  if (Q_opt)
    std::cout << " ± " << ff->GetParError(5); // Q
  else
    std::cout << setw(7) << "";
  std::cout << '\t' << fixed << setprecision(2) << ff->GetChisquare() << "/"
            << setprecision(0) << ff->GetNDF(); // Chi^2/NDF
  std::cout << '\t' << fixed << setprecision(1) << setw(4)
            << ff->GetProb() * 100 << "%"; // prob
  std::cout << '\t' << fixed << setprecision(1) << ff->GetParameter(2) << " ± "
            << ff->GetParError(2) << '\n'; // b
}

void Lab2ff() {
  // gStyle->SetOptStat(111111);

  set_LEO_style();

  double const R = 1.21;
  double const lambdaC = 4.4e-3; // ns^-1
  double const Q = 0.975;
  double const tau0 = 2197; // ns
  double const tauM = 206;  // ns, tau- in iron

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

  ///////                0    1    2     3       4        5
  fullForm->SetParNames("N", "R", "b", "t0", "lambda_C", "Q");
  ///////                0    1    2     3       4        5

  fullForm->SetParameter(0, 2000);

  // fullForm->FixParameter(1, R);
  fullForm->SetParameter(1, R);
  fullForm->SetParLimits(1, 0, 2);

  // fullForm->SetParLimits(2,0,16000);
  fullForm->SetParameter(2, 50);

  // fullForm->SetParLimits(3,0,16000);
  fullForm->SetParameter(3, tau0);

  // fullForm->FixParameter(4, lambdaC);
  fullForm->SetParameter(4, lambdaC);
  fullForm->SetParLimits(4, 1e-3, 1e-2);

  // fullForm->FixParameter(5, Q);
  fullForm->SetParameter(5, 1);
  fullForm->SetParLimits(5, 0.7, 1);

  TF1 *fullForm2 = new TF1("fullForm2", fullEq2, range_min, range_max, 5);

  ///////                  0      1      2      3       4
  fullForm2->SetParNames("N_2", "R_2", "b_2", "t0_2", "t-_2");
  ///////                  0      1      2      3       4

  fullForm2->SetParameter(0, 2000);

  fullForm2->FixParameter(1, R);
  // fullForm2->SetParameter(1, R);
  // fullForm2->SetParLimits(1, 0, 2);

  // fullForm->SetParLimits(2,0,16000);
  fullForm2->SetParameter(2, 50);

  // fullForm->SetParLimits(3,0,16000);
  fullForm2->SetParameter(3, tau0);

  // fullForm->FixParameter(4, 200);
  fullForm2->SetParameter(4, 200);
  fullForm2->SetParLimits(4, 0, 400);

  TF1 *fullForm3 = new TF1("fullForm3", fullEq2, range_min, range_max, 5);

  ///////                  0      1      2      3       4
  fullForm3->SetParNames("N_3", "R_3", "b_3", "t0_3", "t-_3");
  ///////                  0      1      2      3       4

  fullForm3->SetParameter(0, 2000);
  // fullForm->FixParameter(1, R);
  fullForm3->SetParameter(1, R);
  fullForm3->SetParLimits(1, 0, 2);

  // fullForm->SetParLimits(2,0,16000);
  fullForm3->SetParameter(2, 50);

  fullForm3->SetParameter(3, tau0);

  fullForm3->FixParameter(4, tauM);

  TF1 *altFullForm = new TF1("altFullForm", altFullEq, range_min, range_max, 6);

  ///////                   0     1      2       3         4         5
  altFullForm->SetParNames("N+", "N-", "b_a", "t0_a", "lambdaC_a", "Q_a");
  ///////                   0     1      2       3         4         5

  altFullForm->SetParameter(0, 2000);

  altFullForm->SetParameter(1, 2000);

  // altFullForm->SetParLimits(2,0,16000);
  altFullForm->SetParameter(2, 50);

  // altFullForm->SetParLimits(3,0,16000);
  altFullForm->SetParameter(3, tau0);

  // altFullForm->SetParLimits(4,0,16000);
  altFullForm->SetParameter(4, lambdaC);
  altFullForm->SetParLimits(4, 1e-3, 1e-2);

  // altFullForm->FixParameter(5, Q);
  altFullForm->FixParameter(5, 1);
  // altFullForm->SetParLimits(5, 0.7, 1);

  // Reduced equation, fitting only for tau0 and b.
  TF1 *redForm = new TF1("redForm", redEq, range_min, range_max, 3);

  ///////                0      1      2
  redForm->SetParNames("N_r", "b_r", "t0_r");
  ///////                0      1      2

  redForm->SetParameter(0, 4000);
  // redForm->SetParLimits(1,0,1000);
  redForm->SetParameter(1, 50);
  // redForm->FixParameter(2, 4);
  redForm->SetParameter(2, 2000);

  fullForm->SetLineColor(kRed);
  fullForm2->SetLineColor(kBlue);
  fullForm3->SetLineColor(kOrange);
  redForm->SetLineColor(kCyan);
  altFullForm->SetLineColor(kGreen);
  altFullForm->SetLineStyle(kDashed);
  redForm->SetLineColor(kCyan);

  // tp_h->SetFillColor(kPink);
  tp_h->SetMarkerStyle(kFullCircle);
  // tp_h->SetFillColor(kBlue);
  tp_h->SetFillColor(kWhite);
  tp_h->SetMarkerColor(kBlack);
  tp_h->SetLineWidth(1);
  tp_h->SetLineColor(kBlack);
  // TFitResultPtr fullFit = tp_h->Fit(fullForm, "NSR", "");
  TFitResultPtr fullFit2 = tp_h->Fit(fullForm2, "NSR+", "");
  TFitResultPtr fullFit3 = tp_h->Fit(fullForm3, "NSR+", "");

  std::cout << "\033[1m" << '\t' << "τ_0" << "\t\t" << "τ-" << "\t\t" << "R"
            << "\t\t" << "Λ_C" << "\t\t\t" << "Q" << "\t\t" << "χ^2/NDF"
            << "\t\t" << "prob.\033[0m" << '\n';
  fitNPrint(tp_h, fullForm, 0, 0, 0); // FIT0
  fitNPrint(tp_h, fullForm, 1, 0, 0); // FIT1
  fitNPrint(tp_h, fullForm, 0, 1, 0); // FIT2
  fitNPrint(tp_h, fullForm, 0, 0, 1); // FIT3
  fitNPrint(tp_h, fullForm, 1, 1, 0); // FIT4
  fitNPrint(tp_h, fullForm, 1, 0, 1); // FIT5
  fitNPrint(tp_h, fullForm, 0, 1, 1); // FIT6
  fitNPrint(tp_h, fullForm, 1, 1, 1); // FIT7

  TFitResultPtr altFullFit = tp_h->Fit(altFullForm, "NSR+", "");
  std::cout << '\n'
            << '\t' << "alt_tau- is "
            << 1. / (altFullForm->GetParameter(4) +
                     (Q / altFullForm->GetParameter(3)))
            << '\n';
  TFitResultPtr redFit = tp_h->Fit(redForm, "NSR+", "", 750, range_max);

  TF1 *testForm = new TF1("testForm", test_exp, range_min, range_max, 2);
  testForm->SetLineColor(kViolet);
  testForm->SetParNames("N_t", "tau");
  testForm->SetParameter(0, 4000);
  testForm->SetParameter(1, 2000);
  TFitResultPtr testFit = tp_h->Fit(testForm, "NSR+", "", 1000, 8000);

  TCanvas *canvas1 = new TCanvas("canvas1", "Final Fitting", 720, 720);
  canvas1->SetLogy();
  tp_h->Draw("pe1"); // "p" points, "e1" error bars with end bars
  fullForm->Draw("same");
  fullForm2->Draw("same");
  fullForm3->Draw("same");
  // altFullForm->Draw("same");
  // fullForm3->GetParameter("");
  canvas1->BuildLegend(.5, .75, .9, .93);

  gPad->Update();

  std::cout << '\n' << "Chi^2:" << '\n';
  std::cout << '\t' << "Full fit: " << fullForm->GetChisquare() << '\n';
  std::cout << '\t' << "Reduced fit: " << redForm->GetChisquare() << '\n';

  double const temp =
      (1. / fullForm2->GetParameter(4)) - (Q / fullForm2->GetParameter(3));
  std::cout << std::setprecision(4) << std::fixed
            << "Fullfit2 Λ_C value: " << std::setw(8) << temp << " (should be "
            << lambdaC << ")" << '\n';
  std::cout << "Fullfit2 τ_C value: " << 1. / temp << " (should be "
            << 1. / lambdaC << ")" << '\n';
  std::cout << fullForm2->GetParameter(3) << '\n';
  canvas1->Print("graphs/Lab2ff/Lab2ff.pdf");
  std::cout << "\033[1;32mLEO_INFO: Files saved!\033[22m If errors appear "
               "then first create a \"graphs/Lab2ff\" folder in your "
               "directory.\033[0m"
            << '\n';
}