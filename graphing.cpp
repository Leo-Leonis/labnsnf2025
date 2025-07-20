// This file reads in the calibration files in order to get the points with
// stardard deviation by linear fitting and to graph them outputting the pdf
// file

#include "TCanvas.h"
#include "TF1.h"
#include "TFunction.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <vector>

Double_t lin_regr(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t a = par[0]; // intercept
  Double_t b = par[1]; // slope

  return a + b * xx;
}

void run(Bool_t do_print = 0) {

  gStyle->SetOptFit(0);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleFont(62, "");
  gStyle->SetTitleSize(.06, "");

  TF1 *f1 = new TF1("f1", lin_regr, 0, 11000, 2);
  TF1 *f2 = new TF1("f2", lin_regr, 0, 11000, 2);
  TF1 *f3 = new TF1("f3", lin_regr, 0, 11000, 2);

  std::vector<TF1 *> f_vector;
  f_vector.push_back(f1);
  f_vector.push_back(f2);
  f_vector.push_back(f3);

  Double_t a = 0; // intercept
  Double_t b = 1; // slope

  auto h_config_lambda = [a, b](TF1 *f) {
    f->SetParameter(0, a);
    f->SetParameter(1, b);

    f->SetParName(0, "intercept");
    f->SetParName(1, "slope");

    // f->SetParLimits(0, -100, 100);
    // f->SetParLimits(1, 0.9943, 1);
  };

  std::for_each(f_vector.begin(), f_vector.end(), h_config_lambda);

  TString ts1("1");
  TString ts2("2");
  TString ts3("3");

  std::vector<TString *> ts_v;
  ts_v.push_back(&ts1);
  ts_v.push_back(&ts2);
  ts_v.push_back(&ts3);

  TGraphErrors *g1 =
      new TGraphErrors("data/cali_tot/cali_tot_1.txt", "%lg %lg %lg", "");
  TGraphErrors *g2 =
      new TGraphErrors("data/cali_tot/cali_tot_2.txt", "%lg %lg %lg", "");
  TGraphErrors *g3 =
      new TGraphErrors("data/cali_tot/cali_tot_3.txt", "%lg %lg %lg", "");

  std::vector<TGraphErrors *> g_v;
  g_v.push_back(g1);
  g_v.push_back(g2);
  g_v.push_back(g3);

  TCanvas *c1 = new TCanvas("c1", "FPGA 1 Calibration", 720, 720);
  TCanvas *c2 = new TCanvas("c2", "FPGA 2 Calibration", 720, 720);
  TCanvas *c3 = new TCanvas("c3", "FPGA 3 Calibration", 720, 720);

  std::vector<TCanvas *> c_v;
  c_v.push_back(c1);
  c_v.push_back(c2);
  c_v.push_back(c3);

  TF1 *res1;
  TF1 *res2;
  TF1 *res3;

  std::vector<TF1 *> res_v;
  res_v.push_back(res1);
  res_v.push_back(res2);
  res_v.push_back(res3);

  TLegend *l1 = new TLegend(.65, 0.12, 0.95, 0.28);
  TLegend *l2 = new TLegend(.65, 0.12, 0.95, 0.28);
  TLegend *l3 = new TLegend(.65, 0.12, 0.95, 0.28);

  std::vector<TLegend *> l_v;
  l_v.push_back(l1);
  l_v.push_back(l2);
  l_v.push_back(l3);

  TH1F *hist;

  for (int i = 0; i < 3; i++) {
    c_v[i]->cd();

    g_v[i]->SetTitle("FPGA " + *ts_v[i] +
                     " calibration;t_{osc} (ns);t_{FPGA} (ns)");
    g_v[i]->SetMarkerStyle(20);
    g_v[i]->SetMarkerColor(kBlue);
    g_v[i]->SetLineWidth(2);
    g_v[i]->SetMarkerSize(1.5);

    g_v[i]->Fit("f" + *ts_v[i], "R");

    c_v[i]->SetGridx();
    c_v[i]->SetGridy();
    g_v[i]->Draw("APE");

    res_v[i] = g_v[i]->GetFunction("f" + *ts_v[i]);
    std::cout << "case " + *ts_v[i] + ":" << '\n'
              << "          Slope: " << res_v[i]->GetParameter(1) << " ± "
              << res_v[i]->GetParError(1) << '\n'
              << "      Intercept: " << res_v[i]->GetParameter(0) << " ± "
              << res_v[i]->GetParError(0) << '\n';

    l_v[i]->AddEntry(g_v[i], "data points", "pe");
    l_v[i]->AddEntry(f1, "linear fit");
    l_v[i]->SetTextSize(0.05);
    l_v[i]->Draw("same");

    hist = g_v[i]->GetHistogram();

    if (hist) {
      // Change the axis title size ticks label sizes
      hist->GetXaxis()->SetTitleSize(0.06);
      hist->GetYaxis()->SetTitleSize(0.06);

      hist->GetXaxis()->SetTitleOffset(0.85);
      hist->GetYaxis()->SetTitleOffset(1.65);

      hist->GetXaxis()->SetLabelSize(0.05);
      hist->GetYaxis()->SetLabelSize(0.05);

      hist->SetTitleSize(0.10, "");
      // hist->SetTitleFont(22);

      // Force a canvas update to see the effects
      c1->Modified();
      c1->Update();
    }
  }

  if (!do_print) {
    std::cout << "\033[1;31mLEO_WARNING: do_print = 0, no outfile will be "
                 "generated.\033[0m"
              << '\n';
  } else {
    c1->Print("cali_1.pdf");
    c2->Print("cali_2.pdf");
    c3->Print("cali_3.pdf");
  }

  return;
}