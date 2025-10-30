// #include "RooCategory.h"
// #include "RooDataSet.h"
// #include "RooFitResult.h"
// #include "RooHist.h"
// #include "RooPlot.h"
// #include "RooRealVar.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRatioPlot.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>
#include <vector>

// (LEO) Compares two ints from p1 and p2 and returns the smaller one
// bool Isp1Smaller(int pl1, int pl2) {
//   if (pl1 < pl2)
//     return true;
//   else
//     return false;
// }

// (LEO) Static Casts as a DOUBLE the input integer to simplify writing
double scd(const int a) { return static_cast<double>(a); }

// (LEO) takes in 2 values and returns the smaller one calibrated as a double
// double cali_2(const int a, const int b, const int id) {
//   double const p1_cal[2] = {40.6139, 4.003800};
//   double const p2_cal[2] = {43.0851, 4.003570};
//   double const p3_cal[2] = {47.2083, 4.001930};

//   if (a < b) {
//   }
// }

void Lab2hist(const int filepath_option = 0, const bool do_debug = false) {

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

  if (do_debug == 0)
    std::cout << "LEO_INFO: debug mode (=stop at 20 events) not activated"
              << '\n';

  std::string filepath_string;
  switch (filepath_option) {
  case 0:
    std::cout
        << "LEO_ERROR: Please execute the file by typing "
           "\"Lab2hist(<file_id>)\" with file_id = 1 for Makar, 2 for Leo."
        << '\n';
    return;
    break;

  case 1:
    filepath_string = "5b_data_conv.txt";
    break;

  case 2:
    filepath_string = "data/main_data_dec/5b_data_conv_old.txt";
    break;

  default:
    std::cout
        << "LEO_ERROR: invalid filepath_option. Please execute the file by "
           "typing \"Lab2hist(<file_id>)\" with file_id = 1 for Makar, 2 for "
           "Leo."
        << '\n';
    return;
    break;
  }

  // cal[0] = intercept, cal[1] = slope
  double const p1_cal[2] = {40.6139, 4.003800};
  double const p2_cal[2] = {43.0851, 4.003570};
  double const p3_cal[2] = {47.2083, 4.001930};

  std::array<TString, 9> ev_str;
  ev_str[0] = "\"no no no\"";
  ev_str[1] = "\"yes no no\"";
  ev_str[2] = "\"no yes no\"";
  ev_str[3] = "\"no no yes\"";
  ev_str[4] = "\"yes yes no\"";
  ev_str[5] = "\"yes no yes\"";
  ev_str[6] = "\"no yes yes\"";
  ev_str[7] = "\"yes yes yes\"";
  ev_str[8] = "\"yes yes no (FP)\"";

  // collection of all histograms
  std::vector<TH1D *> all_h;
  all_h.reserve(9);
  for (size_t i = 0; i != 9; i++) {
    all_h.push_back(new TH1D(ev_str[i] + " evs",
                             ev_str[i] + " evs;Stop time (ns);Entries", n_bins,
                             min_x, max_x));
  }

  // TODO: VETTORI DI INSTAGRAMI DI DIVISIONE.

  // sum of all histograms
  TH1D *total_h =
      new TH1D("total_h", "all events time distribution;Stop time (ns);Entries",
               n_bins, min_x, max_x);

  // histogram that contains all true positive events
  TH1D *tp_h = new TH1D(
      "tp_h", "all true positive time distribution;Stop time (ns);Entries",
      n_bins, min_x, max_x);

  // distribution of time difference in "yes yes no" events
  TH1D *yyn_diff_h = new TH1D("yyn_diff_h",
                              "t_{PL1} - t_{PL2} in \"yes yes no\" events;time "
                              "difference (ns/4);Entries",
                              40, -20, 20); // bin = 40 is maximum around

  int ev_n, p1, p2, p3;
  double t;
  int n_detections = 0; // number of detections
  // the cases are ("FP" = false positive, "TN" true negative):
  //       0. no no no (TN)
  //       1. yes no no (FP)
  //       2. no yes no (true positive)
  //       3. no no yes (true positive)
  //       4. yes yes no (true positive, |p1-p2| < max_diff)
  //       5. yes no yes (FP)
  //       6. no yes yes (FP)
  //       7. yes yes yes (FP)
  //       8. yes yes no (FP, with |p1-p2| > max_diff)
  int ev_count[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int pn_count[3] = {0, 0, 0}; // indipendent pn counts

  // integer that detects any kind of different case other than triple-fff event
  // bool any_update = false;

  std::ifstream file(filepath_string.c_str());

  double value; // placeholder value

  while (!file.eof()) {
    // read in the values
    file >> ev_n >> t >> p1 >> p2 >> p3;

    if (p1 + p2 + p3 == 12285) { // 0. no no no
      ev_count[0]++;
      continue; // if no stop then skip
    }

    n_detections++;

    if (p1 != max_pn) { // "pn != max_pn" is "yes", otherwise "no"
      // yes xxx xxx
      pn_count[0]++;
      if (p2 != max_pn) {
        // yes yes xxx
        pn_count[1]++;
        if (p3 != max_pn) {
          pn_count[2]++;
          ev_count[7]++; // 7. yes yes yes
          // std::cout << '\t' << ev_n << '\t' << p1 << '\t' << p2 << '\t' << p3
          // << '\n';

          if (p1 <= p2) {
            if (p1 <= p3)
              value = scd(p1) * p1_cal[1] + p1_cal[0];
            else
              value = scd(p3) * p3_cal[1] + p3_cal[0];
          } else {
            if (p2 <= p3)
              value = scd(p2) * p2_cal[1] + p2_cal[0];
            else
              value = scd(p3) * p3_cal[1] + p3_cal[0];
          }

          all_h[7]->Fill(value);
        } else {
          yyn_diff_h->Fill(p1 - p2);

          if (p1 <= p2)
            value = scd(p1) * p1_cal[1] + p1_cal[0];
          else
            value = scd(p2) * p2_cal[1] + p2_cal[0];

          if (TMath::Abs(p1 - p2) <= max_diff) {
            ev_count[4]++; // 4. yes yes no (TP)
            // tp_h->Fill(value);
            all_h[4]->Fill(value);
          } else {
            ev_count[8]++; // 8. yes yes no (FP)
            all_h[8]->Fill(value);
          }
        }
      } else {
        // yes no xxx
        if (p3 != max_pn) {
          pn_count[2]++;
          ev_count[5]++; // 5. yes no yes

          if (p1 <= p3)
            all_h[5]->Fill(scd(p1) * p1_cal[1] + p1_cal[0]);
          else
            all_h[5]->Fill(scd(p3) * p3_cal[1] + p3_cal[0]);
        } else {
          ev_count[1]++; // 1. yes no no
          // value = scd(p1) * p1_cal[1] + p1_cal[0];
          all_h[1]->Fill(scd(p1) * p1_cal[1] + p1_cal[0]);
          // tp_h->Fill(value);
        }
      }
    } else {
      // no xxx xxx
      if (p2 != max_pn) {
        // no yes xxx
        pn_count[1]++;
        if (p3 != max_pn) {
          pn_count[2]++;
          ev_count[6]++; // 6. no yes yes

          if (p2 <= p3)
            all_h[6]->Fill(scd(p2) * p2_cal[1] + p2_cal[0]);
          else
            all_h[6]->Fill(scd(p3) * p3_cal[1] + p3_cal[0]);

        } else {
          ev_count[2]++; // 2. no yes no
          // value = scd(p2) * p2_cal[1] + p2_cal[0];

          all_h[2]->Fill(scd(p2) * p2_cal[1] + p2_cal[0]);
          // tp_h->Fill(value);
        }
      } else {
        // no no xxx
        if (p3 != max_pn) {
          pn_count[2]++;
          ev_count[3]++; // 3. no no yes
          // value = scd(p3) * p3_cal[1] + p3_cal[0];

          // tp_h->Fill(value);
          all_h[3]->Fill(scd(p3) * p3_cal[1] + p3_cal[0]);
        } else { // no no no (check case if all correct)
          if (p1 + p2 + p3 != max_pn * 3) {
            std::cout << "LEO_ERROR: event classification error: p1, p2 "
                         "and p3 are "
                      << p1 << ", " << p2 << ", " << p3 << " at ev. " << ev_n
                      << '\n';
          }
        }
      }
    }

    if (do_debug == true) { // log the info with any update
      std::cout << p1 << '\t' << p2 << '\t' << p3 << '\n';
      // std::cout << n_detections << '\n';

      if (n_detections == debug_detections)
        break; // stop detection

      // any_update = 0; // reset any_update
    }
  }

  // adding no no no histogram ("0" is the underflow bin)
  all_h[0]->SetBinContent(0, ev_count[0]);

  // histogram adding for tp_h and total_h
  tp_h->Add(all_h[1]); // 1. yes no no (TP)
  tp_h->Add(all_h[2]); // 2. no yes no (TP)
  tp_h->Add(all_h[3]); // 3. no no yes (TP)
  tp_h->Add(all_h[4]); // 4. yes yes no (TP)

  for (TH1D *hist : all_h) {
    total_h->Add(hist);
  }

  auto sum_tp = tp_h->GetEntries();
  auto ev_cast = scd(ev_n);
  // auto sum_tp_cast = scd(ev_count[1] + ev_count[2] + ev_count[3] +
  // ev_count[4]);
  // std::cout << "DEBUG: " << sum_tp << "and " << sum_tp << '\n';
  // auto sum_fp_cast =
  //     scd(ev_count[1] + ev_count[5] + ev_count[6] + ev_count[7] +
  //     ev_count[8]);
  auto sum_fp = n_detections - sum_tp;
  auto n_det_cast = scd(n_detections);

  std::cout << "TOTAL EVENTS: " << ev_n << '\n';

  std::cout << "TRUE NEGATIVE:" << '\n'
            << '\t' << "no no no " << '\t' << ev_count[0] << '\t' << "("
            << scd(ev_count[0]) * 100 / ev_cast << "%)" << '\n';

  std::cout << "NON-TRIPLE-FFF (NTF) EVENTS: " << n_detections << " ("
            << n_det_cast * 100 / ev_cast << "%)" << '\n';

  std::cout << "TRUE POSITIVES (" << sum_tp * 100 / ev_cast << "% of total, "
            << sum_tp * 100 / n_det_cast << "% of NTFs):" << '\n'
            << '\t' << "no yes no " << '\t' << ev_count[2] << '\t' << "("
            << scd(ev_count[2]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[2]) * 100 / n_det_cast << "% of NTFs)" << '\n'
            << '\t' << "no no yes " << '\t' << ev_count[3] << '\t' << "("
            << scd(ev_count[3]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[3]) * 100 / n_det_cast << "% of NTFs)" << '\n'
            << '\t' << "yes yes no " << '\t' << ev_count[4] << '\t' << "("
            << scd(ev_count[4]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[4]) * 100 / n_det_cast << "% of NTFs)" << '\n';

  std::cout << "FALSE POSITIVES (" << sum_fp * 100 / ev_cast << "% of total, "
            << sum_fp * 100 / n_det_cast << "% of NTFs):" << '\n'
            << '\t' << "yes no no " << '\t' << ev_count[1] << '\t' << "("
            << scd(ev_count[1]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[1]) * 100 / n_det_cast << "% of NTFs)" << '\n'
            << '\t' << "yes no yes " << '\t' << ev_count[5] << '\t' << "("
            << scd(ev_count[5]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[5]) * 100 / n_det_cast << "% of NTFs)" << '\n'
            << '\t' << "no yes yes " << '\t' << ev_count[6] << '\t' << "("
            << scd(ev_count[6]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[6]) * 100 / n_det_cast << "% of NTFs)" << '\n'
            << '\t' << "yes yes yes " << '\t' << ev_count[7] << '\t' << "("
            << scd(ev_count[7]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[7]) * 100 / n_det_cast << "% of NTFs)" << '\n'
            << '\t' << "yes yes no" << '\t' << ev_count[8] << '\t' << "("
            << scd(ev_count[8]) * 100 / ev_cast << "% of total, "
            << scd(ev_count[8]) * 100 / n_det_cast
            << "% of NTFs, with |time_diff|>" << max_diff * 4 << "ns)" << '\n';

  std::cout << "RAW pn CASES:" << '\n'
            << '\t' << "p1: " << pn_count[0] << '\n'
            << '\t' << "p2: " << pn_count[1] << '\n'
            << '\t' << "p3: " << pn_count[2] << '\n';

  // tp_h->Scale(1. / scd(ev_count[2] + ev_count[3] + ev_count[4]), "");
  // // // for (size_t i = 1; i != 5; i++) {
  // // //   all_h[i]->Scale(1. / ev_count[i], "width");
  // // // }
  // all_h[1]->Scale(1. / scd(ev_count[1]), "width");
  // all_h[2]->Scale(1. / scd(ev_count[2]), "width");
  // all_h[3]->Scale(1. / scd(ev_count[3]), "width");
  // all_h[4]->Scale(1. / scd(ev_count[4]), "width");
  // tp_h->SetFillColor(kBlue);
  // tp_h->SetLineColor(kBlue);
  // tp_h->SetMarkerColor(kBlue);
  // tp_h->SetMarkerStyle(20);
  for (TH1D *hist : all_h) {
    hist->SetMarkerStyle(kFullCircle);
    // hist->SetLineColorAlpha(0, 0.);
  }

  THStack *all_sh =
      new THStack("all_sh", "All events histogram; Stop time (ns); Entries");
  // all_sh->Add(tp_h);
  for (TH1D *hist : all_h) {
    all_sh->Add(hist);
  }

  THStack *all_div_sh = new THStack(
      "all_div_sh", "relative presence histogram; Stop time (ns); Entries");
  for (TH1D *hist : all_h) {
    all_div_sh->Add(hist);
  }

  TCanvas *canvas1 = new TCanvas("canvas1", "Joint histogram", 1280, 720);
  canvas1->Divide(2, 1);
  canvas1->cd(1);

  for (size_t i = 0; i != all_h.size(); i++) {
    // all_h[i]->Scale(1. / all_h[i]->Integral());
  }

  gStyle->SetPalette(kRainBow); // "kRainBow" is not colourblind friendly!
  gPad->SetGrid();
  gPad->SetLogy();
  all_sh->Draw("nostack pmc plc hist p");
  gPad->BuildLegend();

  std::cout << tp_h->Integral(min_x, max_x) << '\n';

  canvas1->cd(2);
  for (TH1D *hist : all_h) {
    hist->Divide(total_h);
  }

  gPad->SetLogy(0);
  all_sh->Draw("pfc plc pmc");
  gPad->BuildLegend();

  // auto ratio_h = new TRatioPlot(tp_h, nny_h);
  // canvas1->SetTicks(0, 1);
  // ratio_h->Draw();
  // ratio_h->GetLowerRefYaxis()->SetRangeUser(0., 2.);
  // ratio_h->GetLowerRefGraph()->GetYaxis()->SetRange(0, 2);
  // ratio_h->GetLowerRefYaxis()->SetTitle("ratio");
  // ratio_h->GetLowYaxis()->SetNdivisions(505);
  // tp_h->Draw();
  // nny_h->Draw("same");
  // add proper legend for histograms
  // ratio_h->GetUpperPad()->cd();

  TCanvas *canvas2 =
      new TCanvas("canvas2", "yes yes no time difference", 720, 720);
  yyn_diff_h->SetFillColor(kBlue);
  yyn_diff_h->Draw();
  // canvas2->SetLogy(0);

  TCanvas *canvas3 =
      new TCanvas("canvas3", "true positive histogram", 720, 720);
  tp_h->SetFillColor(kBlue);
  tp_h->Draw();
  gPad->SetLogy(1);

  // writing everything to file
  TFile *res_f = new TFile("result.root", "RECREATE");
  all_sh->Write();
  tp_h->Write();
  yyn_diff_h->Write();
  res_f->Close();

  gPad->SetLogy(0);
  gPad->SetGrid();
}