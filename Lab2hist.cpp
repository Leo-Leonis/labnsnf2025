#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRatioPlot.h"
#include "TStyle.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

/// @brief (LEO) Imports a TGraphs and reads the Y values of each point, then
/// calculates and prints the value of the weighted mean and its error, returns
/// only the weighted mean
/// @param g The TGraph pointer from which the values are taken in
/// @return Weighted mean
double getWeightedMeanY(TGraph *g) {
  // number of bins of g
  const int nPoints = g->GetN();
  // sum of weights (w = 1/err^2 = err^-2)
  double_t sum_weights = 0;
  // ratio numerator
  double_t num = 0;

  for (int i = 0; i != nPoints; i++) {
    // current bin content
    double_t const i_bin_content = g->GetPointY(i);
    // current bin error
    double_t const i_bin_error = g->GetErrorY(i + 1);
    // std::cout << i_bin_content << " and " << i_bin_error << '\n'; // debug
    if (i_bin_content == 0. & i_bin_error == 0.) {
      // std::cout << "skipped" << '\n'; // debug
      continue;
    } else {
      // current weight
      double_t const i_weight = TMath::Power(i_bin_error, -2);
      sum_weights += i_weight;
      num += i_bin_content * i_weight;
      // std::cout << '\t' << sum_weights << " and " << num << '\n'; // debug
    }
  }
  // the weighted mean
  double_t const ratio = num / sum_weights;
  // the w. mean error
  double_t const ratio_err = TMath::Power(sum_weights, -0.5);

  std::cout << "weighted mean: " << ratio << " +/- " << ratio_err << '\n';

  return ratio;
}

// (LEO) short hand for static_cast<double>(a)
double scd(const int a) { return static_cast<double>(a); }

// (LEO) Sets the default style of the graphs produced
void set_style() {
  // gStyle->SetOptStat(111110);
  gStyle->SetOptStat(0);

  // pad
  gStyle->SetPadLeftMargin(0.13);
  // gStyle->SetPadRightMargin(0.1);
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
  gStyle->SetTitleOffset(1.3, "Y"); // ROOT bug: "x" does y-axis and viceversa
}

// cal[0] = intercept, cal[1] = slope
double const p1_cal[2] = {40.6139, 4.003800};
double const p2_cal[2] = {43.0851, 4.003570};
double const p3_cal[2] = {47.2083, 4.001930};

// cal[0] = intercept, cal[1] = slope
double const p1_cal_new[2] = {-10.1525, 0.249765};
double const p2_cal_new[2] = {-10.7402, 0.249769};
double const p3_cal_new[2] = {-11.8005, 0.249881};

/// @brief (LEO) returns the OLD calibrated values (t_osc dependent, FPGA-counts
/// independent)
/// @param pi the value to be calibrated
/// @param id (1 = PL1, 2 = PL2, 3 = PL3, other will cause abort)
/// @return the calibrated value
double cal_val_old(double const pi, int const id) {
  switch (id) {
  case 1:
    return scd(pi) * p1_cal[1] + p1_cal[0];
    break;
  case 2:
    return scd(pi) * p2_cal[1] + p2_cal[0];
    break;
  case 3:
    return scd(pi) * p3_cal[1] + p3_cal[0];
    break;
  default:
    std::cout
        << "\033[1;31mLEO_ERROR: id input of cal_val_old(x, id) is out of "
           "bounds (y is "
        << id << ", possible values are 1, 2 and 3.)\033[0m" << '\n';
    assert(0);
    return 0.;
    break;
  }
}

/// @brief (LEO) returns the NEW calibrated values (t_osc independent,
/// FPGA-counts dependent)
/// @param pi the value to be calibrated
/// @param id (1 = PL1, 2 = PL2, 3 = PL3, other will cause abort)
/// @return the calibrated value
double cal_val_new(int const pi, int const id) {
  switch (id) {
  case 1:
    return (scd(pi) - p1_cal_new[0]) / p1_cal_new[1];
    break;
  case 2:
    return (scd(pi) - p2_cal_new[0]) / p2_cal_new[1];
    break;
  case 3:
    return (scd(pi) - p3_cal_new[0]) / p3_cal_new[1];
    break;
  default:
    std::cout
        << "\033[1;31mLEO_ERROR: id input of cal_val_new(x, id) is out of "
           "bounds (y is "
        << id << ", possible values are 1, 2 and 3.)\033[0m" << '\n';
    assert(0);
    return 0.;
    break;
  }
}

void Lab2hist(const int filepath_option = 0, const bool do_print = false,
              const bool do_debug = false) {

  set_style();

  const int n_bins = 50;
  // (only in debug mode) decides the number of non triple-fff events being
  // recorded and in all these cases the p_n values are printed in the terminal
  const int debug_detections = 1000;
  // maximum value of each pn count
  const int max_pn = 4095;
  // max difference in absolute value between coincident counts (events in which
  // there are more than 2 positive stop signals) in nanoseconds
  const double max_diff = 6.5;
  const double mean_diff = -2.677; // mean from yyn_diff_th histogram
  const int min_x = 0;
  const int max_x = 16500;

  std::string filepath_string;
  switch (filepath_option) {
  case 0:
    std::cout << "\033[1;31mLEO_ERROR: Please execute the file by typing "
                 "\"Lab2hist(<file_id>)\" with file_id = 1 for Makar, 2 for "
                 "Leo.\033[0m"
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
        << "\033[31;1mLEO_ERROR: invalid filepath_option. Please execute the "
           "file by typing \"Lab2hist(<file_id>)\" with file_id = 1 for "
           "Makar, "
           "2 for "
           "Leo.\033[0m"
        << '\n';
    return;
    break;
  }

  if (do_debug == 0)
    std::cout
        << "\033[33;1mLEO_INFO: NO debug mode (=stop at 20 events).\033[22m "
           "Activate debug mode by executing Lab2hist(x,x,1) \033[0m"
        << '\n';

  // array of the strings that contain the names of the histograms
  std::array<TString, 9> ev_str;
  ev_str[0] = "\"no no no\"";
  ev_str[1] = "\"yes no no\"";
  ev_str[2] = "\"no yes no\"";
  ev_str[3] = "\"no no yes\"";
  ev_str[4] = "\"yes yes no (acc.)\"";
  ev_str[5] = "\"yes no yes\"";
  ev_str[6] = "\"no yes yes\"";
  ev_str[7] = "\"yes yes yes\"";
  ev_str[8] = "\"yes yes no (rej.)\"";

  // collection of all histograms
  std::vector<TH1D *> all_h_v;
  all_h_v.reserve(9);
  for (size_t i = 0; i != 9; i++) {
    all_h_v.push_back(
        new TH1D(ev_str[i] + " evs",
                 ev_str[i] + " evs;Stop time (ns);Entries / (330 ns^{-1})",
                 n_bins, min_x, max_x));
  }

  // sum of all histograms
  TH1D *total_h = new TH1D(
      "total_h",
      "all events time histogram;Stop time (ns);Entries / (330 ns^{-1})",
      n_bins, min_x, max_x);

  // histogram of time difference in "yes yes no" events in counts
  TH1D *yyn_diff_ch = new TH1D("yyn_diff_ch",
                               "t_{PL1} - t_{PL2} in \"yes yes no\" evs;Time "
                               "difference (FPGA counts);Entries",
                               41, -20.5, 20.5); // bin = 40 is maximum

  // histogram of time difference in "yes yes no" events in effective time
  TH1D *yyn_diff_th = new TH1D("yyn_diff_th",
                               "t_{PL1} - t_{PL2} in \"yes yes no\" evs;Time "
                               "difference (ns);Entries",
                               31, -62, 62); // bin = 40 is maximum around

  int ev_n, p1, p2, p3; // event number
  double t;             // time of event
  int n_detections = 0; // number of detections
  // the cases are ("FP" = false positive, "TN" true negative):
  //       0. no no no (TN)
  //       1. yes no no (FP)
  //       2. no yes no (TP)
  //       3. no no yes (TP)
  //       4. yes yes no (TP, |p1-p2| < max_diff)
  //       5. yes no yes (FP)
  //       6. no yes yes (FP)
  //       7. yes yes yes (FP)
  //       8. yes yes no (FP, with |p1-p2| > max_diff)
  int ev_count[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int pn_count[3] = {0, 0, 0}; // indipendent pn counts

  // integer that detects any kind of different case other than triple-fff event
  // bool any_update = false;

  // input file
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

          // fill the smallest value
          if (p1 <= p2) {
            if (p1 <= p3)
              value = cal_val_new(p1, 1);
            else
              value = cal_val_new(p3, 3);
          } else {
            if (p2 <= p3)
              value = cal_val_new(p2, 2);
            else
              value = cal_val_new(p3, 3);
          }

          all_h_v[7]->Fill(value);
        } else { // "yes yes no" case, we need to identify the time difference
          yyn_diff_ch->Fill(p1 - p2);
          double const p1_calibrated = cal_val_new(p1, 1);
          double const p2_calibrated = cal_val_new(p2, 2);
          double const p1p2diff = p1_calibrated - p2_calibrated;
          yyn_diff_th->Fill(p1p2diff);

          value = p2_calibrated; // p2 will be the reference value to be filled

          // yes yes no ev selection criteria
          if (TMath::Abs(p1p2diff - mean_diff) <= max_diff) {
            ev_count[4]++; // 4. yes yes no (TP)
            // tp_h->Fill(value);
            all_h_v[4]->Fill(value);
          } else {
            ev_count[8]++; // 8. yes yes no (FP)
            all_h_v[8]->Fill(value);
          }
        }
      } else {
        // yes no xxx
        if (p3 != max_pn) {
          pn_count[2]++;
          ev_count[5]++; // 5. yes no yes

          if (p1 <= p3)
            all_h_v[5]->Fill(cal_val_new(p1, 1));
          else
            all_h_v[5]->Fill(cal_val_new(p3, 3));
        } else {
          ev_count[1]++; // 1. yes no no
          // value = cal_val_new(p1,1);
          all_h_v[1]->Fill(cal_val_new(p1, 1));
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
            all_h_v[6]->Fill(cal_val_new(p2, 2));
          else
            all_h_v[6]->Fill(cal_val_new(p3, 3));

        } else {
          ev_count[2]++; // 2. no yes no
          // value = cal_val_new(p2,2);

          all_h_v[2]->Fill(cal_val_new(p2, 2));
          // tp_h->Fill(value);
        }
      } else {
        // no no xxx
        if (p3 != max_pn) {
          pn_count[2]++;
          ev_count[3]++; // 3. no no yes
          // value = cal_val_new(p3, 3);

          // tp_h->Fill(value);
          all_h_v[3]->Fill(cal_val_new(p3, 3));
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
      std::cout << value << '\n';

      if (n_detections == debug_detections)
        break; // stop detection

      // any_update = 0; // reset any_update
    }
  }

  // adding no no no histogram ("0" is the underflow bin)
  all_h_v[0]->SetBinContent(0, ev_count[0]);

  // histogram that contains all accepted NFTs (w/ ynn) (legacy)
  TH1D *tp_h = new TH1D(
      "tp_h",
      "Accepted NTFs VS all NTFs;Stop time (ns);Entries / (330 ns^{-1})",
      n_bins, min_x, max_x);
  tp_h->Add(all_h_v[1]); // 1. yes no no
  tp_h->Add(all_h_v[2]); // 2. no yes no
  tp_h->Add(all_h_v[3]); // 3. no no yes
  tp_h->Add(all_h_v[4]); // 4. yes yes no

  // histogram that contains all acc. NTFs (w/o ynn)
  TH1D *tp_h_new = new TH1D(
      "tp_h_new",
      "Accepted NTFs VS all NTFs;Stop time (ns);Entries / (330 ns^{-1})",
      n_bins, min_x, max_x);
  tp_h_new->Add(all_h_v[2]); // 2. no yes no
  tp_h_new->Add(all_h_v[3]); // 3. no no yes
  tp_h_new->Add(all_h_v[4]); // 4. yes yes no

  // add all histograms to total_h
  for (TH1D *hist : all_h_v) {
    total_h->Add(hist);
  }

  auto const sum_tp = tp_h_new->GetEntries(); // sum of events of acc. NTFs
  auto const ev_cast = scd(ev_n);             // event number but in double
  // auto const sum_tp_cast = scd(ev_count[1] + ev_count[2] + ev_count[3] +
  // ev_count[4]);
  // std::cout << "DEBUG: " << sum_tp << "and " << sum_tp << '\n';
  // auto const sum_fp_cast =
  //     scd(ev_count[1] + ev_count[5] + ev_count[6] + ev_count[7] +
  //     ev_count[8]);
  auto const sum_fp = n_detections - sum_tp; // sum of events of rej. NTFs
  auto const n_det_cast = scd(n_detections); // number of NTFs

  std::cout << "TOTAL EVENTS: " << ev_n << '\n';

  std::cout << "TRIPLE-FFF EVENTS:" << '\n'
            << '\t' << "no no no " << '\t' << ev_count[0] << '\t' << "("
            << scd(ev_count[0]) * 100 / ev_cast << "%)" << '\n';

  std::cout << "NON-TRIPLE-FFF (NTF) EVENTS: " << n_detections << " ("
            << n_det_cast * 100 / ev_cast << "%)" << '\n';

  std::cout << "ACCEPTED NTFs (" << sum_tp * 100 / ev_cast << "% of total, "
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

  std::cout << "REJECTED NTFs (" << sum_fp * 100 / ev_cast << "% of total, "
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
            << "% of NTFs, with |time_diff|>" << max_diff << "ns)" << '\n';

  std::cout << "RAW pn CASES:" << '\n'
            << '\t' << "p1: " << pn_count[0] << '\n'
            << '\t' << "p2: " << pn_count[1] << '\n'
            << '\t' << "p3: " << pn_count[2] << '\n';

  std::cout << "sum_tp = " << sum_tp << " and sum_fp = " << sum_fp << '\n';

  // graphical options for the histograms
  for (TH1D *hist : all_h_v) {
    hist->SetMarkerStyle(kFullCircle);
    // hist->SetMarkerSize(1.4);
    // hist->SetLineColorAlpha(0, 0.);
  }

  // the stacked histogram
  THStack *all_sh = new THStack(
      "all_sh",
      "All events histogram; Stop time (ns); Entries / (330 ns^{-1})");
  // all_sh->Add(tp_h);
  for (size_t i = 1; i != all_h_v.size(); i++) {
    all_sh->Add(all_h_v[i]);
  }

  /////////////////////// canvas 1 ///////////////////////////////////
  // this will draw the double canvas (the colorful one)

  TCanvas *canvas1 = new TCanvas("canvas1", "Joint histogram", 0, 0, 1440, 720);
  canvas1->Divide(2, 1, .001); // divide in in 2 horizonal parts w/ 0.001 margin
  canvas1->cd(1); // go to the first subdivision (the left one in this case)

  gStyle->SetPalette(kBird); // "kRainBow" is not colourblind friendly!
  gPad->SetLogy();
  all_sh->Draw("nostack pmc plc pfc hist p1");
  gPad->BuildLegend(.53 + 0.09, .62 - 0.05, .9 + 0.09, .93 - 0.05, "", "f");
  // canvas1->SetFillColor(kRed); // debug

  canvas1->cd(2); // go to the second subdivision (the right one in this case)
  // collection of all ratio histograms
  std::vector<TH1D *> all_div_h_v;
  all_div_h_v.reserve(all_h_v.size());

  // stacked histogram of all ratio histograms
  THStack *all_div_sh = new THStack(
      "all_div_sh",
      "Relative presence histogram; Stop time (ns); Ratio / (330 ns^{-1})");

  // create clone histograms, divide them by total, then add them to the stacked
  // histogram (sh), but first managing the first (i=0) histogram out of sh
  // because it's the no-no-no evs
  all_div_h_v.push_back(new TH1D(*all_h_v[0]));
  all_div_h_v[0]->Divide(total_h);
  for (size_t i = 1; i != all_h_v.size(); i++) {
    all_div_h_v.push_back(new TH1D(*all_h_v[i]));
    all_div_h_v[i]->Divide(total_h);
    all_div_sh->Add(all_div_h_v[i]);
  }
  // all_div_sh->GetHistogram()->GetYaxis()->SetLabelOffset(.01);
  gPad->SetLogy(0);
  all_div_sh->Draw("pfc plc pmc");
  // gPad->BuildLegend();
  all_div_sh->GetYaxis()->SetLabelOffset(0.01);

  /////////////////////// canvas 2 ///////////////////////////////////

  TCanvas *canvas2 = new TCanvas("canvas2", "\"yes yes no\" time difference", 0,
                                 720, 720, 720);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.07);
  // yyn_diff_ch->SetTitleOffset(1.5, "X");
  yyn_diff_th->SetFillColor(kBlue);
  yyn_diff_th->Draw();
  std::cout << "yyn_diff_th->GetMean() = " << yyn_diff_th->GetMean() << " +/- "
            << yyn_diff_th->GetMeanError() << '\n'
            << "yyn_diff_th->GetStdDev() = " << yyn_diff_th->GetStdDev()
            << " +/- " << yyn_diff_th->GetStdDevError() << '\n';
  canvas2->SetLogy();

  /////////////////////// canvas 3 ///////////////////////////////////

  TCanvas *canvas3 =
      new TCanvas("canvas3", "true positive histogram", 720, 720, 720, 720);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.13);
  total_h->SetFillColor(kRed);
  total_h->SetLineWidth(0);
  total_h->SetTitle("True and false positives stop time histogram");
  // total_h->GetYaxis()->SetTitleOffset(.9);
  total_h->Draw();
  // tp_h->SetFillColor(kBlue);
  // tp_h->SetLineWidth(0);
  // tp_h->Draw("same");
  tp_h_new->SetFillColor(kBlue);
  tp_h_new->SetLineWidth(0);
  tp_h_new->Draw("same");

  TLegend *leg3 = new TLegend(.7, .8, .9, .93);
  leg3->AddEntry(tp_h_new, "true positives", "f");
  // leg3->AddEntry(tp_h, "true positives (old)", "f");
  leg3->AddEntry(total_h, "false positives", "f");
  // NB: here total_h is painted before and tp_h painted after, so there will be
  // tp_h and the difference between total_h and tp_h, i.e. a histogram stacked
  // on top of tp_h which contains the false positives
  leg3->Draw("same");

  gPad->RedrawAxis(); // redraw because "same" printed on top
  gPad->RedrawAxis("G");

  /////////////////////// canvas 4 ///////////////////////////////////

  TCanvas *canvas4 = new TCanvas("canvas4", "Ratio", 1440, 720, 720, 720);
  auto ratio_h = new TRatioPlot(tp_h_new, total_h);
  gPad->SetLogy();
  gPad->SetBottomMargin(1.05);
  ratio_h->SetLeftMargin(0.13);
  total_h->SetMarkerColor(kRed);
  total_h->SetMarkerStyle(kWhite);
  total_h->SetLineWidth(1);
  total_h->SetLineColor(kRed);
  total_h->SetMarkerStyle(kFullCircle);

  ratio_h->SetH2DrawOpt("e1"); // drawing option for total_h ("H2")
  ratio_h->GetLowYaxis()->SetNdivisions(505);
  ratio_h->SetSeparationMargin(0.03); // margin between low and top pad
  ratio_h->Draw(); // hides the horizontal dashed lines in lower plot

  auto ratio_g = ratio_h->GetLowerRefGraph(); // ratio graph
  std::cout << "Ratio ";
  // auto w_mean = getWeightedMeanY(ratio_g);
  auto mean_line = std::vector<double>{getWeightedMeanY(ratio_g)};
  ratio_g->SetMarkerStyle(kFullCircle);
  ratio_h->SetGridlines(mean_line);
  ratio_h->Paint(); // updatew grid lines

  ratio_h->GetLowerRefGraph()->SetMinimum(0.3); // lower (ratio) range
  ratio_h->GetLowerRefGraph()->SetMaximum(1.1);
  // ratio_h->GetUpperRefYaxis()->SetTitleOffset(1);   // upper y axis
  ratio_h->GetLowerRefYaxis()->SetTitle("Ratio");   // lower y axis
  ratio_h->GetLowerRefXaxis()->SetTitleOffset(0.9); // lower x axis

  ratio_h->GetUpperPad()->cd();
  TLegend *leg4 = new TLegend(.65, .67, .9, .9);
  leg4->AddEntry(tp_h_new, "Acc. NTFs", "f");
  leg4->AddEntry(total_h, "All NTFs", "pe1");
  leg4->Draw("same");

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  // generate .root file and save graphs if do_print is 1
  if (!do_print) {
    std::cout
        << "\033[1;31mLEO_INFO: graphs and .root files won't be saved!\033[22m "
           "To save the files set do_print to 1 by executing with "
           "Lab2hist(x,1).\033[0m"
        << '\n';
  } else {
    // saving all histograms in one .root file
    TFile *res_f = new TFile("results/result_6.5_adj.root", "RECREATE");
    all_sh->Write();
    total_h->Write();
    tp_h->Write();
    tp_h_new->Write();
    yyn_diff_ch->Write();
    yyn_diff_th->Write();
    all_h_v[4]->Write();
    res_f->Close();

    // saving all canvases produced in different pdf's
    canvas1->Print("graphs/Lab2hist/all_h.pdf");
    canvas2->Print("graphs/Lab2hist/yyn_time_diff.pdf");
    canvas3->Print("graphs/Lab2hist/tp_total_h.pdf");
    canvas4->Print("graphs/Lab2hist/ratio.pdf");

    std::cout << "\033[1;32mLEO_INFO: Files saved!\033[22m If errors appear "
                 "then first create a \"graphs/Lab2hist\" folder in your "
                 "directory.\033[0m"
              << '\n';
    ;
  }
}