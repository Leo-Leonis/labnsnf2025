// #include "RooCategory.h"
// #include "RooDataSet.h"
// #include "RooFitResult.h"
// #include "RooHist.h"
// #include "RooPlot.h"
// #include "RooRealVar.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>

// (LEO) Compares two ints from p1 and p2 and returns the smaller one
// int SmallerOfTwo(int a, int b) {
//   if (a < b)
//     return a;
//   else
//     return b;
// }

void Lab2hist2(int filepath_option = 0, bool do_debug = false) {

  gStyle->SetOptStat(111111);

  // (only in debug mode) decides the number of non triple-fff events being
  // recorded and in all these cases the p_n values are printed in the terminal
  const int debug_detections = 1000;
  // maximum value of each pn count
  const int max_pn = 4095;
  // max difference between coincident counts (events in which there are more
  // than 2 positive stop signals)
  const int max_diff = 10;
  const int n_bins = 150;

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
  double const p1_cal[2] = {10.1525, 0.249765};
  double const p2_cal[2] = {10.7402, 0.249769};
  double const p3_cal[2] = {11.8005, 0.249881};

  // main histogram
  auto *joint_hist =
      new TH1D("Joint histogram", "All TP stops;Stop time (ns);Entries", n_bins,
               0, 16500);
  joint_hist->SetFillColor(kBlue);
  auto *p3_hist = new TH1D("p3 stops", "PL3 stops;Stop time (ns);Entries",
                           n_bins, 0, 16500);
  p3_hist->SetFillColor(kRed);

  int ev, p1, p2, p3;
  double t;
  int n_detections = 0; // number of detections
  int ev_type_count[8] = {0, 0, 0, 0,
                          0, 0, 0, 0}; // all different cases except triple-fff

  /* the cases are: ("FP" = false positive, "TN" true negative)
        0. no no no (TN)
        1. yes no no (FP)
        2. no yes no (true positive)
        3. no no yes (true positive)
        4. yes yes no (true positive)
        5. yes no yes (FP)
        6. no yes yes (FP)
        7. yes yes yes (FP)
                                8. yes yes no, but p1 < p2 (FP)
  */

  // integer that detects any kind of different case other than triple fff
  // bool any_update = false;

  std::ifstream file(filepath_string.c_str());

  while (!file.eof()) {
    // read in the values
    file >> ev >> t >> p1 >> p2 >> p3;

    if (p1 + p2 + p3 == 12285) {
      ev_type_count[0]++;
      continue; // if no stop then skip
    }

    if (p1 != max_pn) { // "pn != max_pn" is "yes", otherwise "no"
      if (p2 != max_pn) {
        if (p3 != max_pn) {
          ev_type_count[7]++; // yes yes yes
        } else {
          if ((p1 >= p2) && ((p1 - p2) < max_diff)) {
            ev_type_count[4]++; // yes yes no (TP)
            joint_hist->Fill((static_cast<double>(p2) + p2_cal[0]) / p2_cal[1]);
          }
        }
      } else {
        if (p3 != max_pn) {
          ev_type_count[5]++; // yes no yes
        } else {
          ev_type_count[1]++; // yes no no
        }
      }
    } else {
      if (p2 != max_pn) {
        if (p3 != max_pn) {
          ev_type_count[6]++; // no yes yes
        } else {
          ev_type_count[2]++; // no yes no (TP)
          joint_hist->Fill((static_cast<double>(p2) + p2_cal[0]) / p2_cal[1]);
        }
      } else {
        if (p3 != max_pn) {
          ev_type_count[3]++; // no no yes (TP)
          joint_hist->Fill((static_cast<double>(p3) + p3_cal[0]) / p3_cal[1]);
          p3_hist->Fill((static_cast<double>(p3) + p3_cal[0]) / p3_cal[1]);
        } else { // no no no (check case if all correct)
          if (p1 + p2 + p3 == max_pn * 3) {
            std::cout << "LEO_ERROR: event classification error: p1, p2 "
                         "and p3 are "
                      << p1 << ", " << p2 << ", " << p3 << " at ev. " << ev
                      << '\n';
          }
        }
      }
    }

#if 0
    if (p1 != max_pn) {
      joint_hist->Fill((static_cast<double>(p1) + p1_cal[0]) / p1_cal[1]);
      any_update = true;
      ev_type_count[0]++;
    } else if (p2 != max_pn) {
      joint_hist->Fill((static_cast<double>(p2) + p2_cal[0]) / p2_cal[1]);
      any_update = true;
      ev_type_count[1]++;
    } else if (p3 != max_pn) {
      joint_hist->Fill((static_cast<double>(p3) + p3_cal[0]) / p3_cal[1]);
      any_update = true;
      ev_type_count[2]++;
    }
#endif

    if (do_debug == true) { // log the info with any update
      std::cout << p1 << '\t' << p2 << '\t' << p3 << '\n';
      n_detections++;
      // std::cout << n_detections << '\n';

      if (n_detections == debug_detections)
        break; // stop detection

      // any_update = 0; // reset any_update
    }
  }

  std::cout << '\t' << "no stops: " << ev_type_count[0] << '\n'
            << '\t' << "p1 stops: " << ev_type_count[1] << '\n'
            << '\t' << "p2 stops: " << ev_type_count[2] << '\n'
            << '\t' << "p3 stops: " << ev_type_count[3] << '\n'
            << '\t' << "yes yes no: " << ev_type_count[4] << '\n';

  TCanvas *canvas = new TCanvas("canvas", "Joint histogram", 1280, 720);
  canvas->cd(1);
  gPad->SetLogy();
  gPad->SetGrid();
  gStyle->SetPalette(1);
  joint_hist->Draw();
  p3_hist->Draw("same");
  // add proper legend for histograms
  gPad->BuildLegend();
}