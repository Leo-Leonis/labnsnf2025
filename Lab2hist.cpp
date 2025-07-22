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

  const int n_bins = 250;
  // (only in debug mode) decides the number of non triple-fff events being
  // recorded and in all these cases the p_n values are printed in the terminal
  const int debug_detections = 1000;
  // maximum value of each pn count
  const int max_pn = 4095;
  // max difference in absolute value between coincident counts (events in which
  // there are more than 2 positive stop signals)
  const int max_diff = 12;

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

  // main histogram
  TH1D *joint_hist =
      new TH1D("Joint histogram", "All TP stops;Stop time (ns);Entries", n_bins,
               0, 16500);
  // stop 3 "PL3" histogram
  TH1D *p3_hist = new TH1D("PL3 stops", "PL3 stops;Stop time (ns);Entries",
                           n_bins, 0, 16500);
  TH1D *yyn_diff_h = new TH1D("yyn_diff_h",
                              "t_{PL1} - t_{PL2} in \"yes yes no\" events;time "
                              "difference (ns/4);Entries",
                              40, -20, 20);

  int ev, p1, p2, p3;
  double t;
  int n_detections = 0;                               // number of detections
  int ev_type_count[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // all different cases

  /* the cases are: ("FP" = false positive, "TN" true negative)
        0. no no no (TN)
        1. yes no no (FP)
        2. no yes no (true positive)
        3. no no yes (true positive)
        4. yes yes no (true positive, p2 <= p1)
        5. yes no yes (FP)
        6. no yes yes (FP)
        7. yes yes yes (FP)
        8. yes yes no, but p1 < p2 (FP)
  */

  // integer that detects any kind of different case other than triple-fff event
  // bool any_update = false;

  std::ifstream file(filepath_string.c_str());

  while (!file.eof()) {
    // read in the values
    file >> ev >> t >> p1 >> p2 >> p3;

    if (p1 + p2 + p3 == 12285) {
      ev_type_count[0]++;
      continue; // if no stop then skip
    }

    n_detections++;

    if (p1 != max_pn) { // "pn != max_pn" is "yes", otherwise "no"
      if (p2 != max_pn) {
        if (p3 != max_pn) {
          ev_type_count[7]++; // yes yes yes
        } else {
          yyn_diff_h->Fill(/* TMath::Abs */ (p1 - p2));
          if (/* (p1 >= p2) && */ (TMath::Abs(p1 - p2) <= max_diff)) {
            ev_type_count[4]++; // yes yes no (TP)
            joint_hist->Fill(static_cast<double>(p2) * p2_cal[1] + p2_cal[0]);
          } else /* if (p1 < p2) */ {
            ev_type_count[8]++; // yes yes no (FP)
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
          joint_hist->Fill(static_cast<double>(p2) * p2_cal[1] + p2_cal[0]);
        }
      } else {
        if (p3 != max_pn) {
          ev_type_count[3]++; // no no yes (TP)
          joint_hist->Fill(static_cast<double>(p3) * p3_cal[1] + p3_cal[0]);
          p3_hist->Fill(static_cast<double>(p3) * p3_cal[1] + p3_cal[0]);
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
      joint_hist->Fill(static_cast<double>(p1) * p1_cal[1] + p1_cal[0]);
      any_update = true;
      ev_type_count[0]++;
    } else if (p2 != max_pn) {
      joint_hist->Fill(static_cast<double>(p2) * p2_cal[1] + p2_cal[0]);
      any_update = true;
      ev_type_count[1]++;
    } else if (p3 != max_pn) {
      joint_hist->Fill(static_cast<double>(p3) * p3_cal[1] + p3_cal[0]);
      any_update = true;
      ev_type_count[2]++;
    }
#endif

    if (do_debug == true) { // log the info with any update
      std::cout << p1 << '\t' << p2 << '\t' << p3 << '\n';
      // std::cout << n_detections << '\n';

      if (n_detections == debug_detections)
        break; // stop detection

      // any_update = 0; // reset any_update
    }
  }

  auto ev_cast = static_cast<double>(ev);
  auto sum_tp_cast = static_cast<double>(ev_type_count[2] + ev_type_count[3] +
                                         ev_type_count[4]);
  auto sum_fp_cast = static_cast<double>(ev_type_count[1] + ev_type_count[5] +
                                         ev_type_count[6] + ev_type_count[7] +
                                         ev_type_count[8]);
  auto n_det_cast = static_cast<double>(n_detections);

  std::cout << "TOTAL EVENTS: " << ev << '\n';

  std::cout << "TRUE NEGATIVE:" << '\n'
            << '\t' << "no no no " << '\t' << ev_type_count[0] << '\t' << "("
            << static_cast<double>(ev_type_count[0]) * 100 / ev_cast << "%)"
            << '\n';

  std::cout << "NON-TRIPLE FFF (NTF) EVENTS: " << n_detections << " ("
            << n_det_cast * 100 / ev_cast << "%)" << '\n';

  std::cout << "TRUE POSITIVES (" << sum_tp_cast * 100 / ev_cast
            << "% of total, " << sum_tp_cast * 100 / n_det_cast
            << "% of NTFs):" << '\n'
            << '\t' << "no yes no " << '\t' << ev_type_count[2] << '\t' << "("
            << static_cast<double>(ev_type_count[2]) * 100 / ev_cast
            << "% of total, "
            << static_cast<double>(ev_type_count[2]) * 100 / n_det_cast
            << "% of NTFs)" << '\n'
            << '\t' << "no no yes " << '\t' << ev_type_count[3] << '\t' << "("
            << static_cast<double>(ev_type_count[3]) * 100 / ev_cast
            << "% of total, "
            << static_cast<double>(ev_type_count[3]) * 100 / n_det_cast
            << "% of NTFs)" << '\n'
            << '\t' << "yes yes no " << '\t' << ev_type_count[4] << '\t' << "("
            << static_cast<double>(ev_type_count[4]) * 100 / ev_cast
            << "% of total, "
            << static_cast<double>(ev_type_count[4]) * 100 / n_det_cast
            << "% of NTFs)" << '\n';

  std::cout
      << "FALSE POSITIVES (" << sum_fp_cast * 100 / ev_cast << "% of total, "
      << sum_fp_cast * 100 / n_det_cast << "% of NTFs):" << '\n'
      << '\t' << "yes no no " << '\t' << ev_type_count[1] << '\t' << "("
      << static_cast<double>(ev_type_count[1]) * 100 / ev_cast << "% of total, "
      << static_cast<double>(ev_type_count[1]) * 100 / n_det_cast
      << "% of NTFs)" << '\n'
      << '\t' << "yes no yes " << '\t' << ev_type_count[5] << '\t' << "("
      << static_cast<double>(ev_type_count[5]) * 100 / ev_cast << "% of total, "
      << static_cast<double>(ev_type_count[5]) * 100 / n_det_cast
      << "% of NTFs)" << '\n'
      << '\t' << "no yes yes " << '\t' << ev_type_count[6] << '\t' << "("
      << static_cast<double>(ev_type_count[6]) * 100 / ev_cast << "% of total, "
      << static_cast<double>(ev_type_count[6]) * 100 / n_det_cast
      << "% of NTFs)" << '\n'
      << '\t' << "yes yes yes " << '\t' << ev_type_count[7] << '\t' << "("
      << static_cast<double>(ev_type_count[7]) * 100 / ev_cast << "% of total, "
      << static_cast<double>(ev_type_count[7]) * 100 / n_det_cast
      << "% of NTFs)" << '\n'
      << '\t' << "yes yes no" << '\t' << ev_type_count[8] << '\t' << "("
      << static_cast<double>(ev_type_count[8]) * 100 / ev_cast << "% of total, "
      << static_cast<double>(ev_type_count[8]) * 100 / n_det_cast
      << "% of NTFs, with |time_diff|>" << max_diff * 4 << "ns)" << '\n';

  // joint_hist->Scale(
  //     1. / (ev_type_count[2] + ev_type_count[3] + ev_type_count[4]), "width");
  // p3_hist->Scale(1. / ev_type_count[3], "width");
  joint_hist->SetFillColor(kBlue);
  joint_hist->SetLineColor(kBlue);
  joint_hist->SetMarkerColor(kBlue);
  joint_hist->SetMarkerStyle(20);
  p3_hist->SetFillColor(kRed);
  p3_hist->SetLineColor(kRed);
  p3_hist->SetMarkerStyle(22);
  p3_hist->SetMarkerColor(kRed);

  TH1D *div_h = new TH1D(*joint_hist);
  div_h->Divide(p3_hist);
  div_h->SetLineColor(kBlack);
  div_h->SetMarkerColor(kBlack);
  div_h->SetLineColor(kBlack);
  div_h->SetMarkerStyle(21);

  TCanvas *canvas1 = new TCanvas("canvas1", "Joint histogram", 1280, 720);
  canvas1->cd(1);
  gStyle->SetPalette(1);

  joint_hist->Draw();
  p3_hist->Draw("same");
  // add proper legend for histograms
  gPad->BuildLegend();
  gPad->SetLogy();
  gPad->SetGrid();

  TCanvas *canvas2 =
      new TCanvas("canvas2", "yes yes no time difference", 1280, 720);
  yyn_diff_h->Draw();
  // canvas2->SetLogy(0);

  gPad->SetLogy();
  gPad->SetGrid();
}