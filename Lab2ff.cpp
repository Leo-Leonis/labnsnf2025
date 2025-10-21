#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRatioPlot.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "RooRealVar.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace RooFit;

double scd(const int a) { return static_cast<double>(a); }

void Lab2ff(){
  std::string filepath_string="/labnsnf2025/5b_data_conv.txt";
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
  double const p1_cal[2] = {40.6139, 4.003800};
  double const p2_cal[2] = {43.0851, 4.003570};
  double const p3_cal[2] = {47.2083, 4.001930};
  double const R=1.21;

  //Graphing/fitting variables
  RooRealVar tau_0("tau_0", "mean of gaussians",  5279, 5220, 5320);
  RooRealVar tau_mu_minus("tau_mu_minus", "mean of gaussians",  5279, 5220, 5320);
  RooRealVar b("b", "Background",  5279, 5220, 5320);
  std::array<TString, 3> ev_str;
  ev_str[0] = "\"no yes no\"";
  ev_str[1] = "\"no no yes\"";
  ev_str[2] = "\"yes yes no\"";

  std::vector<TH1D *> all_h;
  all_h.reserve(3);
  for (size_t i = 0; i != 3; i++) {
    all_h.push_back(new TH1D(ev_str[i] + " evs",
                             ev_str[i] + " evs;Stop time (ns);Entries", n_bins,
                             min_x, max_x));
  }


  // histogram that contains all true positive events
  TH1D *tp_h = new TH1D(
      "tp_h", "all true positive time distribution;Stop time (ns);Entries",
      n_bins, min_x, max_x);

  // distribution of time difference in "yes yes no" events, not used yet to rule out entries
  TH1D *yyn_diff_h = new TH1D("yyn_diff_h",
                              "t_{PL1} - t_{PL2} in \"yes yes no\" events;time "
                              "difference (ns/4);Entries",
                              40, -20, 20); // bin = 40 is maximum around

  int ev_n, p1, p2, p3;
  double t;
  int n_detections = 0; // number of detections
  //       2. no yes no (true positive)
  //       3. no no yes (true positive)
  //       4. yes yes no (true positive, |p1-p2| < max_diff)
  int tp_count = 0;
  int pn_count[3] = {0, 0, 0}; // independent pn counts

  // integer that detects any kind of different case other than triple-fff event
  // bool any_update = false;

  std::ifstream file(filepath_string.c_str());

  double value; // placeholder value

   while (!file.eof()) {
    // read in the values
    file >> ev_n >> t >> p1 >> p2 >> p3;

    //Deleted all conditions not for (even conditionally) TP, as per list
    if (p1 + p2 + p3 == 12285) { // 0. no no no no no wait wait wait wait
      continue; // if no stop then skip
    }

    n_detections++;
    if (p1 == max_pn) {
      // no xxx xxx
      if ((p2 == max_pn)&&(p3 != max_pn)) {
        // no yes no
          tp_count++; 
          value = scd(p2) * p2_cal[1] + p2_cal[0];

          all_h[0]->Fill(scd(p2) * p2_cal[1] + p2_cal[0]);
          tp_h->Fill(value);
      } else if((p2 == max_pn)&&(p3 != max_pn)){
        // no no yes
          pn_count[1]++;
          tp_count++; 
          value = scd(p3) * p3_cal[1] + p3_cal[0];

          tp_h->Fill(value);
          all_h[3]->Fill(scd(p3) * p3_cal[1] + p3_cal[0]);
      }
        
      }  else if((p2!=max_pn)&&(p3==max_pn)){       
        // yes yes no
        pn_count[2]++;

          yyn_diff_h->Fill(p1 - p2);

          if (p1 <= p2)
            value = scd(p1) * p1_cal[1] + p1_cal[0];
          else
            value = scd(p2) * p2_cal[1] + p2_cal[0];
      
          if (TMath::Abs(p1 - p2) <= max_diff) {
            tp_count++; // 4. yes yes no (TP)
            tp_h->Fill(value);
            all_h[2]->Fill(value);
          }


    }
  }


  for (TH1D *hist : all_h) {
    tp_h->Add(hist);
    hist->SetMarkerStyle(kFullCircle);
  }
  TCanvas *canvas1 = new TCanvas("canvas1", "Final Fitting", 1280, 720);
  gStyle->SetPalette(kRainBow); // "kRainBow" is not colourblind friendly!
  gPad->SetGrid();
  gPad->SetLogy();
  tp_h->Draw("Final Fantasy");
  gPad->BuildLegend();
}