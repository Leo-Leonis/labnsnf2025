// #include "RooCategory.h"
// #include "RooDataSet.h"
// #include "RooFitResult.h"
// #include "RooHist.h"
// #include "RooPlot.h"
// #include "RooRealVar.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>

void Lab2hist2(int filepath_option = 0, int do_debug = false) {

  // (only in debug mode) decides the number of  non triple-fff events being
  // recorded and in all these cases the p_n values are printed in the terminal
  const int debug_detections = 20;

  const int max_pn_value = 4095;

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

  // param[0] = intercept, param[1] = slope
  double const p1_param[2] = {40.6223, 0.999062};
  double const p2_param[2] = {40.6223, 0.999077};
  double const p3_param[2] = {40.6223, 0.999524};

  // main histogram
  auto *joint_hist = new TH1D(
      "Joint histogram", "Detected events;Stop time (ns);Entries", 50, 0, 5000);
  joint_hist->SetFillColor(kBlue);

  int ev, p1, p2, p3;
  double t;
  int n_detections = 0; // number of detections
  int n_events[3] = {0, 0, 0};
  // integer that detects any kind of different case other than triple fff
  int any_update;

  std::ifstream file(filepath_string.c_str());
  while (!file.eof()) {

    // read in the values
    file >> ev >> t >> p1 >> p2 >> p3;
    // TODO: maybe add error argument
    if (p1 < max_pn_value) {
      joint_hist->Fill((p1 + p1_param[0]) * 4 / p1_param[1]);
      any_update++;
      n_events[0]++;
    } else if (p2 < max_pn_value) {
      joint_hist->Fill((p2 + p2_param[0]) * 4 / p2_param[1]);
      any_update++;
      n_events[1]++;
    } else if (p3 < max_pn_value) {
      joint_hist->Fill((p3 + p3_param[0]) * 4 / p3_param[1]);
      any_update++;
      n_events[2]++;
    }

    if ((do_debug == 1) && (any_update != 0)) { // log the info with any update
      std::cout << p1 << '\t' << p2 << '\t' << p3 << '\n';
      n_detections++;
      // std::cout << n_detections << '\n';

      if (n_detections == debug_detections)
        break; // stop detection

      any_update = 0; // reset any_update
    }
  }

  std::cout << '\t' << "p1 stops: " << n_events[0] << '\n'
            << '\t' << "p2 stops: " << n_events[1] << '\n'
            << '\t' << "p3 stops: " << n_events[2] << '\n';

  TCanvas *canvas = new TCanvas("canvas", "Joint histogram", 1280, 720);
  canvas->cd(1);
  gStyle->SetPalette(1);
  joint_hist->Draw();
  // frame1->GetYaxis()->SetLimits(0, 100);
  // frame1->GetYaxis()->SetTitle("Events/ 3 GeV");
  // frame1->GetXaxis()->SetTitle("m_{#font[12]{4l}} (GeV)");
  // add proper legend for histograms
  gPad->BuildLegend();
}