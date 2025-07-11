#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include <fstream>

void Lab2hist2(int filepath_option = 0) {

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
        << "LEO_ERROR: invalid file_id. Please execute the file by typing "
           "\"Lab2hist(<file_id>)\" with file_id = 1 for Makar, 2 for Leo."
        << '\n';
    return;
    break;
  }

  // main histogram
  auto *hist = new TH1D("Data", "Detected events", 50, 0, 5000);
  double ev, t, p1, p2, p3;
  std::ifstream file(filepath_string.c_str());
  int line_n = 0;

  while (!file.eof()) {
    file >> ev >> t >> p1 >> p2 >> p3;
    // maybe add error argument
    if (p1 < 4095) {
      hist->Fill((p1 + 40.6223) * 4 / 0.999062);
    } else if (p2 < 4095) {
      hist->Fill((p2 + 42.9588) * 4 / 0.999077);
    } else if (p3 < 4095) {
      hist->Fill((p3 + 47.2021) * 4 / 0.999524);
    }
  }

  TCanvas *can = new TCanvas("canvas", "Joint histogram", 1200, 1200);
  can->cd(1);
  gStyle->SetPalette(1);
  hist->Draw();
  // frame1->GetYaxis()->SetLimits(0, 100);
  // frame1->GetYaxis()->SetTitle("Events/ 3 GeV");
  // frame1->GetXaxis()->SetTitle("m_{#font[12]{4l}} (GeV)");
  // add proper legend for histograms
  gPad->BuildLegend();
}