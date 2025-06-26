// This file reads in the calibration files and makes histograms, providing the
// mean and the RMS of each of the 3 histograms

#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

int analyse(int file_id = 0) {

  gStyle->SetOptStat(111111);

  TString a; // file id case
  int n_events;
  int min_bin_value;
  int max_bin_value;

  switch (file_id) {
  case 1:
    a = '1';
    n_events = 1359;
    min_bin_value = 124;
    max_bin_value = 168;
    break;
  case 2:
    a = '2';
    n_events = 1004;
    min_bin_value = 2128;
    max_bin_value = 2172;
    break;
  case 3:
    a = '3';
    n_events = 1004;
    min_bin_value = 4124;
    max_bin_value = 4168;
    break;
  case 4:
    a = '4';
    n_events = 1066;
    min_bin_value = 6124;
    max_bin_value = 6168;
    break;
  case 5:
    a = '5';
    n_events = 1005;
    min_bin_value = 8124;
    max_bin_value = 8168;
    break;
  case 6:
    a = '6';
    n_events = 1006;
    min_bin_value = 10116;
    max_bin_value = 10160;
    break;

  default:
    break;
  }

  if (file_id < 1 || file_id > 6) {
    std::cout << "LEO_ERROR: invalid file_id (input= " << file_id
              << "), please select from 1 to 6." << '\n';
    return 0;
  }

  std::ifstream ifile("cali_dec/5b_cali_" + a + "_conv.txt");

  int n_bin = 11;

  TH1D *h1 = new TH1D("h1", "FPGA calibration; time (ns); Entries", n_bin,
                      min_bin_value, max_bin_value);
  TH1D *h2 = new TH1D("h2", "FPGAz calibration; time (ns); Entries", n_bin,
                      min_bin_value, max_bin_value);
  TH1D *h3 = new TH1D("h3", "FPGA calibration; time (ns); Entries", n_bin,
                      min_bin_value, max_bin_value);

  std::string line_string;
  int value; // value to be used to fill the histogram
  float seconds;

  for (int line_n = 0; line_n < n_events; line_n++) {
    /* if (!getline(ifile, line_string)) {
      std::cout << line_string << '\n';
      std::cout << "LEO_ERROR: couldn't read line " << line_n
      << " of file 5b_cali_" << a << "_conv.txt" << '\n';
      return 0;
      } */
    ifile >> value; // skip event number
    // std::cout << value << '\t';
    ifile >> seconds; // skip seconds
    // std::cout << seconds << '\t';
    ifile >> value; // acquire first value
    // std::cout << value * 4 << '\t';
    h1->Fill(value * 4);
    ifile >> value; // acquire second value
    // std::cout << value * 4 << '\t';
    h2->Fill(value * 4);
    // std::cout << value * 4 << '\n';
    ifile >> value; // acquire third value
    h3->Fill(value * 4);
  }

  TCanvas *canvas1 = new TCanvas("canvas1", "FPGA 1", 720, 720);
  h1->Draw("h");
  TCanvas *canvas2 = new TCanvas("canvas2", "FPGA 2", 720, 720);
  h2->Draw("h");
  TCanvas *canvas3 = new TCanvas("canvas3", "FPGA 3", 720, 720);
  h3->Draw("h");

  std::cout << '\n' << h1->GetMean() << '\t' << h1->GetRMS() << '\n';
  std::cout << '\n' << h2->GetMean() << '\t' << h2->GetRMS() << '\n';
  std::cout << '\n' << h3->GetMean() << '\t' << h3->GetRMS() << '\n';

  return 0;
}