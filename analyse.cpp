// This file reads in the calibration files and makes histograms, providing the
// mean and the RMS of each of the 3 histograms

// #include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"
// #include "TStyle.h"
// #include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int analyse(int do_debug = 0 /* int file_id = -1 */) {

  // gStyle->SetOptStat(111111);

  int const n_bin = 11;
  std::string const hist_title_str = "FPGA calibration; FPGA count; Entries";

  TString const file_id_char_v[6] = {'1', '2', '3',
                                     '4', '5', '6'}; // file id case
  int const n_events_v[6] = {1359, 1004, 1004, 1066, 1005, 1006};
  int const min_bin_value_v[6] = {31, 532, 1031, 1531, 2031, 2529};
  int const max_bin_value_v[6] = {42, 543, 1042, 1542, 2042, 2540};
  int const x_value[6] = {200, 2200, 4200, 6200, 8200, 10200};

  std::ifstream ifile;

  std::array<std::array<TH1F *, 3>, 6> h_master;

  for (size_t i = 0; i != 3; i++) {
    for (size_t j = 0; j != 6; j++) {

      // std::cout << i << " " << j << '\n'; // DEBUG
      h_master[j][i] =
          new TH1F("FPGA" + file_id_char_v[i] + "." + file_id_char_v[j],
                   hist_title_str.c_str(), n_bin, min_bin_value_v[j],
                   max_bin_value_v[j]);
    }
  }

  // std::cout << "debug lol!" << '\n';

  int value; // value to be used to fill the histogram
  float seconds;

  for (size_t j = 0; j != 6; j++) {
    // std::cout << "debug lol 1" << '\n';
    TString filename =
        "data/cali_dec/5b_cali_" + file_id_char_v[j] + "_conv.txt";
    // ifile.open("data/cali_dec/5b_cali_" + file_id_char_v[j] + "_conv.txt");
    ifile.open(filename);
    // std::cout << "debug lol 2" << '\n';
    if (ifile.is_open()) {
      for (int line_n = 0; line_n < n_events_v[j]; line_n++) {
        // std::cout << "debug lol!" << '\n';
        ifile >> value >> seconds; // skip event number and seconds
        // std::cout << value << '\t';
        // std::cout << seconds << '\t';
        ifile >> value; // acquire first value
        // std::cout << value << '\t';
        h_master[j][0]->Fill(value);
        ifile >> value; // acquire second value
        // std::cout << value << '\t';
        h_master[j][1]->Fill(value);
        ifile >> value; // acquire third value
        // std::cout << value << '\n';
        h_master[j][2]->Fill(value);
      }
    } else {
      std::cout << "LEO_ERROR: file reading unsuccessful, file "
                   "data/cali_dec/5b_cali_" +
                       file_id_char_v[j] + "_conv.txt not found"
                << '\n';
    }
    ifile.close();
  }

  for (std::array<TH1F *, 3> h_v : h_master) {
    for (TH1F *hist : h_v) {
      std::cout << std::setprecision(10) << hist->GetMean() << '\t'
                << hist->GetRMS() << '\t';
    }
    std::cout << '\n';
  }

  std::ofstream ofile;
  for (size_t i = 0; i != 3; i++) {
    ofile.open("data/cali_tot/cali_tot_" + file_id_char_v[i] + ".txt",
               std::ofstream::out | std::ofstream::trunc);
    if (ofile.is_open()) {
      for (size_t j = 0; j != 6; j++) {
        ofile << x_value[j] << '\t' << std::setprecision(10)
              << h_master[j][i]->GetMean() << '\t' << h_master[j][i]->GetRMS()
              << '\n';
      }
    } else {
      std::cout << "LEO_ERROR: could NOT open file \"data/cali_tot/cali_tot_" +
                       file_id_char_v[i] + ".txt\" "
                << '\n';
    }
    ofile.close();
  }

  return 0;
}