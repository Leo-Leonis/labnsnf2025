// This file reads-in a data file with a name given below (to be modified), and
// converts the numbers in the 3rd, 4th and 5th column (considering "_" as the
// first column) from hex base to decimal base and outputs everything in the
// output file (with a name given below), without the "_" column.
// This file does not require ROOT.

// TO COMPILE: (assumes that g++ is installed, idk about other c++ compilers...)
// g++ -Wall -Wextra convertHex2Dec.cpp -o convertHex2Dec

// TO EXECUTE:
// ./convertHex2Dec

#include <fstream>
#include <iomanip> // for std::fixed and std::setprecision
#include <iostream>
#include <string>

int main() {
  // MODIFY THESE NAMES BELOW
  std::string i_filename("main_data_hex/5b_data.txt");      // input file name
  std::string o_filename("main_data_dec/5b_data_conv.txt"); // output file name
  std::ifstream ifile(i_filename.c_str());                  // input file
  std::ofstream ofile;                                      // output file
  ofile.open(o_filename, std::ofstream::out |
                             std::ofstream::trunc); // open and clear the file

  if (!ifile.good()) { // verify that the file is opened successfully
    std::cout << "LEO_ERROR: couldn't open file " << i_filename << ".";
  }

  std::string underscore_txt; // underscore character (as a control)
  int line_n = 0;             // line number
  int col_n = 0;              // column number
  int i_event;                // number of the current event
  double seconds;             // total seconds
  int hex_value;              // value of hex number in decimal base

  // double slope[3] = {0.999062, 0.999077, 0.999524};
  // double intercept[3] = {-40.6223, -42.9588, -47.2021};

  while (1) {

    line_n++;                    // add line count
    col_n = 1;                   // restart column count
    ifile >> underscore_txt;     // read string and put it in underscore_txt
    if (underscore_txt != "_") { // if no "_" was read then stop the program
      std::cout << "LEO_ERROR: no \"_\" on line " << line_n << " and column "
                << col_n << '\n';
      return 0;
    }
    col_n++; // add column count

    ifile >> std::dec >> i_event; // read in decimal base
    // std::cout  << "i_event = " << i_event << "line_n = " << line_n << '\n';
    if (i_event != line_n) { // i_event must be equal to line_n
      std::cout << "LEO_ERROR: line_n and i_event do not match (at line "
                << line_n << " and column " << col_n << ", line_n: " << line_n
                << ", i_event: " << i_event << ")" << '\n';
      return 0; // ^^^AT THE END OF FILE THIS MESSAGE WILL BE OUTPUTTED^^^
    }
    col_n++;
    ofile << i_event << '\t';
    if (i_event % 1000 == 0) // print out result for every 1000th event
      std::cout << i_event << '\t';

    ifile >> seconds;
    col_n++;
    ofile << std::fixed << std::setprecision(3) << seconds; // 3 decimal digits
    if (i_event % 1000 == 0)
      std::cout << std::fixed << std::setprecision(3) << seconds;

    for (int i = 0; i != 3; i++) {
      ifile >> std::hex >> hex_value; // "std::hex" to read values in hex base
      col_n++;
      ofile << '\t' << hex_value; // outputs in decimal base (automatically)
      if (i_event % 1000 == 0)
        std::cout << '\t' << hex_value;
    }

    ofile << '\n';
    if (i_event % 1000 == 0) {
      std::cout << '\n';
    }

    /* if (line_n == 20)
      break; // debug test */
  }

  ofile.close();
  return 0;
}