 #include "RooRealVar.h"
 #include <fstream>
 #include "RooDataSet.h"
 #include "TCanvas.h"
 #include "RooPlot.h"
 #include "TAxis.h"
 #include "RooHist.h"
 #include "TCanvas.h"
 #include "TLegend.h"
 #include "TStyle.h"
 #include "RooFitResult.h"
 #include "TH1D.h"
 #include "RooCategory.h"

using namespace RooFit;
using namespace RooStats;
void Lab2hist2(){

auto *data = new TH1D("Data", "Detected events", 50, 0, 5000);
double ev,t,p1,p2,p3;
std::ifstream file1("5b_data_conv.txt");
while (!file1.eof()) {
file1 >> ev >> t >> p1 >> p2 >> p3;
//maybe add error argument
if(p1<4095){
data->Fill((p1+40.6223)/0.999062);
}
else if(p2<4095){
data->Fill((p2+42.9588)/0.999077); 
}
else if(p3<4095){
data->Fill((p3+47.2021)/0.999524);  
}
}

TCanvas *can = new TCanvas("can", "Joint histogram", 1200, 1200);
can->cd(1);
gStyle->SetPalette(1);
data->Draw();
//frame1->GetYaxis()->SetLimits(0, 100);
//frame1->GetYaxis()->SetTitle("Events/ 3 GeV");
//frame1->GetXaxis()->SetTitle("m_{#font[12]{4l}} (GeV)");
//add proper legend for histograms
gPad->BuildLegend();


}