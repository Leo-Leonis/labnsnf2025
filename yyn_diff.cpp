
#include "RooBreitWigner.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"

#include <string>

void yyn_diff() {

  int const nbins = 100; // number of bins
  std::string const filename = "result.root";

  RooWorkspace w("spaziodilavoro");

  TFile *b0_f = new TFile(filename.c_str(), "READ");
  if (b0_f->IsZombie()) {
    std::cout << "LEO_ERROR: could not open file \"" << filename << "\"."
              << '\n';
    return;
  }

  w.factory("BreitWigner::bw(m[-20.5,20.5], bw_m[-1,1], bw_w[0.01,4])");
  w.factory("Gaussian::g(m, g_m[-1,1], g_s[0.01,2])");
  auto m = w.var("m"); // yyn variable
  m->SetTitle("time difference (FPGA counts)");
  auto bw = w.pdf("bw");
  bw->SetTitle("Breit-Wigner distribution");
  auto g = w.pdf("g");
  g->SetTitle("Guassian distribution");

  // TH1F *yyn_diff_h = static_cast<TH1F *>(b0_f.Get("massaB0"));
  TH1F *yyn_diff_h = (TH1F *)b0_f->Get("yyn_diff_h");
  if (!yyn_diff_h) {
    std::cout << "LEO_ERROR: could not find histogram \"yyn_diff_h\" in file"
              << filename << "." << '\n';
    b0_f->ls();
    return;
  }
  m->setBins(yyn_diff_h->GetNbinsX());

  RooDataHist ds("yyn_diff_h", "t_{PL1} - t_{PL2} in \"yes yes no\" events", *m,
                 RooFit::Import(*yyn_diff_h));

  bw->fitTo(ds);
  g->fitTo(ds);

  TCanvas *c1 =
      new TCanvas("canvas1", "b0 invariant mass distribution", 720, 720);
  auto mframe =
      m->frame(RooFit::Title("t_{PL1} - t_{PL2} in \"yes yes no\" events"));
  {
    using namespace RooFit;
    ds.plotOn(mframe, Name("b0_graph"), MarkerStyle(kFullCircle));
    bw->plotOn(mframe, Name("bw_graph"));
    g->plotOn(mframe, LineColor(kRed), Name("g_graph"));
  }
  mframe->Draw();

  TLegend *legend = new TLegend(.7, .7, .9, .9);
  legend->AddEntry("b0_graph", "data", "pe");
  legend->AddEntry("bw_graph", "Breit-Wigner distribution", "l");
  legend->AddEntry("g_graph", "Gaussian distribution", "l");
  legend->Draw("same");

  c1->Print("graphs/yyn_diff.pdf");

  return;
}