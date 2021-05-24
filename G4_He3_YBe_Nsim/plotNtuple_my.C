// ROOT macro file for plotting example He3 ntuple
// 
// Can be run from ROOT session:
// root[0] .x plotNtuple.C

void plotNtuple_my(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile f("He3_2inch_1MeV.root");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  c1->Divide(2,2);
  
  // Get ntuple
  TNtuple* ntuple = (TNtuple*)f.Get("He3");

  // Draw Eabs histogram in the pad 1
  c1->cd(1);
  gPad->SetLogy(1);
  ntuple->Draw("x");
  
  // Draw Labs histogram in the pad 2
  c1->cd(2);
  gPad->SetLogy(1);
  ntuple->Draw("y");
  
  // Draw Egap histogram in the pad 3
  // with logaritmic scale for y  ?? how to do this?
  c1->cd(3);
  gPad->SetLogy(1);
  ntuple->Draw("z");
  
  // Draw Lgap histogram in the pad 4
  // with logaritmic scale for y  ?? how to do this?
  c1->cd(4);
  //gPad->SetLogy(1);
  ntuple->Draw("x:y");
}  
