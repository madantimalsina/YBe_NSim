// ROOT macro file for plotting example He3 ntuple
// 
// Can be run from ROOT session:
// root[0] .x plotNtuple.C

{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile f("He3_10MeV.root");
/*
  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  c1->Divide(2,2);
  
  // Get ntuple
  TNtuple* ntuple = (TNtuple*)f.Get("He3");

  // Draw Eabs histogram in the pad 1
  c1->cd(1);
  ntuple->Draw("Eabs");
  
  // Draw Labs histogram in the pad 2
  c1->cd(2);
  ntuple->Draw("ETube");
  
  // Draw Egap histogram in the pad 3
  // with logaritmic scale for y  ?? how to do this?
  c1->cd(3);
  gPad->SetLogy(1);
  ntuple->Draw("Eabs");
  
  // Draw Lgap histogram in the pad 4
  // with logaritmic scale for y  ?? how to do this?
  c1->cd(4);
  gPad->SetLogy(1);
  ntuple->Draw("ETube");
  c1->cd(4);*/

  TCanvas* c2 = new TCanvas("c2", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Double_t x;
  Int_t y = 0;
  Int_t z = 0;
  TNtuple* ntuple = (TNtuple*)f.Get("He3");
  ntuple->Draw("ETube");
  TH1F *myh1 = (TH1F *) htemp->Clone();
  myh1->SetName("ETube");
  //myh1->SetXTitle("Energy [MeV]");
  //myh1->GetXaxis()->SetRangeUser(0, 12);
  myh1->Draw();
  //x = myh1->GetBinContent(1);
 /* 
  for (int i = 30000; i > 0; i--) {
      x = myh1->GetBinContent(i+1);
      cout << x << endl;

    }
   */ for (int i = 1; i < 30000; ++i) {
      x = myh1->GetBinContent(i+1);
      if (x > 0) {cout << x << endl;
                  y++;
                  z += x;}
    }
    
    cout << "Total = "<< y << endl;
    cout << "NEW Total = "<< z << endl;
}  
