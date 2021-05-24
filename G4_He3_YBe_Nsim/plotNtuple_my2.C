// ROOT macro file for plotting example He3 ntuple
// 
// Can be run from ROOT session:
// root[0] .x plotNtuple.C

void plotNtuple_my2 () {
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile *f = TFile::Open("He3.root");

  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y1 = 0;
  TNtuple* ntuple = (TNtuple*)f->Get("He3");
  ntuple->Draw("Eabs");
  TTree* tree = (TTree*)f->Get("He3");
  uint nentries = tree->GetEntries();

  Double_t Eabs;
  tree->SetBranchAddress("Eabs",&Eabs);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs > 0)
    { cout<<"Center Tube: "<<Eabs<<endl;
    y1++;}
  }


  TCanvas* c2 = new TCanvas("c2", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y2 = 0;
  ntuple->Draw("Eabs1");
  Double_t Eabs1;
  tree->SetBranchAddress("Eabs1",&Eabs1);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs1 > 0)
    { cout<<"Side Tube1: "<<Eabs1<<endl;
    y2++;}
  }


  TCanvas* c3 = new TCanvas("c3", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Double_t x;
  Int_t y3 = 0;
  ntuple->Draw("Eabs2");
  Double_t Eabs2;
  tree->SetBranchAddress("Eabs2",&Eabs2);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs2 > 0)
    { cout<<"Side Tube2: "<<Eabs2<<endl;
    y3++;}
  }
    cout << "Center Tube: " << y1 << endl;
    cout << "Side Tube1: " << y2 << endl;
    cout << "Side Tube2: " << y3 << endl; 
    cout << "Total = "<< y1 + y2 + y3 << endl;
  }
