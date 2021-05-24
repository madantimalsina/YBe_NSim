// ROOT macro file for plotting example He3 ntuple
// 
// Can be run from ROOT session:
// root[0] .x plotNtuple_12MeV.C

void plotNtuple_5MeV () {
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile *f = TFile::Open("He3_2p5inch_5MeV.root");

  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y1 = 0;
  Int_t z1 = 0;
  TNtuple* ntuple = (TNtuple*)f->Get("He3");
  ntuple->Draw("Eabs");
  TTree* tree = (TTree*)f->Get("He3");
  uint nentries = tree->GetEntries();
 
  // Central Tube
  Double_t Eabs;
  tree->SetBranchAddress("Eabs",&Eabs);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs > 0)
    { cout<<"Center Tube: "<<Eabs<<endl;
    y1++;}

    if(Eabs > 0 && Eabs < 0.1)
    { cout<<"nHits (0 < Edep < 0.1 MeV) .\n Center Tube: "<< Eabs <<endl;
    z1++;}
  }
  // Side Tube 1 (Left)
  TCanvas* c2 = new TCanvas("c2", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y2 = 0;
  Int_t z2 = 0;
  ntuple->Draw("Eabs1");
  Double_t Eabs1;
  tree->SetBranchAddress("Eabs1",&Eabs1);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs1 > 0)
    { cout<<"Side Tube1: "<<Eabs1<<endl;
    y2++;}

    if(Eabs1 > 0 && Eabs1 < 0.1)
    { cout<<"nHits (0 < Edep < 0.1 MeV) .\n Center Tube: "<< Eabs1 <<endl;
    z2++;}
  }

  // side Tube 2 (Right)
  TCanvas* c3 = new TCanvas("c3", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Double_t x;
  Int_t y3 = 0;
  Int_t z3 = 0;
  ntuple->Draw("Eabs2");
  Double_t Eabs2;
  tree->SetBranchAddress("Eabs2",&Eabs2);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs2 > 0)
    { cout<<"Side Tube2: "<<Eabs2<<endl;
    y3++;}

    if(Eabs2 > 0 && Eabs2 < 0.1)
    { cout<<"nHits (0 < Edep < 0.1 MeV) .\n Center Tube: "<< Eabs2 <<endl;
    z3++;}
  } 


  // Inside He3 Tube
  // Central Tube
  TCanvas* c1_T = new TCanvas("c1_T", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y1_T = 0;
  Int_t z1_T = 0;
  ntuple->Draw("Eabs_inT");
  Double_t Eabs_inT;
  tree->SetBranchAddress("Eabs_inT",&Eabs_inT);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs_inT > 0)
    { cout<<"Center Tube (Inside He3 Radius): "<<Eabs_inT<<endl;
    y1_T++;}

    if(Eabs_inT > 0 && Eabs_inT < 0.1)
    { cout<<"nHits (0 < Edep < 0.1 MeV) .\n Center Tube (Inside He3 Tube): "<< Eabs_inT <<endl;
    z1_T++;}
  }

  // Side Tube 1 (Left)
  TCanvas* c2_T = new TCanvas("c2_T", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y2_T = 0;
  Int_t z2_T = 0;
  ntuple->Draw("Eabs1_inT");
  Double_t Eabs1_inT;
  tree->SetBranchAddress("Eabs1_inT",&Eabs1_inT);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs1_inT > 0)
    { cout<<"Side Tube 1 (Inside He3 Radius): "<<Eabs1_inT<<endl;
    y2_T++;}

    if(Eabs1_inT > 0 && Eabs1_inT < 0.1)
    { cout<<"nHits (0 < Edep < 0.1 MeV) .\n Side Tube 1 (Inside He3 Tube): "<< Eabs1_inT <<endl;
    z2_T++;}
  }



  // Side Tube 2 (Right)
  TCanvas* c3_T = new TCanvas("c3_T", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y3_T = 0;
  Int_t z3_T = 0;
  ntuple->Draw("Eabs2_inT");
  Double_t Eabs2_inT;
  tree->SetBranchAddress("Eabs2_inT",&Eabs2_inT);
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs2_inT > 0)
    { cout<<"Side Tube 2 (Inside He3 Radius): "<<Eabs2_inT<<endl;
    y3_T++;}

    if(Eabs2_inT > 0 && Eabs2_inT < 0.1)
    { cout<<"nHits (0 < Edep < 0.1 MeV) .\n Side Tube 2 (Inside He3 Tube): "<< Eabs2_inT <<endl;
    z3_T++;}
  }

    cout << endl;
    cout << "Center Tube: " << y1 << endl;
    cout << "Side Tube1: " << y2 << endl;
    cout << "Side Tube2: " << y3 << endl; 
    cout << "Total = "<< y1 + y2 + y3 << endl;

    cout << "\n" << endl;
    cout << "nHits (0 < Edep < 0.1 MeV)" << endl;
    cout << "Center Tube: " << z1 << endl;
    cout << "Side Tube1: " << z2 << endl;
    cout << "Side Tube2: " << z3 << endl; 
    cout << "Total = "<< z1 + z2 + z3 << endl;

    cout << "\n" << endl;
    cout << "NOW INSIDE HE3 TUBES:::::.\n" << endl;
    cout << "Center Tube (Inside): " << y1_T << endl;
    cout << "Side Tube1 (Inside): " << y2_T << endl;
    cout << "Side Tube2 (Inside): " << y3_T << endl; 
    cout << "Total (Inside) = "<< y1_T + y2_T + y3_T << endl;

    cout << "\n" << endl;
    cout << "nHits (0 < Edep < 0.1 MeV) Inside" << endl;
    cout << "Center Tube (Inside): " << z1_T << endl;
    cout << "Side Tube1 (Inside): " << z2_T << endl;
    cout << "Side Tube2 (Inside): " << z3_T << endl; 
    cout << "Total (Inside) = "<< z1_T + z2_T + z3_T << endl;

  }
