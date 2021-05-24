// ROOT macro file for plotting example He3 ntuple
// 
// Can be run from ROOT session:
// root[0] .x plotNtuple_12MeV.C

void plotNtuple_YBe() {
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile *f = TFile::Open("He3_1p5inch_YBe_1m_Tungsten_SAir.root");
  // Get ntuple
  TNtuple* ntuple = (TNtuple*)f->Get("He3");
  TTree* tree = (TTree*)f->Get("He3");
  uint nentries = tree->GetEntries();

  // Inside He3 Tube
  // Central Tube
  TCanvas* c1_T = new TCanvas("c1_T", "", 20, 20, 800, 800);
  gPad->SetLogy(1);
  Int_t y1_T = 0;
  Int_t z1_T = 0;
  Double_t Eabs_inT;
  tree->SetBranchAddress("Eabs_inT",&Eabs_inT);
  ntuple->Draw("Eabs_inT");
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    if(Eabs_inT > 0)
    { cout<<"Center Tube: "<<Eabs_inT<<endl;
    y1_T++;}

    if(Eabs_inT >= 0.15)
    { cout<<"Center Tube (Edep >= 0.15): "<< Eabs_inT <<endl;
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
    { cout<<"Left:Side Tube 1 : "<<Eabs1_inT<<endl;
    y2_T++;}

    if(Eabs1_inT >= 0.15)
    { cout<<"Left:Side Tube 1 (Edep >= 0.15): "<< Eabs1_inT <<endl;
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
    { cout<<"Rigth:Side Tube 2: "<<Eabs2_inT<<endl;
    y3_T++;}

    if(Eabs2_inT >= 0.15)
    { cout<<"Right:Side Tube 2 (Edep >= 0.15): "<< Eabs2_inT <<endl;
    z3_T++;}
  }

    cout << "\n" << endl;
    cout << "NOW INSIDE HE3 TUBES:::::.\n" << endl;
    cout << "Center Tube: " << y1_T << endl;
    cout << "Left:Side Tube1: " << y2_T << endl;
    cout << "Right:Side Tube2: " << y3_T << endl; 
    cout << "Total = "<< y1_T + y2_T + y3_T << endl;

    cout << "\n" << endl;
    cout << "nHits (Edep >= 0.15)" << endl;
    cout << "Center Tube: " << z1_T << endl;
    cout << "Left:Side Tube1: " << z2_T << endl;
    cout << "Right:Side Tube2: " << z3_T << endl; 
    cout << "Total = "<< z1_T + z2_T + z3_T << endl;
    cout << "Efficiency = " << (z1_T + z2_T + z3_T)/1000000.0 << endl;

  }
