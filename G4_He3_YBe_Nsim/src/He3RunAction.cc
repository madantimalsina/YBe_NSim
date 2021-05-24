//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file He3RunAction.cc
/// \brief Implementation of the He3RunAction class

#include "He3RunAction.hh"
#include "He3Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3RunAction::He3RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in He3Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("My_He3");
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH1("Eabs","Edep in center tube", 100, 0., 1*MeV);
  analysisManager->CreateH1("Egap","Edep in Moderator (UHMWPE)", 100, 0., 1*MeV);
  //analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  //analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);
  analysisManager->CreateH1("Eabs1","Edep in side tube1", 100, 0., 1*MeV);
  analysisManager->CreateH1("Eabs2","Edep in side tube2", 100, 0., 1*MeV);
  analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  analysisManager->CreateH1("Labs1","trackL in absorber1", 100, 0., 1*m);
  analysisManager->CreateH1("Labs2","trackL in absorber2", 100, 0., 1*m);

  analysisManager->CreateH1("Eabs_inT","Edep in center tube w/ p", 100, 0., 1*MeV);
  analysisManager->CreateH1("Eabs1_inT","Edep in side tube1 w/ p", 100, 0., 1*MeV);
  analysisManager->CreateH1("Eabs2_inT","Edep in side tube2 w/ p", 100, 0., 1*MeV);

  analysisManager->CreateH2("Center_xy","xy central tube", 250, -500., 500., 1000, -500., 500.);
  analysisManager->CreateH2("Center_yz","yz central tube", 1000, -500., 500., 1000, -500., 500.);
  analysisManager->CreateH2("Center_zx","zx central tube", 1000, -500., 500., 250, -500., 500.);

  analysisManager->CreateH2("Left_xy","xy side tube 1 (left)", 250, -500., 500., 1000, -500., 500.);
  analysisManager->CreateH2("Left_yz","yz side tube 1 (left)", 1000, -500., 500., 1000, -500., 500.);
  analysisManager->CreateH2("Left_zx","zx side tube 1 (left)", 1000, -500., 500., 250, -500., 500.);

  analysisManager->CreateH2("Right_xy","xy side tube 2 (right)", 250, -500., 500., 1000, -500., 500.);
  analysisManager->CreateH2("Right_yz","yz side tube 2 (right)", 1000, -500., 500., 1000, -500., 500.);
  analysisManager->CreateH2("Right_zx","zx side tube 2 (right)", 1000, -500., 500., 250, -500., 500.);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("He3", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Eabs");
  analysisManager->CreateNtupleDColumn("Egap");
  //analysisManager->CreateNtupleDColumn("Labs");
  //analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->CreateNtupleDColumn("Eabs1");
  analysisManager->CreateNtupleDColumn("Eabs2");
  analysisManager->CreateNtupleDColumn("Labs");
  analysisManager->CreateNtupleDColumn("Labs1");
  analysisManager->CreateNtupleDColumn("Labs2");

  analysisManager->CreateNtupleDColumn("Eabs_inT");
  analysisManager->CreateNtupleDColumn("Eabs1_inT");
  analysisManager->CreateNtupleDColumn("Eabs2_inT");

  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3RunAction::~He3RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  //auto analysisManager = G4AnalysisManager::Instance();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Open an output file
  //
  //G4String fileName = "He3";
  //analysisManager->OpenFile(fileName);

  // Open an output file control from mac file
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl; 
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl; 
    }
    
    G4cout << " EAbs (center tube) : mean = " 
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    
    G4cout << " EGap (UHMWPE) : mean = " 
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " EAbs1 (Side tube 1) : mean = " 
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;

    G4cout << " EAbs2 (side tube 2) : mean = " 
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Energy") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Energy") << G4endl;
  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
