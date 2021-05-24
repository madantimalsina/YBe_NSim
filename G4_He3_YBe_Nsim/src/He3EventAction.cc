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
/// \file He3aEventAction.cc
/// \brief Implementation of the He3EventAction class

#include "He3EventAction.hh"
#include "He3RunAction.hh"
#include "He3Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3EventAction::He3EventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   fEnergyAbs1(0.),
   fEnergyAbs2(0.),
   fTrackLAbs1(0.),
   fTrackLAbs2(0.),
   fpositionHe3{0.},
   fpositionHe3c{0.},
   fpositionHe3s1{0.},
   fpositionHe3s2{0.}
   //fpositionX(0.),
   //fpositionY(0.),
   //fpositionZ(0.)

{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3EventAction::~He3EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3EventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
  fEnergyAbs1 = 0.;
  fEnergyAbs2 = 0.;
  fTrackLAbs1 = 0.;
  fTrackLAbs2 = 0.;
  //fpositionX = 0.;
  //fpositionY = 0.;
  //fpositionZ = 0.;
  fpositionHe3[0] = 0.0;
  fpositionHe3[1] = 0.0;
  fpositionHe3[2] = 0.0;

  fpositionHe3c[0] = 0.0;
  fpositionHe3c[1] = 0.0;
  fpositionHe3c[2] = 0.0;
  //fpositionHe3[3] = 0.0;

  fpositionHe3s1[0] = 0.0;
  fpositionHe3s1[1] = 0.0;
  fpositionHe3s1[2] = 0.0;

  fpositionHe3s2[0] = 0.0;
  fpositionHe3s2[1] = 0.0;
  fpositionHe3s2[2] = 0.0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3EventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  if (fEnergyAbs > 0) analysisManager->FillH1(0, fEnergyAbs);
  if (fEnergyGap > 0) analysisManager->FillH1(1, fEnergyGap);
  //analysisManager->FillH1(2, fTrackLAbs);
  //analysisManager->FillH1(3, fTrackLGap);
  if (fEnergyAbs1 > 0) analysisManager->FillH1(2, fEnergyAbs1);
  if (fEnergyAbs2 > 0) analysisManager->FillH1(3, fEnergyAbs2);

  if (fTrackLAbs > 0) analysisManager->FillH1(4, fTrackLAbs);
  if (fTrackLAbs1 > 0) analysisManager->FillH1(5, fTrackLAbs1);
  if (fTrackLAbs2 > 0) analysisManager->FillH1(6, fTrackLAbs2);

  // Central Tube with cuts
  // 48.05 is the z-height for 2.5" moderator thickness ( 2.5"/2 + 0.5" + 3.6 mm)

  if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
       //std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + (((fpositionHe3[2] - (58.7375)) * (fpositionHe3[2] - (58.7375)) ))) <= 11.8872
       std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
     )
      {  if (fEnergyAbs > 0)  analysisManager->FillH1(7, fEnergyAbs); } 

  // Side Tube 1 (left) with cuts 
  if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
       std::sqrt( ((fpositionHe3[1] + 127.0) * (fpositionHe3[1] + 127.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
     )
     {  if (fEnergyAbs1 > 0) analysisManager->FillH1(8, fEnergyAbs1); }

  // Side Tube 2 (right) with cuts
  if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
       std::sqrt( ((fpositionHe3[1] - 127.0) * (fpositionHe3[1] - 127.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
     )
      { if (fEnergyAbs2 > 0) analysisManager->FillH1(9, fEnergyAbs2); }


     //filling 2D histogram --> Central tube

      if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
           std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
         )
      {
        analysisManager->FillH2(0, fpositionHe3[0], fpositionHe3[1]); 
        analysisManager->FillH2(1, fpositionHe3[1], fpositionHe3[2]); 
        analysisManager->FillH2(2, fpositionHe3[2], fpositionHe3[0]);
      }  

      //filling 2D histogram --> Side tube 1 (left)

      if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
           std::sqrt( ((fpositionHe3[1] + 127.0) * (fpositionHe3[1] + 127.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
         )
      {
        analysisManager->FillH2(3, fpositionHe3[0], fpositionHe3[1]); 
        analysisManager->FillH2(4, fpositionHe3[1], fpositionHe3[2]); 
        analysisManager->FillH2(5, fpositionHe3[2], fpositionHe3[0]);
      }



      //filling 2D histogram --> Side tube 2 (Right)

      if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
           std::sqrt( ((fpositionHe3[1] - 127.0) * (fpositionHe3[1] - 127.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
         )
      {
        analysisManager->FillH2(6, fpositionHe3[0], fpositionHe3[1]); 
        analysisManager->FillH2(7, fpositionHe3[1], fpositionHe3[2]); 
        analysisManager->FillH2(8, fpositionHe3[2], fpositionHe3[0]);
      }

/*      if ( fpositionHe3s2[0] > -165.1 && fpositionHe3s2[0] < 165.1 &&
          fpositionHe3s2[1] > - (0.5 - 0.032)*25.4 && fpositionHe3s2[1] < (0.5 - 0.032)*25.4
         )
      {
        analysisManager->FillH2(6, fpositionHe3s2[0], fpositionHe3s2[1] );  
      }

      if (fpositionHe3s2[1] > - (0.5 - 0.032)*25.4 && fpositionHe3s2[1] < (0.5 - 0.032)*25.4
           && fpositionHe3s2[2] > 78. && fpositionHe3s2[2] < 104.)
      {
        analysisManager->FillH2(7, fpositionHe3s2[1], fpositionHe3s2[2]);
      }

      if ( fpositionHe3s2[0] > -165.1 && fpositionHe3s2[0] < 165.1
           && fpositionHe3s2[2] > 78. && fpositionHe3s2[2] < 104.
         )
      {
        analysisManager->FillH2(8, fpositionHe3s2[2], fpositionHe3s2[0]);
      }*/

        
/*  analysisManager->FillH2(3, fpositionHe3s1[0], fpositionHe3s1[1]);  
    analysisManager->FillH2(4, fpositionHe3s1[1], fpositionHe3s1[2]); 
    analysisManager->FillH2(5, fpositionHe3s1[2], fpositionHe3s1[0]);

    analysisManager->FillH2(6, fpositionHe3s2[0], fpositionHe3s2[1]);  
    analysisManager->FillH2(7, fpositionHe3s2[1], fpositionHe3s2[2]); 
    analysisManager->FillH2(8, fpositionHe3s2[2], fpositionHe3s2[0]);*/


  
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  analysisManager->FillNtupleDColumn(1, fEnergyGap);
  //analysisManager->FillNtupleDColumn(2, fTrackLAbs);
  //analysisManager->FillNtupleDColumn(3, fTrackLGap);
  analysisManager->FillNtupleDColumn(2, fEnergyAbs1);
  analysisManager->FillNtupleDColumn(3, fEnergyAbs2);
  analysisManager->FillNtupleDColumn(4, fTrackLAbs);
  analysisManager->FillNtupleDColumn(5, fTrackLAbs1);
  analysisManager->FillNtupleDColumn(6, fTrackLAbs2);

  // Central tubes with cuts 
  if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
       std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
     )
      {  if (fEnergyAbs > 0)  analysisManager->FillNtupleDColumn(7, fEnergyAbs); } 

  // Side Tube 1 (left) with cuts 
  if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
       std::sqrt( ((fpositionHe3[1] + 127.0) * (fpositionHe3[1] + 127.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
     )
     {  if (fEnergyAbs1 > 0) analysisManager->FillNtupleDColumn(8, fEnergyAbs1); }

  // Side Tube 2 (right) with cuts
  if ( fpositionHe3[0] > (-165.1 + (0.15*25.4))  && fpositionHe3[0] < (165.1 + (0.15*25.4)) &&
       std::sqrt( ((fpositionHe3[1] - 127.0) * (fpositionHe3[1] - 127.0)) + ((fpositionHe3[2] - 48.05) * (fpositionHe3[2] - 48.05)) ) <= 11.8872
     )
      { if (fEnergyAbs2 > 0) analysisManager->FillNtupleDColumn(9, fEnergyAbs2); }
   

  analysisManager->FillNtupleDColumn(10, fpositionHe3[0]);
  analysisManager->FillNtupleDColumn(11, fpositionHe3[1]);
  analysisManager->FillNtupleDColumn(12, fpositionHe3[2]);

  analysisManager->AddNtupleRow();  
  
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    G4cout
       << "   Center Tube: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyAbs,"Energy")
       << G4endl
       << "   Side Tube 1: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyAbs1,"Energy")
       << G4endl
       << "   Side Tube 2: total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyAbs2,"Energy")
       << G4endl
       << "   Moderator (UHMWPE): total energy: " << std::setw(7)
                                        << G4BestUnit(fEnergyGap,"Energy")
       << G4endl                                 
       << "   Center Tube: total Length: " << std::setw(7)
                                        << G4BestUnit(fTrackLAbs,"Length")  
       << G4endl                                 
       << "   Side Tube 1 (Left): total Length: " << std::setw(7)
                                        << G4BestUnit(fTrackLAbs1,"Length")
       << G4endl                                 
       << "   Side Tube 2 (Right): total Length: " << std::setw(7)
                                        << G4BestUnit(fTrackLAbs2,"Length") 
       << G4endl                                 
       << "   position Z: " << std::setw(7)
                            << G4BestUnit(fpositionHe3[0],"Length")
       << G4endl                                 
       << "   position Y: " << std::setw(7)
                            << G4BestUnit(fpositionHe3[1],"Length")
        << G4endl        
        << "   position X: " << std::setw(7)
                            << G4BestUnit(fpositionHe3[2],"Length")    
                             << G4BestUnit(fTrackLAbs1,"Length")                                       
                              << G4BestUnit(fTrackLAbs1,"Length")                                                                                                       
       << G4endl;
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
