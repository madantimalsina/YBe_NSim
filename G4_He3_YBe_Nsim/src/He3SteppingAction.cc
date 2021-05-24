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
/// \file He3SteppingAction.cc
/// \brief Implementation of the He3SteppingAction class

#include "He3SteppingAction.hh"
#include "He3EventAction.hh"
#include "He3DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3SteppingAction::He3SteppingAction(
                      const He3DetectorConstruction* detectorConstruction,
                      He3EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3SteppingAction::~He3SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3SteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  G4Track* track = step->GetTrack();
  //G4ThreeVector position = track->GetPosition();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();

    //G4ThreeVector position = step->GetPosition();
    //G4double time = track->GetGlobalTime();

  }
      
  if ( volume == fDetConstruction->GetAbsorberPV() ) {
    G4ThreeVector position = track->GetPosition();
    fEventAction->AddAbs(edep,stepLength, position);
  }
    if ( volume == fDetConstruction->GetAbsorberPV1() ) {
    G4ThreeVector position = track->GetPosition();
    fEventAction->AddAbs1(edep,stepLength, position);
  }
    if ( volume == fDetConstruction->GetAbsorberPV2() ) {
    G4ThreeVector position = track->GetPosition();
    fEventAction->AddAbs2(edep,stepLength, position);
  }
  
  if ( volume == fDetConstruction->GetGapPV() ) {
    //G4ThreeVector position = track->GetPosition();
    fEventAction->AddGap(edep,stepLength);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
