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
/// \file He3ActionInitialization.hh
/// \brief Definition of the He3ActionInitialization class

#ifndef He3ActionInitialization_h
#define He3ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class He3DetectorConstruction;

/// Action initialization class.
///

class He3ActionInitialization : public G4VUserActionInitialization
{
  public:
    He3ActionInitialization(He3DetectorConstruction*);
    virtual ~He3ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    He3DetectorConstruction* fDetConstruction;
};

#endif

    
