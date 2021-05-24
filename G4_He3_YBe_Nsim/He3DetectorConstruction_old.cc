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
/// \file He3DetectorConstruction.cc
/// \brief Implementation of the He3DetectorConstruction class

#include "He3DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4SubtractionSolid.hh"
#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

 #include "G4MultiUnion.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* He3DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3DetectorConstruction::He3DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fGapPV(nullptr),
   fAbsorberPV1(nullptr),
   fAbsorberPV2(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3DetectorConstruction::~He3DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* He3DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3DetectorConstruction::DefineMaterials()
{ 
  // Define materials

  G4double A;  // atomic mass (mass of a mole)
  G4double Z;  // atomic number (mean number of protons)
  G4double d;  // density
  G4double fractionMass; // Fraction by Mass (Weight %)
  G4double abundance;
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  //nistManager->FindOrBuildMaterial("G4_TEFLON");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE"); 
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_W");

//*******************************************************************

  //General Material defination
  G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A = 1.01*g/mole);
  G4Element* elBe = new G4Element("Beryllium","Be", Z=4, A=9.1218*g/mole);
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A = 12.011*g/mole);
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A = 16.00*g/mole);
  G4Element* elCr  = new G4Element("Chromium","Cr",Z = 24.,A = 52.00*g/mole);
  G4Element* elFe = new G4Element("Iron","Fe", Z=26., A = 55.85*g/mole);
  G4Element* elNi = new G4Element("Nickel","Ni", Z=28., A = 58.69*g/mole);
  G4Element* elW = new G4Element("Tungsten","W",  Z=74., A = 183.84*g/mole);

 //Molybdenum definition
    G4Isotope* Mo92 = new G4Isotope("Mo92", 42, 92, 91.906809 * g / mole);
    G4Isotope* Mo94 = new G4Isotope("Mo94", 42, 94, 93.9050853 * g / mole);
    G4Isotope* Mo95 = new G4Isotope("Mo95", 42, 95, 94.9058411 * g / mole);
    G4Isotope* Mo96 = new G4Isotope("Mo96", 42, 96, 95.9046785 * g / mole);
    G4Isotope* Mo97 = new G4Isotope("Mo97", 42, 97, 96.9060205 * g / mole);
    G4Isotope* Mo98 = new G4Isotope("Mo98", 42, 98, 97.9054073 * g / mole);
    G4Isotope* Mo100 = new G4Isotope("Mo100", 42, 100, 99.907477 * g / mole);

  G4Element* elMo = new G4Element("Natural Mo", "elMo", 7);
    elMo->AddIsotope(Mo92, abundance=14.84*perCent);
    elMo->AddIsotope(Mo94, abundance=9.25*perCent);
    elMo->AddIsotope(Mo95, abundance=15.92*perCent);
    elMo->AddIsotope(Mo96, abundance=16.68*perCent);
    elMo->AddIsotope(Mo97, abundance=9.55*perCent);
    elMo->AddIsotope(Mo98, abundance=24.13*perCent);
    elMo->AddIsotope(Mo100, abundance=9.63*perCent);

  //YBe Source Assembly Materials
  //MT-185 alloys of Tungsten used to make YBe Pig (LZYBePig)
  // Midwest Tungsten Service -> 97.0% (W), 2.1% (Ni), 0.9% (Fe)   
  // Mi Tech Tungsten Metals -> 97.05% (W), 2.04% (Ni), 0.91% (Fe)  
  d = 18.55*g/cm3;
  G4Material* YBeTungsten_MT185 = new G4Material("Tungsten_MT185",d,3);
  YBeTungsten_MT185->AddElement(elW, fractionMass=97.05*perCent);
  YBeTungsten_MT185->AddElement(elNi, fractionMass=2.04*perCent);
  YBeTungsten_MT185->AddElement(elFe, fractionMass=0.91*perCent);
   //YBeTungsten_MT185->AddElement(elW,fractionMass=0.97);
   //YBeTungsten_MT185->AddElement(elNi,fractionMass=0.021);
   //YBeTungsten_MT185->AddElement(elFe,fractionMass=0.009);


  // LZ YBe Source (Beryllium metal)
  d = 1.85*g/cm3;
  G4Material* Beryllium = new G4Material("Beryllium", d, 1);
  Beryllium->AddElement(elBe, fractionMass=100.0*perCent);

  // LZYBeDisk (Made up of SS-316)
  d = 7.99*g/cm3;
  G4Material* SS316 = new G4Material("SS316",d,4);
  SS316->AddElement(elFe, fractionMass=68.5*perCent);
  SS316->AddElement(elCr, fractionMass=17.0*perCent);
  SS316->AddElement(elNi, fractionMass=12.0*perCent);
  SS316->AddElement(elMo, fractionMass=2.5*perCent);
  

//************************************************************************
  // plexiglass, lucite 
  d = 1.19*g/cm3;
  G4Material* matplexiglass = new G4Material("Plexiglass",d,3);
  matplexiglass->AddElement(elH, fractionMass=8.0*perCent);
  matplexiglass->AddElement(elC, fractionMass=60.0*perCent);
  matplexiglass->AddElement(elO, fractionMass=32.0*perCent);

  //Volume of the He3 tubes
  G4double volume_He3 = (13*2.54*pi*1.18872*1.18872);

  // Build He3 gas
  G4int protons=2, neutrons=1, nucleons=protons+neutrons;
  G4double elements;
  G4double atomicMass_He3 = 3.016*g/mole; //molar mass
  G4Isotope* isoHe3 = new G4Isotope("He3", protons, nucleons, atomicMass_He3);
  G4Element* elHe3 = new G4Element("Helium3", "He3", 1);
  elHe3->AddIsotope(isoHe3, 100*perCent);
  G4double pressure_He3 = 11.02*atmosphere;
  G4double temperature = 293.15*kelvin;
  G4double molar_constant = Avogadro*k_Boltzmann;
  //G4double density = (atomicMass*pressure)/(temperature*molar_constant);
  G4double density_He3 = (atomicMass_He3*pressure_He3)/(temperature*molar_constant);

  G4Material* Helium3 = new G4Material("Helium3", density_He3, elements=1, kStateGas, temperature, pressure_He3);
  Helium3->AddElement(elHe3, fractionMass=100*perCent);

  // Argon 
  G4double atomicMass_Ar = 39.948*g/mole;
  G4double pressure_Ar = 0.58*atmosphere;
  G4double density_Ar = (atomicMass_Ar*pressure_Ar)/(temperature*molar_constant);
  G4Element* elAr = new G4Element("Argon", "Ar", Z=18., atomicMass_Ar);
  G4Material* Argon = new G4Material("Argon"  , density_Ar, 1,  kStateGas, temperature, pressure_Ar);
  Argon->AddElement(elAr, 1);

  // 95% He3 + 4.95% Ar + 0.05% CH4,
  G4double atomicMass_He3Ar = ((0.95 *3.016) + (0.05* 39.948))*g/mole;
  G4double pressure_He3Ar = 11.60*atmosphere;
  G4double density_He3Ar = (atomicMass_He3Ar*pressure_He3Ar)/(temperature*molar_constant);

  G4Material* He3Ar = new G4Material("He3Ar"  , density_He3Ar, 2,  kStateGas, temperature, pressure_He3Ar);
  He3Ar->AddMaterial( Helium3, (0.95 *3.016) / ((0.95 *3.016) + (0.05* 39.948)) );
  He3Ar->AddMaterial( Argon, (0.05* 39.948) / ((0.95 *3.016) + (0.05* 39.948)) );

  // UHMW (Ultra High Molecular Weight Polyethylene)
  d = 0.94*g/cm3;
  nistManager->BuildMaterialWithNewDensity("UHMWPE","G4_POLYETHYLENE",d);


  // UHMW for Thermal Scattering of neutron
  G4Element *H = new G4Element("TS_H_of_Polyethylene", "H", 1., 1.0079*g/mole);
  G4Material *POLYMAT = new G4Material("POLYMAT", 0.94*g/cm3, 2, kStateSolid, 293.15*kelvin);
  POLYMAT->AddElement(H, fractionMass=14.3711*perCent);
  POLYMAT->AddElement(elC, fractionMass=85.6289*perCent);


  // POLYPROPYLENE
  d = 913.43685*mg/cm3; // Geant4 density: 900.000 mg/cm3
  nistManager->BuildMaterialWithNewDensity("POLYPROPYLENE","G4_POLYPROPYLENE",d);

  // POLYPROPYLENE for Thermal Scattering of neutron
  G4Material *POLYMAT2 = new G4Material("POLYMAT2", 913.43685*mg/cm3, 2, kStateSolid, 293.15*kelvin);
  POLYMAT2->AddElement(H, fractionMass=14.3711*perCent);
  POLYMAT2->AddElement(elC, fractionMass=85.6289*perCent);


  G4cout<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl; 
  G4cout << "MY CHECK....!!!!  " << "pressure_He3  =  " << pressure_He3/atmosphere << G4endl;
  G4cout << "MY CHECK....!!!!  " << "pressure_Ar  =  " << pressure_Ar/atmosphere << G4endl;
  G4cout << "MY CHECK....!!!!  " << "pressure_He3Ar  =  " << pressure_He3Ar/atmosphere << G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_He3  =  " << density_He3/(mg/cm3) << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_Ar  =  " << density_Ar/(mg/cm3) << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_He3Ar  =  " << density_He3Ar/(mg/cm3) << G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl; 
  G4cout << "MY CHECK....!!!!  " << "volume of the tube  =  " << volume_He3 << " cm3"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_He3  =  " << (density_He3/(mg/cm3))*volume_He3 << " mg" << G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_Ar  =  " << (density_Ar/(mg/cm3))*volume_He3 << " mg"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_He3Ar  =  " << (density_He3Ar/(mg/cm3))*volume_He3 << " mg"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* He3DetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;

  // Outer layer of He3 tube
  G4double outerRadius_THe3 = 25.4/2*mm;
  G4double thickness_AlTube = (0.032*25.4)*mm;


  // Parameters: Inner layer of He3
  // Tubes filled with He3 gas
  G4double outerRadius_gHe3 = outerRadius_THe3 - thickness_AlTube;
  G4double innerRadius_gHe3 = 0.*mm;

  // This is the effective length (13") of He-3 tubes
  // divide by 2 because of G4 style (-33.02/2 to +33.02/2)
  G4double hz_gHe3 = 33.02/2*cm; 
  //G4double startAngle = 0.*deg;
  //G4double spanningAngle = 360.*deg;

  // Gap between UHMWPE and outer layer of He3 tube
  // Polypropylene U-Shaped Channel of Dimension,   
  // 2-1/2"(2.5" --> 6.35 Cm) Base,   
  // 2-1/4"(2.25" --> 5.715 Cm ) Legs X 3/16" Wall Thick
  // gapPEHe3T = Legs - Diameter of He3 Tubes
  //G4double gapPEHe3T = 5.715*cm - 2*outerRadius_THe3; // Legs 2.25" (5.715 cm)
  
  // Rotation of the Tubes 
  // (both outer and inner filled with He3 gas)
  G4RotationMatrix* rotD1 = new G4RotationMatrix();
  G4RotationMatrix* rotD2 = new G4RotationMatrix();
  G4RotationMatrix* rotD3 = new G4RotationMatrix();
  rotD1->rotateZ(90.*deg);
  rotD2->rotateX(90.*deg);
  rotD3->rotateY(90.*deg);


  G4int nofLayers = 1;
  G4double gapThickness = 2.5*2.540*cm;  // Thickness of moderator (Varies: 1", 2", 2.5", 5" etc)
  G4double calorSizeY  = 30.48*cm;   // 12" 
  G4double calorSizeX  = 30.48*cm;  // 12"
  //G4double calorSizeX  = 31.115*cm;  // 12.25"
  G4double calorSizeXY = 62.0*cm;

  auto layerThickness2 = 100.0*cm;
  auto calorThickness2 = nofLayers * layerThickness2;
  auto worldSizeXY = 5.0 * calorSizeX; 
  auto worldSizeZ  = 5.0 * (layerThickness2);
  
  // Parameters: source holder YBe source
  G4double innerRadius_sHold = 0.*cm;
  G4double outerRadius_sHold = 10.*cm;
  G4double hz_YBesHold = (20./2.)*cm;
  G4double outerRadius_YBesHoldCap = (10.52/2.0)*cm;
  G4double hz_YBesHoldCap = (4./2.)*cm;
  G4double Position_YBesHold = gapThickness/2 + hz_YBesHold; // just add thickness here for table height etc.
  G4double Position_YBesHoldCap = Position_YBesHold+ hz_YBesHold + hz_YBesHoldCap;
  G4double outerRadius_YBeBeO = (2.54/2.)*cm;
  G4double hz_YBeBeO = (6.46/2.)*cm;
  G4double Position_YBeBeO =  Position_YBesHold + 2.0*cm + 3.23*cm; //Bottom of BeO start at 12 cm from base of YBe Pig
  G4double outerRadius_YDisk = (2.35)*mm;
  G4double hz_YDisk = (4.6/2.)*mm;
  G4double Position_YDisk = Position_YBeBeO + 1.23*cm;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
  //auto absorberMaterial = G4Material::GetMaterial("Helium3");
  auto absorberMaterial = G4Material::GetMaterial("He3Ar");
  //auto gapMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  //auto gapMaterial = G4Material::GetMaterial("UHMWPE");
  auto gapMaterial = G4Material::GetMaterial("POLYMAT");
  //auto gapMaterial2 = G4Material::GetMaterial("POLYPROPYLENE");
  auto gapMaterial2 = G4Material::GetMaterial("POLYMAT2");
  //auto holdMaterial = G4Material::GetMaterial("Plexiglass");
  auto gapMaterialT = G4Material::GetMaterial("G4_Al");
  //auto holdMaterial_Tungsten = G4Material::GetMaterial("G4_W");
  auto holdMaterial_Tungsten_MT185 = G4Material::GetMaterial("Tungsten_MT185");
  auto BeMaterial = G4Material::GetMaterial("Beryllium");
  auto YMaterial = G4Material::GetMaterial("SS316");

  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("He3DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
  // Volume defination:
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness2/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                                 
  // Layer
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness2/2.); // its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness2);  // witdth of replica

  //
  // Source holder for YBe neutron source
  //YBe source Pig Cap 
  //             
  G4Tubs* sourceTubeYBeCap = 
    new G4Tubs("SourceCap", innerRadius_sHold, outerRadius_YBesHoldCap, hz_YBesHoldCap, startAngle, spanningAngle);

  G4LogicalVolume* logicRing3 =                         
    new G4LogicalVolume(sourceTubeYBeCap,            //its solid
                        holdMaterial_Tungsten_MT185, //its material
                        "SourceCap"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -Position_YBesHoldCap), // its position (gapThickness/2 + hz_YBesHold)
                 logicRing3,            // its logical volume                         
                 "SourceCap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  //
  //YBe Pig
  //               
  G4Tubs* sourceTube = 
    new G4Tubs("SourcePig", innerRadius_sHold, outerRadius_sHold, hz_YBesHold, startAngle, spanningAngle);
     
     // outerRadius_YBeBeO+0.3937*mm, 0.7874/2 added because of 1.03 inch of source inner can diameter 
     // hz_YBeBeO+2.63*mm, 5.25/2 added because of 2.75 inch of source inner can height 

  G4Tubs* sourceTube_BeOVol = 
    new G4Tubs("Source_BeoVol", innerRadius_sHold, outerRadius_YBeBeO+0.3937*mm, hz_YBeBeO+2.63*mm, startAngle, spanningAngle); 
  
     // substraction of Beo overlap volume from Tungsten Pig
  G4VSolid* subtract_BeoVol = new G4SubtractionSolid("SourcePig-Source_BeoVol", sourceTube, sourceTube_BeOVol, 0, G4ThreeVector(0., 0., -5.18*cm));

  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(subtract_BeoVol,           //its solid
                        holdMaterial_Tungsten_MT185,//its material
                        "SourcePig"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -Position_YBesHold), // its position (gapThickness/2 + hz_YBesHold)
                 logicRing,            // its logical volume                         
                 "SourcePig",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Source BeO
  //
  G4Tubs* sourceTubeBeO = 
    new G4Tubs("SourceBeo", innerRadius_sHold, outerRadius_YBeBeO, hz_YBeBeO, startAngle, spanningAngle);

  G4Tubs* sourceTube_YDiskVol = 
    new G4Tubs("Source_YDiskVol", innerRadius_sHold, outerRadius_YDisk+0.01*mm, hz_YDisk+0.04*mm, startAngle, spanningAngle);

    // substraction of overlap Y-88 Disk volume from BeO
  G4VSolid* subtract_YDiskVol = new G4SubtractionSolid("SourceBeo-Source_YDiskVol", sourceTubeBeO, sourceTube_YDiskVol, 0, G4ThreeVector(0., 0., -1.23357*cm));
  
  G4LogicalVolume* logicRing1 =                         
    new G4LogicalVolume(subtract_YDiskVol,           //its solid
                        BeMaterial,                  //its material
                        "SourceBeo"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(Position_YBeBeO)), // its position (gapThickness/2 + hz_YBesHold)
                 logicRing1,            // its logical volume                         
                 "SourceBeo",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 



  //
  // source YDisk
  //
     G4Tubs* sourceTubeYDisk = 
    new G4Tubs("SourceYDisk", innerRadius_sHold, outerRadius_YDisk, hz_YDisk, startAngle, spanningAngle);

  G4LogicalVolume* logicRing2 =                         
    new G4LogicalVolume(sourceTubeYDisk,           //its solid
                        YMaterial,                 //its material
                        "SourceYDisk"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(Position_YDisk)), // its position (gapThickness/2 + hz_YBesHold)
                 logicRing2,            // its logical volume                         
                 "SourceYDisk",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

//***************************************************************************************************//                               
  // Gap
  // Here Gap material is UHMWPE as neutron moderator
  //
  auto gapS 
    = new G4Box("Gap",             // its name
                 calorSizeX/2, calorSizeY/2, gapThickness/2); // its size
                         
  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "Gap");           // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0.), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
//polypropylene U-shaped channel for He3 tubes support
// below the tube
//(3/16)*25.4*mm thicknesss of polypropylene U-shaped channel for He3 tubes support
// (2.25 + 3/16")*25.4*mm is total length of legs of polypropylene U-shaped channel for He3 tubes support
// (2.50 - (2*3/16))*25.4*mm is the base of the U-shapped channel   
// A constant value 59.53125*mm --> (2.25 + 3/(2*16)")*25.4*mm legs of polypropylene U-shaped channel for He3 tubes support
// A constant value 120.65*mm --> (6 - 2.125/2 -3/16) distance from centre to centre of u channel on the sides
  
  G4double gapThicknessB2 = 4.7625*mm; //(0.1875)*25.4*mm polypropylene U-shaped channel for He3 tubes support
  G4double gapThicknessB2z = 61.9125*mm; //(2.25 + 3/16")*25.4*mm legs of polypropylene U-shaped channel for He3 tubes support
  
  auto gapS2b 
    = new G4Box("Gap2b",calorSizeY/2, (2.125*2.54*cm)/2, gapThicknessB2/2);
                         
  auto gapLV2b
    = new G4LogicalVolume(
                 gapS2b,             // its solid
                 gapMaterial2,      // its material
                 "Gap2b");           // its name
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (gapThickness/2 + 59.53125)*mm), // 59.53125 --> //(2.25 + 3/(2*16)")*25.4*mm legs of polypropylene U-shaped channel for He3 tubes support
                 gapLV2b,            // its logical volume                         
                 "Gap2b",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

// making a copy 
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 120.65*mm, (gapThickness/2 + 59.53125)*mm),
                 gapLV2b,            // its logical volume                         
                 "Gap2b_S1",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -120.65*mm, (gapThickness/2 + 59.53125)*mm),
                 gapLV2b,            // its logical volume                         
                 "Gap2b_S2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps  

//inner side wall of the U-channel (double thick wall)
// +/- 26.9875 + 2.38125 mm --> (2.125/2 + 3/32)" is the side for the central U rail
//+/- 150.01875 mm --> (6-(3/32))" is the side for the end side of the ednd U-rail
// +/- 91.28125 mm --> (6 -2.125 -3/16-3/32)" is the side for the near side of the side U-rail  
  auto gapS2s 
    = new G4Box("Gap2s",             // its name
                 calorSizeY/2, gapThicknessB2/2, gapThicknessB2z/2); // its size
                         
  auto gapLV2s
    = new G4LogicalVolume(
                 gapS2s,             // its solid
                 gapMaterial2,      // its material
                 "Gap2s");           // its name
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 29.36875*mm, gapThickness/2. + gapThicknessB2z/2.),
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s_Ca",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -29.36875*mm, gapThickness/2. + gapThicknessB2z/2.), 
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s_Cb",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 91.28125*mm, gapThickness/2. + gapThicknessB2z/2.), 
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s_S1a",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps  

      fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -91.28125*mm, gapThickness/2. + gapThicknessB2z/2.), 
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s_S1b",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 150.01875*mm, gapThickness/2. + gapThicknessB2z/2.),
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s2_s2a",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -150.01875*mm, gapThickness/2. + gapThicknessB2z/2.),
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s2_s2b",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  


//For side (front and back) U-rail to support He3 tubes
G4double Support_He3_x = (4.7625)/2.*mm; //(3/16")
G4double Support_He3_y = (12.0*25.4)/2.*mm;
G4double Support_He3_z = (2.25*25.4)/2.*mm;
// for base
G4double Support_He3b_x = (2.25*25.4)/2.*mm;
G4double Support_He3b_y = (12.0*25.4)/2.*mm;
G4double Support_He3b_z = (4.7625)/2.*mm; //(3/16")

  //The hole diameter should then be (25.4 + 2.6) mm = 28.0 mm
  //That makes the top gap from top of the tubes to top of L-rail (1.0 + 2.6) mm = 3.6 mm
  // z-position is 57.15/2 - 28/2 + 1 = 15.575
  //And in case you're confused about the math: 
  //12.0" total width of trimmed L-rail / 2  =  6.0"  =  w/2 in Madan's note + 1.0" distance from outer edge of L-rail
  //So center of holes are at +/-2.5" and +/-5.0"
  //Note: central hole at 0.00 mm, 2-holes at +/- 63.50 mm ant other 2-holes at +/- 127.00 mm (SolidWorks drawing He3TubeSupportPlate_v2.pdf)

  auto Support_He3 = new G4Box("S_He3",Support_He3_x, Support_He3_y, Support_He3_z);

  //The hole diameter should then be (25.4 + 2.6) mm = 28.0 mm
  // I just assume 14.25 mm (actually 14 mm) to remove some of the overlap 
  G4Tubs* Support_He3T_hole =
    new G4Tubs("He3T_hole", innerRadius_sHold, 14.0*mm, 31.75*mm, startAngle, spanningAngle); // 31.75 --> (2.5/2" base of U channel)

    // Combining 3- holes
      G4RotationMatrix rotm  = G4RotationMatrix();
      G4ThreeVector position1 = G4ThreeVector(0., 0., -12.275);
      //G4ThreeVector position2 = G4ThreeVector(0., (155.58-28.575), -15.575); // (155.58-28.575) --> 1" from 12" moderator and 1.125 from 12.25" u-rail
      G4ThreeVector position2 = G4ThreeVector(0., (152.40-25.4), -12.275); // (152.40-25.4) --> 1" from 12" moderator
      G4ThreeVector position3 = G4ThreeVector(0., -(152.40-25.4), -12.275);
      G4ThreeVector position4 = G4ThreeVector(0., 63.50, -12.275);
      G4ThreeVector position5 = G4ThreeVector(0., -63.50, -12.275);
      G4Transform3D tr1 = G4Transform3D(rotm,position1);
      G4Transform3D tr2 = G4Transform3D(rotm,position2);
      G4Transform3D tr3 = G4Transform3D(rotm,position3);
      G4Transform3D tr4 = G4Transform3D(rotm,position4);
      G4Transform3D tr5 = G4Transform3D(rotm,position5);
  G4MultiUnion* union_holes = new G4MultiUnion("holes_Union");
      union_holes->AddNode(*Support_He3T_hole,tr1);
      union_holes->AddNode(*Support_He3T_hole,tr2);
      union_holes->AddNode(*Support_He3T_hole,tr3);
      union_holes->AddNode(*Support_He3T_hole,tr4);
      union_holes->AddNode(*Support_He3T_hole,tr5);
      union_holes->Voxelize();

    
  // substraction of overlap volume
   G4VSolid* subtract_He3T_hole = new G4SubtractionSolid("S_He3-He3T_holes", Support_He3, union_holes, rotD3, G4ThreeVector(0., 0., -12.275));
  //G4VSolid* subtract_He3T_hole = new G4SubtractionSolid("S_He3-He3T_hole", Support_He3, Support_He3T_hole, rotD3, G4ThreeVector(0., 0., -12.275));
  //G4VSolid* subtract_He3T_hole_s1 = new G4SubtractionSolid("S_He3-He3T_hole1", Support_He3, Support_He3T_hole, rotD3, G4ThreeVector(0., 155.58, -12.275));
                      
  auto Support_He3LV = new G4LogicalVolume(subtract_He3T_hole, gapMaterial2,"S_He3-He3T_holes");  
  //auto Support_He3LV_s1 = new G4LogicalVolume(subtract_He3T_hole_s1, gapMaterial,"S_He3-He3T_hole"); 
                                             
  fGapPV2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(152.4, 0., gapThickness/2 + Support_He3_z), // its position
                 Support_He3LV,            // its logical volume                         
                 "S_He3-He3T_holes_f1",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  fGapPV2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(-152.4, 0., gapThickness/2 + Support_He3_z), // its position
                 Support_He3LV,            // its logical volume                         
                 "S_He3-He3T_holes_b1",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  fGapPV2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(((6*25.4) +(2*Support_He3b_x-2.38125)), 0., gapThickness/2 + Support_He3_z ), // its position
                 Support_He3LV,            // its logical volume                         
                 "S_He3-He3T_holes_f2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 2,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  fGapPV2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(-((6*25.4) +(2*Support_He3b_x-2.38125)), 0., gapThickness/2 + Support_He3_z), // its position
                 Support_He3LV,            // its logical volume                         
                 "S_He3-He3T_holes_b2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 3,                // copy number
                 fCheckOverlaps);  // checking overlaps 



//For Base U-rail to support He3 tubes

  auto Support_He3b = new G4Box("S_He3b",Support_He3b_x, Support_He3b_y, Support_He3b_z);                      
  auto Support_He3LVb = new G4LogicalVolume(Support_He3b, gapMaterial2,"S_He3b");                                 
  fGapPV2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(6*25.4 + Support_He3_z, 0., gapThickness/2. + 2*Support_He3_z + 2.38125*mm), // its position
                 Support_He3LVb,            // its logical volume                         
                 "S_He3b1",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  fGapPV2 = new G4PVPlacement(
                 0,                // no rotation
                 //G4ThreeVector(-(6*25.4+Support_He3_z), 0., gapThickness/2. +59.53125*mm), // its position
                 G4ThreeVector(-(6*25.4+Support_He3_z), 0., gapThickness/2 + 2*Support_He3_z + 2.38125*mm), // its position
                 Support_He3LVb,            // its logical volume                         
                 "S_He3b2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  //
  // For He3 tube outerlayer
  // Tube 1 (center tube) outerlayer
  // centre tube is off by 0.635cm(0.25")
  // Length of the tube (not symmetric)
  // 1.6"-1.25" = 0.35" active length extending out of moderator at NEAR END
  // 1.6"-0.95" = 0.65" active length extending out of moderator at FAR END
  // So the X-center is pushed by the 0.15"   

  //The hole diameter should then be (25.4 + 2.6) mm = 28.0 mm
  //That makes the top gap from top of the tubes to top of L-rail (1.0 + 2.6) mm = 3.6 mm

  G4double Z_He3Tube = (gapThickness/2.) + outerRadius_THe3 + 3.6*mm; 
  G4double Z_He3gas = (gapThickness/2.) + outerRadius_THe3 + 3.6*mm;
  
  // Almunium layer 
  auto gapsT1 
    = new G4Tubs("T1", outerRadius_gHe3, outerRadius_THe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto gapLVT1
    = new G4LogicalVolume(
                 gapsT1,        // its solid
                 gapMaterialT, // its material
                 "T1");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 //G4ThreeVector((0.15 * 2.54)*cm, (0.5625*2.54)*cm, 58.7375*mm),
                 G4ThreeVector((0.15 * 2.54)*cm, 0, Z_He3Tube),
                 gapLVT1,       // its logical volume                         
                 "T1",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // Tube 2 outerlayer
  auto gapsT2 
    = new G4Tubs("T2", outerRadius_gHe3, outerRadius_THe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto gapLVT2
    = new G4LogicalVolume(
                 gapsT2,        // its solid
                 gapMaterialT, // its material
                 "T2");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, -127.0*mm, Z_He3Tube),
                 gapLVT2,       // its logical volume                         
                 "T2",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps  
  
  // Tube 3 outerlayer
  auto gapsT3 
    = new G4Tubs("T3", outerRadius_gHe3, outerRadius_THe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto gapLVT3
    = new G4LogicalVolume(
                 gapsT3,        // its solid
                 gapMaterialT, // its material
                 "T3");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, 127.0*mm, Z_He3Tube), // its position
                 gapLVT3,       // its logical volume                         
                 "T3",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Absorber
  // Absorber is the He3 gas
  // For He3 gas Tube 1 (Center tube) 
  // Centre tube is off by 0.635cm (0.25")
  // Length of the tube (not symmetric)
  // 1.6"-1.25" = 0.35" active length extending out of moderator at NEAR END
  // 1.6"-0.95" = 0.65" active length extending out of moderator at FAR END
  // So the X-center is pushed by the 0.15"  
  auto absorberS 
    = new G4Tubs("Abso", innerRadius_gHe3, outerRadius_gHe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "Abso");          // its name
                                   
  fAbsorberPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, 0, Z_He3gas), // its position Centre tube is off by 0.635cm (0.25")
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // For He3 gas Tube 2
  auto absorberS2 
    = new G4Tubs("Abso2", innerRadius_gHe3, outerRadius_gHe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto absorberLV2
    = new G4LogicalVolume(
                 absorberS2,        // its solid
                 absorberMaterial, // its material
                 "Abso2");          // its name
                                   
  fAbsorberPV1
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, -127*mm, Z_He3gas), // its position
                 absorberLV2,       // its logical volume                         
                 "Abso2",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // For He3 gas Tube 3
  auto absorberS3 
    = new G4Tubs("Abso3", innerRadius_gHe3, outerRadius_gHe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto absorberLV3
    = new G4LogicalVolume(
                 absorberS3,        // its solid
                 absorberMaterial, // its material
                 "Abso3");          // its name
                                   
  fAbsorberPV2
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, 127*mm, Z_He3gas), // its position
                 absorberLV3,       // its logical volume                         
                 "Abso3",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // print parameters
  ///////
  G4cout<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "Gap Thickness (Moderator Thickness)  =  " << gapThickness << " mm"<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Position of Y-88 source disk  =  " << Position_YDisk << " mm"<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Position of center of BeO volume  =  " << Position_YBeBeO << " mm"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl; 
  G4cout << "MY CHECK....!!!!  " << "X-position of all He3 tubes (const. /w M.T.)  =  " << 0.15*2.54 << " mm "<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Y-position of 2 side He3 tubes (const. /w M.T.)  =  +/-" << 5.0*25.4 << " mm "<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Y-position of central He3 tube (const. /w M.T.)  =  " << 0.0*25.4 << " mm "<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Z-position of all He3 tubes (change /w M.T.)  =  " << Z_He3Tube << " mm "<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  ///////
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofLayers << " layers of: [ YBe source assembly + "
    << gapThickness/mm << " mm of " << gapMaterial->GetName() << " + "
    << 2*outerRadius_THe3/mm << " mm of " << absorberMaterial->GetName() << " Tubes" << G4endl
    << " and Polypropylene U-shaped channel to support He-3 tubes] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //

  //G4VisAttributes* grey_solid = new G4VisAttributes(G4Colour::Grey());
  //G4VisAttributes* grey_clear = new G4VisAttributes(G4Colour(G4Colour::Grey().GetRed(), G4Colour::Grey().GetGreen(), G4Colour::Grey().GetBlue(), 0.4));
  G4VisAttributes* blue_clear = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.4));
  G4VisAttributes* aqua_clear = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.4));
  G4VisAttributes* brown_clear = new G4VisAttributes(G4Colour(0.45, 0.25, 0.0, 0.7));
  G4VisAttributes* yellow_clear = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));
  //G4VisAttributes* m_clear = new G4VisAttributes(G4Colour(0.45, 0.25, 0.0, 0.1));
  
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  calorLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  layerLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  //layerLV->SetVisAttributes(G4Colour(0.0,0.0,0.0)); //Black color

  
  //SBoxLV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
  logicRing->SetVisAttributes(brown_clear);
  logicRing1->SetVisAttributes(yellow_clear);
  logicRing2->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  logicRing3->SetVisAttributes(brown_clear);

  absorberLV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  absorberLV2->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  absorberLV3->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  gapLVT1->SetVisAttributes(blue_clear);
  gapLVT2->SetVisAttributes(blue_clear);
  gapLVT3->SetVisAttributes(blue_clear);
  //gapLVT1_a->SetVisAttributes(m_clear);

  //gapLVT1->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  //gapLVT2->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  //gapLVT3->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

  gapLV->SetVisAttributes(aqua_clear);
  //gapLV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

  //auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,2.0,3.0));
  //simpleBoxVisAtt->SetVisibility(true);
  //calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
