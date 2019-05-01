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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPWorld(0), fLWorld(0), fWorldMater(0), fDetectorMessenger(0)
{
  fWorldSizeX = 2*m;
  fWorldSizeY = 2*m;
  fWorldSizeZ = 30*m;
  DefineMaterials();
  
  // liquid argon target
  fTargetLength = 150.*cm; 
  fTargetRadius = 2.*cm; 
  
  // liquid argon container
  fLArcontainerLength = 150.*cm; 
  fLArcontainerRadius = 50.*cm; 
  
  // gas insulator
  fGasinsulatorLength = 160.*cm; 
  fGasinsulatorRadius = 55.*cm; 
  
  // gas container
  fGascontainerLength = 160.*cm; 
  fGascontainerRadius = 57.*cm; 
  
  // kapton window to be defined
  
  
  // neutron collimator to be defined
  fCollimatorShieldThickness = 10.*cm;
  fCollimatorHollowLength = 90.*cm;   
  fCollimatorHollowRadius = fTargetRadius;
  fCollimatorSolidLength = fCollimatorHollowLength + fCollimatorShieldThickness;    
  fCollimatorSolidRadius = fCollimatorHollowRadius + fCollimatorShieldThickness;    

  
  // toy neutron detector
  fDetectorLength = 10.*cm;
  fDetectorRadius = 3.*cm;
  fDetectorPositionZ = -10.*m;
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)

  G4int Z, A, a, ncomponents, natoms;
  G4double fractionmass, abundance;

  // pressurized water
  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Material* H2O = 
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                         kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  // heavy water
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);  
  G4Material* D2O = new G4Material("HeavyWater", 1.11*g/cm3, ncomponents=2,
                        kStateLiquid, 293.15*kelvin, 1*atmosphere);
  D2O->AddElement(D, natoms=2);
  D2O->AddElement(O, natoms=1);
  
  // world material
  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);    
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2, kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);

  // vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  G4Material* Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);
  
  // graphite
  G4Isotope* C12 = new G4Isotope("C12", 6, 12);  
  G4Element* C   = new G4Element("TS_C_of_Graphite","C", ncomponents=1);
  C->AddIsotope(C12, 100.*perCent);
  G4Material* graphite = new G4Material("graphite", 2.27*g/cm3, ncomponents=1,
                         kStateSolid, 293*kelvin, 1*atmosphere);
  graphite->AddElement(C, natoms=1);  
  
  // Define elements for all materials not found in the NIST database  
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Li = man->FindOrBuildElement("Li");
  G4Element* B = man->FindOrBuildElement("B");
  G4Element* F = man->FindOrBuildElement("F");
  G4Element* Na = man->FindOrBuildElement("Na");
  G4Element* Mg = man->FindOrBuildElement("Mg");
  G4Element* Al = man->FindOrBuildElement("Al");
  G4Element* Si = man->FindOrBuildElement("Si");
  G4Element* K = man->FindOrBuildElement("K");
  G4Element* Ca = man->FindOrBuildElement("Ca");
  G4Element* Ti = man->FindOrBuildElement("Ti");
  G4Element* Cr = man->FindOrBuildElement("Cr");
  G4Element* Mn = man->FindOrBuildElement("Mn");
  G4Element* Fe = man->FindOrBuildElement("Fe");
  G4Element* Ni = man->FindOrBuildElement("Ni");
  G4Element* Sb = man->FindOrBuildElement("Sb");
  G4Element* Xe = man->FindOrBuildElement("Xe");
  G4Element* Cs = man->FindOrBuildElement("Cs");
  G4Element* Bi = man->FindOrBuildElement("Bi");
  

  
  // stainless steel
  G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
      StainlessSteel->AddElement(C, fractionmass=0.015); 
      StainlessSteel->AddElement(Si, fractionmass=0.008);
      StainlessSteel->AddElement(Cr, fractionmass=0.18);
      StainlessSteel->AddElement(Mn, fractionmass=0.01);
      StainlessSteel->AddElement(Fe, fractionmass=0.697);
      StainlessSteel->AddElement(Ni, fractionmass=0.09);

	
  // MgF2
  G4Material* MgF2 = new G4Material("MgF2", 3.15*g/cm3, ncomponents=2, kStateSolid);
      MgF2->AddElement(Mg, natoms=1);
      MgF2->AddElement(F, natoms=2);

  
  // TiF3
  G4Material* TiF3 = new G4Material("TiF3", 3.4*g/cm3, ncomponents=2, kStateSolid);
      TiF3->AddElement(Ti, natoms=1);
      TiF3->AddElement(F, natoms=3);

	
  // Fe-56 isotope
  G4Isotope* iso_Fe = new G4Isotope("iso_Fe", Z=26, A=56, a=55.9349363*g/mole);
  G4Element* ele_Fe = new G4Element("ele_Fe", "Fe", ncomponents=1);
  ele_Fe->AddIsotope(iso_Fe,abundance=100.*perCent);
  G4Material* mat_Fe=new G4Material("mat_Fe",7.874*g/cm3, ncomponents = 1);
	mat_Fe->AddElement(ele_Fe, fractionmass = 1 );

	
  // Li-6 isotope
  G4Isotope* iso_Li = new G4Isotope("iso_Li", Z=3, A=6, a=6.015122795*g/mole);
  G4Element* ele_Li = new G4Element("ele_Li", "Li", ncomponents=1);
  ele_Li->AddIsotope(iso_Li,abundance=100.*perCent);
  G4Material* mat_Li=new G4Material("mat_Li",0.534*g/cm3, ncomponents = 1);
    mat_Li->AddElement(ele_Li, fractionmass = 1 );

  //Lithium Polyethylene
  G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* LiPoly = new G4Material("LiPoly", 1.06*g/cm3, ncomponents=2);
    LiPoly->AddElement (Li, 7.54*perCent);
    LiPoly->AddMaterial (polyethylene, 92.46*perCent);
  	
  // world mater
  fWorldMater = Air20;
  
  // gas container
  fGascontainerMater = StainlessSteel;
  
  // gas insulator
  fGasinsulatorMater = Air20;
  
  // LAr container
  fLArcontainerMater = LiPoly;
  	
  // liquid argon target
  fTargetMater = man->FindOrBuildMaterial("G4_lAr");
  //fTargetMater = Vacuum;
  
  // neutron collimator
  fCollimatorMater = LiPoly;
  
  // neutron detector
  fDetectorMater = H2O;

  
 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // world
  G4Box*
  sWorld = new G4Box("Container", 
                   fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2); 

  fLWorld = new G4LogicalVolume(sWorld, 
                             fWorldMater,
                             fWorldMater->GetName()); 

  fPWorld = new G4PVPlacement(0, 
                            G4ThreeVector(),            
                            fLWorld,                    
                            fWorldMater->GetName(),     
                            0,                          
                            false,                      
                            0);                         
  
  // gas container
  G4Tubs* 
  sGascontainer = new G4Tubs("Gascontainer_s",                                               
                        0, fGascontainerRadius, 0.5*fGascontainerLength, 0.,CLHEP::twopi);   

  fLogicGascontainer = new G4LogicalVolume(sGascontainer,     
                                      fGascontainerMater,     
                                      "Gascontainer_l");      

  fPhysiGascontainer = new G4PVPlacement(0,                   
                            				G4ThreeVector(),          
                                    fLogicGascontainer ,      
                                    "Gascontainer_p",         
                                    fLWorld,                  
                                    false,                    
                                    0);                                                      
  // gas insulator
  G4Tubs* 
  sGasinsulator = new G4Tubs("Gasinsulator_s",                                               
                        0, fGasinsulatorRadius, 0.5*fGasinsulatorLength, 0.,CLHEP::twopi);   

  fLogicGasinsulator = new G4LogicalVolume(sGasinsulator,     
                                      fGasinsulatorMater,     
                                      "Gasinsulator_l");      

  fPhysiGasinsulator = new G4PVPlacement(0,                   
                            				G4ThreeVector(),          
                                    fLogicGasinsulator ,      
                                    "Gasinsulator_p",         
                                    fLogicGascontainer,       
                                    false,                    
                                    0);                       
                                    
  // LAr container
  G4Tubs* 
  sLArcontainer = new G4Tubs("LArcontainer_s",                                               
                        0, fLArcontainerRadius, 0.5*fLArcontainerLength, 0.,CLHEP::twopi);   

  fLogicLArcontainer = new G4LogicalVolume(sLArcontainer,   
                                      fLArcontainerMater,   
                                      "LArcontainer_l");    

  fPhysiLArcontainer = new G4PVPlacement(0,                 
                            				G4ThreeVector(),        
                                    fLogicLArcontainer ,    
                                    "LArcontainer_p",       
                                    fLogicGasinsulator,     
                                    false,                  
                                    0);                                                      
  
  // liquid argon target
  G4Tubs* 
  sTarget = new G4Tubs("Target_s",                                               
                        0, fTargetRadius, 0.5*fTargetLength, 0.,CLHEP::twopi);   

  fLogicTarget = new G4LogicalVolume(sTarget,                                    
                                      fTargetMater,                              
                                      "Target_l");           

  fPhysiTarget = new G4PVPlacement(0,                        
                            				G4ThreeVector(),         
                                    fLogicTarget ,           
                                    "Target_p",              
                                    fLogicLArcontainer,      
                                    false,                   
                                    0);                      
  
  // neutron collimator
  G4Tubs*
  sCollimatorSolid = new G4Tubs("CollimatorSolid",                        
                   0, fCollimatorSolidRadius, fCollimatorSolidLength/2, 0.,CLHEP::twopi ); 
  G4Tubs*
  sCollimatorHollow = new G4Tubs("CollimatorHollow",                        
                   0, fCollimatorHollowRadius, fCollimatorHollowLength/2, 0.,CLHEP::twopi );  
  G4ThreeVector zTransCollimator(0, 0, fCollimatorShieldThickness/2); 
  G4SubtractionSolid* sCollimator =
  new G4SubtractionSolid("CollimatorSolid-CollimatorHollow", sCollimatorSolid, sCollimatorHollow, 0, zTransCollimator); 
  
  fLogicCollimator = new G4LogicalVolume(sCollimator,                                    
                                      fCollimatorMater,                              
                                      "sCollimator_l");           

  fPhysiCollimator = new G4PVPlacement(0,                        
                            				G4ThreeVector(0., 0., fDetectorPositionZ - fDetectorLength/2 - fCollimatorShieldThickness + fCollimatorSolidLength/2),         
                                    fLogicCollimator ,           
                                    "sCollimator_p",              
                                    fLWorld,      
                                    false,                   
                                    0);                      
  
  // neutron detector 
  G4Tubs* 
  sDetector = new G4Tubs("Detector_s",                                               
                        0, fDetectorRadius, 0.5*fDetectorLength, 0.,CLHEP::twopi);   

  fLogicDetector = new G4LogicalVolume(sDetector,                                    
                                      fDetectorMater,                              
                                      "Detector_l");           

  fPhysiDetector = new G4PVPlacement(0,                        
                            				G4ThreeVector(0., 0., fDetectorPositionZ),         
                                    fLogicDetector ,           
                                    "Detector_p",              
                                    fLWorld,      
                                    false,                   
                                    0);              
  
  // set VisAttributes
  fLWorld->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  G4VisAttributes* VisAttGascontainer= new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  fLogicGascontainer->SetVisAttributes(VisAttGascontainer);  
  G4VisAttributes* VisAttGasinsulator= new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  fLogicGasinsulator->SetVisAttributes(VisAttGasinsulator);  
  G4VisAttributes* VisAttLArcontainer= new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  fLogicLArcontainer->SetVisAttributes(VisAttLArcontainer);  
  G4VisAttributes* VisAttTarget= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  fLogicTarget->SetVisAttributes(VisAttTarget);   
  G4VisAttributes* VisAttDetector= new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  fLogicDetector->SetVisAttributes(VisAttDetector);                 
  
    
  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fWorldMater->GetName() 
         << "\n \n" << fWorldMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fWorldMater != pttoMaterial) {
      fWorldMater = pttoMaterial;
      if(fLWorld) { fLWorld->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fWorldSizeZ = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

