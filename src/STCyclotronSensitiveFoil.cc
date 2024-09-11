#include "STCyclotronRun.hh"
#include "STCyclotronSensitiveFoil.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

STCyclotronSensitiveFoil::STCyclotronSensitiveFoil(const G4String& name,
					           STCyclotronDetectorConstruction* det)
  : G4VSensitiveDetector(name), 
    fDet(det)
{ 
  fTempTrack = 0;
  fTempTrack1 = 0;
  fTempEnergy = 0.;
  fTempVector = G4ThreeVector(0.,0.,0.);  
  fRun =0; 
}

STCyclotronSensitiveFoil::~STCyclotronSensitiveFoil()
{ 
delete fRun;
}

G4bool STCyclotronSensitiveFoil::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  fRun  = static_cast<STCyclotronRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4Track* fTrack = aStep->GetTrack();

  auto analysisManager = G4AnalysisManager::Instance();
  
  //Step/track information
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ThreeVector momentumDirection = aStep->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector vectorPosition = aStep->GetPreStepPoint()->GetPosition();
  G4String name = fTrack->GetDefinition()->GetParticleName();
  
  //Collect general information concerning all of the particles
  fRun->EnergyDepositionFoil(edep);

  //Collect information about protons
  if(name == "proton" || name == "deuteron"){

    if(fTrack->GetTrackID()!=fTempTrack && (momentumDirection.getZ()>0.) &&
       vectorPosition.getX()< 7.5 &&  
       vectorPosition.getX()>-7.5 &&
       vectorPosition.getY()< 7.5 &&  
       vectorPosition.getY()>-7.5){
    
      analysisManager->FillH2(1,vectorPosition.getX(),vectorPosition.getY());
      analysisManager->FillH1(1,energy);
    
      fTempTrack = fTrack->GetTrackID();
    
    }
    
    if(fTempTrack1 == 0){
      fTempTrack1 = fTrack->GetTrackID();
    }

    if(fTrack->GetTrackID()!=fTempTrack1 && (momentumDirection.getZ()>0.) &&
       vectorPosition.getX()< 7.5 &&  
       vectorPosition.getX()>-7.5 &&
       vectorPosition.getY()< 7.5 &&  
       vectorPosition.getY()>-7.5 ){

      analysisManager->FillH2(5,fTempVector.getX(),fTempVector.getY());
      analysisManager->FillH1(3,fTempEnergy);

      fTempTrack1 = fTrack->GetTrackID();

    }
 
    fTempVector =  aStep->GetPostStepPoint()->GetPosition();//vectorPosition;
    fTempEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();//energy;
  }
  }
  
  
