#ifndef STCyclotronSensitiveFoil_h
#define STCyclotronSensitiveFoil_h 1

#include "G4VSensitiveDetector.hh"
#include "STCyclotronDetectorConstruction.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class STCyclotronRun;

class STCyclotronSensitiveFoil : public G4VSensitiveDetector
{

public:

  STCyclotronSensitiveFoil(const G4String& name,
			   STCyclotronDetectorConstruction* det);
  virtual ~STCyclotronSensitiveFoil();
  
  //methods from base class
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
 
private:
  STCyclotronDetectorConstruction* fDet;
  STCyclotronRun* fRun;
  G4int fTempTrack;
  G4int fTempTrack1;
  G4ThreeVector fTempVector;
  G4double fTempEnergy;

};

#endif
