#ifndef G4SteppingAction_h
#define G4SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4EventAction;
class G4DetectorConstruction;
class G4EventAction;
class G4LogicalVolume;

/// Stepping action class
///

class G4SteppingAction : public G4UserSteppingAction
{
  public:
    G4SteppingAction(const G4DetectorConstruction* detectorConstruction, G4EventAction* eventAction);
    virtual ~G4SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    const G4DetectorConstruction* fDetConstruction;
    G4EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
