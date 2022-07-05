#include "G4ActionInitialization.hh"
#include "G4PrimaryGeneratorAction.hh"
#include "G4RunAction.hh"
#include "G4EventAction.hh"
#include "G4SteppingAction.hh"
#include "G4DetectorConstruction.hh"

class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ActionInitialization::G4ActionInitialization(G4DetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
   fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ActionInitialization::~G4ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ActionInitialization::BuildForMaster() const
{
  G4RunAction* runAction = new G4RunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ActionInitialization::Build() const
{
  SetUserAction(new G4PrimaryGeneratorAction);

  G4RunAction* runAction = new G4RunAction;
  SetUserAction(runAction);

  G4EventAction* eventAction = new G4EventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new G4SteppingAction(fDetConstruction,eventAction));
}
