#include "G4SteppingAction.hh"
#include "G4EventAction.hh"
#include "G4DetectorConstruction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
//#include "G4Definition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SteppingAction::G4SteppingAction(const G4DetectorConstruction* detectorConstruction, G4EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fDetConstruction(detectorConstruction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SteppingAction::~G4SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingAction::UserSteppingAction(const G4Step* step)
{

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume == fScoringVolume) {
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);
  }

  if (volume == fDetConstruction->GetMCP1LV()) {


    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);

    if (step->GetTrack()->GetParticleDefinition() != Instance) {

      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->AddfMCP1_x(position[0]/CLHEP::mm);
      fEventAction->AddfMCP1_y(position[1]/CLHEP::mm);
      fEventAction->AddfMCP1_z(position[2]/CLHEP::mm);
      fEventAction->AddfMCP1_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);

    }
  }

  if (volume == fDetConstruction->GetMCP2LV()) {

    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);

    if (step->GetTrack()->GetParticleDefinition() != Instance) {

      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->AddfMCP2_x(position[0]/CLHEP::mm);
      fEventAction->AddfMCP2_y(position[1]/CLHEP::mm);
      fEventAction->AddfMCP2_z(position[2]/CLHEP::mm);
      fEventAction->AddfMCP2_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);

    }
  }

  if (volume == fDetConstruction->GetDSSSD1LV()) {

    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    if (step->GetTrack()->GetParticleDefinition() != Instance) {

      fEventAction->AddfDSSSD1_energy(step->GetTotalEnergyDeposit()/CLHEP::MeV);
      fEventAction->AddfDSSSD1_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);

    }
  }

  if (volume == fDetConstruction->GetDSSSD2LV()) {

    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    if (step->GetTrack()->GetParticleDefinition() != Instance) {

      fEventAction->AddfDSSSD2_energy(step->GetTotalEnergyDeposit()/CLHEP::MeV);
      fEventAction->AddfDSSSD2_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);

    }
  }

  if (volume == fDetConstruction->GetscintLV()) {
    fEventAction->Addfscint_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);
    fEventAction->Addfscint_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
  }

  if (volume == fDetConstruction->GetgeLV()) {

    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);

    if (step->GetTrack()->GetParticleDefinition() == Instance) {



      //if (energy>1.8*CLHEP::MeV) {

        G4double energy = step->GetPostStepPoint()->GetKineticEnergy();
        //if ((energy>2.599*CLHEP::MeV)&(energy<2.630*CLHEP::MeV) {
        G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
        fEventAction->Addfge_x(position[0]/CLHEP::mm);
        fEventAction->Addfge_y(position[1]/CLHEP::mm);
        fEventAction->Addfge_z(position[2]/CLHEP::mm);
        fEventAction->Addfge_energy(energy/CLHEP::MeV);
        fEventAction->Addfge_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);
      //}
    }
  }

  if (volume == fDetConstruction->Getafter_degLV()) {

      G4double energy = step->GetPostStepPoint()->GetKineticEnergy();
      fEventAction->Addfafter_deg(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);

  }



/*

  if (volume == fDetConstruction->GettarLV()) {
    G4double dep1 = step->GetPostStepPoint()->GetKineticEnergy();
    G4double dep2 = step->GetTotalEnergyDeposit();
    fEventAction->AddfEnergytar(dep1);
    fEventAction->AddfEnergytar_dep(dep2);
  }

if (volume == fDetConstruction->GettarinLV()) {
    G4double in = step->GetPostStepPoint()->GetKineticEnergy();
    fEventAction->AddfEnergy_tar_in(in);
  }

if (volume == fDetConstruction->GettaroutLV()) {
    G4double out = step->GetPostStepPoint()->GetKineticEnergy();
    fEventAction->AddfEnergy_tar_out(out);
  }

  if (volume == fDetConstruction->GetDSSDinLV()) {

      G4String name = "e-";
      // search in particle table]
      G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* Instance = pTable->FindParticle(name);

      if (step->GetTrack()->GetParticleDefinition() == Instance) {
        G4double in = step->GetPostStepPoint()->GetKineticEnergy();
        fEventAction->AddfEnergy_DSSD_in(in);
      }

    }

  if (volume == fDetConstruction->GetDSSDoutLV()) {
      G4double out = step->GetPostStepPoint()->GetKineticEnergy();
      fEventAction->AddfEnergy_DSSD_out(out);

      G4String name = "gamma";
      // search in particle table]
      G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* Instance = pTable->FindParticle(name);

       //;
      if (step->GetTrack()->GetParticleDefinition() == Instance) {
        fEventAction->GetAllGammaTime(step->GetPreStepPoint()->GetLocalTime());
        fEventAction->GetDegGammaTime(step->GetPreStepPoint()->GetGlobalTime());
      }

    }

    if (volume == fDetConstruction->GetbeamLV()) {

        //G4ParticleDefinition particle = step->GetTrack()->GetParticleDefinition();


          G4double in = step->GetPostStepPoint()->GetKineticEnergy();
          fEventAction->AddfEnergybeam(in);

        //G4double out = step->GetPostStepPoint()->GetKineticEnergy();
        //fEventAction->AddfEnergybeam(out);


        G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();

        G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();

        G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();



        fEventAction->Addfbeam_pos_x(pos[0]);
        fEventAction->Addfbeam_pos_y(pos[1]);
        fEventAction->Addfbeam_mass(mass);
        fEventAction->Addfbeam_charge(charge);
      }

      if (volume == fDetConstruction->GetdegbeamLV()) {
          G4double out = step->GetPostStepPoint()->GetKineticEnergy();
          fEventAction->AddfEnergydeg_beam(out);
        }

*/

}
