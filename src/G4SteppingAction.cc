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

  //G4cout<<volume->GetName()<<G4endl;

  if (volume == fDetConstruction->GetMCP1LV()) {

    //G4cout<<"MCP1"<<G4endl;
    G4String name = "gamma";
    G4String name2 = "e-";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);

    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2) {
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->AddfMCP1_pos(position);
      fEventAction->AddfMCP1_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);
    }
    }
  }

  if (volume == fDetConstruction->GetMCP2LV()) {

    //G4cout<<"MCP2"<<G4endl;
    G4String name = "gamma";
G4String name2 = "e-";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);

    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2) {

      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {

      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->AddfMCP2_pos(position);
      fEventAction->AddfMCP2_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);
    }
    }
  }



  if (volume == fDetConstruction->GetDSSSD1LV()) {
          G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2) {

      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {
      fEventAction->AddfDSSSD1_energy(step->GetTotalEnergyDeposit()/CLHEP::MeV);
      fEventAction->AddfDSSSD1_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);
    }

    }
  }

  if (volume == fDetConstruction->GetDSSSD2LV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {

      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();

      if (mass>5) {
      //if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->AddfDSSSD2_pos(position);
      //G4cout<<"DSSSD2 hit"<<" x "<<position.getX()<<" y "<<position.getY()<<" z "<<position.getZ()<<G4endl;

      fEventAction->Addfmass(mass);
      fEventAction->Addfcharge(charge);
      fEventAction->AddfDSSSD2_energy(step->GetTotalEnergyDeposit()/CLHEP::MeV);
      fEventAction->AddfDSSSD2_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);
    }
    //}
    }
  }


  if (volume == fDetConstruction->GetdegLV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {


      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();

      if (mass>5) {
      if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addfdeg_pos(position);
      fEventAction->Addfdeg_mass(mass);
      fEventAction->Addfdeg_charge(charge);
      fEventAction->Addfdeg_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
      }
    }
    }
  }

  if (volume == fDetConstruction->Getz1LV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {

      //G4cout<<volume<<G4endl;
      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();

      if (mass>5) {
      if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addfz1_pos(position);
      fEventAction->Addfz1_mass(mass);
      fEventAction->Addfz1_charge(charge);
      fEventAction->Addfz1_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
      }
    }
    }
  }


  if (volume == fDetConstruction->Getz2LV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {
      //G4cout<<volume<<G4endl;
      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {
      if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addfz2_pos(position);
      fEventAction->Addfz2_mass(mass);
      fEventAction->Addfz2_charge(charge);
      fEventAction->Addfz2_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
      }
    }
    }
  }

  if (volume == fDetConstruction->Getz3LV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {
      //G4cout<<volume<<G4endl;
      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {
      if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addfz3_pos(position);
      fEventAction->Addfz3_mass(mass);
      fEventAction->Addfz3_charge(charge);
      fEventAction->Addfz3_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
      }
    }
    }
  }

  if (volume == fDetConstruction->Getz4LV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {
      //G4cout<<volume<<G4endl;
      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {
      if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addfz4_pos(position);
      fEventAction->Addfz4_mass(mass);
      fEventAction->Addfz4_charge(charge);
      fEventAction->Addfz4_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
      }
    }
    }
  }

  if (volume == fDetConstruction->Getz5LV()) {

      G4String name2 = "e-";
    G4String name = "gamma";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {
      //G4cout<<volume<<G4endl;
      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {
      if (step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addfz5_pos(position);
      fEventAction->Addfz5_mass(mass);
      fEventAction->Addfz5_charge(charge);
      fEventAction->Addfz5_energy(step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
      }
    }
    }
  }






  G4String name = "gamma";
  G4String name2 = "e-";
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Instance = pTable->FindParticle(name);
  G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);

  G4double howmuch = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;

  //G4cout<<volume->GetName()<<" "<<howmuch<<G4endl;

//  G4cout<<" x "<<step->GetDeltaPosition().getX()<<
//  " y "<<step->GetDeltaPosition().getY()<<" z "<<step->GetDeltaPosition().getZ()<<G4endl;




if (step->GetTrack()->GetParticleDefinition() == Instance || step->GetTrack()->GetParticleDefinition()==Instance2) {


for (G4int i=1; i<10; i++)
{

  G4double E_remain = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;

  sprintf(Ge_crystal_name_a, " Ge_crystal_0%da_log", i);
  sprintf(Ge_crystal_name_tens_a, " Ge_crystal_%da_log", i);

  sprintf(Ge_crystal_name_b, " Ge_crystal_0%db_log", i);
  sprintf(Ge_crystal_name_tens_b, " Ge_crystal_%db_log", i);

  sprintf(Ge_crystal_name_c, " Ge_crystal_0%dc_log", i);
  sprintf(Ge_crystal_name_tens_c, " Ge_crystal_%dc_log", i);

  if( volume->GetName() == Ge_crystal_name_a ||
   volume->GetName() == Ge_crystal_name_tens_a ||
   volume->GetName() == Ge_crystal_name_b ||
    volume->GetName() == Ge_crystal_name_tens_b ||
    volume->GetName() == Ge_crystal_name_c ||
     volume->GetName() == Ge_crystal_name_tens_c ) {

       G4double TotEdep = step->GetTotalEnergyDeposit()/CLHEP::MeV;
       G4double en = 0;



  //G4cout<<" type "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<" parentID "<<step->GetTrack()->GetParentID()<<" trackID "<<step->GetTrack()->GetTrackID()<<" dep "<<TotEdep<<" energy "<<step->GetPreStepPoint()->GetKineticEnergy()<<" local "<<step->GetTrack()->GetLocalTime()<<" global "<<step->GetTrack()->GetGlobalTime()<<G4endl;

  //if( volume->GetName() == Ge_crystal_name_a || volume->GetName() == Ge_crystal_name_tens_a && step->GetTrack()->GetParentID()==1 && step->GetTrack()->GetTrackID() == 1) {
  if( volume->GetName() == Ge_crystal_name_a || volume->GetName() == Ge_crystal_name_tens_a) {
    if (step->GetTrack()->GetParticleDefinition() == Instance && step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
      if (fEventAction->fge_array_pos[i][0]>G4ThreeVector(0,0,0)&&fEventAction->fge_energy[i][0]==0&&fEventAction->fge_energy[i][1]==0&&fEventAction->fge_energy[i][2]==0) {
      fEventAction->fge_array_pos[i][0] = position;
      //G4cout<<" type "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<" parentID "<<step->GetTrack()->GetParentID()<<" trackID "<<step->GetTrack()->GetTrackID()<<" energy "<<TotEdep<<G4endl;
      //G4cout<<" A "<<" x "<<position.getX()<<" y "<<position.getY()<<" z "<<position.getZ()<<G4endl;
      }
    }

      //if (step->GetPreStepPoint()->GetKineticEnergy()>60/1000) {fEventAction->fge_energy[i][0] += TotEdep;}
      fEventAction->fge_energy[i][0] += TotEdep;


     }

  if( volume->GetName() == Ge_crystal_name_b || volume->GetName() == Ge_crystal_name_tens_b)  {
    if (step->GetTrack()->GetParticleDefinition() == Instance && step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
      if (fEventAction->fge_array_pos[i][1]>G4ThreeVector(0,0,0)&&fEventAction->fge_energy[i][1]==0&&fEventAction->fge_energy[i][1]==0&&fEventAction->fge_energy[i][2]==0) {
      fEventAction->fge_array_pos[i][1] = position;
      //G4cout<<" type "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<" parentID "<<step->GetTrack()->GetParentID()<<" trackID "<<step->GetTrack()->GetTrackID()<<" energy "<<TotEdep<<G4endl;
      //G4cout<<" B "<<" x "<<position.getX()<<" y "<<position.getY()<<" z "<<position.getZ()<<G4endl;
      }
    }

      //if (step->GetPreStepPoint()->GetKineticEnergy()>60/1000) {fEventAction->fge_energy[i][0] += TotEdep;}
      fEventAction->fge_energy[i][1] += TotEdep;
     }

  if( volume->GetName() == Ge_crystal_name_c || volume->GetName() == Ge_crystal_name_tens_c) {
    if (step->GetTrack()->GetParticleDefinition() == Instance && step->IsFirstStepInVolume()) {
      G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
      if (fEventAction->fge_array_pos[i][2]>G4ThreeVector(0,0,0)&&fEventAction->fge_energy[i][2]==0&&fEventAction->fge_energy[i][1]==0&&fEventAction->fge_energy[i][2]==0) {
      fEventAction->fge_array_pos[i][2] = position;
      //G4cout<<" C "<<" type "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<" parentID "<<step->GetTrack()->GetParentID()<<" trackID "<<step->GetTrack()->GetTrackID()<<" energy "<<TotEdep<<G4endl;
      //G4cout<<" x "<<position.getX()<<" y "<<position.getY()<<" z "<<position.getZ()<<G4endl;
      }

    }

    //if (step->GetPreStepPoint()->GetKineticEnergy()>60/1000) {fEventAction->fge_energy[i][0] += TotEdep;}
    fEventAction->fge_energy[i][2] += TotEdep;
     }

     //G4cout<<volume->GetName()<<" dep = "<<TotEdep/1.345<<" record = "<<en/1.345<<" remaining = "<<E_remain/1.345<<" "<<step->GetTrack()->GetParentID()<<" "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()<<G4endl;

}

}

}



/*

char crystalNameDEGAS[50];
for(int j = 0; j<5; j++){
for (int j1=0; j1<3; j1++){
  sprintf(crystalNameDEGAS,"det%i_crystal%i",j,j1 );
  //G4cout<<volume->GetName()<<" "<<crystalNameDEGAS<<G4endl;

  if (volume->GetName() == crystalNameDEGAS) {

    G4String name = "gamma";
    G4String name2 = "e-";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);

    if (step->GetTrack()->GetParticleDefinition() == Instance || step->GetTrack()->GetParticleDefinition()==Instance2) {

    //if (step->GetTrack()->GetParentID()==2) {break;}

    G4double TotEdep = step->GetTotalEnergyDeposit()/CLHEP::MeV;
    G4double en = 0;
    G4double E_remain = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;

    fEventAction->fge_energy[j][j1] += TotEdep;

    //if (step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV<0.1*CLHEP::MeV) {

    for(int j = 0; j<5; j++){
      for (int j1=0; j1<3; j1++){
        en = en + fEventAction->fge_energy[j][j1];
      }
    }

    //if (step->GetTrack()->GetParticleDefinition() == Instance && step->IsFirstStepInVolume() ) {
    if (step->GetTrack()->GetParticleDefinition() == Instance && step->IsFirstStepInVolume() && E_remain>0.5) {
      G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
      fEventAction->Addfge_pos(position);
      //G4cout<<j<<" "<<j1<<" "<<volume->GetName()<<" "<<TotEdep*1000<<" "<<en*1000<<" "<<E_remain*1000<<" "<<step->GetTrack()->GetTrackID()<<" gamma"<<" wisely "<<G4endl;
    }


    if (step->GetTrack()->GetParticleDefinition() == Instance) {
      G4cout<<j<<" "<<j1<<" "<<volume->GetName()<<" "<<TotEdep*1000<<" "<<en*1000<<" "<<E_remain*1000<<" "<<step->GetTrack()->GetTrackID()<<" gamma"<<G4endl;
    }

    if (step->GetTrack()->GetParticleDefinition() == Instance2) {
      G4cout<<j<<" "<<j1<<" "<<volume->GetName()<<" "<<TotEdep*1000<<" "<<en*1000<<" "<<E_remain*1000<<" "<<step->GetTrack()->GetTrackID()<<" e-"<<G4endl;
    }


    }
  }

  }
}

*/


  if (volume == fDetConstruction->GetgeLV()) {

    G4String name = "gamma";
    G4String name2 = "e-";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);

      //G4duble TotEdep = step->GetTotalEnergyDeposit()/CLHEP::MeV;
      //fEventAction->fge_energy += TotEdep;
      if (step->GetTrack()->GetParticleDefinition() == Instance || step->GetTrack()->GetParticleDefinition()==Instance2) {
        G4double energy = step->GetPostStepPoint()->GetKineticEnergy();

        //if (energy>1*CLHEP::MeV) {G4cout<<"coulex?";}

        //if ((energy>2.599*CLHEP::MeV)&(energy<2.630*CLHEP::MeV) {
        G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
        //fEventAction->Addfge_pos(position);
        fEventAction->Addfge_energy(energy/CLHEP::MeV);
        fEventAction->Addfge_time(step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns);

      }
  }


  if (volume == fDetConstruction->GettarLV()) {
    G4String name = "gamma";
    G4String name2 = "e-";
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* Instance = pTable->FindParticle(name);
    G4ParticleDefinition* Instance2 = pTable->FindParticle(name2);
    if (step->GetTrack()->GetParticleDefinition() != Instance &&step->GetTrack()->GetParticleDefinition() != Instance2  ) {
      //G4cout<<volume<<G4endl;
      G4double charge = step->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
      G4double mass = step->GetTrack()->GetParticleDefinition()->GetAtomicMass();
      if (mass>5) {

      G4double energy = step->GetPostStepPoint()->GetKineticEnergy();
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      fEventAction->Addftar_pos(position);
      fEventAction->Addftar_energy(energy/CLHEP::MeV);

    }
  }
}
}
