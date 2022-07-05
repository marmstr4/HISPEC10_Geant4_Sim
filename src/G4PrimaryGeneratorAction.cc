#include "G4PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ChargedGeantino.hh"

#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4hMultipleScattering.hh"
#include "G4ProcessManager.hh"
#include "G4IonConstructor.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PrimaryGeneratorAction::G4PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  G4String tickle = "chargedgeantino";

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(tickle);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-20*cm));
  fParticleGun->SetParticleEnergy(250*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PrimaryGeneratorAction::~G4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  G4double ionCharge = 0.*eplus;
  G4double excitEnergy = 0.*MeV;

  G4double select = 100*G4UniformRand();

  //if (select<99.99999999) {
  if (select<99.99999999) {
    G4int mass = 208;
    G4int Z = 82;
        G4double Energy = 7.5*MeV, dEnergy = 1.0*MeV;
        //G4double Energy = 250*MeV, dEnergy = 5*MeV;
        //G4double Energy = 136*MeV, dEnergy = 1*MeV;
        //G4double Energy = 136*MeV, dEnergy = 1*MeV;
        //G4double Energy = 13*MeV, dEnergy = 8*MeV;

    if (particle == G4ChargedGeantino::ChargedGeantino()) {
      G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,mass,excitEnergy);
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(ionCharge);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(mass*Energy,mass*dEnergy));

    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    //G4double dx0 = 11.0*mm, dy0 = 20*mm, dz0 = 0*cm;
    G4double dx0 = 1.0*mm, dy0 = 1.0*mm, dz0 = 0*cm;

    fParticleGun->SetParticlePosition(G4ThreeVector(
      G4RandGauss::shoot(x0,dx0),
      G4RandGauss::shoot(y0,dy0),
      -13*cm));
      //-15*cm));

    };


/*
    if (select>99.0 && select<99.3) {
    G4int mass = 61;
    G4int Z = 27;
    G4double Energy = 16.2*MeV, dEnergy = 7.8*MeV;

  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,mass,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(mass*Energy,mass*dEnergy));
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    G4double dx0 = 10.0*mm, dy0 = 20.0*mm, dz0 = 0*cm;
    fParticleGun->SetParticlePosition(G4ThreeVector(
      G4RandGauss::shoot(x0,dx0),
      G4RandGauss::shoot(y0,dy0),
      -20*cm));
  };

  if (select>99.3 && select<99.5) {
    G4int mass = 63;
    G4int Z = 28;
    G4double Energy = 5.3*MeV, dEnergy = 3.1*MeV;

  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,mass,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(mass*Energy,mass*dEnergy));
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    G4double dx0 = 17*mm, dy0 = 22*mm, dz0 = 0*cm;
    fParticleGun->SetParticlePosition(G4ThreeVector(
      G4RandGauss::shoot(x0,dx0),
      G4RandGauss::shoot(y0,dy0),
      -20*cm));
  };

  if (select>99.5 && select<99.7) {
    G4int mass = 62;
    G4int Z = 27;
    G4double Energy = 21.6*MeV, dEnergy = 8.6*MeV;

    if (particle == G4ChargedGeantino::ChargedGeantino()) {
      G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,mass,excitEnergy);
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(ionCharge);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(mass*Energy,mass*dEnergy));
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    G4double dx0 = 10.6*mm, dy0 = 20*mm, dz0 = 0*cm;
    fParticleGun->SetParticlePosition(G4ThreeVector(
      G4RandGauss::shoot(x0,dx0),
      G4RandGauss::shoot(y0,dy0),
      -20*cm));
  };

  if (select>99.7 && select<99.87) {
    G4int mass = 59;
    G4int Z = 26;
    G4double Energy = 27*MeV, dEnergy = 1.16*MeV;

    if (particle == G4ChargedGeantino::ChargedGeantino()) {
      G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,mass,excitEnergy);
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(ionCharge);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(mass*Energy,mass*dEnergy));
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    G4double dx0 = 10.6*mm, dy0 = 20*mm, dz0 = 0*cm;
    fParticleGun->SetParticlePosition(G4ThreeVector(
      G4RandGauss::shoot(x0,dx0),
      G4RandGauss::shoot(y0,dy0),
      -20*cm));
  };

  if (select>99.87&&select<99.88) {
    G4int mass = 64;
    G4int Z = 26;
    G4double Energy = 31*MeV, dEnergy = 12.6*MeV;


    if (particle == G4ChargedGeantino::ChargedGeantino()) {
      G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,mass,excitEnergy);
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(ionCharge);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(mass*Energy,mass*dEnergy));
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    G4double dx0 = 10.4*mm, dy0 = 20*mm, dz0 = 0*cm;
    fParticleGun->SetParticlePosition(G4ThreeVector(
      G4RandGauss::shoot(x0,dx0),
      G4RandGauss::shoot(y0,dy0),
      -20*cm));
  };
*/
fParticleGun->GeneratePrimaryVertex(anEvent);

}
