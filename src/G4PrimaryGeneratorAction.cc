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
#include <stdio.h>

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
  G4ParticleDefinition* particle = particleTable->FindParticle(tickle);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,53*cm,-200*mm));
  fParticleGun->SetParticleEnergy(0.001*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  std::ifstream fin("beam_sample.dat");
  std::ifstream fin2("random_eurica0.dat");
  std::string line;
  std::string line2;

  while (std::getline(fin,line)) {
    beam_sample.push_back(line);
  }

  while (std::getline(fin2,line2)) {
  gamma_sample.push_back(line2);
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PrimaryGeneratorAction::~G4PrimaryGeneratorAction()
{
  delete fParticleGun;
}


void G4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{


  std::string event_sample = beam_sample.at(rand()%(41589));

  std::vector<double> properties;
  std::string dummy;

  for (int i = 0; i < event_sample.length(); i++){

    dummy+=event_sample.at(i);

    if (event_sample.at(i)==' ' || i == event_sample.length()-1 ) {
      properties.push_back(std::stod(dummy));
      dummy.clear();
    }

  }

  G4double rand4 = G4UniformRand();

  //G4cout<<rand4<<G4endl;

  //if (rand4>0.951) {

    energy = properties.at(0);
    mass = round(properties.at(1));
    z = round(properties.at(2));
    x = properties.at(3)*10; //convert from cm to mm
    y = properties.at(4)*10; //convert from cm to mm
    A = properties.at(5);
    B = properties.at(6);
    //x = G4RandGauss::shoot(0,9); //convert from cm to mm
    //y = G4RandGauss::shoot(0,6.2); //convert from cm to mm
    //A = G4RandGauss::shoot(0,2*22.1); //convert from z = ; to rad
    //B = G4RandGauss::shoot(0,2*21.8); //convert from z = ; to rad
    /*
    if (mass==0 && z==0 || mass==4 || z==7) {

    energy = 9999;
    mass = 65;
    z = 28;
    x = 0; //convert from cm to mm
    y = 0; //convert from cm to mm
    A = 0.00002; //convert from z = ; to rad
    B = 0.00002; //convert from z = ; to rad

    }
  }

  if (rand4<0.961) {

    mass = 64;
    z = 28;
    energy = mass*G4RandGauss::shoot(4.6,2.263);
    x = G4RandGauss::shoot(0,9); //convert from cm to mm
    y = G4RandGauss::shoot(0,6.2); //convert from cm to mm
    A = G4RandGauss::shoot(0,22.1); //convert from z = ; to rad
    B = G4RandGauss::shoot(0,21.8); //convert from z = ; to rad
  }
  */

  //G4double A = G4RandGauss::shoot(0,31.49);
  //G4double B = G4RandGauss::shoot(0,28.98);

  G4double ux = std::sin(A/1000)*std::cos(B/1000);
  G4double uy = std::sin(A/1000)*std::sin(B/1000);
  G4double uz = std::cos(A/1000);

  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  G4double select = G4UniformRand();
  G4double ionCharge = 0*eplus;
  G4double excitEnergy = 0.*MeV;

  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(z,mass,excitEnergy);
  //G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(42,100,excitEnergy);
  fParticleGun->SetParticleDefinition(ion);
  fParticleGun->SetParticleCharge(ionCharge);

  G4ThreeVector directions(ux,uy,uz);
  //G4ThreeVector directions(0,0,1);
  directions = directions.unit();

  fParticleGun->SetParticleMomentumDirection(directions);

  fParticleGun->SetParticleEnergy(energy);
  //fParticleGun->SetParticleEnergy(1000);
  G4double MCP_Sep = 150*cm;

  //if (G4UniformRand()>0.5) {x=-1*x;}
  //if (G4UniformRand()>0.5) {y=-1*y;}

  fParticleGun->SetParticlePosition(G4ThreeVector(x,y,-159*mm));
  //fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-159*mm));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~gamma_part

  /*

  G4ParticleTable* particleTable2;
  particleTable2 = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* gamma = particleTable2->FindParticle("gamma");

  G4double x = 0*cm;
  G4double y = 0*cm;
  G4double A = 0;
  G4double B = 0;

  G4ThreeVector directions2 = G4ThreeVector(0,0,1);

  fParticleGun->SetParticleDefinition(gamma);
  fParticleGun->SetParticleMomentumDirection(directions2);
  fParticleGun->SetParticleEnergy(1.345*CLHEP::MeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-149*mm));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  */
  /*
  G4ParticleTable* particleTable2;
  particleTable2 = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* gamma = particleTable2->FindParticle("gamma");

  G4int select2 = 1;
  fParticleGun->SetNumberOfParticles(2);
  for (int i = 0; i<select2; i++) {
    G4double gam_energy = std::stod(gamma_sample.at(rand()%(999997)));
    //G4double gam_energy = std::stod(gamma_sample.at(rand()%(999999)));
    //G4double gam_energy = 1.345*MeV;

    G4double theta = CLHEP::twopi*G4UniformRand(), phi = acos(2*G4UniformRand()-1);
    G4double ux2 = std::sin(phi)*std::cos(theta);
    G4double uy2 = std::sin(phi)*std::sin(theta);
    G4double uz2 = std::cos(phi);

    G4ThreeVector directions2 = G4ThreeVector(ux2,uy2,uz2);

    fParticleGun->SetParticleDefinition(gamma);
    fParticleGun->SetParticleMomentumDirection(directions2);
    fParticleGun->SetParticleEnergy(gam_energy);
    fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-159*mm));
    fParticleGun->GeneratePrimaryVertex(anEvent);

  }
  */
}

void select_ion(G4double rand) {



}
