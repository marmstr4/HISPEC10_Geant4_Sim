#include "detector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *R0hist)
{
  G4Track *track = aStep->GetTrack();

  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  G4ThreeVector posParticle = preStepPoint->GetPosition();


}
