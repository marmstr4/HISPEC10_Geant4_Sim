#include "G4RunAction.hh"
#include "G4PrimaryGeneratorAction.hh"
#include "G4DetectorConstruction.hh"
// #include "G4Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RunAction::G4RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{
  // add new units for dose
  //
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);

  // Register accumulable to the accumulable manager
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->RegisterAccumulable(fEdep);

  G4RunManager::GetRunManager()->SetPrintProgress(1);
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  analysisManager->SetVerboseLevel(1);
  //analysisManager->SetNtupleMerging(true);
  analysisManager->SetFirstHistoId(1);

  analysisManager->CreateNtuple("HISPEC", "HISPEC RECONSTRUCTION");
  analysisManager->CreateNtupleDColumn("MCP1_x");
  analysisManager->CreateNtupleDColumn("MCP1_y");
  analysisManager->CreateNtupleDColumn("MCP1_z");
  analysisManager->CreateNtupleDColumn("MCP1_timing");
  analysisManager->CreateNtupleDColumn("MCP2_x");
  analysisManager->CreateNtupleDColumn("MCP2_y");
  analysisManager->CreateNtupleDColumn("MCP2_z");
  analysisManager->CreateNtupleDColumn("MCP2_timing");
  analysisManager->CreateNtupleDColumn("MCP_TOF");
  analysisManager->CreateNtupleDColumn("DSSSD1_energy");
  analysisManager->CreateNtupleDColumn("DSSSD1_timing");
  analysisManager->CreateNtupleDColumn("DSSSD2_energy");
  analysisManager->CreateNtupleDColumn("DSSSD2_timing");
  analysisManager->CreateNtupleDColumn("Scintillator_timing");
  analysisManager->CreateNtupleDColumn("Scintillator_energy");
  analysisManager->CreateNtupleDColumn("Ge_x");
  analysisManager->CreateNtupleDColumn("Ge_y");
  analysisManager->CreateNtupleDColumn("Ge_z");
  analysisManager->CreateNtupleDColumn("Ge_timing");
  analysisManager->CreateNtupleDColumn("Ge_energy");
  analysisManager->CreateNtupleDColumn("Ge_energy_doppler");
  analysisManager->CreateNtupleDColumn("Theta");
  analysisManager->CreateNtupleDColumn("Beta");
  analysisManager->CreateNtupleDColumn("gamma_angle");
  analysisManager->CreateNtupleDColumn("ion_angle");

  analysisManager->FinishNtuple();

  //analysisManager->CreateH1("MCP_TOF","MCP_TOF",1000,0,100*ns);
  //analysisManager->CreateH1("MCP_Velocity","MCP_Velocity",100,-1,1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RunAction::~G4RunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RunAction::BeginOfRunAction(const G4Run*)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Reset();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4String fileName = "Edep";
  analysisManager->OpenFile(fileName);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RunAction::EndOfRunAction(const G4Run* run)
{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;


  const G4DetectorConstruction* detectorConstruction
   = static_cast<const G4DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());


  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const G4PrimaryGeneratorAction* generatorAction
   = static_cast<const G4PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
}
