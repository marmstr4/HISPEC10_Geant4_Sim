
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
  analysisManager->CreateNtupleDColumn("Ge_x");
  analysisManager->CreateNtupleDColumn("Ge_y");
  analysisManager->CreateNtupleDColumn("Ge_z");
  analysisManager->CreateNtupleDColumn("Ge_timing");
  analysisManager->CreateNtupleDColumn("Ge_energy");
  analysisManager->CreateNtupleDColumn("Ge_energy_doppler");
  analysisManager->CreateNtupleDColumn("beta_TOF");
  analysisManager->CreateNtupleDColumn("angle");
  analysisManager->CreateNtupleDColumn("mass");
  analysisManager->CreateNtupleDColumn("charge");
  analysisManager->CreateNtupleDColumn("DSSSD2_x");
  analysisManager->CreateNtupleDColumn("DSSSD2_y");
  analysisManager->CreateNtupleDColumn("DSSSD2_z");
  analysisManager->CreateNtupleDColumn("coulex_energy");
  analysisManager->CreateNtupleDColumn("shifted_energy");
  analysisManager->CreateNtupleDColumn("perfect_energy");
  analysisManager->CreateNtupleDColumn("coulex_beta");
  analysisManager->CreateNtupleDColumn("coulex_angle");
  analysisManager->CreateNtupleDColumn("beta_DSSSD");
  analysisManager->CreateNtupleDColumn("ftar_x");
  analysisManager->CreateNtupleDColumn("ftar_y");
  analysisManager->CreateNtupleDColumn("ftar_z");
  analysisManager->CreateNtupleDColumn("ftar_energy");
  analysisManager->CreateNtupleDColumn("beta_tar");
  analysisManager->CreateNtupleDColumn("angle_diff");
  analysisManager->CreateNtupleDColumn("beta_diff");
  analysisManager->CreateNtupleDColumn("factor_diff");
  analysisManager->CreateNtupleDColumn("MCPvsreal_tar_x_diff");// mid target position vs MCP reconstruction
  analysisManager->CreateNtupleDColumn("MCPvsreal_tar_y_diff");
  analysisManager->CreateNtupleDColumn("MCPvsreal_tar_z_diff");
  analysisManager->CreateNtupleDColumn("coulexvsMCP_tar_x_diff"); //MCP reconstrution vs real coulex position
  analysisManager->CreateNtupleDColumn("coulexvsMCP_tar_y_diff");
  analysisManager->CreateNtupleDColumn("coulexvsMCP_tar_z_diff");
  analysisManager->CreateNtupleDColumn("coulex_ion_direction_vs_ion_scatter_dir");//real coulex direction vs dir between ge and real position
  analysisManager->CreateNtupleDColumn("coulex_direction_vs_reco_gamma_dir");//dir between ge and tar vs dir between ge and real position
  analysisManager->CreateNtupleDColumn("reco_coulex_angle_vs_read_out");
  analysisManager->CreateNtupleDColumn("scattering");
  analysisManager->CreateNtupleDColumn("scattering_lab");
  analysisManager->CreateNtupleDColumn("scattering_energy");
  analysisManager->CreateNtupleDColumn("polar_lab");
  analysisManager->CreateNtupleDColumn("polar_com");
  analysisManager->CreateNtupleDColumn("lab_ion_angle");
  analysisManager->CreateNtupleDColumn("target_energy");
  analysisManager->CreateNtupleDColumn("safe_coulex");
  analysisManager->CreateNtupleDColumn("crystal");
  analysisManager->CreateNtupleDColumn("dead_crystal");
  analysisManager->CreateNtupleDColumn("deg_x");
  analysisManager->CreateNtupleDColumn("deg_y");
  analysisManager->CreateNtupleDColumn("deg_z");
  analysisManager->CreateNtupleDColumn("deg_energy");
  analysisManager->CreateNtupleDColumn("deg_mass");
  analysisManager->CreateNtupleDColumn("deg_charge");
  analysisManager->CreateNtupleDColumn("z1_x");
  analysisManager->CreateNtupleDColumn("z1_y");
  analysisManager->CreateNtupleDColumn("z1_z");
  analysisManager->CreateNtupleDColumn("z1_energy");
  analysisManager->CreateNtupleDColumn("z1_mass");
  analysisManager->CreateNtupleDColumn("z1_charge");
  analysisManager->CreateNtupleDColumn("z2_x");
  analysisManager->CreateNtupleDColumn("z2_y");
  analysisManager->CreateNtupleDColumn("z2_z");
  analysisManager->CreateNtupleDColumn("z2_energy");
  analysisManager->CreateNtupleDColumn("z2_mass");
  analysisManager->CreateNtupleDColumn("z2_charge");
  analysisManager->CreateNtupleDColumn("z3_x");
  analysisManager->CreateNtupleDColumn("z3_y");
  analysisManager->CreateNtupleDColumn("z3_z");
  analysisManager->CreateNtupleDColumn("z3_energy");
  analysisManager->CreateNtupleDColumn("z3_mass");
  analysisManager->CreateNtupleDColumn("z3_charge");
  analysisManager->CreateNtupleDColumn("z4_x");
  analysisManager->CreateNtupleDColumn("z4_y");
  analysisManager->CreateNtupleDColumn("z4_z");
  analysisManager->CreateNtupleDColumn("z4_energy");
  analysisManager->CreateNtupleDColumn("z4_mass");
  analysisManager->CreateNtupleDColumn("z4_charge");
  analysisManager->CreateNtupleDColumn("z5_x");
  analysisManager->CreateNtupleDColumn("z5_y");
  analysisManager->CreateNtupleDColumn("z5_z");
  analysisManager->CreateNtupleDColumn("z5_energy");
  analysisManager->CreateNtupleDColumn("z5_mass");
  analysisManager->CreateNtupleDColumn("z5_charge");

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
