#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include "G4DetectorConstruction.hh"
#include "G4ActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VMultipleScattering.hh"
#include "G4DecayPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4StoppingPhysics.hh"
#include "PhysicsList.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UserLimits.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLivermorePhysics.hh"
int main()
{

  // construct the default run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();

  // set mandatory initialization classes
  // Detector construction
  G4DetectorConstruction* detConstruction = new G4DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  PhysicsList* physics = new PhysicsList;




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  G4int ver = 1;
  //G4VModularPhysicsList *physics = new G4VModularPhysicsList();

  //auto *physics = new QGSP_BERT;

  //physics->RegisterPhysics( new G4EmStandardPhysicsWVI(ver) );
  //physics->RegisterPhysics( new G4EmStandardPhysicsWVI(ver) );
  //physics->RegisterPhysics( new G4EmStandardPhysics_option3(ver) );
  //physics->RegisterPhysics(new QGSP_BERT);
  //physics->RegisterPhysics(new G4EmLivermorePhysics());

  //physics->RegisterPhysics( new G4EmStandardPhysics_option4(ver) );
  /*
  physics->RegisterPhysics( new G4EmExtraPhysics(ver) );
  physics->RegisterPhysics( new G4DecayPhysics(ver) );
  physics->RegisterPhysics( new G4HadronElasticPhysics(ver) );
  physics->RegisterPhysics( new G4HadronPhysicsQGSP_BERT(ver));
  physics->RegisterPhysics( new G4StoppingPhysics(ver));
  physics->RegisterPhysics( new G4IonPhysics(ver));
  physics->RegisterPhysics(new G4IonBinaryCascadePhysics());
  physics->RegisterPhysics( new G4NeutronTrackingCut(ver));
  physics->RegisterPhysics( new G4StepLimiterPhysics());
  */



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //G4VUserPhysicsList *phys = new PhysicsList();





  //runManager->SetUserInitialization(phys);
  runManager->SetUserInitialization(physics);
  // User action initialization
  G4ActionInitialization* actionInitialization
     = new G4ActionInitialization(detConstruction);
  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();


  //G4UImanager* UImanager = G4UImanager::GetUIpointer();
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // initialize G4 kernel
  runManager->Initialize();


  G4UIsession* ui = new G4UIterminal;
  UImanager->ApplyCommand("/control/execute init_vis.mac");
  ui->SessionStart();
  delete ui;

  // job termination
  delete visManager;
  delete runManager;

}
