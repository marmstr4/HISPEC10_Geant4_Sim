#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector>
#include "G4CrossSectionInelastic.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4HadronicInteraction.hh"
#include "G4PreCompoundModel.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4FTFBuilder.hh"
#include "G4QMDReaction.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronInelasticProcess.hh"

class G4VPhysicsConstructor;
class G4ProductionCuts;

class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();

  void ConstructParticle();

  void AddReaction();
  void AddDecay();
  void ConstructProcess();
  void AddRadioactiveDecay();

  void AddCoulomb();

  void AddProcess(const G4String&,
		  G4ParticleDefinition*,
		  G4BinaryLightIonReaction*,
		  G4QMDReaction*,
		  G4HadronicInteraction*,
		  G4VCrossSectionDataSet*);

      G4double eminQMD;
      G4double emaxQMD;
      G4double overlap;

private:

  // hide assignment operator
  PhysicsList & operator=(const PhysicsList &right);
  PhysicsList(const PhysicsList&);
  G4VPhysicsConstructor* fEmPhysicsList;
    G4String               fEmName;
  G4VPhysicsConstructor*  emPhysicsList;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
