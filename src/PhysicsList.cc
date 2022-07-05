#include "PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmProcessOptions.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

//ions
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"
#include "G4hMultipleScattering.hh"
#include "G4ICRU49NuclearStoppingModel.hh"
//#include "Coulomb.hh"
#include "G4NuclideTable.hh"
#include "G4NuclearLevelData.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4Radioactivation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4VAtomDeexcitation.hh"
#include "globals.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Decay.hh"
#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"
#include "G4NuclideTable.hh"
#include "G4RadioactiveDecay.hh"

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
#include "G4InelasticCoulombScattering.hh"


PhysicsList::PhysicsList() : G4VModularPhysicsList()
{

  fEmName = G4String("emstandard_opt4");
  fEmPhysicsList = new G4EmStandardPhysics_option4();
  G4LossTableManager::Instance();
  defaultCutValue = 1.*CLHEP::km;
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(10*eV, 1*GeV);
  SetDefaultCutValue(1*mm);

  G4DeexPrecoParameters* deex =
    G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetCorrelatedGamma(false);
  deex->SetStoreAllLevels(true);
  deex->SetIsomerProduction(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                /std::log(2.));

}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();

}

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  fEmPhysicsList->ConstructProcess();

  AddReaction();

  AddDecay();

  AddRadioactiveDecay();

  AddCoulomb();

}


void PhysicsList::AddCoulomb()
{
  //Coulomb* Coulex = new Coulomb();

  G4InelasticCoulombScattering* Coulex = new G4InelasticCoulombScattering();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "GenericIon") {
    pmanager->AddDiscreteProcess(Coulex);
  }
  }

}



void PhysicsList::AddReaction()
{
   //eminQMD  = 4999.*GeV;
   //eminQMD  = 10.*MeV;
   eminQMD  = 4999.*GeV;
  emaxQMD  = 5000.*GeV;
   overlap  = 10*MeV;
//  SetPhysicsType(bIons);
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  param->SetDeexChannelsType(fCombined);

  G4HadronicInteraction* p =
  G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4PreCompoundModel* thePreCompound = static_cast<G4PreCompoundModel*>(p);
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel; }

  G4BinaryLightIonReaction* theIonBC = new G4BinaryLightIonReaction(thePreCompound);
  theIonBC->SetMaxEnergy(eminQMD + overlap);

  G4double emax = G4HadronicParameters::Instance()->GetMaxEnergy();
  emaxQMD = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
  G4HadronicInteraction* theFTFP = nullptr;
  if(emax > emaxQMD) {
    G4FTFBuilder theFTFPBuilder("FTFP",thePreCompound);
    theFTFP = theFTFPBuilder.GetModel();
    theFTFP->SetMinEnergy(emaxQMD - overlap);
    theFTFP->SetMaxEnergy(emax);
  }

  G4QMDReaction* theQMD = new G4QMDReaction();
  theQMD->SetMinEnergy(eminQMD);
  theQMD->SetMaxEnergy(emaxQMD);

  G4VCrossSectionDataSet* theNuclNuclData =
    new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc() );

  AddProcess("dInelastic", G4Deuteron::Deuteron(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("tInelastic", G4Triton::Triton(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("He3Inelastic", G4He3::He3(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("alphaInelastic", G4Alpha::Alpha(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("ionInelastic", G4GenericIon::GenericIon(), theIonBC, theQMD, theFTFP, theNuclNuclData);

  }

  void PhysicsList::AddProcess(const G4String& name,
  			         G4ParticleDefinition* p,
  				 G4BinaryLightIonReaction* BIC,
  				 G4QMDReaction* QMD,
  				 G4HadronicInteraction* FTFP,
  				 G4VCrossSectionDataSet* theNuclNuclData)
  {

    G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
    G4ProcessManager* pManager = p->GetProcessManager();
    pManager->AddDiscreteProcess(hadi);

    hadi->AddDataSet(theNuclNuclData);

    hadi->RegisterMe(BIC);
    hadi->RegisterMe(QMD);
    if(FTFP) { hadi->RegisterMe(FTFP); }


  }


void PhysicsList::AddDecay()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // Decay Process
  //
  G4Decay* fDecayProcess = new G4Decay();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    if (fDecayProcess->IsApplicable(*particle))
      ph->RegisterProcess(fDecayProcess, particle);
  }
}

void PhysicsList::AddRadioactiveDecay()
{
  G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();

  radioactiveDecay->SetARM(true);                //Atomic Rearangement

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

  // mandatory for G4NuclideTable
  //
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
}
