#ifndef APCAD_PhysicsList_h
#define APCAD_PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicsConstructor;
class G4ProductionCuts;

class APCAD_PhysicsList: public G4VModularPhysicsList
{
public:
  APCAD_PhysicsList();
  virtual ~APCAD_PhysicsList();

  void ConstructParticle();

  void ConstructProcess();

private:

  // hide assignment operator
  APCAD_PhysicsList & operator=(const APCAD_PhysicsList &right);
  APCAD_PhysicsList(const APCAD_PhysicsList&);

  G4VPhysicsConstructor*  emPhysicsList;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
