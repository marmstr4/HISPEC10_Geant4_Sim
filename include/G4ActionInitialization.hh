#ifndef G4ActionInitialization_h
#define G4ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class G4DetectorConstruction;

/// Action initialization class.

class G4ActionInitialization : public G4VUserActionInitialization
{
  public:
    G4ActionInitialization(G4DetectorConstruction*);
    virtual ~G4ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    G4DetectorConstruction* fDetConstruction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
