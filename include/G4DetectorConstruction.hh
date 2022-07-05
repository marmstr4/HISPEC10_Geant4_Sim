#ifndef G4DetectorConstruction_h
#define G4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "detector.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class G4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    G4DetectorConstruction();
    virtual ~G4DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();


    const G4LogicalVolume* GetMCP1LV() const;
    const G4LogicalVolume* GetMCP2LV() const;
    const G4LogicalVolume* GetDSSSD1LV() const;
    const G4LogicalVolume* GetDSSSD2LV() const;
    const G4LogicalVolume* GetscintLV() const;
    const G4LogicalVolume* GetgeLV() const;
    const G4LogicalVolume* Getafter_degLV() const;

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  private:

    virtual void ConstructSDandField();
    G4LogicalVolume* fMCP1;
    G4LogicalVolume* fMCP2;
    G4LogicalVolume* fDSSSD1;
    G4LogicalVolume* fDSSSD2;
    G4LogicalVolume* fscint;
    G4LogicalVolume* fge;
    G4LogicalVolume* fafter_deg;

  protected:
    G4LogicalVolume*  fScoringVolume;
};

inline const G4LogicalVolume* G4DetectorConstruction::Getafter_degLV() const {
  return fafter_deg;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetMCP1LV() const {
  return fMCP1;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetMCP2LV() const {
  return fMCP2;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetDSSSD1LV() const {
  return fDSSSD1;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetDSSSD2LV() const {
  return fDSSSD2;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetscintLV() const {
  return fscint;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetgeLV() const {
  return fge;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
