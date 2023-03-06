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

    G4RotationMatrix* inclination = new G4RotationMatrix();
    G4double MCP_Sep = 150*CLHEP::cm;
    G4double foilsizeXY = 19*CLHEP::cm;
    G4double foil_pos = -9.9*CLHEP::cm;
    G4double film_pos = -10*CLHEP::cm;

    G4LogicalVolume* logicWorld;

    const G4LogicalVolume* GetMCP1LV() const;
    const G4LogicalVolume* GetMCP2LV() const;
    const G4LogicalVolume* GetDSSSD1LV() const;
    const G4LogicalVolume* GetDSSSD2LV() const;

    const G4LogicalVolume* GetgeLV() const;
    const G4LogicalVolume* GettarLV() const;
    const G4LogicalVolume* GetdegLV() const;

    const G4LogicalVolume* Getz1LV() const;
    const G4LogicalVolume* Getz2LV() const;
    const G4LogicalVolume* Getz3LV() const;
    const G4LogicalVolume* Getz4LV() const;
    const G4LogicalVolume* Getz5LV() const;

    G4LogicalVolume* alu1_log;
    G4LogicalVolume* alu2_log;
    G4LogicalVolume* mylar1_log;
    G4LogicalVolume* mylar2_log;
    G4LogicalVolume* ge_log;
    G4LogicalVolume* tar_log;
    G4LogicalVolume* DSSD1_log;
    G4LogicalVolume* DSSD2_log;
    G4LogicalVolume* deg_log;
    G4LogicalVolume* z1_log;
    G4LogicalVolume* z2_log;
    G4LogicalVolume* z3_log;
    G4LogicalVolume* z4_log;
    G4LogicalVolume* z5_log;

    G4LogicalVolume* fLogicLaBr;
    G4LogicalVolume* fLogic_for_color;

    G4VPhysicalVolume* fPhysiLaBr;
    G4VPhysicalVolume* fPhysiGe;
    G4VPhysicalVolume* physWorld;
    std::vector<G4VPhysicalVolume*> fLogicGe;

    virtual void ConstructDEGAS();
    //virtual void ConstructDEGASAll();
    virtual void ConstructHISPEC();
    virtual void ConstructCubicGe();
    virtual void ConstructWires();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    G4Material* vacuum;
    G4Material* gold_mat;
    G4Material* alu_mat;
    G4Material* mylar_mat;
    G4Material* wire_mat;
    G4Material* DSSD_mat;
    G4Material* ge_mat;
    G4Material* plastic_mat;
    G4Material* steel_mat;
    G4Material* lanthanum_mat;
    G4Material* lead_mat;

    G4bool checkOverlaps = true;

    G4LogicalVolume* chamber_log;
    G4LogicalVolume* chamber2_log;
    G4LogicalVolume* flange_log;
    G4LogicalVolume* flange2_log;

  private:

    G4LogicalVolume* fMCP1;
    G4LogicalVolume* fMCP2;
    G4LogicalVolume* fDSSSD1;
    G4LogicalVolume* fDSSSD2;
    G4LogicalVolume* fscint;
    G4LogicalVolume* fge;
    //G4LogicalVolume* fafter_deg;
    G4LogicalVolume* ftar_log;
    G4LogicalVolume* fdeg;
    G4LogicalVolume* fz1_log;
    G4LogicalVolume* fz2_log;
    G4LogicalVolume* fz3_log;
    G4LogicalVolume* fz4_log;
    G4LogicalVolume* fz5_log;



  protected:
    G4LogicalVolume*  fScoringVolume;
};

//inline const G4LogicalVolume* G4DetectorConstruction::Getafter_degLV() const {
  //return fafter_deg;
//}

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

inline const G4LogicalVolume* G4DetectorConstruction::GetdegLV() const {
  return fdeg;
}

inline const G4LogicalVolume* G4DetectorConstruction::GetgeLV() const {
  return fge;
}

inline const G4LogicalVolume* G4DetectorConstruction::GettarLV() const {
  return ftar_log;
}

inline const G4LogicalVolume* G4DetectorConstruction::Getz1LV() const {
  return fz1_log;
}

inline const G4LogicalVolume* G4DetectorConstruction::Getz2LV() const {
  return fz2_log;
}

inline const G4LogicalVolume* G4DetectorConstruction::Getz3LV() const {
  return fz3_log;
}

inline const G4LogicalVolume* G4DetectorConstruction::Getz4LV() const {
  return fz4_log;
}

inline const G4LogicalVolume* G4DetectorConstruction::Getz5LV() const {
  return fz5_log;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
