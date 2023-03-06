#include "G4DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4AutoDelete.hh"
#include "G4SDManager.hh"
#include "G4Colour.hh"
#include "G4VSolid.hh"
#include "G4WireSD.hh"
#include "G4PVReplica.hh"
#include <map>
#include <cmath>
#include <string>
#include "G4UserLimits.hh"
#include "G4GDMLParser.hh"
#include "G4UnionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVReplica.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DetectorConstruction::G4DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fMCP1(0),
  fMCP2(0),
  fDSSSD1(0),
  fDSSSD2(0),
  fge(0),
  ftar_log(0),
  fdeg(0),
  fz1_log(0),
  fz2_log(0),
  fz3_log(0),
  fz4_log(0),
  fz5_log(0)

  {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DetectorConstruction::~G4DetectorConstruction()
{
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4DetectorConstruction::Construct()
{
  //Setup Materials for detectors
  G4NistManager* nist = G4NistManager::Instance();
  G4double atomicNumber = 1;
  G4double massofMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 2.e-18*pascal;
  /*
  gold_mat = nist->FindOrBuildMaterial("G4_Au");
  alu_mat = nist->FindOrBuildMaterial("G4_Al");
  mylar_mat = nist->FindOrBuildMaterial("G4_MYLAR");
  wire_mat = nist->FindOrBuildMaterial("G4_W");
  DSSD_mat = nist->FindOrBuildMaterial("G4_Si");
  ge_mat = nist->FindOrBuildMaterial("G4_Ge");
  plastic_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  //steel_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  steel_mat = nist->FindOrBuildMaterial("G4_Al");
  lanthanum_mat = nist->FindOrBuildMaterial("G4_La");
  lead_mat = nist->FindOrBuildMaterial("G4_Pb");
*/
  gold_mat = nist->FindOrBuildMaterial("G4_Au");
  alu_mat = nist->FindOrBuildMaterial("G4_Al");
  mylar_mat = nist->FindOrBuildMaterial("G4_MYLAR");
  wire_mat = nist->FindOrBuildMaterial("G4_W");
  DSSD_mat = nist->FindOrBuildMaterial("G4_Si");
  ge_mat = nist->FindOrBuildMaterial("G4_Ge");
  plastic_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  //steel_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  steel_mat = nist->FindOrBuildMaterial("G4_Fe");
  lanthanum_mat = nist->FindOrBuildMaterial("G4_La");
  lead_mat = nist->FindOrBuildMaterial("G4_Pb");

  //Construct World
  G4double world_sizeXY = 100*cm, world_sizeZ = 500*cm;
  vacuum = new G4Material("vaccum",atomicNumber,massofMole, density,kStateGas,temperature,pressure);
  G4Material* world_mat = vacuum;
  G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");


  ConstructHISPEC();
  //ConstructDEGAS();
  //ConstructWires();

physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World",  0, false,  0,checkOverlaps);
logicWorld->SetVisAttributes(G4VisAttributes::Invisible);


const G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
std::vector<G4VPhysicalVolume*>::const_iterator pvciter;
for( pvciter = pvs->begin(); pvciter != pvs->end(); pvciter++ ){
//      G4cout << (*pvciter)->GetName() << G4endl;
    G4ThreeVector det_pos = (*pvciter)->GetTranslation();
    G4cout<<(*pvciter)->GetLogicalVolume()->GetName()<<" = G4ThreeVector("<<det_pos.getX()<<","<<det_pos.getY()<<","<<det_pos.getZ()<<");"<<G4endl;

  }


return physWorld;

}

void G4DetectorConstruction::ConstructDEGAS() {

  G4GDMLParser fParser;
  fParser.Read("DEGAS_V8_GeAllLaBrAll_SnCeOfFiLaRiSi_1_b_15.gdml");

  G4VPhysicalVolume* fWorldPhysVol1 = fParser.GetWorldVolume();

  G4cout << "Parsing through physical volumes" << G4endl;
  const G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume*>::const_iterator pvciter;
  for( pvciter = pvs->begin(); pvciter != pvs->end(); pvciter++ ) {

    char Ge_crystal_name_a[50];
    char Ge_crystal_name_tens_a[50];
    char Ge_crystal_name_b[50];
    char Ge_crystal_name_tens_b[50];
    char Ge_crystal_name_c[50];
    char Ge_crystal_name_tens_c[50];

    char Ge_crystal_solid[50];
    char Ge_crystal_rot[50];
    char Ge_crystal_trans[50];
    char Ge_crystal_mat[50];
    char Ge_crystal_log[50];
    char Ge_crystal_phys[50];

    char Ge_DZ_name_a[50];
    char Ge_DZ_name_tens_a[50];
    char Ge_DZ_name_b[50];
    char Ge_DZ_name_tens_b[50];
    char Ge_DZ_name_c[50];
    char Ge_DZ_name_tens_c[50];

    char Ge_DZ_solid[50];
    char Ge_DZ_rot[50];
    char Ge_DZ_trans[50];
    char Ge_DZ_mat[50];
    char Ge_DZ_log[50];
    char Ge_DZ_phys[50];

    char Ge_Cap_name_a[50];
    char Ge_Cap_name_tens_a[50];
    char Ge_Cap_name_b[50];
    char Ge_Cap_name_tens_b[50];
    char Ge_Cap_name_c[50];
    char Ge_Cap_name_tens_c[50];

    char Ge_Cap_solid[50];
    char Ge_Cap_rot[50];
    char Ge_Cap_trans[50];
    char Ge_Cap_mat[50];
    char Ge_Cap_log[50];
    char Ge_Cap_phys[50];

    char Ge_Cap_Lid_name_a[50];
    char Ge_Cap_Lid_name_tens_a[50];
    char Ge_Cap_Lid_name_b[50];
    char Ge_Cap_Lid_name_tens_b[50];
    char Ge_Cap_Lid_name_c[50];
    char Ge_Cap_Lid_name_tens_c[50];

    char Ge_Cap_Lid_solid[50];
    char Ge_Cap_Lid_rot[50];
    char Ge_Cap_Lid_trans[50];
    char Ge_Cap_Lid_mat[50];
    char Ge_Cap_Lid_log[50];
    char Ge_Cap_Lid_phys[50];

    char Ge_End_Cap_name_a[50];
    char Ge_End_Cap_name_tens_a[50];

    char Ge_End_Cap_name_b[50];
    char Ge_End_Cap_name_tens_b[50];

    char Ge_End_Cap_name_c[50];
    char Ge_End_Cap_name_tens_c[50];

    char Ge_End_Cap_solid[50];
    char Ge_End_Cap_rot[50];
    char Ge_End_Cap_trans[50];
    char Ge_End_Cap_mat[50];
    char Ge_End_Cap_log[50];
    char Ge_End_Cap_phys[50];

    char Ge_Head_Lid_name_a[50];
    char Ge_Head_Lid_name_tens_a[50];

    char Ge_Head_Lid_name_b[50];
    char Ge_Head_Lid_name_tens_b[50];

    char Ge_Head_Lid_name_c[50];
    char Ge_Head_Lid_name_tens_c[50];

    char Ge_Head_Lid_solid[50];
    char Ge_Head_Lid_rot[50];
    char Ge_Head_Lid_trans[50];
    char Ge_Head_Lid_mat[50];
    char Ge_Head_Lid_log[50];
    char Ge_Head_Lid_phys[50];

    char Ge_Cold_Frame_name_a[50];
    char Ge_Cold_Frame_name_tens_a[50];

    char Ge_Cold_Frame_name_b[50];
    char Ge_Cold_Frame_name_tens_b[50];

    char Ge_Cold_Frame_name_c[50];
    char Ge_Cold_Frame_name_tens_c[50];

    char Ge_Cold_Frame_solid[50];
    char Ge_Cold_Frame_rot[50];
    char Ge_Cold_Frame_trans[50];
    char Ge_Cold_Frame_mat[50];
    char Ge_Cold_Frame_log[50];
    char Ge_Cold_Frame_phys[50];

    for (G4int i=1; i<13; i++)
  	{

    sprintf(Ge_crystal_name_a, "box_Ge_0%da", i);
	  sprintf(Ge_crystal_name_tens_a, "box_Ge_%da", i);

    sprintf(Ge_crystal_name_b, "box_Ge_0%db", i);
	  sprintf(Ge_crystal_name_tens_b, "box_Ge_%db", i);

    sprintf(Ge_crystal_name_c, "box_Ge_0%dc", i);
	  sprintf(Ge_crystal_name_tens_c, "box_Ge_%dc", i);

    sprintf(Ge_DZ_name_a, "box_GeDZ_0%da", i);
	  sprintf(Ge_DZ_name_tens_a, "box_GeDZ_%da", i);

    sprintf(Ge_DZ_name_b, "box_GeDZ_0%db", i);
	  sprintf(Ge_DZ_name_tens_b, "box_GeDZ_%db", i);

    sprintf(Ge_DZ_name_c, "box_GeDZ_0%dc", i);
	  sprintf(Ge_DZ_name_tens_c, "box_GeDZ_%dc", i);

    sprintf(Ge_Cap_name_a, "box_GeCap_0%da", i);
	  sprintf(Ge_Cap_name_tens_a, "box_GeCap_%da", i);

    sprintf(Ge_Cap_name_b, "box_GeCap_0%db", i);
	  sprintf(Ge_Cap_name_tens_b, "box_GeCap_%db", i);

    sprintf(Ge_Cap_name_c, "box_GeCap_0%dc", i);
	  sprintf(Ge_Cap_name_tens_c, "box_GeCap_%dc", i);

    sprintf(Ge_Cap_Lid_name_a, "box_GeCapLid_0%da", i);
	  sprintf(Ge_Cap_Lid_name_tens_a, "box_GeCapLid_%da", i);

    sprintf(Ge_Cap_Lid_name_b, "box_GeCapLid_0%db", i);
	  sprintf(Ge_Cap_Lid_name_tens_b, "box_GeCapLid_%db", i);

    sprintf(Ge_Cap_Lid_name_c, "box_GeCapLid_0%dc", i);
	  sprintf(Ge_Cap_Lid_name_tens_c, "box_GeCapLid_%dc", i);

    sprintf(Ge_End_Cap_name_a, "box_GeEndCap_0%da", i);
	  sprintf(Ge_End_Cap_name_tens_a, "box_GeEndCap_%da", i);

    sprintf(Ge_End_Cap_name_b, "box_GeEndCap_0%db", i);
	  sprintf(Ge_End_Cap_name_tens_b, "box_GeEndCap_%db", i);

    sprintf(Ge_End_Cap_name_c, "box_GeEndCap_0%dc", i);
	  sprintf(Ge_End_Cap_name_tens_c, "box_GeEndCap_%dc", i);

    sprintf(Ge_Head_Lid_name_a, "box_GeColdFrame_0%da", i);
	  sprintf(Ge_Head_Lid_name_tens_a, "box_GeColdFrame_%da", i);

    sprintf(Ge_Head_Lid_name_b, "box_GeColdFrame_0%db", i);
	  sprintf(Ge_Head_Lid_name_tens_b, "box_GeColdFrame_%db", i);

    sprintf(Ge_Head_Lid_name_c, "box_GeColdFrame_0%dc", i);
	  sprintf(Ge_Head_Lid_name_tens_c, "box_GeColdFrame_%dc", i);

    sprintf(Ge_Cold_Frame_name_a, "box_GeHeadLid_0%da", i);
	  sprintf(Ge_Cold_Frame_name_tens_a, "box_GeHeadLid_%da", i);

    sprintf(Ge_Cold_Frame_name_b, "box_GeHeadLid_0%db", i);
	  sprintf(Ge_Cold_Frame_name_tens_b, "box_GeHeadLid_%db", i);

    sprintf(Ge_Cold_Frame_name_c, "box_GeHeadLid_0%dc", i);
	  sprintf(Ge_Cold_Frame_name_tens_c, "box_GeHeadLid_%dc", i);



    if ((*pvciter)->GetName() == Ge_crystal_name_a )  {

    G4VSolid* Ge_crystal_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_crystal_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_crystal_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_crystal_mat = (*pvciter)->GetLogicalVolume()->GetMaterial();

    if (i==1) {
      G4LogicalVolume* Ge_crystal_01a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_01a_log");
      Ge_crystal_01a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_01a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_01a_log, " Ge_crystal_01a_phys", logicWorld, false, 0,false);
    }
    if (i==2) {
      G4LogicalVolume* Ge_crystal_02a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_02a_log");
      Ge_crystal_02a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_02a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_02a_log, " Ge_crystal_02a_phys", logicWorld, false, 0,false);
    }
    if (i==3) {
      G4LogicalVolume* Ge_crystal_03a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_03a_log");
      Ge_crystal_03a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_03a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_03a_log, " Ge_crystal_03a_phys", logicWorld, false, 0,false);
    }
    if (i==4) {
      G4LogicalVolume* Ge_crystal_04a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_04a_log");
      Ge_crystal_04a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_04a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_04a_log, " Ge_crystal_04a_phys", logicWorld, false, 0,false);
    }
    if (i==5) {
      G4LogicalVolume* Ge_crystal_05a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_05a_log");
      Ge_crystal_05a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_05a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_05a_log, " Ge_crystal_05a_phys", logicWorld, false, 0,false);
    }
    if (i==6) {
      G4LogicalVolume* Ge_crystal_06a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_06a_log");
      Ge_crystal_06a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_06a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_06a_log, " Ge_crystal_06a_phys", logicWorld, false, 0,false);
    }
    if (i==7) {
      G4LogicalVolume* Ge_crystal_07a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_07a_log");
      Ge_crystal_07a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_07a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_07a_log, " Ge_crystal_07a_phys", logicWorld, false, 0,false);
    }
    if (i==8) {
      G4LogicalVolume* Ge_crystal_08a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_08a_log");
      Ge_crystal_08a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_08a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_08a_log, " Ge_crystal_08a_phys", logicWorld, false, 0,false);
    }
    if (i==9) {
      G4LogicalVolume* Ge_crystal_09a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_09a_log");
      Ge_crystal_09a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_09a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_09a_log, " Ge_crystal_09a_phys", logicWorld, false, 0,false);
    }

    }

    if ((*pvciter)->GetName() == Ge_crystal_name_tens_a )  {

    G4VSolid* Ge_crystal_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_crystal_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_crystal_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_crystal_mat = (*pvciter)->GetLogicalVolume()->GetMaterial();

    if (i==10) {
      G4LogicalVolume* Ge_crystal_10a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_10a_log");
      Ge_crystal_10a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_10a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_10a_log, " Ge_crystal_10a_phys", logicWorld, false, 0,false);
    }

    if (i==11) {
      G4LogicalVolume* Ge_crystal_11a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_11a_log");
      Ge_crystal_11a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_11a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_11a_log, " Ge_crystal_11a_phys", logicWorld, false, 0,false);
    }

    if (i==12) {
      G4LogicalVolume* Ge_crystal_12a_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_12a_log");
      Ge_crystal_12a_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_12a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_12a_log, " Ge_crystal_12a_phys", logicWorld, false, 0,false);
    }

    }


    if ((*pvciter)->GetName() == Ge_crystal_name_b )  {

    G4VSolid* Ge_crystal_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_crystal_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_crystal_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_crystal_mat = (*pvciter)->GetLogicalVolume()->GetMaterial();

    if (i==1) {
      G4LogicalVolume* Ge_crystal_01b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_01b_log");
      Ge_crystal_01b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_01b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_01b_log, " Ge_crystal_01b_phys", logicWorld, false, 0,false);
    }
    if (i==2) {
      G4LogicalVolume* Ge_crystal_02b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_02b_log");
      Ge_crystal_02b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_02b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_02b_log, " Ge_crystal_02b_phys", logicWorld, false, 0,false);
    }
    if (i==3) {
      G4LogicalVolume* Ge_crystal_03b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_03b_log");
      Ge_crystal_03b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_03b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_03b_log, " Ge_crystal_03b_phys", logicWorld, false, 0,false);
    }
    if (i==4) {
      G4LogicalVolume* Ge_crystal_04b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_04b_log");
      Ge_crystal_04b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_04b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_04b_log, " Ge_crystal_04b_phys", logicWorld, false, 0,false);
    }
    if (i==5) {
      G4LogicalVolume* Ge_crystal_05b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_05b_log");
      Ge_crystal_05b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_05b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_05b_log, " Ge_crystal_05b_phys", logicWorld, false, 0,false);
    }
    if (i==6) {
      G4LogicalVolume* Ge_crystal_06b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_06b_log");
      Ge_crystal_06b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_06b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_06b_log, " Ge_crystal_06b_phys", logicWorld, false, 0,false);
    }
    if (i==7) {
      G4LogicalVolume* Ge_crystal_07b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_07b_log");
      Ge_crystal_07b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_07b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_07b_log, " Ge_crystal_07b_phys", logicWorld, false, 0,false);
    }
    if (i==8) {
      G4LogicalVolume* Ge_crystal_08b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_08b_log");
      Ge_crystal_08b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_08b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_08b_log, " Ge_crystal_08b_phys", logicWorld, false, 0,false);
    }
    if (i==9) {
      G4LogicalVolume* Ge_crystal_09b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_09b_log");
      Ge_crystal_09b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_09b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_09b_log, " Ge_crystal_09b_phys", logicWorld, false, 0,false);
    }

    }

    if ((*pvciter)->GetName() == Ge_crystal_name_tens_b )  {

    G4VSolid* Ge_crystal_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_crystal_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_crystal_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_crystal_mat = (*pvciter)->GetLogicalVolume()->GetMaterial();

    if (i==10) {
      G4LogicalVolume* Ge_crystal_10b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_10b_log");
      Ge_crystal_10b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_10b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_10b_log, " Ge_crystal_10b_phys", logicWorld, false, 0,false);
    }

    if (i==11) {
      G4LogicalVolume* Ge_crystal_11b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_11b_log");
      Ge_crystal_11b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_11b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_11b_log, " Ge_crystal_11b_phys", logicWorld, false, 0,false);
    }

    if (i==12) {
      G4LogicalVolume* Ge_crystal_12b_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_12b_log");
      Ge_crystal_12b_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_12b_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_12b_log, " Ge_crystal_12b_phys", logicWorld, false, 0,false);
    }

    }

    if ((*pvciter)->GetName() == Ge_crystal_name_c )  {

    G4VSolid* Ge_crystal_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_crystal_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_crystal_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_crystal_mat = (*pvciter)->GetLogicalVolume()->GetMaterial();

    if (i==1) {
      G4LogicalVolume* Ge_crystal_01c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_01c_log");
      Ge_crystal_01c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_01c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_01c_log, " Ge_crystal_01c_phys", logicWorld, false, 0,false);
    }
    if (i==2) {
      G4LogicalVolume* Ge_crystal_02c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_02c_log");
      Ge_crystal_02c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_02c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_02c_log, " Ge_crystal_02c_phys", logicWorld, false, 0,false);
    }
    if (i==3) {
      G4LogicalVolume* Ge_crystal_03c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_03c_log");
      Ge_crystal_03c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_03c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_03c_log, " Ge_crystal_03c_phys", logicWorld, false, 0,false);
    }
    if (i==4) {
      G4LogicalVolume* Ge_crystal_04c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_04c_log");
      Ge_crystal_04c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_04c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_04c_log, " Ge_crystal_04c_phys", logicWorld, false, 0,false);
    }
    if (i==5) {
      G4LogicalVolume* Ge_crystal_05c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_05c_log");
      Ge_crystal_05c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_05c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_05c_log, " Ge_crystal_05c_phys", logicWorld, false, 0,false);
    }
    if (i==6) {
      G4LogicalVolume* Ge_crystal_06c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_06c_log");
      Ge_crystal_06c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_06c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_06c_log, " Ge_crystal_06c_phys", logicWorld, false, 0,false);
    }
    if (i==7) {
      G4LogicalVolume* Ge_crystal_07c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_07c_log");
      Ge_crystal_07c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_07c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_07c_log, " Ge_crystal_07c_phys", logicWorld, false, 0,false);
    }
    if (i==8) {
      G4LogicalVolume* Ge_crystal_08c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_08c_log");
      Ge_crystal_08c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_08c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_08c_log, " Ge_crystal_08c_phys", logicWorld, false, 0,false);
    }
    if (i==9) {
      G4LogicalVolume* Ge_crystal_09c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_09c_log");
      Ge_crystal_09c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_09a_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_09c_log, " Ge_crystal_09a_phys", logicWorld, false, 0,false);
    }

  }

    if ((*pvciter)->GetName() == Ge_crystal_name_tens_c )  {

    G4VSolid* Ge_crystal_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_crystal_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_crystal_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_crystal_mat = (*pvciter)->GetLogicalVolume()->GetMaterial();

    if (i==10) {
      G4LogicalVolume* Ge_crystal_10c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_10c_log");
      Ge_crystal_10c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_10c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_10c_log, " Ge_crystal_10c_phys", logicWorld, false, 0,false);
    }

    if (i==11) {
      G4LogicalVolume* Ge_crystal_11c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_11c_log");
      Ge_crystal_11c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_11c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_11c_log, " Ge_crystal_11c_phys", logicWorld, false, 0,false);
    }

    if (i==12) {
      G4LogicalVolume* Ge_crystal_12c_log = new G4LogicalVolume(Ge_crystal_solid, Ge_crystal_mat, " Ge_crystal_12c_log");
      Ge_crystal_12c_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
      G4VPhysicalVolume*  Ge_crystal_12c_phys = new G4PVPlacement(Ge_crystal_rot,Ge_crystal_trans+G4ThreeVector(0,0,MCP_Sep+24*cm), Ge_crystal_12c_log, " Ge_crystal_12c_phys", logicWorld, false, 0,false);
    }

    }



    if ((*pvciter)->GetName() == Ge_DZ_name_a ||
  (*pvciter)->GetName() == Ge_DZ_name_b ||
(*pvciter)->GetName() == Ge_DZ_name_c ||
(*pvciter)->GetName() == Ge_DZ_name_tens_a ||
(*pvciter)->GetName() == Ge_DZ_name_tens_b ||
(*pvciter)->GetName() == Ge_DZ_name_tens_c )  {

    sprintf(Ge_DZ_solid, "box_Ge_0%da_solid", i);
    sprintf(Ge_DZ_rot, "box_Ge_0%da_rot", i);
    sprintf(Ge_DZ_trans, "box_Ge_0%da_trans", i);
    sprintf(Ge_DZ_mat, "box_Ge_0%da_mat", i);
    sprintf(Ge_DZ_log, "box_Ge_0%da__vol", i);
    sprintf(Ge_DZ_phys, "box_Ge_0%da", i);

    G4VSolid* Ge_DZ_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_DZ_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_DZ_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_DZ_mat = vacuum;
    G4LogicalVolume* Ge_DZ_log = new G4LogicalVolume(Ge_DZ_solid, Ge_DZ_mat, "Ge_DZ_log");
    Ge_DZ_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
    G4VPhysicalVolume* Ge_DZ_phys = new G4PVPlacement(Ge_DZ_rot,Ge_DZ_trans+G4ThreeVector(0,0,MCP_Sep+24*cm),Ge_DZ_log, "Ge_DZ_phys", logicWorld, false, 0,false);

    }

    if ((*pvciter)->GetName() == Ge_Cap_name_a ||
  (*pvciter)->GetName() == Ge_Cap_name_b ||
(*pvciter)->GetName() == Ge_Cap_name_c ||
(*pvciter)->GetName() == Ge_Cap_name_tens_a ||
(*pvciter)->GetName() == Ge_Cap_name_tens_b ||
(*pvciter)->GetName() == Ge_Cap_name_tens_c )  {

    sprintf(Ge_Cap_solid, "box_Ge_0%da_solid", i);
    sprintf(Ge_Cap_rot, "box_Ge_0%da_rot", i);
    sprintf(Ge_Cap_trans, "box_Ge_0%da_trans", i);
    sprintf(Ge_Cap_mat, "box_Ge_0%da_mat", i);
    sprintf(Ge_Cap_log, "box_Ge_0%da__vol", i);
    sprintf(Ge_Cap_phys, "box_Ge_0%da", i);

    G4VSolid* Ge_Cap_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_Cap_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_Cap_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_Cap_mat = vacuum;
    G4LogicalVolume* Ge_Cap_log = new G4LogicalVolume(Ge_Cap_solid, Ge_Cap_mat, "Ge_Cap_log");
    Ge_Cap_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
    G4VPhysicalVolume* Ge_Cap_phys = new G4PVPlacement(Ge_Cap_rot,Ge_Cap_trans+G4ThreeVector(0,0,MCP_Sep+24*cm),Ge_Cap_log, "Ge_Cap_phys", logicWorld, false, 0,false);

    }

    if ((*pvciter)->GetName() == Ge_End_Cap_name_a ||
  (*pvciter)->GetName() == Ge_End_Cap_name_b ||
(*pvciter)->GetName() == Ge_End_Cap_name_c ||
(*pvciter)->GetName() == Ge_End_Cap_name_tens_a ||
(*pvciter)->GetName() == Ge_End_Cap_name_tens_b ||
(*pvciter)->GetName() == Ge_End_Cap_name_tens_c )  {

    sprintf(Ge_End_Cap_solid, "box_Ge_0%da_solid", i);
    sprintf(Ge_End_Cap_rot, "box_Ge_0%da_rot", i);
    sprintf(Ge_End_Cap_trans, "box_Ge_0%da_trans", i);
    sprintf(Ge_End_Cap_mat, "box_Ge_0%da_mat", i);
    sprintf(Ge_End_Cap_log, "box_Ge_0%da__vol", i);
    sprintf(Ge_End_Cap_phys, "box_Ge_0%da", i);

    G4VSolid* Ge_End_Cap_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_End_Cap_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_End_Cap_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_End_Cap_mat = vacuum;
    G4LogicalVolume* Ge_End_Cap_log = new G4LogicalVolume(Ge_End_Cap_solid, Ge_End_Cap_mat, "Ge_End_Cap_log");
    Ge_End_Cap_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
    G4VPhysicalVolume* Ge_End_Cap_phys = new G4PVPlacement(Ge_End_Cap_rot,Ge_End_Cap_trans+G4ThreeVector(0,0,MCP_Sep+24*cm),Ge_End_Cap_log, "Ge_End_Cap_phys", logicWorld, false, 0,false);

    }

    if ((*pvciter)->GetName() == Ge_Cap_Lid_name_a ||
  (*pvciter)->GetName() == Ge_Cap_Lid_name_b ||
(*pvciter)->GetName() == Ge_Cap_Lid_name_c ||
(*pvciter)->GetName() == Ge_Cap_Lid_name_tens_a ||
(*pvciter)->GetName() == Ge_Cap_Lid_name_tens_b ||
(*pvciter)->GetName() == Ge_Cap_Lid_name_tens_c )  {

    sprintf(Ge_Cap_Lid_solid, "box_Ge_0%da_solid", i);
    sprintf(Ge_Cap_Lid_rot, "box_Ge_0%da_rot", i);
    sprintf(Ge_Cap_Lid_trans, "box_Ge_0%da_trans", i);
    sprintf(Ge_Cap_Lid_mat, "box_Ge_0%da_mat", i);
    sprintf(Ge_Cap_Lid_log, "box_Ge_0%da__vol", i);
    sprintf(Ge_Cap_Lid_phys, "box_Ge_0%da", i);

    G4VSolid* Ge_Cap_Lid_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_Cap_Lid_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_Cap_Lid_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_Cap_Lid_mat = vacuum;
    G4LogicalVolume* Ge_Cap_Lid_log = new G4LogicalVolume(Ge_Cap_Lid_solid, Ge_Cap_Lid_mat, "Ge_Cap_Lid_log");
    Ge_Cap_Lid_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
    G4VPhysicalVolume* Ge_Cap_Lid_phys = new G4PVPlacement(Ge_Cap_Lid_rot,Ge_Cap_Lid_trans+G4ThreeVector(0,0,MCP_Sep+24*cm),Ge_Cap_Lid_log, "Ge_Cap_Lid_phys", logicWorld, false, 0,false);

    }

    if ((*pvciter)->GetName() == Ge_Cold_Frame_name_a ||
  (*pvciter)->GetName() == Ge_Cold_Frame_name_b ||
(*pvciter)->GetName() == Ge_Cold_Frame_name_c ||
(*pvciter)->GetName() == Ge_Cold_Frame_name_tens_a ||
(*pvciter)->GetName() == Ge_Cold_Frame_name_tens_b ||
(*pvciter)->GetName() == Ge_Cold_Frame_name_tens_c )  {

    sprintf(Ge_Cold_Frame_solid, "box_Ge_0%da_solid", i);
    sprintf(Ge_Cold_Frame_rot, "box_Ge_0%da_rot", i);
    sprintf(Ge_Cold_Frame_trans, "box_Ge_0%da_trans", i);
    sprintf(Ge_Cold_Frame_mat, "box_Ge_0%da_mat", i);
    sprintf(Ge_Cold_Frame_log, "box_Ge_0%da__vol", i);
    sprintf(Ge_Cold_Frame_phys, "box_Ge_0%da", i);

    G4VSolid* Ge_Cold_Frame_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_Cold_Frame_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_Cold_Frame_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_Cold_Frame_mat = vacuum;
    G4LogicalVolume* Ge_Cold_Frame_log = new G4LogicalVolume(Ge_Cold_Frame_solid, Ge_Cold_Frame_mat, "Ge_Cold_Frame_log");
    Ge_Cold_Frame_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
    G4VPhysicalVolume* Ge_Cold_Frame_phys = new G4PVPlacement(Ge_Cold_Frame_rot,Ge_Cold_Frame_trans+G4ThreeVector(0,0,MCP_Sep+24*cm),Ge_Cold_Frame_log, "Ge_Cold_Frame_phys", logicWorld, false, 0,false);

    }

    if ((*pvciter)->GetName() == Ge_Head_Lid_name_a ||
  (*pvciter)->GetName() == Ge_Head_Lid_name_b ||
(*pvciter)->GetName() == Ge_Head_Lid_name_c ||
(*pvciter)->GetName() == Ge_Head_Lid_name_tens_a ||
(*pvciter)->GetName() == Ge_Head_Lid_name_tens_b ||
(*pvciter)->GetName() == Ge_Head_Lid_name_tens_c )  {

    sprintf(Ge_Head_Lid_solid, "box_Ge_0%da_solid", i);
    sprintf(Ge_Head_Lid_rot, "box_Ge_0%da_rot", i);
    sprintf(Ge_Head_Lid_trans, "box_Ge_0%da_trans", i);
    sprintf(Ge_Head_Lid_mat, "box_Ge_0%da_mat", i);
    sprintf(Ge_Head_Lid_log, "box_Ge_0%da__vol", i);
    sprintf(Ge_Head_Lid_phys, "box_Ge_0%da", i);

    G4VSolid* Ge_Head_Lid_solid = (*pvciter)->GetLogicalVolume()->GetSolid();
    G4RotationMatrix* Ge_Head_Lid_rot = (*pvciter)->GetRotation();
    G4ThreeVector Ge_Head_Lid_trans = (*pvciter)->GetTranslation();
    G4Material* Ge_Head_Lid_mat = vacuum;
    G4LogicalVolume* Ge_Head_Lid_log = new G4LogicalVolume(Ge_Head_Lid_solid, Ge_Head_Lid_mat, "Ge_Head_Lid_log");
    Ge_Head_Lid_log->SetVisAttributes(G4Colour(0.0,0.67,1.0));
    G4VPhysicalVolume* Ge_Head_Lid_phys = new G4PVPlacement(Ge_Head_Lid_rot,Ge_Head_Lid_trans+G4ThreeVector(0,0,MCP_Sep+24*cm),Ge_Head_Lid_log, "Ge_Head_Lid_phys", logicWorld, false, 0,false);

    }

  }

  }

}

/*
void G4DetectorConstruction::ConstructDEGAS() {

  G4GDMLParser fParser;
  fParser.Read("DEGAS_V8_GeAllLaBrAll_SnCeOfFiLaRiSi_1_b_15.gdml");

  G4VPhysicalVolume* fWorldPhysVol1 = fParser.GetWorldVolume();
  G4LogicalVolume* degas_log = fWorldPhysVol1->GetLogicalVolume();
  G4VPhysicalVolume* degas_phys = new G4PVPlacement(0,G4ThreeVector(0,0,MCP_Sep+24*cm),degas_log, "degas_phys", logicWorld, false, 0,checkOverlaps);

  G4cout << "Parsing through physical volumes" << G4endl;
  const G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume*>::const_iterator pvciter;
  for( pvciter = pvs->begin(); pvciter != pvs->end(); pvciter++ ) {

      G4cout << (*pvciter)->GetName()<<G4endl;

      char box_LaBr_name_a[50];
      char box_LaBr_name_tens_a[50];

      char box_LaBr_name_b[50];
      char box_LaBr_name_tens_b[50];

      char box_LaBr_name_c[50];
      char box_LaBr_name_tens_c[50];

      char box_LaBrLightGuide_name_a[50];
      char box_LaBrLightGuide_name_tens_a[50];

      char box_LaBrLightGuide_name_b[50];
      char box_LaBrLightGuide_name_tens_b[50];

      char box_LaBrLightGuide_name_c[50];
      char box_LaBrLightGuide_name_tens_c[50];

      char box_LaBrPMT_name_a[50];
      char box_LaBrPMT_name_tens_a[50];

      char box_LaBrPMT_name_b[50];
      char box_LaBrPMT_name_tens_b[50];

      char box_LaBrPMT_name_c[50];
      char box_LaBrPMT_name_tens_c[50];

      char box_LaBrPMTTube_name_a[50];
      char box_LaBrPMTTube_name_tens_a[50];

      char box_LaBrPMTTube_name_b[50];
      char box_LaBrPMTTube_name_tens_b[50];

      char box_LaBrPMTTube_name_c[50];
      char box_LaBrPMTTube_name_tens_c[50];

      char box_LaBrPbShield_name_a[50];
      char box_LaBrPbShield_name_tens_a[50];

      char box_LaBrPbShield_name_b[50];
      char box_LaBrPbShield_name_tens_b[50];

      char box_LaBrPbShield_name_c[50];
      char box_LaBrPbShield_name_tens_c[50];

      char box_LaBrCasing_name_a[50];
      char box_LaBrCasing_name_tens_a[50];

      char box_LaBrCasing_name_b[50];
      char box_LaBrCasing_name_tens_b[50];

      char box_LaBrCasing_name_c[50];
      char box_LaBrCasing_name_tens_c[50];

      char box_LaBrCap_name_a[50];
      char box_LaBrCap_name_tens_a[50];

      char box_LaBrCap_name_b[50];
      char box_LaBrCap_name_tens_b[50];

      char box_LaBrCap_name_c[50];
      char box_LaBrCap_name_tens_c[50];

      char box_LaBrShortRingBack_name_a[50];
      char box_LaBrShortRingBack_name_tens_a[50];

      char box_LaBrShortRingBack_name_b[50];
      char box_LaBrShortRingBack_name_tens_b[50];

      char box_LaBrShortRingBack_name_c[50];
      char box_LaBrShortRingBack_name_tens_c[50];

      char box_LaBrShortRingEnd_name_a[50];
      char box_LaBrShortRingEnd_name_tens_a[50];

      char box_LaBrShortRingEnd_name_b[50];
      char box_LaBrShortRingEnd_name_tens_b[50];

      char box_LaBrShortRingEnd_name_c[50];
      char box_LaBrShortRingEnd_name_tens_c[50];

      char box_LaBrShortRing_name_a[50];
      char box_LaBrShortRing_name_tens_a[50];

      char box_LaBrShortRing_name_b[50];
      char box_LaBrShortRing_name_tens_b[50];

      char box_LaBrShortRing_name_c[50];
      char box_LaBrShortRing_name_tens_c[50];

      char snout_Head[50];
      char snout_Tube[50];
      char AIDA_bplast[50];
      char AIDA_si[50];

    	for (G4int i=0; i<13; i++) {

        sprintf(snout_Head, "snoutHead_%d", i);
        sprintf(snout_Tube, "snoutTube_%d", i);
        sprintf(AIDA_bplast, "AIDA_bPlast_%d", i);
        sprintf(AIDA_si, "AIDA_si_%d", i);

        sprintf(box_LaBrLightGuide_name_a, "box_LaBrLightGuide_0%da", i);
    	  sprintf(box_LaBrLightGuide_name_tens_a, "box_LaBrLightGuide_%da", i);

        sprintf(box_LaBrLightGuide_name_b, "box_LaBrLightGuide_0%db", i);
    	  sprintf(box_LaBrLightGuide_name_tens_b, "box_LaBrLightGuide_%db", i);

        sprintf(box_LaBrLightGuide_name_c, "box_LaBrLightGuide_0%dc", i);
    	  sprintf(box_LaBrLightGuide_name_tens_c, "box_LaBrLightGuide_%dc", i);

        sprintf(box_LaBrPMT_name_a, "box_LaBrPMT_0%da", i);
    	  sprintf(box_LaBrPMT_name_tens_a, "box_LaBrPMT_%da", i);

        sprintf(box_LaBrPMT_name_b, "box_LaBrPMT_0%db", i);
    	  sprintf(box_LaBrPMT_name_tens_b, "box_LaBrPMT_%db", i);

        sprintf(box_LaBrPMT_name_c, "box_LaBrPMT_0%dc", i);
    	  sprintf(box_LaBrPMT_name_tens_c, "box_LaBrPMT_%dc", i);

        sprintf(box_LaBrPMTTube_name_a, "box_LaBrPMTTube_0%da", i);
    	  sprintf(box_LaBrPMTTube_name_tens_a, "box_LaBrPMTTube_%da", i);

        sprintf(box_LaBrPMTTube_name_b, "box_LaBrPMTTube_0%db", i);
    	  sprintf(box_LaBrPMTTube_name_tens_b, "box_LaBrPMTTube_%db", i);

        sprintf(box_LaBrPMTTube_name_c, "box_LaBrPMTTube_0%dc", i);
    	  sprintf(box_LaBrPMTTube_name_tens_c, "box_LaBrPMTTube_%dc", i);

        sprintf(box_LaBrPbShield_name_a, "box_LaBrPbShield_0%da", i);
    	  sprintf(box_LaBrPbShield_name_tens_a, "box_LaBrPbShield_%da", i);

        sprintf(box_LaBrPbShield_name_b, "box_LaBrPbShield_0%db", i);
    	  sprintf(box_LaBrPbShield_name_tens_b, "box_LaBrPbShield_%db", i);

        sprintf(box_LaBrPbShield_name_c, "box_LaBrPbShield_0%dc", i);
    	  sprintf(box_LaBrPbShield_name_tens_c, "box_LaBrPbShield_%dc", i);

        sprintf(box_LaBrCasing_name_a, "box_LaBrCasing_0%da", i);
    	  sprintf(box_LaBrCasing_name_tens_a, "box_LaBrCasing_%da", i);

        sprintf(box_LaBrCasing_name_b, "box_LaBrCasing_0%db", i);
    	  sprintf(box_LaBrCasing_name_tens_b, "box_LaBrCasing_%db", i);

        sprintf(box_LaBrCasing_name_c, "box_LaBrCasing_0%dc", i);
    	  sprintf(box_LaBrCasing_name_tens_c, "box_LaBrCasing_%dc", i);

        sprintf(box_LaBrCap_name_a, "box_LaBrCap_0%da", i);
    	  sprintf(box_LaBrCap_name_tens_a, "box_LaBrCap_%da", i);

        sprintf(box_LaBrCap_name_b, "box_LaBrCap_0%db", i);
    	  sprintf(box_LaBrCap_name_tens_b, "box_LaBrCap_%db", i);

        sprintf(box_LaBrCap_name_c, "box_LaBrCap_0%dc", i);
    	  sprintf(box_LaBrCap_name_tens_c, "box_LaBrCap_%dc", i);

        sprintf(box_LaBr_name_a, "box_LaBr_0%da", i);
    	  sprintf(box_LaBr_name_tens_a, "box_LaBr_%da", i);

        sprintf(box_LaBr_name_b, "box_LaBr_0%db", i);
    	  sprintf(box_LaBr_name_tens_b, "box_LaBr_%db", i);

        sprintf(box_LaBr_name_c, "box_LaBr_0%dc", i);
    	  sprintf(box_LaBr_name_tens_c, "box_LaBr_%dc", i);

        sprintf(box_LaBrShortRingBack_name_a, "box_LaBrShortRingBack_0%da", i);
    	  sprintf(box_LaBrShortRingBack_name_tens_a, "box_LaBrShortRingBack_%da", i);

        sprintf(box_LaBrShortRingBack_name_b, "box_LaBrShortRingBack_0%db", i);
    	  sprintf(box_LaBrShortRingBack_name_tens_b, "box_LaBrShortRingBack_%db", i);

        sprintf(box_LaBrShortRingBack_name_c, "box_LaBrShortRingBack_0%dc", i);
    	  sprintf(box_LaBrShortRingBack_name_tens_c, "box_LaBrShortRingBack_%dc", i);

        sprintf(box_LaBrShortRingEnd_name_a, "box_LaBrShortRingEnd_0%da", i);
    	  sprintf(box_LaBrShortRingEnd_name_tens_a, "box_LaBrShortRingEnd_%da", i);

        sprintf(box_LaBrShortRingEnd_name_b, "box_LaBrShortRingEnd_0%db", i);
    	  sprintf(box_LaBrShortRingEnd_name_tens_b, "box_LaBrShortRingEnd_%db", i);

        sprintf(box_LaBrShortRingEnd_name_c, "box_LaBrShortRingEnd_0%dc", i);
    	  sprintf(box_LaBrShortRingBack_name_tens_c, "box_LaBrShortRingEnd_%dc", i);

        sprintf(box_LaBrShortRing_name_a, "box_LaBrShortRing_0%da", i);
    	  sprintf(box_LaBrShortRing_name_tens_a, "box_LaBrShortRing_%da", i);

        sprintf(box_LaBrShortRing_name_b, "box_LaBrShortRing_0%db", i);
    	  sprintf(box_LaBrShortRing_name_tens_b, "box_LaBrShortRing_%db", i);

        sprintf(box_LaBrShortRing_name_c, "box_LaBrShortRing_0%dc", i);
    	  sprintf(box_LaBrShortRingBack_name_tens_c, "box_LaBrShortRing_%dc", i);

        //G4cout<<(*pvciter)->GetName()<<G4endl;

      if(
      (*pvciter)->GetName() == box_LaBrLightGuide_name_a ||
      (*pvciter)->GetName() == box_LaBrLightGuide_name_b ||
      (*pvciter)->GetName() == box_LaBrLightGuide_name_c ||
      (*pvciter)->GetName() == box_LaBrLightGuide_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrLightGuide_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrLightGuide_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrPMT_name_a ||
      (*pvciter)->GetName() == box_LaBrPMT_name_b ||
      (*pvciter)->GetName() == box_LaBrPMT_name_c ||
      (*pvciter)->GetName() == box_LaBrPMT_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrPMT_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrPMT_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrPbShield_name_a ||
      (*pvciter)->GetName() == box_LaBrPbShield_name_b ||
      (*pvciter)->GetName() == box_LaBrPbShield_name_c ||
      (*pvciter)->GetName() == box_LaBrPbShield_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrPbShield_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrPbShield_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrCasing_name_a ||
      (*pvciter)->GetName() == box_LaBrCasing_name_b ||
      (*pvciter)->GetName() == box_LaBrCasing_name_c ||
      (*pvciter)->GetName() == box_LaBrCasing_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrCasing_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrCasing_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrCap_name_a ||
      (*pvciter)->GetName() == box_LaBrCap_name_b ||
      (*pvciter)->GetName() == box_LaBrCap_name_c ||
      (*pvciter)->GetName() == box_LaBrCap_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrCap_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrCap_name_tens_c ||
      (*pvciter)->GetName() == box_LaBr_name_a ||
      (*pvciter)->GetName() == box_LaBr_name_b ||
      (*pvciter)->GetName() == box_LaBr_name_c ||
      (*pvciter)->GetName() == box_LaBr_name_tens_a ||
      (*pvciter)->GetName() == box_LaBr_name_tens_b ||
      (*pvciter)->GetName() == box_LaBr_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrShortRingBack_name_a ||
      (*pvciter)->GetName() == box_LaBrShortRingBack_name_b ||
      (*pvciter)->GetName() == box_LaBrShortRingBack_name_c ||
      (*pvciter)->GetName() == box_LaBrShortRingBack_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrShortRingBack_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrShortRingBack_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrShortRingEnd_name_a ||
      (*pvciter)->GetName() == box_LaBrShortRingEnd_name_b ||
      (*pvciter)->GetName() == box_LaBrShortRingEnd_name_c ||
      (*pvciter)->GetName() == box_LaBrShortRingEnd_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrShortRingEnd_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_tens_c ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_a ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_b ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_c ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_tens_a ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_tens_b ||
      (*pvciter)->GetName() == box_LaBrShortRing_name_tens_c ||
      (*pvciter)->GetName() == snout_Head ||
      (*pvciter)->GetName() == snout_Tube ||
      (*pvciter)->GetName() == AIDA_bplast ||
      (*pvciter)->GetName() == AIDA_si  ||
      (*pvciter)->GetName() == "snoutFrame" ||
      (*pvciter)->GetName() == "snoutTube" ) {

        G4cout<<" removed "<<(*pvciter)->GetName()<<G4endl;
        //degas_log->GetDaughter((*pvciter))->GetName();
        degas_log->RemoveDaughter((*pvciter));
        //logicWorld->RemoveDaughter((*pvciter));

        }

      }

  }


}

*/
void G4DetectorConstruction::ConstructHISPEC() {


  //Detector positions
  G4ThreeVector deg_pos = G4ThreeVector(0,0,-158*mm);
  G4ThreeVector after_deg_pos = G4ThreeVector(0,0,-14.99*cm);
  G4ThreeVector alu1_pos = G4ThreeVector(0, 0, film_pos);
  G4ThreeVector mylar1_pos = G4ThreeVector(0, 0, foil_pos);
  G4ThreeVector alu2_pos = G4ThreeVector(0, 0, film_pos+MCP_Sep);
  G4ThreeVector mylar2_pos = G4ThreeVector(0, 0, foil_pos+MCP_Sep);
  G4ThreeVector wires1_pos = G4ThreeVector(0, 0, foil_pos-1*mm);
  G4ThreeVector wires2_pos = G4ThreeVector(0, 0, foil_pos+MCP_Sep-1*mm);
  G4ThreeVector DSSD1_pos = G4ThreeVector(0, 0, MCP_Sep+26*cm);
  G4ThreeVector DSSD2_pos = G4ThreeVector(0, 0, MCP_Sep+29*cm);
  G4ThreeVector ge_pos = G4ThreeVector(0,0,MCP_Sep+40*cm);
  G4ThreeVector tar_pos = G4ThreeVector(0,0,MCP_Sep+24*cm);
  G4ThreeVector DSSSD_sphere_pos = G4ThreeVector(0,0,MCP_Sep+24*cm);

  //Detector Geometry
  G4Box* solid_tar = new G4Box("solid_tar", 0.5*8*cm, 0.5*8*cm, 0.5*1*um);
  G4Box* solid_alu1 = new G4Box("solid_alu1", foilsizeXY, foilsizeXY, 0.5*200*nm);
  G4Box* solid_mylar1 = new G4Box("solid_mylar1", foilsizeXY, foilsizeXY, 0.5*10*nm);
  G4Box* solid_DSSD1 = new G4Box("solid_DSSD1", 0.5*8*cm, 0.5*8*cm, 0.5*20*um);
  G4Box* solid_DSSD2 = new G4Box("solid_DSSD2", 0.5*8*cm, 0.5*8*cm, 0.5*300*um);
  G4Sphere* solid_ge = new G4Sphere("solid_ge", 50*cm, 50.1*cm, 0.,twopi, 0., pi);
  G4VSolid* solid_chamber = new G4Tubs("solid_chamber",200*mm,205*mm,(MCP_Sep+24*cm)/2,0,twopi);
  G4VSolid* solid_chamber2 = new G4Tubs("solid_chamber",100*mm,105*mm,(30*cm)/2+2.5*cm,0,twopi);
  G4VSolid* solid_flange = new G4Tubs("solid_flange",0,105*mm,1*cm,0,twopi);
  G4VSolid* solid_flange2 = new G4Tubs("solid_flange2",105*mm,125*mm,1*cm,0,twopi);
  G4VSolid* solid_wall = new G4Tubs("solid_wall",131*mm,50*cm,0.5*20*cm,0,twopi);
  G4VSolid* solid_deg = new G4Box("solid_deg",200*cm,200*cm,5*mm);
  G4VSolid* solid_chamber3 = new G4Box("solid_chamber3",100*cm,100*cm,1*mm);

  //G4Sphere* solid_ge = new G4Sphere("solid_ge", 50*cm, 50.1*cm, 0.,twopi, 0., pi);
  //G4VSolid* DSSSD_sphere = new G4Sphere("DSSSD_sphere",,0,twopi,0,pi);
  G4VSolid* DSSSD_sphere = new G4Sphere("DSSSD_sphere",100*mm-300*um,100*mm,0,twopi,0,pi);
  G4VSolid* DSSSD_cylinder = new G4Tubs("DSSSD_sphere_cylinder",0,45.25*mm,15*cm,0,twopi);
  G4SubtractionSolid* solid_DSSSD_sphere = new G4SubtractionSolid("solid_DSSSD_sphere",DSSSD_sphere,DSSSD_cylinder);

  //Detector Logical volumes
  G4LogicalVolume* tar_log = new G4LogicalVolume(solid_tar,vacuum,"tar_log");
  G4LogicalVolume* alu1_log = new G4LogicalVolume(solid_alu1,vacuum,"alu1_log");
  G4LogicalVolume* mylar1_log = new G4LogicalVolume(solid_mylar1,vacuum,"mylar1_log");
  G4LogicalVolume* alu2_log = new G4LogicalVolume(solid_alu1,vacuum,"alu2_log");
  G4LogicalVolume* mylar2_log = new G4LogicalVolume(solid_mylar1,vacuum,"mylar2_log");
  G4LogicalVolume* DSSD1_log = new G4LogicalVolume(solid_DSSD1,vacuum,"DSSD1_log");
  G4LogicalVolume* DSSD2_log = new G4LogicalVolume(solid_DSSD2,vacuum,"DSSD2_log");
  G4LogicalVolume* ge_log = new G4LogicalVolume(solid_ge,vacuum,"ge_log");
  G4LogicalVolume* deg_log = new G4LogicalVolume(solid_chamber3, vacuum, "deg_log");
  //G4LogicalVolume* wall_log = new G4LogicalVolume(solid_wall,vacuum,"wall_log");
  chamber_log = new G4LogicalVolume(solid_chamber,vacuum , "chamber_log");
  chamber2_log = new G4LogicalVolume(solid_chamber2, vacuum, "chamber_log2");
  flange_log = new G4LogicalVolume(solid_flange, vacuum, "flange_log");
  flange2_log = new G4LogicalVolume(solid_flange2, vacuum, "flange2_log");

  G4LogicalVolume* DSSSD_sphere_log = new G4LogicalVolume(solid_DSSSD_sphere,vacuum,"DSSSD_sphere_log");
  G4LogicalVolume* z1_log = new G4LogicalVolume(solid_chamber3, vacuum, "z1_log");
  G4LogicalVolume* z2_log = new G4LogicalVolume(solid_chamber3, vacuum, "z2_log");
  G4LogicalVolume* z3_log = new G4LogicalVolume(solid_chamber3, vacuum, "z3_log");
  G4LogicalVolume* z4_log = new G4LogicalVolume(solid_chamber3, vacuum, "z4_log");
  G4LogicalVolume* z5_log = new G4LogicalVolume(solid_chamber3, vacuum, "z5_log");
  /*
  //Detector Logical volumes
  G4LogicalVolume* tar_log = new G4LogicalVolume(solid_tar,gold_mat,"tar_log");
  G4LogicalVolume* alu1_log = new G4LogicalVolume(solid_alu1,alu_mat,"alu1_log");
  G4LogicalVolume* mylar1_log = new G4LogicalVolume(solid_mylar1,vacuum,"mylar1_log");
  G4LogicalVolume* alu2_log = new G4LogicalVolume(solid_alu1,alu_mat,"alu2_log");
  G4LogicalVolume* mylar2_log = new G4LogicalVolume(solid_mylar1,vacuum,"mylar2_log");
  G4LogicalVolume* DSSD1_log = new G4LogicalVolume(solid_DSSD1,DSSD_mat,"DSSD1_log");
  G4LogicalVolume* DSSD2_log = new G4LogicalVolume(solid_DSSD2,DSSD_mat,"DSSD2_log");
  G4LogicalVolume* ge_log = new G4LogicalVolume(solid_ge,vacuum,"ge_log");
  G4LogicalVolume* wall_log = new G4LogicalVolume(solid_wall,lead_mat,"wall_log");
  chamber_log = new G4LogicalVolume(solid_chamber, steel_mat, "chamber_log");
  chamber2_log = new G4LogicalVovacuumlume(solid_chamber2, steel_mat, "chamber_log2");
  flange_log = new G4LogicalVolume(solid_flange, steel_mat, "flange_log");
  flange2_log = new G4LogicalVolume(solid_flange2, steel_mat, "flange2_log");
  */

  //Detector Physical volumes
  inclination->rotateX(315.*deg);
  G4VPhysicalVolume* ge_phys = new G4PVPlacement(0,ge_pos,ge_log,"ge_phys",logicWorld,false,0,checkOverlaps);
  G4VPhysicalVolume* tar_phys = new G4PVPlacement(0,tar_pos,tar_log,"tar_phys",logicWorld,false,9,checkOverlaps);
  G4VPhysicalVolume* DSSD1_phys = new G4PVPlacement(0, DSSD1_pos,DSSD1_log, "DSSD1_phys", logicWorld, false,  0, checkOverlaps);
  G4VPhysicalVolume* DSSD2_phys = new G4PVPlacement(0, DSSD2_pos,DSSD2_log, "DSSD2_phys", logicWorld, false,  0, checkOverlaps);
  G4VPhysicalVolume* alu1_phys = new G4PVPlacement(inclination,   alu1_pos,  alu1_log,  "alu1_phys", logicWorld, false, 0, checkOverlaps);        //overlaps checking
  G4VPhysicalVolume* mylar1_phys = new G4PVPlacement(inclination, mylar1_pos, mylar1_log, "mylar1_phys", logicWorld, false, 0,checkOverlaps);        //overlaps checking
  G4VPhysicalVolume* alu2_phys = new G4PVPlacement(inclination,  alu2_pos,   alu2_log, "alu2_phys", logicWorld, false, 0, checkOverlaps);        //overlaps checking
  G4VPhysicalVolume* mylar2_phys = new G4PVPlacement(inclination, mylar2_pos,   mylar2_log,  "mylar2_phys", logicWorld,   false, 0, checkOverlaps);        //overlaps checking
  G4VPhysicalVolume* chamber_phys = new G4PVPlacement(0,G4ThreeVector(0,0,(MCP_Sep)/2-10*cm),chamber_log,"chamber_phys",logicWorld,false,0,checkOverlaps);
  G4VPhysicalVolume* chamber_phys2 = new G4PVPlacement(0,G4ThreeVector(0,0,MCP_Sep+20*cm),chamber2_log,"chamber_phys2",logicWorld,false,0,checkOverlaps);
  G4VPhysicalVolume* flange_phys = new G4PVPlacement(0,G4ThreeVector(0,0,MCP_Sep+20*cm+(30*cm)/2+2.5*cm),flange_log,"flange_phys",logicWorld,false,0,checkOverlaps);
  G4VPhysicalVolume* flange2_phys = new G4PVPlacement(0,G4ThreeVector(0,0,(MCP_Sep)/2-10*cm+(MCP_Sep+24*cm)/2+0.5*cm),flange2_log,"flange2_phys",logicWorld,false,0,checkOverlaps);
  //G4VPhysicalVolume* wall_phys = new G4PVPlacement(0,G4ThreeVector(0,0,(MCP_Sep)/2-10*cm+(MCP_Sep+24*cm)/2+0.5*cm+5*cm),wall_log,"wall_phys",logicWorld,false,0,false);
  G4ThreeVector add = G4ThreeVector(0,0,-1*cm);
  G4VPhysicalVolume* DSSSD_sphere_phys = new G4PVPlacement(0,DSSSD_sphere_pos,DSSSD_sphere_log,"DSSSD_sphere_phys",logicWorld,false,0,false);
  G4VPhysicalVolume* deg_phys = new G4PVPlacement(0,deg_pos,deg_log,"deg_phys", logicWorld, false, 0, false);
  G4VPhysicalVolume* z1_phys = new G4PVPlacement(0,mylar1_pos+add,z1_log,"z1_phys", logicWorld, false, 0, false);
  G4VPhysicalVolume* z2_phys = new G4PVPlacement(0,mylar2_pos+add,z2_log,"z2_phys", logicWorld, false, 0, false);
  G4VPhysicalVolume* z3_phys = new G4PVPlacement(0,tar_pos+add,z3_log,"z3_phys", logicWorld, false, 0, false);
  G4VPhysicalVolume* z4_phys = new G4PVPlacement(0,DSSD1_pos+add,z4_log,"z4_phys", logicWorld, false, 0, false);
  G4VPhysicalVolume* z5_phys = new G4PVPlacement(0,DSSD2_pos+add,z5_log,"z5_phys", logicWorld, false, 0, false);

  alu1_log->SetVisAttributes(G4Colour(1,0,0));
  alu2_log->SetVisAttributes(G4Colour(1,0,0));
  mylar1_log->SetVisAttributes(G4Colour(1,0,0));
  mylar2_log->SetVisAttributes(G4Colour(1,0,0));
  DSSD1_log->SetVisAttributes(G4Colour(0,1,0));
  DSSD2_log->SetVisAttributes(G4Colour(0,1,0));
  //tar_log->SetVisAttributes(G4Colour(0.756,0.624,0.196));
  tar_log->SetVisAttributes(G4Colour(1,0,0));

  /*
  alu1_log->SetVisAttributes(G4VisAttributes::Invisible);
  alu2_log->SetVisAttributes(G4VisAttributes::Invisible);
  mylar1_log->SetVisAttributes(G4VisAttributes::Invisible);
  mylar2_log->SetVisAttributes(G4VisAttributes::Invisible);
  DSSD1_log->SetVisAttributes(G4VisAttributes::Invisible);
  DSSD2_log->SetVisAttributes(G4VisAttributes::Invisible);
  tar_log->SetVisAttributes(G4VisAttributes::Invisible);
  */
  ge_log->SetVisAttributes(G4VisAttributes::Invisible);
  //DSSSD_sphere_log->SetVisAttributes(G4Colour(1,0,0));
  //Feed logical volumes to event action
  fMCP1 = mylar1_log;
  fMCP2 = mylar2_log;
  fDSSSD1 = DSSD1_log;
  fdeg = deg_log;
  fDSSSD2 = DSSD2_log;
  //fDSSSD2 = DSSSD_sphere_log;
  fge = ge_log;
  ftar_log = tar_log;
  fz1_log = z1_log;
  fz2_log = z2_log;
  fz3_log = z3_log;
  fz4_log = z4_log;
  fz5_log = z5_log;
}

void G4DetectorConstruction::ConstructCubicGe() {
  G4double minRadiusGer1 = 0*cm;        	 // Inner radius
  G4double maxRadiusGer1 = 3.5*cm;	// outer radius
  G4double HalfLengthGer1 = 3.9*cm;     // half length
  G4double startPhiGer1 = 0*deg;
  G4double deltaPhiGer1 = 360*deg;
  //Dimention of Borehole in Ge crystal
  G4double minRadiusHole1 = 0*cm;
  G4double maxRadiusHole1 = 0.6*cm;
  G4double HalfLengthHole1 = 3.2*cm;
  G4double startPhiHole1 = 0*deg;
  G4double deltaPhiHole1 = 360*deg;

   G4double w1dx2 = 3.5/2*cm;
   G4double w1dx1 = 3.5/2*cm;
   G4double w1dy2 = 0.55*cm;
   G4double w1dy1 = 0.001*cm;
   G4double w1dz =  3.8*cm;

  G4ThreeVector positionHole1 = G4ThreeVector(0.0*cm,0.0*cm,-0.71*cm);

  G4Tubs* solidG1   = new G4Tubs("G1",minRadiusGer1,maxRadiusGer1,HalfLengthGer1,startPhiGer1,deltaPhiGer1);
  G4Trd*  solidWed1 = new G4Trd("Wed1",w1dx1,w1dx2,w1dy1,w1dy2,w1dz);
  G4Tubs* solidHole = new G4Tubs("Hole",minRadiusHole1,maxRadiusHole1,HalfLengthHole1,startPhiHole1,deltaPhiHole1);

  G4ThreeVector positionWed1 = G4ThreeVector(-0*cm,(35),0.11*cm); //
  G4SubtractionSolid* G1minusW1 = new G4SubtractionSolid("Ger1",solidG1,solidWed1,0,positionWed1); //

  G4ThreeVector positionWed2 = G4ThreeVector(-0*cm,-35,0.11*cm); //
  G4SubtractionSolid* G1minusW1W2 = new G4SubtractionSolid("Ger1",G1minusW1,solidWed1,0,positionWed2); //

  G4ThreeVector positionWed3 = G4ThreeVector(3.50*cm * std::cos(30*deg), 3.50*cm * std::sin(30*deg),0.11*cm); //
  G4RotationMatrix *rotW3  = new G4RotationMatrix;
  rotW3->rotateZ(60*deg);
  G4SubtractionSolid* G1minusW1W2W3 = new G4SubtractionSolid("Ger1",G1minusW1W2,solidWed1,rotW3,positionWed3); //

  G4ThreeVector positionWed4 = G4ThreeVector(-3.50*cm * std::cos(30*deg), -3.50*cm * std::sin(30*deg),0.11*cm); //
  G4RotationMatrix *rotW4  = new G4RotationMatrix;
  rotW4->rotateZ(60*deg);
  G4SubtractionSolid* G1minusW1W2W3W4 = new G4SubtractionSolid("Ger1",G1minusW1W2W3,solidWed1,rotW4,positionWed4); //

  G4ThreeVector positionWed5 = G4ThreeVector(3.50*cm * std::cos(30*deg), -3.50*cm * std::sin(30*deg),0.11*cm); //
  G4RotationMatrix *rotW5  = new G4RotationMatrix;
  rotW5->rotateZ(-60*deg);
  G4SubtractionSolid* G1minusW1W2W3W4W5 = new G4SubtractionSolid("Ger1",G1minusW1W2W3W4,solidWed1,rotW5,positionWed5); //

  G4ThreeVector positionWed6 = G4ThreeVector(-3.50*cm * std::cos(30*deg), 3.50*cm * std::sin(30*deg),0.11*cm); //
  G4RotationMatrix *rotW6  = new G4RotationMatrix;
  rotW6->rotateZ(-60*deg);
  G4SubtractionSolid* G1minusW1W2W3W4W5W6 = new G4SubtractionSolid("Ger1",G1minusW1W2W3W4W5,solidWed1,rotW6,positionWed6); //

  G4SubtractionSolid* G1minusW1W2W3W4W5W6H1 = new G4SubtractionSolid("Ger11",G1minusW1W2W3W4W5W6,solidHole,0,positionHole1);

  // crystal CAPSULE
  G4Tubs* solidG1_capsule   = new G4Tubs("G1_capsule",minRadiusGer1,(maxRadiusGer1+1), (HalfLengthGer1+0.5), startPhiGer1, deltaPhiGer1);
  G4Trd*  solidWed1_capsule = new G4Trd("Wed1",(w1dx1+0.5),(w1dx2+0.5),(w1dy1+0.5),(w1dy2+0.5), (w1dz+0.5) );

  G4ThreeVector positionWed1_capsule = G4ThreeVector(-0*cm, (35 + 1), (0.11+0.1)*cm); //
  G4SubtractionSolid* G1minusW1_capsule = new G4SubtractionSolid("Ger1_capsule",solidG1_capsule, solidWed1_capsule, 0,positionWed1_capsule); //

  G4ThreeVector positionWed2_capsule = G4ThreeVector(-0*cm,-(35+1), (0.11+0.1)*cm); //
  G4SubtractionSolid* G1minusW1W2_capsule = new G4SubtractionSolid("Ger1_capsule",G1minusW1_capsule, solidWed1_capsule, 0, positionWed2_capsule); //

  G4ThreeVector positionWed3_capsule = G4ThreeVector((3.5+0.1)*cm * std::cos(30*deg), (3.5+0.1)*cm * std::sin(30*deg),(0.11+0.1)*cm); //
  G4RotationMatrix *rotW3_capsule  = new G4RotationMatrix;
  rotW3_capsule->rotateZ(60*deg);
  G4SubtractionSolid* G1minusW1W2W3_capsule = new G4SubtractionSolid("Ger1_capsule",G1minusW1W2_capsule, solidWed1_capsule, rotW3_capsule, positionWed3_capsule); //

  G4ThreeVector positionWed4_capsule = G4ThreeVector(-(3.5+0.1)*cm * std::cos(30*deg), -(3.5+0.1)*cm * std::sin(30*deg),(0.11+0.1)*cm); //
  G4RotationMatrix *rotW4_capsule  = new G4RotationMatrix;
  rotW4_capsule->rotateZ(60*deg);
  G4SubtractionSolid* G1minusW1W2W3W4_capsule = new G4SubtractionSolid("Ger1_capsule",G1minusW1W2W3_capsule, solidWed1_capsule, rotW4_capsule, positionWed4_capsule); //

  G4ThreeVector positionWed5_capsule = G4ThreeVector((3.5+0.1)*cm * std::cos(30*deg), -(3.5+0.1)*cm * std::sin(30*deg), (0.11+0.1)*cm); //
  G4RotationMatrix *rotW5_capsule  = new G4RotationMatrix;
  rotW5_capsule->rotateZ(-60*deg);
  G4SubtractionSolid* G1minusW1W2W3W4W5_capsule = new G4SubtractionSolid("Ger1_capsule",G1minusW1W2W3W4_capsule, solidWed1_capsule, rotW5_capsule, positionWed5_capsule); //

  G4ThreeVector positionWed6_capsule = G4ThreeVector(-(3.5+0.1)*cm * std::cos(30*deg), (3.5+0.1)*cm * std::sin(30*deg), (0.11+0.1)*cm); //
  G4RotationMatrix *rotW6_capsule  = new G4RotationMatrix;
  rotW6_capsule->rotateZ(-60*deg);
  G4SubtractionSolid* G1minusW1W2W3W4W5W6_capsule = new G4SubtractionSolid("Ger1_capsule",G1minusW1W2W3W4W5_capsule, solidWed1_capsule, rotW6_capsule, positionWed6_capsule); //

  G4SubtractionSolid* G1minusW1W2W3W4W5W6_capsuleShell = new G4SubtractionSolid("Ger1_capsule",G1minusW1W2W3W4W5W6_capsule, G1minusW1W2W3W4W5W6, 0, G4ThreeVector(0,0,0) );

  // Al housing
  G4double phiStart_poly = 0.0*deg;
  G4double phiTotal_poly= 360 *deg;
  G4double z_Al_casing[] = { 0.0*cm, 10*cm };
  G4double rInner_Al_casing[] = { 0.0*cm, 0.0*cm };
  G4double Outer_Al_casing[] =  { 66/2, 80/2 };
  G4double Outer_Al_casingIn[] =  { (66/2-0.5), (80/2-0.5) };
  G4int numSide = 6;
  G4int numZPlanes = 2;
  G4double z_vacumme_length[] = { 0.0*cm, 9.8*cm };
  G4double Outer_Al_casing2[] =  { 66/2, 66/2 };

  G4Polyhedra *outherMost_Al_solid1 = new G4Polyhedra("outherMost_solid1", phiStart_poly, phiTotal_poly, numSide, numZPlanes, z_Al_casing, rInner_Al_casing, Outer_Al_casing );

  G4Polyhedra *vacumme_after_outherMost_Al_solid = new G4Polyhedra("vacumme_after_outherMost_Al_solid", 0*deg, 360*deg, 6, numZPlanes, z_vacumme_length, rInner_Al_casing, Outer_Al_casingIn);


  //outermost Al housing part1
    G4Box* support_solid_forUnion = new G4Box("support_solid_forUnion",0.00001,0.00001,0.000001);
    G4RotationMatrix *rot1_alHouse  = new G4RotationMatrix;
    rot1_alHouse->rotateX(-4.1*deg);
    rot1_alHouse->rotateY(4.1/2*deg);
    rot1_alHouse->rotateZ(-120*deg);
  G4UnionSolid* outherMost_Al_solidPart1 = new G4UnionSolid("outherMost_Al_solidPart1",support_solid_forUnion, outherMost_Al_solid1, rot1_alHouse, G4ThreeVector(1,0,-4.6*cm) ); //

  //outermost Al housing part2
    G4RotationMatrix *rot2_alHouse  = new G4RotationMatrix;
    rot2_alHouse->rotateX(4.1*deg);
    rot2_alHouse->rotateY(4.1/2*deg);
  G4UnionSolid* outherMost_Al_solidPart2 = new G4UnionSolid("outherMost_Al_solidPart2",outherMost_Al_solidPart1, outherMost_Al_solid1, rot2_alHouse, G4ThreeVector(1,65.93,-4.6*cm) ); //

  //outermost Al housing part3
    G4RotationMatrix *rot3_alHouse  = new G4RotationMatrix;
    rot3_alHouse->rotateY(-4.1*deg);
    rot3_alHouse->rotateY(-4.1/4*deg);
    rot3_alHouse->rotateZ(120*deg);
  G4UnionSolid* outherMost_Al_solidPart123 = new G4UnionSolid("outherMost_Al_solidPart3",outherMost_Al_solidPart2, outherMost_Al_solid1, rot3_alHouse, G4ThreeVector(58, 32.6, -4.6*cm) ); //

  // Vacumme inside Outermost Al casing
  // part1
  G4UnionSolid* vacumme_outherMost_Al_solidPart1 = new G4UnionSolid("vacumme_outherMost_Al_solidPart1",support_solid_forUnion, vacumme_after_outherMost_Al_solid, rot1_alHouse, G4ThreeVector(1,0,-4.6*cm) ); //
  // part2
  G4UnionSolid* vacumme_outherMost_Al_solidPart2 = new G4UnionSolid("vacumme_outherMost_Al_solidPart2",vacumme_outherMost_Al_solidPart1, vacumme_after_outherMost_Al_solid, rot2_alHouse, G4ThreeVector(1,65.93,-4.6*cm) ); //
  // part3
  G4UnionSolid* vacumme_outherMost_Al_solidPart123 = new G4UnionSolid("vacumme_outherMost_Al_solidPart3",vacumme_outherMost_Al_solidPart2, vacumme_after_outherMost_Al_solid, rot3_alHouse, G4ThreeVector(58, 32.6, -4.6*cm) ); //

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4int numRZ3 = 8;
  G4double z3[]      = {0.0*cm,  7.76*cm, 7.76001*cm, (7.76+34.35)*cm, (7.76+34.3500001)*cm ,(34.3500001 + 63.76)*cm,(34.3500002 +63.76)*cm,101*cm };
  G4double rInner3[] = {0.0*cm,  0.00*cm,       0*cm , 0*cm   ,            0.0*cm,   		0.0*cm,                   0.0*cm,         0.0*cm  };
  G4double rOuter3[] = {4.7/2*cm, 4.7/2*cm, 5.7/2*cm,    5.7/2*cm ,       11.8/2*cm, 		11.8/2*cm  ,              5.5/2*cm, 5.5/2*cm};
  G4Polycone* Al_Back_pipe_part1_solid = new G4Polycone("Al_Back_pipe_part1_solid", 0*deg, 360.0*deg, numRZ3, z3, rInner3, rOuter3);

  G4LogicalVolume* Al_Back_part_logi = new G4LogicalVolume(Al_Back_pipe_part1_solid, alu_mat, "Al_back_part");

  char  crystal_no[30], detector_no[30], crystal_Al_casing[30], Al_casing_of_crystal[50], crystal[50], crystal1[50], outermost_Al_casing[50], DetCapsule1[30], DetCapsule2[30], DetCapsule3[30], vacumme_outermost_Al_casing[40];
  G4LogicalVolume* crystal_logi[5][3];
  G4LogicalVolume* capsule_logi[5][3];
  G4LogicalVolume* outherMost_Al_logi[5];
  G4LogicalVolume* vacumme_outherMost_Al_logi[5];

  G4VisAttributes *red_color = new G4VisAttributes(true,G4Colour(1.0,0.0,0.0,1.));
  G4VisAttributes *green_colour = new G4VisAttributes(true,G4Colour(0.0,1.0,0.0));
  G4VisAttributes *blue_colour = new G4VisAttributes(true,G4Colour(0,0.5,1.0));

  for (int id = 0; id<5; id++){
    for(int id1= 0; id1<3; id1++){

      sprintf(Al_casing_of_crystal,"Al_casing_of_crysta%i_%i",id,id1);
    sprintf(crystal,"det%i_crystal%i",id,id1);
    sprintf(crystal_no,"DetVac%i",id);

  // All logical volume with proper name:
      crystal_logi[id][id1] = new G4LogicalVolume(G1minusW1W2W3W4W5W6H1, ge_mat, crystal);
    }

  crystal_logi[id][0]->SetVisAttributes(red_color);
  crystal_logi[id][1]->SetVisAttributes(green_colour);
  crystal_logi[id][2]->SetVisAttributes(blue_colour);
  sprintf(outermost_Al_casing,"Outermost_Al_casing%i",id);
  outherMost_Al_logi[id] = new G4LogicalVolume(outherMost_Al_solidPart123, alu_mat, outermost_Al_casing);
  sprintf(vacumme_outermost_Al_casing,"vacumme_Outermost_Al_casing%i",id);
  vacumme_outherMost_Al_logi[id] = new G4LogicalVolume(vacumme_outherMost_Al_solidPart123, vacuum, vacumme_outermost_Al_casing);

  // placement of crystal inside vacumme
  sprintf(DetCapsule1,"Det%i_crys0",id);
  sprintf(DetCapsule2,"Det%i_crys1",id);
  sprintf(DetCapsule3,"Det%i_crys2",id);

  // Ge1
  G4RotationMatrix *rot1  = new G4RotationMatrix;
  rot1->rotateY(-180*deg);
  rot1->rotateX(4.1*deg);
  rot1->rotateY(4.1/2*deg);
  new G4PVPlacement(rot1,G4ThreeVector(-0.2,-1.6,0), crystal_logi[id][0],DetCapsule1,  vacumme_outherMost_Al_logi[id],false,(3*id+0),checkOverlaps);                          //overlaps checking

  // Ge2
  G4RotationMatrix *rot2  = new G4RotationMatrix;
  rot2->rotateY(-180*deg);
  rot2->rotateX(-4.1*deg);
  rot2->rotateY(4.1*deg);
  new G4PVPlacement(rot2,G4ThreeVector(0,67,0), crystal_logi[id][1],DetCapsule2, vacumme_outherMost_Al_logi[id],false,(3*id+0),checkOverlaps);                          //overlaps checking

  // Ge3
  G4RotationMatrix *rot3  = new G4RotationMatrix;
  rot3->rotateY(-180*deg);
  rot3->rotateY(-4.1*deg);
  new G4PVPlacement(rot3,G4ThreeVector(59,33.0,0), crystal_logi[id][2], DetCapsule3, vacumme_outherMost_Al_logi[id], false, (3*id+1), checkOverlaps);                          // overlaps checking


  // placement of vacumme inside outermost Al housing
  new G4PVPlacement(0, G4ThreeVector(0,0,1),  vacumme_outherMost_Al_logi[id], crystal_no, outherMost_Al_logi[id],    false,  id,   checkOverlaps);
  }

  G4RotationMatrix *rot_det1  = new G4RotationMatrix;

  rot_det1->rotateX(90*deg);
  new G4PVPlacement(rot_det1,G4ThreeVector(-1.75*cm,31.5*cm-3.5*cm,2.5*cm+MCP_Sep+24*cm+1*cm), outherMost_Al_logi[0], "AlHouse0",logicWorld, false, 0, checkOverlaps);
  G4RotationMatrix *rot_det2  = new G4RotationMatrix;

  rot_det2->rotateX(-90*deg);
  new G4PVPlacement(rot_det2,G4ThreeVector(-2.75*cm,-25*cm-3.5*cm,-2.5*cm+MCP_Sep+23*cm), outherMost_Al_logi[1], "AlHouse1",logicWorld, false, 1, checkOverlaps);
  G4RotationMatrix *rot_det3  = new G4RotationMatrix;

  rot_det3->rotateY(-90*deg);
  new G4PVPlacement(rot_det3,G4ThreeVector(29*cm-3.5*cm,-3.5*cm,2.5*cm+MCP_Sep+24*cm), outherMost_Al_logi[2], "AlHouse2",logicWorld, false, 2, checkOverlaps);
  G4RotationMatrix *rot_det4  = new G4RotationMatrix;

  rot_det4->rotateY(90*deg);
  new G4PVPlacement(rot_det4,G4ThreeVector(-21*cm-3.5*cm,-3.5*cm,-2.5*cm+MCP_Sep+24*cm), outherMost_Al_logi[3], "AlHouse3",logicWorld, false, 3, checkOverlaps);
  G4RotationMatrix *rot_det5  = new G4RotationMatrix;//

  new G4PVPlacement(rot_det5,G4ThreeVector(-3.5*cm,-3.5*cm,29*cm+MCP_Sep+24*cm), outherMost_Al_logi[4], "AlHouse4",logicWorld, false, 4, checkOverlaps);
  G4RotationMatrix *rot_det6  = new G4RotationMatrix;
}

void G4DetectorConstruction::ConstructWires()
{
  G4double seperation = 0.1*cm;
  G4double radians = 45*4*atan(1.0)/180.0;
  G4Box* solid_wire_hori = new G4Box("solid_wire_hori", foilsizeXY, 0.1*mm, 0.1*mm);
  //G4Box* solid_wire_hori = new G4Box("solid_wire_hori", foilsizeXY, 50*nm, 50*nm);
  G4double Nwires = 50;

  //First MCP
  G4int wireNo = 0;
  G4ThreeVector wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire01 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire01");
  G4VPhysicalVolume* phys_wire01 = new G4PVPlacement(0,wire_pos,logical_wire01,"phys_wire01",logicWorld,false,0);

  wireNo = 1;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire11 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire11");
  G4VPhysicalVolume* phys_wire11 = new G4PVPlacement(0,wire_pos,logical_wire11,"phys_wire11",logicWorld,false,0);

  wireNo = 2;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire21 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire21");
  G4VPhysicalVolume* phys_wire21 = new G4PVPlacement(0,wire_pos,logical_wire21,"phys_wire21",logicWorld,false,0);

  wireNo = 3;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire31 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire31");
  G4VPhysicalVolume* phys_wire31 = new G4PVPlacement(0,wire_pos,logical_wire31,"phys_wire31",logicWorld,false,0);

  wireNo = 4;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire41 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire41");
  G4VPhysicalVolume* phys_wire41 = new G4PVPlacement(0,wire_pos,logical_wire41,"phys_wire41",logicWorld,false,0);

  wireNo = 5;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire51 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire51");
  G4VPhysicalVolume* phys_wire51 = new G4PVPlacement(0,wire_pos,logical_wire51,"phys_wire51",logicWorld,false,0);

  wireNo = 6;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire61 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire61");
  G4VPhysicalVolume* phys_wire611 = new G4PVPlacement(0,wire_pos,logical_wire61,"phys_wire61",logicWorld,false,0);

  wireNo = 7;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire71 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire71");
  G4VPhysicalVolume* phys_wire71 = new G4PVPlacement(0,wire_pos,logical_wire71,"phys_wire71",logicWorld,false,0);

  wireNo = 8;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire81 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire81");
  G4VPhysicalVolume* phys_wire81 = new G4PVPlacement(0,wire_pos,logical_wire81,"phys_wire81",logicWorld,false,0);

  wireNo = 9;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire91 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire91");
  G4VPhysicalVolume* phys_wire91 = new G4PVPlacement(0,wire_pos,logical_wire91,"phys_wire91",logicWorld,false,0);

  wireNo = 10;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire101 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire101");
  G4VPhysicalVolume* phys_wire101 = new G4PVPlacement(0,wire_pos,logical_wire101,"phys_wire101",logicWorld,false,0);

  wireNo = 11;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire111 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire111");
  G4VPhysicalVolume* phys_wire111 = new G4PVPlacement(0,wire_pos,logical_wire111,"phys_wire111",logicWorld,false,0);

  wireNo = 12;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire121 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire121");
  G4VPhysicalVolume* phys_wire121 = new G4PVPlacement(0,wire_pos,logical_wire121,"phys_wire121",logicWorld,false,0);

  wireNo = 13;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire131 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire131");
  G4VPhysicalVolume* phys_wire131 = new G4PVPlacement(0,wire_pos,logical_wire131,"phys_wire131",logicWorld,false,0);

  wireNo = 14;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire141 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire141");
  G4VPhysicalVolume* phys_wire141 = new G4PVPlacement(0,wire_pos,logical_wire141,"phys_wire141",logicWorld,false,0);

  wireNo = 15;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire151 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire151");
  G4VPhysicalVolume* phys_wire151 = new G4PVPlacement(0,wire_pos,logical_wire151,"phys_wire151",logicWorld,false,0);

  wireNo = 16;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire161 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire161");
  G4VPhysicalVolume* phys_wire161 = new G4PVPlacement(0,wire_pos,logical_wire161,"phys_wire161",logicWorld,false,0);

  wireNo = 17;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire171 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire171");
  G4VPhysicalVolume* phys_wire171 = new G4PVPlacement(0,wire_pos,logical_wire171,"phys_wire171",logicWorld,false,0);

  wireNo = 18;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire181 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire181");
  G4VPhysicalVolume* phys_wire181 = new G4PVPlacement(0,wire_pos,logical_wire181,"phys_wire181",logicWorld,false,0);

  wireNo = 19;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire191 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire191");
  G4VPhysicalVolume* phys_wire191 = new G4PVPlacement(0,wire_pos,logical_wire191,"phys_wire191",logicWorld,false,0);

  wireNo = 20;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire201 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire201");
  G4VPhysicalVolume* phys_wire201 = new G4PVPlacement(0,wire_pos,logical_wire201,"phys_wire201",logicWorld,false,0);

  wireNo = 21;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire211 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire211");
  G4VPhysicalVolume* phys_wire211 = new G4PVPlacement(0,wire_pos,logical_wire211,"phys_wire211",logicWorld,false,0);

  wireNo = 22;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire221 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire221");
  G4VPhysicalVolume* phys_wire221 = new G4PVPlacement(0,wire_pos,logical_wire221,"phys_wire221",logicWorld,false,0);

  wireNo = 23;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire231 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire231");
  G4VPhysicalVolume* phys_wire231 = new G4PVPlacement(0,wire_pos,logical_wire231,"phys_wire231",logicWorld,false,0);

  wireNo = 24;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire241 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire241");
  G4VPhysicalVolume* phys_wire241 = new G4PVPlacement(0,wire_pos,logical_wire241,"phys_wire241",logicWorld,false,0);

  wireNo = 25;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire251 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire251");
  G4VPhysicalVolume* phys_wire251 = new G4PVPlacement(0,wire_pos,logical_wire251,"phys_wire251",logicWorld,false,0);

  wireNo = 26;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire261 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire261");
  G4VPhysicalVolume* phys_wire261 = new G4PVPlacement(0,wire_pos,logical_wire261,"phys_wire261",logicWorld,false,0);

  wireNo = 27;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire271 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire271");
  G4VPhysicalVolume* phys_wire271 = new G4PVPlacement(0,wire_pos,logical_wire271,"phys_wire271",logicWorld,false,0);

  wireNo = 28;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire281 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire281");
  G4VPhysicalVolume* phys_wire281 = new G4PVPlacement(0,wire_pos,logical_wire281,"phys_wire281",logicWorld,false,0);

  wireNo = 29;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire291 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire291");
  G4VPhysicalVolume* phys_wire291 = new G4PVPlacement(0,wire_pos,logical_wire291,"phys_wire291",logicWorld,false,0);

  wireNo = 30;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire301 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire301");
  G4VPhysicalVolume* phys_wire301 = new G4PVPlacement(0,wire_pos,logical_wire301,"phys_wire301",logicWorld,false,0);

  wireNo = 31;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire311 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire311");
  G4VPhysicalVolume* phys_wire311 = new G4PVPlacement(0,wire_pos,logical_wire311,"phys_wire311",logicWorld,false,0);

  wireNo = 32;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire321 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire321");
  G4VPhysicalVolume* phys_wire321 = new G4PVPlacement(0,wire_pos,logical_wire321,"phys_wire321",logicWorld,false,0);

  wireNo = 33;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire331 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire331");
  G4VPhysicalVolume* phys_wire331 = new G4PVPlacement(0,wire_pos,logical_wire331,"phys_wire331",logicWorld,false,0);

  wireNo = 34;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire341 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire341");
  G4VPhysicalVolume* phys_wire341 = new G4PVPlacement(0,wire_pos,logical_wire341,"phys_wire341",logicWorld,false,0);

  wireNo = 35;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire351 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire351");
  G4VPhysicalVolume* phys_wire351 = new G4PVPlacement(0,wire_pos,logical_wire351,"phys_wire351",logicWorld,false,0);

  wireNo = 36;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire361 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire361");
  G4VPhysicalVolume* phys_wire361 = new G4PVPlacement(0,wire_pos,logical_wire361,"phys_wire361",logicWorld,false,0);

  wireNo = 37;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire371 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire371");
  G4VPhysicalVolume* phys_wire371 = new G4PVPlacement(0,wire_pos,logical_wire371,"phys_wire371",logicWorld,false,0);

  wireNo = 38;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire381 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire381");
  G4VPhysicalVolume* phys_wire381 = new G4PVPlacement(0,wire_pos,logical_wire381,"phys_wire381",logicWorld,false,0);

  wireNo = 39;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire391 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire391");
  G4VPhysicalVolume* phys_wire391 = new G4PVPlacement(0,wire_pos,logical_wire391,"phys_wire391",logicWorld,false,0);

  wireNo = 40;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire401 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire401");
  G4VPhysicalVolume* phys_wire401 = new G4PVPlacement(0,wire_pos,logical_wire401,"phys_wire401",logicWorld,false,0);

  wireNo = 41;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire411 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire411");
  G4VPhysicalVolume* phys_wire411 = new G4PVPlacement(0,wire_pos,logical_wire411,"phys_wire411",logicWorld,false,0);

  wireNo = 42;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire421 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire421");
  G4VPhysicalVolume* phys_wire421 = new G4PVPlacement(0,wire_pos,logical_wire421,"phys_wire421",logicWorld,false,0);

  wireNo = 43;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire431 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire431");
  G4VPhysicalVolume* phys_wire431 = new G4PVPlacement(0,wire_pos,logical_wire431,"phys_wire431",logicWorld,false,0);

  wireNo = 44;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire441 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire441");
  G4VPhysicalVolume* phys_wire441 = new G4PVPlacement(0,wire_pos,logical_wire441,"phys_wire441",logicWorld,false,0);

  wireNo = 45;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire451 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire451");
  G4VPhysicalVolume* phys_wire451 = new G4PVPlacement(0,wire_pos,logical_wire451,"phys_wire451",logicWorld,false,0);

  wireNo = 46;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire461 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire461");
  G4VPhysicalVolume* phys_wire461 = new G4PVPlacement(0,wire_pos,logical_wire461,"phys_wire461",logicWorld,false,0);

  wireNo = 47;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire471 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire471");
  G4VPhysicalVolume* phys_wire471 = new G4PVPlacement(0,wire_pos,logical_wire471,"phys_wire471",logicWorld,false,0);

  wireNo = 48;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire481 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire481");
  G4VPhysicalVolume* phys_wire481 = new G4PVPlacement(0,wire_pos,logical_wire481,"phys_wire481",logicWorld,false,0);

  wireNo = 49;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707);
  G4LogicalVolume* logical_wire491 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire491");
  G4VPhysicalVolume* phys_wire491 = new G4PVPlacement(0,wire_pos,logical_wire491,"phys_wire491",logicWorld,false,0);



  //Second MCP
  wireNo = 0;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire02 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire02");
  G4VPhysicalVolume* phys_wire02 = new G4PVPlacement(0,wire_pos,logical_wire02,"phys_wire02",logicWorld,false,0);

  wireNo = 1;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire12 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire12");
  G4VPhysicalVolume* phys_wire12 = new G4PVPlacement(0,wire_pos,logical_wire12,"phys_wire12",logicWorld,false,0);

  wireNo = 2;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire22 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire22");
  G4VPhysicalVolume* phys_wire22 = new G4PVPlacement(0,wire_pos,logical_wire22,"phys_wire22",logicWorld,false,0);

  wireNo = 3;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire32 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire32");
  G4VPhysicalVolume* phys_wire32 = new G4PVPlacement(0,wire_pos,logical_wire32,"phys_wire32",logicWorld,false,0);

  wireNo = 4;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire42 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire42");
  G4VPhysicalVolume* phys_wire42 = new G4PVPlacement(0,wire_pos,logical_wire42,"phys_wire42",logicWorld,false,0);

  wireNo = 5;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire52 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire52");
  G4VPhysicalVolume* phys_wire52 = new G4PVPlacement(0,wire_pos,logical_wire52,"phys_wire52",logicWorld,false,0);

  wireNo = 6;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire62 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire62");
  G4VPhysicalVolume* phys_wire62 = new G4PVPlacement(0,wire_pos,logical_wire62,"phys_wire62",logicWorld,false,0);

  wireNo = 7;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire72 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire72");
  G4VPhysicalVolume* phys_wire72 = new G4PVPlacement(0,wire_pos,logical_wire72,"phys_wire72",logicWorld,false,0);

  wireNo = 8;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire82 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire82");
  G4VPhysicalVolume* phys_wire82 = new G4PVPlacement(0,wire_pos,logical_wire82,"phys_wire82",logicWorld,false,0);

  wireNo = 9;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire92 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire92");
  G4VPhysicalVolume* phys_wire92 = new G4PVPlacement(0,wire_pos,logical_wire92,"phys_wire92",logicWorld,false,0);

  wireNo = 10;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire102 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire102");
  G4VPhysicalVolume* phys_wire102 = new G4PVPlacement(0,wire_pos,logical_wire102,"phys_wire102",logicWorld,false,0);

  wireNo = 11;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire112 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire112");
  G4VPhysicalVolume* phys_wire112 = new G4PVPlacement(0,wire_pos,logical_wire112,"phys_wire112",logicWorld,false,0);

  wireNo = 12;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire122 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire122");
  G4VPhysicalVolume* phys_wire122 = new G4PVPlacement(0,wire_pos,logical_wire122,"phys_wire122",logicWorld,false,0);

  wireNo = 13;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire132 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire132");
  G4VPhysicalVolume* phys_wire132 = new G4PVPlacement(0,wire_pos,logical_wire132,"phys_wire132",logicWorld,false,0);

  wireNo = 14;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire142 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire142");
  G4VPhysicalVolume* phys_wire142 = new G4PVPlacement(0,wire_pos,logical_wire142,"phys_wire142",logicWorld,false,0);

  wireNo = 15;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire152 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire152");
  G4VPhysicalVolume* phys_wire152 = new G4PVPlacement(0,wire_pos,logical_wire152,"phys_wire152",logicWorld,false,0);

  wireNo = 16;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire162 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire162");

  G4VPhysicalVolume* phys_wire162 = new G4PVPlacement(0,wire_pos,logical_wire162,"phys_wire162",logicWorld,false,0);

  wireNo = 17;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire172 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire172");
  G4VPhysicalVolume* phys_wire172 = new G4PVPlacement(0,wire_pos,logical_wire172,"phys_wire172",logicWorld,false,0);

  wireNo = 18;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire182 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire182");
  G4VPhysicalVolume* phys_wire182 = new G4PVPlacement(0,wire_pos,logical_wire182,"phys_wire182",logicWorld,false,0);

  wireNo = 19;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire192 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire192");
  G4VPhysicalVolume* phys_wire192 = new G4PVPlacement(0,wire_pos,logical_wire192,"phys_wire192",logicWorld,false,0);

  wireNo = 20;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire202 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire202");
  G4VPhysicalVolume* phys_wire202 = new G4PVPlacement(0,wire_pos,logical_wire202,"phys_wire202",logicWorld,false,0);

  wireNo = 21;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire212 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire212");
  G4VPhysicalVolume* phys_wire212 = new G4PVPlacement(0,wire_pos,logical_wire212,"phys_wire212",logicWorld,false,0);

  wireNo = 22;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire222 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire222");
  G4VPhysicalVolume* phys_wire222 = new G4PVPlacement(0,wire_pos,logical_wire222,"phys_wire222",logicWorld,false,0);

  wireNo = 23;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire232 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire232");
  G4VPhysicalVolume* phys_wire232 = new G4PVPlacement(0,wire_pos,logical_wire232,"phys_wire232",logicWorld,false,0);

  wireNo = 24;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire242 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire242");
  G4VPhysicalVolume* phys_wire242 = new G4PVPlacement(0,wire_pos,logical_wire242,"phys_wire242",logicWorld,false,0);

  wireNo = 25;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire252 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire252");
  G4VPhysicalVolume* phys_wire252 = new G4PVPlacement(0,wire_pos,logical_wire252,"phys_wire252",logicWorld,false,0);

  wireNo = 26;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire262 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire262");
  G4VPhysicalVolume* phys_wire262 = new G4PVPlacement(0,wire_pos,logical_wire262,"phys_wire262",logicWorld,false,0);

  wireNo = 27;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire272 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire272");
  G4VPhysicalVolume* phys_wire272 = new G4PVPlacement(0,wire_pos,logical_wire272,"phys_wire272",logicWorld,false,0);

  wireNo = 28;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire282 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire282");
  G4VPhysicalVolume* phys_wire282 = new G4PVPlacement(0,wire_pos,logical_wire282,"phys_wire282",logicWorld,false,0);

  wireNo = 29;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire292 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire292");
  G4VPhysicalVolume* phys_wire292 = new G4PVPlacement(0,wire_pos,logical_wire292,"phys_wire292",logicWorld,false,0);

  wireNo = 30;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire302 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire302");
  G4VPhysicalVolume* phys_wire302 = new G4PVPlacement(0,wire_pos,logical_wire302,"phys_wire302",logicWorld,false,0);

  wireNo = 31;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire312 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire312");
  G4VPhysicalVolume* phys_wire312 = new G4PVPlacement(0,wire_pos,logical_wire312,"phys_wire312",logicWorld,false,0);

  wireNo = 32;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire322 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire322");
  G4VPhysicalVolume* phys_wire322 = new G4PVPlacement(0,wire_pos,logical_wire322,"phys_wire322",logicWorld,false,0);

  wireNo = 33;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire332 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire332");
  G4VPhysicalVolume* phys_wire332 = new G4PVPlacement(0,wire_pos,logical_wire332,"phys_wire332",logicWorld,false,0);

  wireNo = 34;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire342 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire342");
  G4VPhysicalVolume* phys_wire342 = new G4PVPlacement(0,wire_pos,logical_wire342,"phys_wire342",logicWorld,false,0);

  wireNo = 35;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire352 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire352");
  G4VPhysicalVolume* phys_wire352 = new G4PVPlacement(0,wire_pos,logical_wire352,"phys_wire352",logicWorld,false,0);

  wireNo = 36;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire362 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire362");
  G4VPhysicalVolume* phys_wire362 = new G4PVPlacement(0,wire_pos,logical_wire362,"phys_wire362",logicWorld,false,0);

  wireNo = 37;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire372 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire372");
  G4VPhysicalVolume* phys_wire372 = new G4PVPlacement(0,wire_pos,logical_wire372,"phys_wire372",logicWorld,false,0);

  wireNo = 38;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire382 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire382");
  G4VPhysicalVolume* phys_wire382 = new G4PVPlacement(0,wire_pos,logical_wire382,"phys_wire382",logicWorld,false,0);

  wireNo = 39;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire392 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire392");
  G4VPhysicalVolume* phys_wire392 = new G4PVPlacement(0,wire_pos,logical_wire392,"phys_wire392",logicWorld,false,0);

  wireNo = 40;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire402 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire402");
  G4VPhysicalVolume* phys_wire402 = new G4PVPlacement(0,wire_pos,logical_wire402,"phys_wire402",logicWorld,false,0);

  wireNo = 41;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire412 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire412");
  G4VPhysicalVolume* phys_wire412 = new G4PVPlacement(0,wire_pos,logical_wire412,"phys_wire412",logicWorld,false,0);

  wireNo = 42;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire422 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire422");
  G4VPhysicalVolume* phys_wire422 = new G4PVPlacement(0,wire_pos,logical_wire422,"phys_wire422",logicWorld,false,0);

  wireNo = 43;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire432 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire432");
  G4VPhysicalVolume* phys_wire432 = new G4PVPlacement(0,wire_pos,logical_wire432,"phys_wire432",logicWorld,false,0);

  wireNo = 44;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire442 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire442");
  G4VPhysicalVolume* phys_wire442 = new G4PVPlacement(0,wire_pos,logical_wire442,"phys_wire442",logicWorld,false,0);

  wireNo = 45;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire452 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire452");
  G4VPhysicalVolume* phys_wire452 = new G4PVPlacement(0,wire_pos,logical_wire452,"phys_wire452",logicWorld,false,0);

  wireNo = 46;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire462 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire462");
  G4VPhysicalVolume* phys_wire462 = new G4PVPlacement(0,wire_pos,logical_wire462,"phys_wire462",logicWorld,false,0);

  wireNo = 47;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire472 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire472");
  G4VPhysicalVolume* phys_wire472 = new G4PVPlacement(0,wire_pos,logical_wire472,"phys_wire472",logicWorld,false,0);

  wireNo = 48;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire482 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire482");
  G4VPhysicalVolume* phys_wire482 = new G4PVPlacement(0,wire_pos,logical_wire482,"phys_wire482",logicWorld,false,0);

  wireNo = 49;
  wire_pos = G4ThreeVector(0,
  2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
  2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+MCP_Sep);
  G4LogicalVolume* logical_wire492 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire492");
  G4VPhysicalVolume* phys_wire492 = new G4PVPlacement(0,wire_pos,logical_wire492,"phys_wire492",logicWorld,false,0);



  //###############################Place Vertical wires here ######################################


  //First MCP
  G4Box* solid_wire_vert = new G4Box("solid_wire_vert", 50*nm, foilsizeXY, 50*nm);

  wireNo= 0;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert01 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert01");
  G4VPhysicalVolume* phys_wire_vert01 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert01,"phys_wire_vert01",logicWorld,false,0);

  wireNo= 1;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert11 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert11");
  G4VPhysicalVolume* phys_wire_vert11 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert11,"phys_wire_vert11",logicWorld,false,0);

  wireNo= 2;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert21 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert21");
  G4VPhysicalVolume* phys_wire_vert21 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert21,"phys_wire_vert21",logicWorld,false,0);

  wireNo= 3;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert31 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert31");
  G4VPhysicalVolume* phys_wire_vert31 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert31,"phys_wire_vert31",logicWorld,false,0);

  wireNo= 4;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert41 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert41");
  G4VPhysicalVolume* phys_wire_vert41 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert41,"phys_wire_vert41",logicWorld,false,0);

  wireNo= 5;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert51 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert51");
  G4VPhysicalVolume* phys_wire_vert51 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert51,"phys_wire_vert51",logicWorld,false,0);

  wireNo= 6;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert61 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert61");
  G4VPhysicalVolume* phys_wire_vert61 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert61,"phys_wire_vert61",logicWorld,false,0);

  wireNo= 7;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert71 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert71");
  G4VPhysicalVolume* phys_wire_vert71 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert71,"phys_wire_vert71",logicWorld,false,0);

  wireNo= 8;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert81 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert81");
  G4VPhysicalVolume* phys_wire_vert81 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert81,"phys_wire_vert81",logicWorld,false,0);

  wireNo= 9;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert91 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert91");
  G4VPhysicalVolume* phys_wire_vert91 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert91,"phys_wire_vert91",logicWorld,false,0);

  wireNo= 10;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert101 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert101");
  G4VPhysicalVolume* phys_wire_vert101 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert101,"phys_wire_vert101",logicWorld,false,0);

  wireNo= 11;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert111 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert111");
  G4VPhysicalVolume* phys_wire_vert111 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert111,"phys_wire_vert111",logicWorld,false,0);

  wireNo= 12;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert121 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert121");
  G4VPhysicalVolume* phys_wire_vert121 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert121,"phys_wire_vert121",logicWorld,false,0);

  wireNo= 13;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert131 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert131");
  G4VPhysicalVolume* phys_wire_vert131 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert131,"phys_wire_vert131",logicWorld,false,0);

  wireNo= 14;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert141 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert141");
  G4VPhysicalVolume* phys_wire_vert141 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert141,"phys_wire_vert141",logicWorld,false,0);

  wireNo= 15;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert151 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert151");
  G4VPhysicalVolume* phys_wire_vert151 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert151,"phys_wire_vert151",logicWorld,false,0);

  wireNo= 16;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert161 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert161");
  G4VPhysicalVolume* phys_wire_vert161 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert161,"phys_wire_vert161",logicWorld,false,0);

  wireNo= 17;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert171 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert171");
  G4VPhysicalVolume* phys_wire_vert171 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert171,"phys_wire_vert171",logicWorld,false,0);

  wireNo= 18;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert181 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert181");
  G4VPhysicalVolume* phys_wire_vert181 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert181,"phys_wire_vert181",logicWorld,false,0);

  wireNo= 19;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert191 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert191");
  G4VPhysicalVolume* phys_wire_vert191 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert191,"phys_wire_vert191",logicWorld,false,0);

  wireNo= 20;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert201 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert201");
  G4VPhysicalVolume* phys_wire_vert201 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert201,"phys_wire_vert201",logicWorld,false,0);

  wireNo= 21;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert211 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert211");
  G4VPhysicalVolume* phys_wire_vert211 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert211,"phys_wire_vert211",logicWorld,false,0);
  wireNo= 22;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert221 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert221");
  G4VPhysicalVolume* phys_wire_vert221 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert221,"phys_wire_vert221",logicWorld,false,0);

  wireNo= 23;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert231 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert231");
  G4VPhysicalVolume* phys_wire_vert231 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert231,"phys_wire_vert231",logicWorld,false,0);

  wireNo= 24;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert241 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert241");
  G4VPhysicalVolume* phys_wire_vert241 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert241,"phys_wire_vert241",logicWorld,false,0);

  wireNo= 25;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert251 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert251");
  G4VPhysicalVolume* phys_wire_vert251 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert251,"phys_wire_vert251",logicWorld,false,0);

  wireNo= 26;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert261 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert261");
  G4VPhysicalVolume* phys_wire_vert261 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert261,"phys_wire_vert261",logicWorld,false,0);

  wireNo= 27;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert271 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert271");
  G4VPhysicalVolume* phys_wire_vert271 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert271,"phys_wire_vert271",logicWorld,false,0);

  wireNo= 28;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert281 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert281");
  G4VPhysicalVolume* phys_wire_vert281 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert281,"phys_wire_vert281",logicWorld,false,0);

  wireNo= 29;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert291 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert291");
  G4VPhysicalVolume* phys_wire_vert291 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert291,"phys_wire_vert291",logicWorld,false,0);

  wireNo= 30;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert301 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert301");
  G4VPhysicalVolume* phys_wire_vert301 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert301,"phys_wire_vert301",logicWorld,false,0);

  wireNo= 31;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert311 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert311");
  G4VPhysicalVolume* phys_wire_vert311 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert311,"phys_wire_vert311",logicWorld,false,0);

  wireNo= 32;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert321 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert321");
  G4VPhysicalVolume* phys_wire_vert321 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert321,"phys_wire_vert321",logicWorld,false,0);

  wireNo= 33;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert331 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert331");
  G4VPhysicalVolume* phys_wire_vert331 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert331,"phys_wire_vert331",logicWorld,false,0);

  wireNo= 34;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert341 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert341");
  G4VPhysicalVolume* phys_wire_vert341 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert341,"phys_wire_vert341",logicWorld,false,0);

  wireNo= 35;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert351 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert351");
  G4VPhysicalVolume* phys_wire_vert351 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert351,"phys_wire_vert351",logicWorld,false,0);

  wireNo= 36;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert361 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert361");
  G4VPhysicalVolume* phys_wire_vert361 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert361,"phys_wire_vert361",logicWorld,false,0);

  wireNo= 37;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert371 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert371");
  G4VPhysicalVolume* phys_wire_vert371 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert371,"phys_wire_vert371",logicWorld,false,0);

  wireNo= 38;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert381 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert381");
  G4VPhysicalVolume* phys_wire_vert381 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert381,"phys_wire_vert381",logicWorld,false,0);

  wireNo= 39;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert391 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert391");
  G4VPhysicalVolume* phys_wire_vert391 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert391,"phys_wire_vert391",logicWorld,false,0);

  wireNo= 40;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert401 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert401");
  G4VPhysicalVolume* phys_wire_vert401 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert401,"phys_wire_vert401",logicWorld,false,0);

  wireNo= 41;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert411 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert411");
  G4VPhysicalVolume* phys_wire_vert411 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert411,"phys_wire_vert411",logicWorld,false,0);

  wireNo= 42;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert421 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert421");
  G4VPhysicalVolume* phys_wire_vert421 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert421,"phys_wire_vert421",logicWorld,false,0);

  wireNo= 43;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert431 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert431");
  G4VPhysicalVolume* phys_wire_vert431 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert431,"phys_wire_vert431",logicWorld,false,0);

  wireNo= 44;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert441 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert441");
  G4VPhysicalVolume* phys_wire_vert441 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert441,"phys_wire_vert441",logicWorld,false,0);

  wireNo= 45;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert451 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert451");
  G4VPhysicalVolume* phys_wire_vert451 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert451,"phys_wire_vert451",logicWorld,false,0);

  wireNo= 46;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert461 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert461");
  G4VPhysicalVolume* phys_wire_vert461 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert461,"phys_wire_vert461",logicWorld,false,0);

  wireNo= 47;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert471 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert471");
  G4VPhysicalVolume* phys_wire_vert471 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert471,"phys_wire_vert471",logicWorld,false,0);

  wireNo= 48;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert481 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert481");
  G4VPhysicalVolume* phys_wire_vert481 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert481,"phys_wire_vert481",logicWorld,false,0);

  wireNo= 49;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos);
  G4LogicalVolume* logical_wire_vert491 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert491");
  G4VPhysicalVolume* phys_wire_vert491 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert491,"phys_wire_vert491",logicWorld,false,0);


  //Second MCP
  wireNo= 0;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert02 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert02");
  G4VPhysicalVolume* phys_wire_vert02 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert02,"phys_wire_vert02",logicWorld,false,0);

  wireNo= 1;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert12 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert12");
  G4VPhysicalVolume* phys_wire_vert12 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert12,"phys_wire_vert12",logicWorld,false,0);

  wireNo= 2;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert22 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert22");
  G4VPhysicalVolume* phys_wire_vert22 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert22,"phys_wire_vert22",logicWorld,false,0);

  wireNo= 3;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert32 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert32");
  G4VPhysicalVolume* phys_wire_vert32 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert32,"phys_wire_vert32",logicWorld,false,0);

  wireNo= 4;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert42 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert42");
  G4VPhysicalVolume* phys_wire_vert42 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert42,"phys_wire_vert42",logicWorld,false,0);

  wireNo= 5;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert52 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert52");
  G4VPhysicalVolume* phys_wire_vert52 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert52,"phys_wire_vert52",logicWorld,false,0);

  wireNo= 6;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert62 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert62");
  G4VPhysicalVolume* phys_wire_vert62 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert62,"phys_wire_vert62",logicWorld,false,0);

  wireNo= 7;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert72 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert72");
  G4VPhysicalVolume* phys_wire_vert72 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert72,"phys_wire_vert72",logicWorld,false,0);

  wireNo= 8;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert82 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert82");
  G4VPhysicalVolume* phys_wire_vert82 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert82,"phys_wire_vert82",logicWorld,false,0);

  wireNo= 9;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert92 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert92");
  G4VPhysicalVolume* phys_wire_vert92 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert92,"phys_wire_vert92",logicWorld,false,0);

  wireNo= 10;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert102 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert102");
  G4VPhysicalVolume* phys_wire_vert102 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert102,"phys_wire_vert102",logicWorld,false,0);

  wireNo= 11;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert112 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert112");
  G4VPhysicalVolume* phys_wire_vert112 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert112,"phys_wire_vert112",logicWorld,false,0);

  wireNo= 12;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert122 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert122");
  G4VPhysicalVolume* phys_wire_vert122 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert122,"phys_wire_vert122",logicWorld,false,0);

  wireNo= 13;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert132 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert132");
  G4VPhysicalVolume* phys_wire_vert132 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert132,"phys_wire_vert132",logicWorld,false,0);

  wireNo= 14;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert142 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert142");
  G4VPhysicalVolume* phys_wire_vert142 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert142,"phys_wire_vert142",logicWorld,false,0);

  wireNo= 15;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert152 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert152");
  G4VPhysicalVolume* phys_wire_vert152 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert152,"phys_wire_vert152",logicWorld,false,0);

  wireNo= 16;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert162 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert162");
  G4VPhysicalVolume* phys_wire_vert162 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert162,"phys_wire_vert162",logicWorld,false,0);

  wireNo= 17;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert172 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert172");
  G4VPhysicalVolume* phys_wire_vert172 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert172,"phys_wire_vert172",logicWorld,false,0);

  wireNo= 18;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert182 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert182");
  G4VPhysicalVolume* phys_wire_vert182 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert182,"phys_wire_vert182",logicWorld,false,0);

  wireNo= 19;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert192 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert192");
  G4VPhysicalVolume* phys_wire_vert192 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert192,"phys_wire_vert192",logicWorld,false,0);

  wireNo= 20;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert202 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert202");
  G4VPhysicalVolume* phys_wire_vert202 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert202,"phys_wire_vert202",logicWorld,false,0);

  wireNo= 21;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert212 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert212");
  G4VPhysicalVolume* phys_wire_vert212 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert212,"phys_wire_vert212",logicWorld,false,0);
  wireNo= 22;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert222 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert222");
  G4VPhysicalVolume* phys_wire_vert222 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert222,"phys_wire_vert222",logicWorld,false,0);

  wireNo= 23;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert232 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert232");
  G4VPhysicalVolume* phys_wire_vert232 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert232,"phys_wire_vert232",logicWorld,false,0);

  wireNo= 24;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert242 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert242");
  G4VPhysicalVolume* phys_wire_vert242 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert242,"phys_wire_vert242",logicWorld,false,0);

  wireNo= 25;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert252 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert252");
  G4VPhysicalVolume* phys_wire_vert252 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert252,"phys_wire_vert252",logicWorld,false,0);

  wireNo= 26;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert262 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert262");
  G4VPhysicalVolume* phys_wire_vert262 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert262,"phys_wire_vert262",logicWorld,false,0);

  wireNo= 27;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert272 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert272");
  G4VPhysicalVolume* phys_wire_vert272 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert272,"phys_wire_vert272",logicWorld,false,0);

  wireNo= 28;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert282 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert282");
  G4VPhysicalVolume* phys_wire_vert282 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert282,"phys_wire_vert282",logicWorld,false,0);

  wireNo= 29;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert292 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert292");
  G4VPhysicalVolume* phys_wire_vert292 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert292,"phys_wire_vert292",logicWorld,false,0);

  wireNo= 30;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert302 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert302");
  G4VPhysicalVolume* phys_wire_vert302 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert302,"phys_wire_vert302",logicWorld,false,0);

  wireNo= 31;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert312 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert312");
  G4VPhysicalVolume* phys_wire_vert312 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert312,"phys_wire_vert312",logicWorld,false,0);

  wireNo= 32;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert322 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert322");
  G4VPhysicalVolume* phys_wire_vert322 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert322,"phys_wire_vert322",logicWorld,false,0);

  wireNo= 33;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert332 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert332");
  G4VPhysicalVolume* phys_wire_vert332 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert332,"phys_wire_vert332",logicWorld,false,0);

  wireNo= 34;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert342 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert342");
  G4VPhysicalVolume* phys_wire_vert342 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert342,"phys_wire_vert342",logicWorld,false,0);

  wireNo= 35;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert352 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert352");
  G4VPhysicalVolume* phys_wire_vert352 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert352,"phys_wire_vert352",logicWorld,false,0);

  wireNo= 36;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert362 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert362");
  G4VPhysicalVolume* phys_wire_vert362 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert362,"phys_wire_vert362",logicWorld,false,0);

  wireNo= 37;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert372 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert372");
  G4VPhysicalVolume* phys_wire_vert372 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert372,"phys_wire_vert372",logicWorld,false,0);

  wireNo= 38;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert382 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert382");
  G4VPhysicalVolume* phys_wire_vert382 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert382,"phys_wire_vert382",logicWorld,false,0);

  wireNo= 39;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert392 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert392");
  G4VPhysicalVolume* phys_wire_vert392 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert392,"phys_wire_vert392",logicWorld,false,0);

  wireNo= 40;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert402 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert402");
  G4VPhysicalVolume* phys_wire_vert402 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert402,"phys_wire_vert402",logicWorld,false,0);

  wireNo= 41;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert412 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert412");
  G4VPhysicalVolume* phys_wire_vert412 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert412,"phys_wire_vert412",logicWorld,false,0);

  wireNo= 42;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert422 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert422");
  G4VPhysicalVolume* phys_wire_vert422 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert422,"phys_wire_vert422",logicWorld,false,0);

  wireNo= 43;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert432 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert432");
  G4VPhysicalVolume* phys_wire_vert432 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert432,"phys_wire_vert432",logicWorld,false,0);

  wireNo= 44;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert442 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert442");
  G4VPhysicalVolume* phys_wire_vert442 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert442,"phys_wire_vert442",logicWorld,false,0);

  wireNo= 45;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert452 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert452");
  G4VPhysicalVolume* phys_wire_vert452 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert452,"phys_wire_vert452",logicWorld,false,0);

  wireNo= 46;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert462 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert462");
  G4VPhysicalVolume* phys_wire_vert462 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert462,"phys_wire_vert462",logicWorld,false,0);

  wireNo= 47;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert472 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert472");
  G4VPhysicalVolume* phys_wire_vert472 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert472,"phys_wire_vert472",logicWorld,false,0);

  wireNo= 48;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert482 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert482");
  G4VPhysicalVolume* phys_wire_vert482 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert482,"phys_wire_vert482",logicWorld,false,0);

  wireNo= 49;
  wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+MCP_Sep);
  G4LogicalVolume* logical_wire_vert492 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert492");
  G4VPhysicalVolume* phys_wire_vert492 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert492,"phys_wire_vert492",logicWorld,false,0);

}
