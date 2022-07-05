

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
#include "G4WireSD.hh"
#include <map>
#include <cmath>
#include <string>
#include "G4UserLimits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DetectorConstruction::G4DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fMCP1(0),
  fMCP2(0),
  fDSSSD1(0),
  fDSSSD2(0),
  fscint(0),
  fge(0),
  fafter_deg(0)

  {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DetectorConstruction::~G4DetectorConstruction()
{
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double world_sizeXY = 200*cm, world_sizeZ = 500*cm;


  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4double atomicNumber = 1;
  G4double massofMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  G4Material* vacuum = new G4Material("vaccum",atomicNumber,massofMole, density,kStateGas,temperature,pressure);
  G4Material* world_mat = vacuum;

  G4double foilsizeXY = 5*cm;
  G4double film_pos = -10*cm;
  G4double foil_pos = -9.9*cm;

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

//Here the target is the delta E telescope 300um
G4Material* tar_mat = nist->FindOrBuildMaterial("G4_W");
G4Material* gold_mat = nist->FindOrBuildMaterial("G4_Au");
G4Material* alu_mat = nist->FindOrBuildMaterial("G4_Al");
G4Material* mylar_mat = nist->FindOrBuildMaterial("G4_MYLAR");
G4Material* wire_mat = nist->FindOrBuildMaterial("G4_W");
G4Material* DSSD_mat = nist->FindOrBuildMaterial("G4_Si");

G4ThreeVector scint_pos = G4ThreeVector(0,0,-18*cm);
G4ThreeVector deg_pos = G4ThreeVector(0,0,-15*cm);
G4ThreeVector after_deg_pos = G4ThreeVector(0,0,-14.99*cm);
G4ThreeVector alu1_pos = G4ThreeVector(0, 0, film_pos);
G4ThreeVector mylar1_pos = G4ThreeVector(0, 0, foil_pos);
G4ThreeVector alu2_pos = G4ThreeVector(0, 0, film_pos+150*cm);
G4ThreeVector mylar2_pos = G4ThreeVector(0, 0, foil_pos+150*cm);
G4ThreeVector DSSD1_pos = G4ThreeVector(0, 0, 190*cm);
G4ThreeVector DSSD2_pos = G4ThreeVector(0, 0, 184*cm);
//G4ThreeVector ge_pos = G4ThreeVector(0, 0, 190*cm);
//G4ThreeVector ge_pos = G4ThreeVector(-22*cm, 0, 190*cm);
G4ThreeVector ge_pos = G4ThreeVector(0,0,190*cm);
G4ThreeVector tar_pos = G4ThreeVector(0,0,174*cm);

G4double half_length = 0.5;

G4Box* solid_dep = new G4Box("solid_dep", 40*cm, 40*cm, 0.001*um);
G4Box* solid_tar = new G4Box("solid_tar", 10*cm, 10*cm, 20*um);
G4Box* solid_alu1 = new G4Box("solid_alu1", foilsizeXY, foilsizeXY, half_length*60*nm);
G4Box* solid_mylar1 = new G4Box("solid_mylar1", foilsizeXY, foilsizeXY, half_length*2*um);
G4Box* solid_wires = new G4Box("solid_wires", foilsizeXY, foilsizeXY, half_length*25*um);
G4Box* solid_DSSD1 = new G4Box("solid_DSSD1", 10*cm, 10*cm, half_length*20*um);
G4Box* solid_DSSD2 = new G4Box("solid_DSSD2", 10*cm, 10*cm, half_length*300*um);
G4Box* solid_deg = new G4Box("solid_deg",15*cm,15*cm,half_length*3100*0.7*um);
//G4Box* solid_ge = new G4Box("solid_ge", 150*cm, 150*cm, 0.001*um);
G4Sphere* solid_ge = new G4Sphere("solid_ge", 34*cm, 34.1*cm, 0.,twopi, 0., pi);

//G4LogicalVolume* tar_log = new G4LogicalVolume(solid_tar,gold_mat,"tar_log");
G4LogicalVolume* deg_log = new G4LogicalVolume(solid_deg,alu_mat,"deg_log");
G4LogicalVolume* after_deg_log = new G4LogicalVolume(solid_deg,alu_mat,"after_deg_log");
G4LogicalVolume* alu1_log = new G4LogicalVolume(solid_alu1,alu_mat,"alu1_log");
G4LogicalVolume* mylar1_log = new G4LogicalVolume(solid_mylar1,mylar_mat,"mylar1_log");
G4LogicalVolume* alu2_log = new G4LogicalVolume(solid_alu1,alu_mat,"alu2_log");
G4LogicalVolume* mylar2_log = new G4LogicalVolume(solid_mylar1,mylar_mat,"mylar2_log");
G4LogicalVolume* DSSD1_log = new G4LogicalVolume(solid_DSSD1,DSSD_mat,"DSSD1_log");
G4LogicalVolume* DSSD2_log = new G4LogicalVolume(solid_DSSD2,DSSD_mat,"DSSD2_log");
G4LogicalVolume* ge_log = new G4LogicalVolume(solid_ge,vacuum,"ge_log");
G4LogicalVolume* scint_log = new G4LogicalVolume(solid_dep,vacuum,"scint_log");

G4VPhysicalVolume* deg_phys = new G4PVPlacement(0,deg_pos,deg_log,"deg_phys",logicWorld,false,0,checkOverlaps);
G4VPhysicalVolume* after_deg_phys = new G4PVPlacement(0,after_deg_pos,after_deg_log,"after_deg_phys",logicWorld,false,0,checkOverlaps);
G4VPhysicalVolume* ge_phys = new G4PVPlacement(0,ge_pos,ge_log,"ge_phys",logicWorld,false,0,checkOverlaps);
G4VPhysicalVolume* scint_phys = new G4PVPlacement(0,scint_pos,scint_log,"scint_phys",logicWorld,false,0,checkOverlaps);
//G4VPhysicalVolume* tar_phys = new G4PVPlacement(0,tar_pos,tar_log,"tar_phys",logicWorld,false,9,checkOverlaps);
G4VPhysicalVolume* DSSD1_phys = new G4PVPlacement(0, DSSD1_pos,DSSD1_log, "DSSD1_phys", logicWorld, false,  0, checkOverlaps);
G4VPhysicalVolume* DSSD2_phys = new G4PVPlacement(0, DSSD2_pos,DSSD2_log, "DSSD2_phys", logicWorld, false,  0, checkOverlaps);

G4RotationMatrix* inclination = new G4RotationMatrix();
inclination->rotateX(315.*deg);

G4VPhysicalVolume* alu1_phys =   new G4PVPlacement(inclination,   alu1_pos,  alu1_log,  "alu1_phys", logicWorld, false, 0, checkOverlaps);        //overlaps checking
G4VPhysicalVolume* mylar1_phys = new G4PVPlacement(inclination, mylar1_pos, mylar1_log, "mylar1_phys", logicWorld, false, 0,checkOverlaps);        //overlaps checking
G4VPhysicalVolume* alu2_phys = new G4PVPlacement(inclination,  alu2_pos,   alu2_log, "alu2_phys", logicWorld, false, 0, checkOverlaps);        //overlaps checking
G4VPhysicalVolume* mylar2_phys = new G4PVPlacement(inclination, mylar2_pos,   mylar2_log,  "mylar2_phys", logicWorld,   false, 0, checkOverlaps);        //overlaps checking

//####################Place horizontal wires here ##########################################

G4double seperation = 0.1*cm;
G4double radians = 45*4*atan(1.0)/180.0;
G4Box* solid_wire_hori = new G4Box("solid_wire_hori", foilsizeXY, 0.001*cm, 0.001*cm);
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
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire02 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire02");
G4VPhysicalVolume* phys_wire02 = new G4PVPlacement(0,wire_pos,logical_wire02,"phys_wire02",logicWorld,false,0);

wireNo = 1;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire12 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire12");
G4VPhysicalVolume* phys_wire12 = new G4PVPlacement(0,wire_pos,logical_wire12,"phys_wire12",logicWorld,false,0);

wireNo = 2;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire22 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire22");
G4VPhysicalVolume* phys_wire22 = new G4PVPlacement(0,wire_pos,logical_wire22,"phys_wire22",logicWorld,false,0);

wireNo = 3;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire32 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire32");
G4VPhysicalVolume* phys_wire32 = new G4PVPlacement(0,wire_pos,logical_wire32,"phys_wire32",logicWorld,false,0);

wireNo = 4;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire42 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire42");
G4VPhysicalVolume* phys_wire42 = new G4PVPlacement(0,wire_pos,logical_wire42,"phys_wire42",logicWorld,false,0);

wireNo = 5;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire52 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire52");
G4VPhysicalVolume* phys_wire52 = new G4PVPlacement(0,wire_pos,logical_wire52,"phys_wire52",logicWorld,false,0);

wireNo = 6;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire62 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire62");
G4VPhysicalVolume* phys_wire62 = new G4PVPlacement(0,wire_pos,logical_wire62,"phys_wire62",logicWorld,false,0);

wireNo = 7;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire72 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire72");
G4VPhysicalVolume* phys_wire72 = new G4PVPlacement(0,wire_pos,logical_wire72,"phys_wire72",logicWorld,false,0);

wireNo = 8;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire82 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire82");
G4VPhysicalVolume* phys_wire82 = new G4PVPlacement(0,wire_pos,logical_wire82,"phys_wire82",logicWorld,false,0);

wireNo = 9;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire92 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire92");
G4VPhysicalVolume* phys_wire92 = new G4PVPlacement(0,wire_pos,logical_wire92,"phys_wire92",logicWorld,false,0);

wireNo = 10;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire102 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire102");
G4VPhysicalVolume* phys_wire102 = new G4PVPlacement(0,wire_pos,logical_wire102,"phys_wire102",logicWorld,false,0);

wireNo = 11;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire112 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire112");
G4VPhysicalVolume* phys_wire112 = new G4PVPlacement(0,wire_pos,logical_wire112,"phys_wire112",logicWorld,false,0);

wireNo = 12;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire122 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire122");
G4VPhysicalVolume* phys_wire122 = new G4PVPlacement(0,wire_pos,logical_wire122,"phys_wire122",logicWorld,false,0);

wireNo = 13;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire132 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire132");
G4VPhysicalVolume* phys_wire132 = new G4PVPlacement(0,wire_pos,logical_wire132,"phys_wire132",logicWorld,false,0);

wireNo = 14;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire142 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire142");
G4VPhysicalVolume* phys_wire142 = new G4PVPlacement(0,wire_pos,logical_wire142,"phys_wire142",logicWorld,false,0);

wireNo = 15;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire152 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire152");
G4VPhysicalVolume* phys_wire152 = new G4PVPlacement(0,wire_pos,logical_wire152,"phys_wire152",logicWorld,false,0);

wireNo = 16;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire162 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire162");

G4VPhysicalVolume* phys_wire162 = new G4PVPlacement(0,wire_pos,logical_wire162,"phys_wire162",logicWorld,false,0);

wireNo = 17;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire172 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire172");
G4VPhysicalVolume* phys_wire172 = new G4PVPlacement(0,wire_pos,logical_wire172,"phys_wire172",logicWorld,false,0);

wireNo = 18;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire182 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire182");
G4VPhysicalVolume* phys_wire182 = new G4PVPlacement(0,wire_pos,logical_wire182,"phys_wire182",logicWorld,false,0);

wireNo = 19;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire192 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire192");
G4VPhysicalVolume* phys_wire192 = new G4PVPlacement(0,wire_pos,logical_wire192,"phys_wire192",logicWorld,false,0);

wireNo = 20;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire202 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire202");
G4VPhysicalVolume* phys_wire202 = new G4PVPlacement(0,wire_pos,logical_wire202,"phys_wire202",logicWorld,false,0);

wireNo = 21;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire212 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire212");
G4VPhysicalVolume* phys_wire212 = new G4PVPlacement(0,wire_pos,logical_wire212,"phys_wire212",logicWorld,false,0);

wireNo = 22;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire222 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire222");
G4VPhysicalVolume* phys_wire222 = new G4PVPlacement(0,wire_pos,logical_wire222,"phys_wire222",logicWorld,false,0);

wireNo = 23;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire232 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire232");
G4VPhysicalVolume* phys_wire232 = new G4PVPlacement(0,wire_pos,logical_wire232,"phys_wire232",logicWorld,false,0);

wireNo = 24;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire242 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire242");
G4VPhysicalVolume* phys_wire242 = new G4PVPlacement(0,wire_pos,logical_wire242,"phys_wire242",logicWorld,false,0);

wireNo = 25;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire252 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire252");
G4VPhysicalVolume* phys_wire252 = new G4PVPlacement(0,wire_pos,logical_wire252,"phys_wire252",logicWorld,false,0);

wireNo = 26;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire262 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire262");
G4VPhysicalVolume* phys_wire262 = new G4PVPlacement(0,wire_pos,logical_wire262,"phys_wire262",logicWorld,false,0);

wireNo = 27;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire272 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire272");
G4VPhysicalVolume* phys_wire272 = new G4PVPlacement(0,wire_pos,logical_wire272,"phys_wire272",logicWorld,false,0);

wireNo = 28;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire282 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire282");
G4VPhysicalVolume* phys_wire282 = new G4PVPlacement(0,wire_pos,logical_wire282,"phys_wire282",logicWorld,false,0);

wireNo = 29;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire292 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire292");
G4VPhysicalVolume* phys_wire292 = new G4PVPlacement(0,wire_pos,logical_wire292,"phys_wire292",logicWorld,false,0);

wireNo = 30;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire302 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire302");
G4VPhysicalVolume* phys_wire302 = new G4PVPlacement(0,wire_pos,logical_wire302,"phys_wire302",logicWorld,false,0);

wireNo = 31;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire312 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire312");
G4VPhysicalVolume* phys_wire312 = new G4PVPlacement(0,wire_pos,logical_wire312,"phys_wire312",logicWorld,false,0);

wireNo = 32;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire322 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire322");
G4VPhysicalVolume* phys_wire322 = new G4PVPlacement(0,wire_pos,logical_wire322,"phys_wire322",logicWorld,false,0);

wireNo = 33;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire332 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire332");
G4VPhysicalVolume* phys_wire332 = new G4PVPlacement(0,wire_pos,logical_wire332,"phys_wire332",logicWorld,false,0);

wireNo = 34;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire342 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire342");
G4VPhysicalVolume* phys_wire342 = new G4PVPlacement(0,wire_pos,logical_wire342,"phys_wire342",logicWorld,false,0);

wireNo = 35;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire352 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire352");
G4VPhysicalVolume* phys_wire352 = new G4PVPlacement(0,wire_pos,logical_wire352,"phys_wire352",logicWorld,false,0);

wireNo = 36;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire362 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire362");
G4VPhysicalVolume* phys_wire362 = new G4PVPlacement(0,wire_pos,logical_wire362,"phys_wire362",logicWorld,false,0);

wireNo = 37;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire372 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire372");
G4VPhysicalVolume* phys_wire372 = new G4PVPlacement(0,wire_pos,logical_wire372,"phys_wire372",logicWorld,false,0);

wireNo = 38;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire382 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire382");
G4VPhysicalVolume* phys_wire382 = new G4PVPlacement(0,wire_pos,logical_wire382,"phys_wire382",logicWorld,false,0);

wireNo = 39;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire392 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire392");
G4VPhysicalVolume* phys_wire392 = new G4PVPlacement(0,wire_pos,logical_wire392,"phys_wire392",logicWorld,false,0);

wireNo = 40;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire402 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire402");
G4VPhysicalVolume* phys_wire402 = new G4PVPlacement(0,wire_pos,logical_wire402,"phys_wire402",logicWorld,false,0);

wireNo = 41;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire412 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire412");
G4VPhysicalVolume* phys_wire412 = new G4PVPlacement(0,wire_pos,logical_wire412,"phys_wire412",logicWorld,false,0);

wireNo = 42;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire422 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire422");
G4VPhysicalVolume* phys_wire422 = new G4PVPlacement(0,wire_pos,logical_wire422,"phys_wire422",logicWorld,false,0);

wireNo = 43;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire432 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire432");
G4VPhysicalVolume* phys_wire432 = new G4PVPlacement(0,wire_pos,logical_wire432,"phys_wire432",logicWorld,false,0);

wireNo = 44;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire442 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire442");
G4VPhysicalVolume* phys_wire442 = new G4PVPlacement(0,wire_pos,logical_wire442,"phys_wire442",logicWorld,false,0);

wireNo = 45;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire452 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire452");
G4VPhysicalVolume* phys_wire452 = new G4PVPlacement(0,wire_pos,logical_wire452,"phys_wire452",logicWorld,false,0);

wireNo = 46;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire462 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire462");
G4VPhysicalVolume* phys_wire462 = new G4PVPlacement(0,wire_pos,logical_wire462,"phys_wire462",logicWorld,false,0);

wireNo = 47;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire472 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire472");
G4VPhysicalVolume* phys_wire472 = new G4PVPlacement(0,wire_pos,logical_wire472,"phys_wire472",logicWorld,false,0);

wireNo = 48;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire482 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire482");
G4VPhysicalVolume* phys_wire482 = new G4PVPlacement(0,wire_pos,logical_wire482,"phys_wire482",logicWorld,false,0);

wireNo = 49;
wire_pos = G4ThreeVector(0,
2*wireNo*0.707*foilsizeXY/Nwires-seperation-foilsizeXY*0.707,
2*wireNo*0.707*foilsizeXY/Nwires+foil_pos+seperation-foilsizeXY*0.707+150*cm);
G4LogicalVolume* logical_wire492 = new G4LogicalVolume(solid_wire_hori,wire_mat,"logical_wire492");
G4VPhysicalVolume* phys_wire492 = new G4PVPlacement(0,wire_pos,logical_wire492,"phys_wire492",logicWorld,false,0);



//###############################Place Vertical wires here ######################################


//First MCP
G4Box* solid_wire_vert = new G4Box("solid_wire_vert", 0.001*cm, foilsizeXY, 0.001*cm);

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
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert02 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert02");
G4VPhysicalVolume* phys_wire_vert02 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert02,"phys_wire_vert02",logicWorld,false,0);

wireNo= 1;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert12 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert12");
G4VPhysicalVolume* phys_wire_vert12 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert12,"phys_wire_vert12",logicWorld,false,0);

wireNo= 2;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert22 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert22");
G4VPhysicalVolume* phys_wire_vert22 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert22,"phys_wire_vert22",logicWorld,false,0);

wireNo= 3;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert32 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert32");
G4VPhysicalVolume* phys_wire_vert32 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert32,"phys_wire_vert32",logicWorld,false,0);

wireNo= 4;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert42 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert42");
G4VPhysicalVolume* phys_wire_vert42 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert42,"phys_wire_vert42",logicWorld,false,0);

wireNo= 5;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert52 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert52");
G4VPhysicalVolume* phys_wire_vert52 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert52,"phys_wire_vert52",logicWorld,false,0);

wireNo= 6;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert62 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert62");
G4VPhysicalVolume* phys_wire_vert62 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert62,"phys_wire_vert62",logicWorld,false,0);

wireNo= 7;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert72 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert72");
G4VPhysicalVolume* phys_wire_vert72 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert72,"phys_wire_vert72",logicWorld,false,0);

wireNo= 8;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert82 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert82");
G4VPhysicalVolume* phys_wire_vert82 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert82,"phys_wire_vert82",logicWorld,false,0);

wireNo= 9;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert92 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert92");
G4VPhysicalVolume* phys_wire_vert92 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert92,"phys_wire_vert92",logicWorld,false,0);

wireNo= 10;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert102 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert102");
G4VPhysicalVolume* phys_wire_vert102 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert102,"phys_wire_vert102",logicWorld,false,0);

wireNo= 11;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert112 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert112");
G4VPhysicalVolume* phys_wire_vert112 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert112,"phys_wire_vert112",logicWorld,false,0);

wireNo= 12;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert122 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert122");
G4VPhysicalVolume* phys_wire_vert122 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert122,"phys_wire_vert122",logicWorld,false,0);

wireNo= 13;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert132 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert132");
G4VPhysicalVolume* phys_wire_vert132 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert132,"phys_wire_vert132",logicWorld,false,0);

wireNo= 14;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert142 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert142");
G4VPhysicalVolume* phys_wire_vert142 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert142,"phys_wire_vert142",logicWorld,false,0);

wireNo= 15;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert152 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert152");
G4VPhysicalVolume* phys_wire_vert152 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert152,"phys_wire_vert152",logicWorld,false,0);

wireNo= 16;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert162 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert162");
G4VPhysicalVolume* phys_wire_vert162 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert162,"phys_wire_vert162",logicWorld,false,0);

wireNo= 17;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert172 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert172");
G4VPhysicalVolume* phys_wire_vert172 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert172,"phys_wire_vert172",logicWorld,false,0);

wireNo= 18;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert182 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert182");
G4VPhysicalVolume* phys_wire_vert182 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert182,"phys_wire_vert182",logicWorld,false,0);

wireNo= 19;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert192 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert192");
G4VPhysicalVolume* phys_wire_vert192 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert192,"phys_wire_vert192",logicWorld,false,0);

wireNo= 20;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert202 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert202");
G4VPhysicalVolume* phys_wire_vert202 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert202,"phys_wire_vert202",logicWorld,false,0);

wireNo= 21;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert212 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert212");
G4VPhysicalVolume* phys_wire_vert212 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert212,"phys_wire_vert212",logicWorld,false,0);
wireNo= 22;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert222 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert222");
G4VPhysicalVolume* phys_wire_vert222 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert222,"phys_wire_vert222",logicWorld,false,0);

wireNo= 23;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert232 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert232");
G4VPhysicalVolume* phys_wire_vert232 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert232,"phys_wire_vert232",logicWorld,false,0);

wireNo= 24;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert242 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert242");
G4VPhysicalVolume* phys_wire_vert242 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert242,"phys_wire_vert242",logicWorld,false,0);

wireNo= 25;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert252 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert252");
G4VPhysicalVolume* phys_wire_vert252 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert252,"phys_wire_vert252",logicWorld,false,0);

wireNo= 26;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert262 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert262");
G4VPhysicalVolume* phys_wire_vert262 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert262,"phys_wire_vert262",logicWorld,false,0);

wireNo= 27;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert272 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert272");
G4VPhysicalVolume* phys_wire_vert272 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert272,"phys_wire_vert272",logicWorld,false,0);

wireNo= 28;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert282 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert282");
G4VPhysicalVolume* phys_wire_vert282 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert282,"phys_wire_vert282",logicWorld,false,0);

wireNo= 29;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert292 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert292");
G4VPhysicalVolume* phys_wire_vert292 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert292,"phys_wire_vert292",logicWorld,false,0);

wireNo= 30;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert302 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert302");
G4VPhysicalVolume* phys_wire_vert302 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert302,"phys_wire_vert302",logicWorld,false,0);

wireNo= 31;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert312 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert312");
G4VPhysicalVolume* phys_wire_vert312 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert312,"phys_wire_vert312",logicWorld,false,0);

wireNo= 32;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert322 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert322");
G4VPhysicalVolume* phys_wire_vert322 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert322,"phys_wire_vert322",logicWorld,false,0);

wireNo= 33;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert332 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert332");
G4VPhysicalVolume* phys_wire_vert332 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert332,"phys_wire_vert332",logicWorld,false,0);

wireNo= 34;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert342 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert342");
G4VPhysicalVolume* phys_wire_vert342 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert342,"phys_wire_vert342",logicWorld,false,0);

wireNo= 35;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert352 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert352");
G4VPhysicalVolume* phys_wire_vert352 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert352,"phys_wire_vert352",logicWorld,false,0);

wireNo= 36;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert362 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert362");
G4VPhysicalVolume* phys_wire_vert362 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert362,"phys_wire_vert362",logicWorld,false,0);

wireNo= 37;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert372 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert372");
G4VPhysicalVolume* phys_wire_vert372 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert372,"phys_wire_vert372",logicWorld,false,0);

wireNo= 38;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert382 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert382");
G4VPhysicalVolume* phys_wire_vert382 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert382,"phys_wire_vert382",logicWorld,false,0);

wireNo= 39;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert392 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert392");
G4VPhysicalVolume* phys_wire_vert392 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert392,"phys_wire_vert392",logicWorld,false,0);

wireNo= 40;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert402 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert402");
G4VPhysicalVolume* phys_wire_vert402 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert402,"phys_wire_vert402",logicWorld,false,0);

wireNo= 41;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert412 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert412");
G4VPhysicalVolume* phys_wire_vert412 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert412,"phys_wire_vert412",logicWorld,false,0);

wireNo= 42;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert422 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert422");
G4VPhysicalVolume* phys_wire_vert422 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert422,"phys_wire_vert422",logicWorld,false,0);

wireNo= 43;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert432 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert432");
G4VPhysicalVolume* phys_wire_vert432 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert432,"phys_wire_vert432",logicWorld,false,0);

wireNo= 44;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert442 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert442");
G4VPhysicalVolume* phys_wire_vert442 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert442,"phys_wire_vert442",logicWorld,false,0);

wireNo= 45;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert452 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert452");
G4VPhysicalVolume* phys_wire_vert452 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert452,"phys_wire_vert452",logicWorld,false,0);

wireNo= 46;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert462 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert462");
G4VPhysicalVolume* phys_wire_vert462 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert462,"phys_wire_vert462",logicWorld,false,0);

wireNo= 47;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert472 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert472");
G4VPhysicalVolume* phys_wire_vert472 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert472,"phys_wire_vert472",logicWorld,false,0);

wireNo= 48;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert482 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert482");
G4VPhysicalVolume* phys_wire_vert482 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert482,"phys_wire_vert482",logicWorld,false,0);

wireNo= 49;
wire_pos = G4ThreeVector(-foilsizeXY+2*wireNo*foilsizeXY/Nwires,-seperation,seperation+foil_pos+150*cm);
G4LogicalVolume* logical_wire_vert492 = new G4LogicalVolume(solid_wire_vert,wire_mat,"logical_wire_vert492");
G4VPhysicalVolume* phys_wire_vert492 = new G4PVPlacement(inclination,wire_pos,logical_wire_vert492,"phys_wire_vert492",logicWorld,false,0);

fMCP1 = mylar1_log;
fMCP2 = mylar2_log;
fDSSSD1 = DSSD1_log;
fDSSSD2 = DSSD2_log;
fscint = scint_log;
fge = ge_log;
fafter_deg = after_deg_log;

return physWorld;

}

void G4DetectorConstruction::ConstructSDandField()
{
  //MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

  //fScoringtar->SetSensitiveDetector(sensDet);

}
