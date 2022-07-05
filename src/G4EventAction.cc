#include "G4EventAction.hh"
#include "G4RunAction.hh"
#include "G4Analysis.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include "Randomize.hh"
#include <cmath>
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EventAction::G4EventAction(G4RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EventAction::~G4EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.0;
  fafter_deg = 0;
  fMCP1_x = 0;
  fMCP1_y = 0;
  fMCP1_z = 0;
  fMCP1_time = 0;
  fMCP2_x = 0;
  fMCP2_y = 0;
  fMCP2_z = 0;
  fMCP2_time = 0;
  fDSSSD1_energy = 0;
  fDSSSD1_time = 0;
  fDSSSD2_energy = 0;
  fDSSSD2_time = 0;
  fscint_time = 0;
  fscint_energy = 0;
  fge_energy = 0;
  fge_time = 0;
  fge_x = 0;
  fge_y = 0;
  fge_z = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EventAction::EndOfEventAction(const G4Event* event)
{

  using namespace std;

  //Reconstruction



  double flip = G4UniformRand();
  G4int mass = 208;
  G4double DSSSD1_z = 1880*mm;
  G4double lightspeed = 299.792; //mm/ns;

  //G4double tof = abs(G4RandGauss::shoot(fMCP2_time,0.1)-G4RandGauss::shoot(fMCP1_time,0.1)); //with resolution

  G4ThreeVector ion_direction(fMCP2_x-fMCP1_x,fMCP2_y-fMCP1_y,fMCP2_z-fMCP1_z);
  G4double tof = abs(G4RandGauss::shoot(fMCP2_time,0.1)-G4RandGauss::shoot(fMCP1_time,0.1));
  G4double ion_veloctiy = ion_direction.mag()/tof;
  G4double beta = ion_veloctiy/lightspeed;

  G4ThreeVector beam_dir(0,0,1);
  G4ThreeVector DSSSD1_position(0,0,DSSSD1_z-30*mm);
  G4ThreeVector gamma_direction(fge_x-DSSSD1_position.getX(),fge_y-DSSSD1_position.getY(),fge_z-DSSSD1_position.getZ());


  G4double ion_angle = ion_direction.angle(beam_dir);
  G4double gamma_angle = gamma_direction.angle(beam_dir);
  G4double angle = gamma_direction.angle(ion_direction);

  G4double factor = sqrt(1.0 - beta*beta)/(1-beta*cos(angle));
  G4double ge_e_dopp = fge_energy*factor;

  fRunAction->AddEdep(fEdep);
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->CreateNtuple("HISPEC", "HISPEC");
  analysisManager->FillNtupleDColumn(0, fMCP1_x); //MCP1_x
  analysisManager->FillNtupleDColumn(1, fMCP1_y); //MCP1_y
  analysisManager->FillNtupleDColumn(2, fMCP1_z); //MCP1_z
  analysisManager->FillNtupleDColumn(3, fMCP1_time); //MCP1_timing
  analysisManager->FillNtupleDColumn(4, fMCP2_x); //MCP2_x
  analysisManager->FillNtupleDColumn(5, fMCP2_y); //MCP2_y
  analysisManager->FillNtupleDColumn(6, fMCP2_z); //MCP2_z
  analysisManager->FillNtupleDColumn(7, fMCP2_time); //MCP2_timing
  analysisManager->FillNtupleDColumn(8, tof); //MCP_TOF
  analysisManager->FillNtupleDColumn(9, fDSSSD1_energy); //DSSSD1_energy
  analysisManager->FillNtupleDColumn(10, fDSSSD1_time); //DSSSD1_timing
  analysisManager->FillNtupleDColumn(11, fDSSSD2_energy); //DSSSD2_energy
  analysisManager->FillNtupleDColumn(12, fDSSSD2_time); //DSSSD2_timing
  analysisManager->FillNtupleDColumn(13, fscint_time); //Scintillator_timing
  analysisManager->FillNtupleDColumn(14, fscint_energy); //Scintillator_energy

  if (gamma_angle>0.05) {
  analysisManager->FillNtupleDColumn(15, fge_x); //Ge_x
  analysisManager->FillNtupleDColumn(16, fge_y); //Ge_y
  analysisManager->FillNtupleDColumn(17, fge_z); //Ge_z
  analysisManager->FillNtupleDColumn(18, fge_time); //Ge_timing
  analysisManager->FillNtupleDColumn(19, fge_energy); //Ge_energy
  analysisManager->FillNtupleDColumn(20, ge_e_dopp); //Ge_energy doppler correction
  analysisManager->FillNtupleDColumn(21, angle); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(22, beta); //Ion beta
  if (fge_energy>2) {analysisManager->FillNtupleDColumn(23, gamma_angle);} //Ion beta
  }
  
  analysisManager->FillNtupleDColumn(24, ion_angle); //Ion beta


  analysisManager->AddNtupleRow();



  //analysisManager->FillH1(0,tof); //ns
  //analysisManager->FillH1(1,(10000000*tof/150)/300000000); //ns per cm

  G4int eventID = event->GetEventID();


  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    //
/*
    G4cout << G4BestUnit(fbeam_charge,"Electric charge") << G4endl ;
    G4cout << G4BestUnit(fbeam_mass,"Mass") << G4endl ;
    */
/*

    G4cout << fEnergydeg_beam<< " should be MeV" << G4endl;
    G4cout << fEnergydeg_beam/keV<< " per keV"<<G4endl;
    G4cout << fEnergydeg_beam*keV<< " times cm" << G4endl;

    G4cout << fbeam_x_pos<< G4endl;
    G4cout << fbeam_x_pos/cm<< G4endl;
    G4cout << fbeam_x_pos*cm<< G4endl;
    */
/*
    G4cout
       << "       Enegy deposited in foil" << fEnergyfoil/keV << " [keV]"
       << "       Enegy deposited in target" << fEnergytar/keV << " [keV]"
       << "       Enegy deposited in wires" << fEnergywires/keV << " [keV]"
//       << "       total track length: " << std::setw(7)
//                                        << G4BestUnit(fTrackLLaBr3,"Length")
    << G4endl
    */
  }
}
