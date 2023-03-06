#include "G4EventAction.hh"
#include "G4RunAction.hh"
#include "G4Analysis.hh"
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
{
  //std::ifstream fin2("random_eurica0.dat");
  //std::ifstream fin2("random.dat");
  //std::string line2;
  //while (std::getline(fin2,line2)) {
    //gamma_sample.push_back(line2);
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EventAction::~G4EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.0;
  fMCP1_pos = G4ThreeVector(999999,999999,999999);
  fMCP1_time = 999999;
  fMCP2_pos = G4ThreeVector(999999,999999,999999);
  fMCP2_time = 999999;
  fDSSSD1_energy = 0;
  fDSSSD1_time = 999999;
  fDSSSD2_energy = 0;
  fDSSSD2_time = 999999;


  for (int i = 0; i <12; i++) {
    for (size_t j = 0; j <3; j++) {
      fge_energy[i][j]=0;
      fge_array_pos[i][j] = G4ThreeVector(999999,999999,999999);
    }
  }

  fDSSSD2_pos = G4ThreeVector(999999,999999,999999);
  fge_time = 999999;
  fge_pos = G4ThreeVector(999999,999999,999999);
  fmass = 999999;
  fcharge = 999999;
  ftar_pos = G4ThreeVector(999999,999999,999999);
  ftar_energy = 999999;
  fge_counter = 999999;

  fdeg_energy = 999999;
  fdeg_mass = 999999;
  fdeg_charge = 999999;
  fdeg_pos =  G4ThreeVector(999999,999999,999999);

  fz1_energy = 999999;
  fz1_mass = 999999;
  fz1_charge = 999999;
  fz1_pos =  G4ThreeVector(999999,999999,999999);

  fz2_energy = 999999;
  fz2_mass = 999999;
  fz2_charge = 999999;
  fz2_pos =  G4ThreeVector(999999,999999,999999);

  fz3_energy = 999999;
  fz3_mass = 999999;
  fz3_charge = 999999;
  fz3_pos =  G4ThreeVector(999999,999999,999999);

  fz4_energy = 999999;
  fz4_mass = 999999;
  fz4_charge = 999999;
  fz4_pos =  G4ThreeVector(999999,999999,999999);

  fz5_energy = 999999;
  fz5_mass = 999999;
  fz5_charge = 999999;
  fz5_pos =  G4ThreeVector(999999,999999,999999);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EventAction::EndOfEventAction(const G4Event* event)
{

  fRunAction->AddEdep(fEdep);
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //if (fge_energy<90000 &&fge_time<9000 &&fmass==64 && fcharge==28 && fDSSSD1_energy<90000 && fDSSSD2_energy<90000) {
  analysisManager->CreateNtuple("HISPEC", "HISPEC");

  using namespace std;

  /*
  fge_array_pos[0][1] = G4ThreeVector(-17.62*mm,282.1*mm,1775*mm);
  fge_array_pos[0][2] = G4ThreeVector(-17-48*mm,272.1*mm,1717*mm);
  fge_array_pos[0][3] = G4ThreeVector(-17.75*mm,282.1*mm,1736*mm);

  fge_array_pos[1][1] = G4ThreeVector(-27.58*mm,-286*mm,1702*mm);
  fge_array_pos[1][2] = G4ThreeVector(-27.45*mm,-284.9*mm,1775*mm);
  fge_array_pos[1][3] = G4ThreeVector(-27.75*mm,-285.5*mm,1700*mm);

  fge_array_pos[2][1] = G4ThreeVector(256.1*mm,-36.6*mm,1765*mm);
  fge_array_pos[2][2] = G4ThreeVector(256.2*mm,32.02*mm,1765*mm);
  fge_array_pos[2][3] = G4ThreeVector(256*mm,-35.15*mm,1706*mm);

  fge_array_pos[3][1] = G4ThreeVector(-245.4*mm,-36.74*mm,1720*mm);
  fge_array_pos[3][2] = G4ThreeVector(-245.4*mm,-35.07*mm,1720*mm);
  fge_array_pos[3][3] = G4ThreeVector(-247.1*mm,-34.93*mm,1775*mm);

  fge_array_pos[4][1] = G4ThreeVector(-35.1*mm+,-36.6*mm,2031*mm);
  fge_array_pos[4][2] = G4ThreeVector(-34.93*mm,-32*mm,2035*mm);
  fge_array_pos[4][3] = G4ThreeVector(-35*mm,-35*mm,2031*mm);
  */

/*
 fge_array_pos2[0][0] = G4ThreeVector(227.264,20.8666,1822.29);
 fge_array_pos2[0][1] = G4ThreeVector(217.257,-39.199,1853.45);
 fge_array_pos2[0][2] = G4ThreeVector(206.243,18.414,1888.75);
 fge_array_pos2[1][0] = G4ThreeVector(176.839,-146.083,1822.29);
 fge_array_pos2[1][1] = G4ThreeVector(126.414,-183.085,1853.45);
 fge_array_pos2[1][2] = G4ThreeVector(160.773,-133.43,1888.75);
 fge_array_pos2[2][0] = G4ThreeVector(20.8666,-227.264,1822.29);
 fge_array_pos2[2][1] = G4ThreeVector(-39.199,-217.257,1853.45);
 fge_array_pos2[2][2] = G4ThreeVector(18.414,-206.243,1888.75);
 fge_array_pos2[3][0] = G4ThreeVector(-146.083,-176.839,1822.29);
 fge_array_pos2[3][1] = G4ThreeVector(-183.085,-126.414,1853.45);
 fge_array_pos2[3][2] = G4ThreeVector(-133.43,-160.773,1888.75);
 fge_array_pos2[4][0] = G4ThreeVector(-227.264,-20.8666,1822.29);
 fge_array_pos2[4][1] = G4ThreeVector(-217.257,39.199,1853.45);
 fge_array_pos2[4][2] = G4ThreeVector(-206.243,-18.414,1888.75);
 fge_array_pos2[5][0] = G4ThreeVector(-176.839,146.083,1822.29);
 fge_array_pos2[5][1] = G4ThreeVector(-126.414,183.085,1853.45);
 fge_array_pos2[5][2] = G4ThreeVector(-160.773,133.43,1888.75);
 fge_array_pos2[6][0] = G4ThreeVector(-20.8666,227.264,1822.29);
 fge_array_pos2[6][1] = G4ThreeVector(39.199,217.257,1853.45);
 fge_array_pos2[6][2] = G4ThreeVector(-18.414,206.243,1888.75);
 fge_array_pos2[7][0] = G4ThreeVector(146.083,176.839,1822.29);
 fge_array_pos2[7][1] = G4ThreeVector(183.085,126.414,1853.45);
 fge_array_pos2[7][2] = G4ThreeVector(133.43,160.773,1888.75);
 fge_array_pos2[8][0] = G4ThreeVector(-77.7449,19.7416,2001.03);
 fge_array_pos2[8][1] = G4ThreeVector(-136.519,35.9584,1967.96);
 fge_array_pos2[8][2] = G4ThreeVector(-88.013,86.3578,1979.49);
 fge_array_pos2[9][0] = G4ThreeVector(19.7416,77.7449,2001.03);
 fge_array_pos2[9][1] = G4ThreeVector(35.9584,136.519,1967.96);
 fge_array_pos2[9][2] = G4ThreeVector(86.3578,88.013,1979.49);
 fge_array_pos2[10][0] = G4ThreeVector(77.7449,-19.7416,2001.03);
 fge_array_pos2[10][1] = G4ThreeVector(136.519,-35.9584,1967.96);
 fge_array_pos2[10][2] = G4ThreeVector(88.013,-86.3578,1979.49);
 fge_array_pos2[11][0] = G4ThreeVector(-19.7416,-77.7449,2001.03);
 fge_array_pos2[11][1] = G4ThreeVector(-35.9584,-136.519,1967.96);
 fge_array_pos2[11][2]= G4ThreeVector(-86.3578,-88.013,1979.49);
 */

 dead_rate[0][0] = 1.4e-5;
 dead_rate[0][1] = 1.7e-5;
 dead_rate[0][2] = 2.7e-5;
 dead_rate[1][0] = 2.6e-5;
 dead_rate[1][1] = 3.6e-5;
 dead_rate[1][2] = 1.6e-5;
 dead_rate[2][0] = 2.0e-5;
 dead_rate[2][1] = 4.2e-5;
 dead_rate[2][2] = 3.1e-5;
 dead_rate[3][0] = 2.7e-5;
 dead_rate[3][1] = 4.6e-5;
 dead_rate[3][2] = 5.4e-5;
 dead_rate[4][0] = 2.8e-5;
 dead_rate[4][1] = 3.7e-5;
 dead_rate[4][2] = 3.8e-5;
 dead_rate[5][0] = 3.2e-5;
 dead_rate[5][1] = 3.5e-5;
 dead_rate[5][2] = 2.8e-5;
 dead_rate[6][0] = 2.0e-5;
 dead_rate[6][1] = 3.8e-5;
 dead_rate[6][2] = 3.8e-5;
 dead_rate[7][0] = 3.4e-5;
 dead_rate[7][1] = 4.9e-6;
 dead_rate[7][2] = 4.8e-5;
 dead_rate[8][0] = 3.7e-5;
 dead_rate[8][1] = 1.01e-5;
 dead_rate[8][2] = 1.03e-5;
 dead_rate[9][0] = 3.5e-5;
 dead_rate[9][1] = 1.23e-5;
 dead_rate[9][2] = 1.39e-5;
 dead_rate[10][0] = 1.29e-4;
 dead_rate[10][1] = 2.0e-4;
 dead_rate[10][2] = 4.73e-4;
 dead_rate[11][0] = 2e-4;
 dead_rate[11][1] = 1.5e-4;
 dead_rate[11][2] = 1.72e-4;

  G4double ge_energy_all = 0;
  G4double max = 0;
  G4double top = 0;
  G4double ener = 0;

  for (int i = 0; i <12; i++) {
    ener = 0;
    for (int j = 0; j <3; j++) {


      /*
        if (G4UniformRand()<dead_rate[i][j] && i<10) {
          G4double gam_energy = std::stod(gamma_sample.at(rand()%(9999999)));
          fge_time = abs(G4RandGauss::shoot(0,5));
          fge_energy[i][j]=gam_energy/1000; //G4cout<<"dead crystal "<<i<<" "<<j<<G4endl;
        }
        */
        /*
        if (fge_energy[i][j]>10/1000) {
          analysisManager->FillNtupleDColumn(57,i+20*j);
          G4cout<<"dead_crystal"<<" "<<i<<" "<<j<<G4endl;
        }
        */
        ener += fge_energy[i][j];

        if (fge_energy[i][j]>max) {
          max = fge_energy[i][j];
          fge_pos = fge_array_pos[i][j];
        }


    }


    if (ener>top) {
      top = ener;
      ge_energy_all = ener;

    }

  }

  fge_time = G4RandGauss::shoot(fge_time,20/2.355);


  ge_energy_all = G4RandGauss::shoot(ge_energy_all,2.35/1000);
  fMCP1_pos = G4ThreeVector(G4RandGauss::shoot(fMCP1_pos.getX(),1),G4RandGauss::shoot(fMCP1_pos.getY(),1),G4RandGauss::shoot(fMCP1_pos.getZ(),1));
  fMCP2_pos = G4ThreeVector(G4RandGauss::shoot(fMCP2_pos.getX(),1),G4RandGauss::shoot(fMCP2_pos.getY(),1),G4RandGauss::shoot(fMCP2_pos.getZ(),1));
  fMCP1_time = G4RandGauss::shoot(fMCP1_time,120/1000);
  fMCP2_time = G4RandGauss::shoot(fMCP2_time,120/1000);

  //Set Constants
  double flip = G4UniformRand();
  G4double MCP_Sep = 150*cm;
  G4double target_z = MCP_Sep+24*cm;
  G4double DSSSD1_z = MCP_Sep+34*cm;
  G4double lightspeed = 299.792; //mm/ns;
  G4double amu = 1.6605e-27;
  G4double e = 1.6602e-19;

  //Obtain betas
  G4double mass_kg = fmass*amu;
  G4double rest_mass_MeV = fmass*mass_kg*2.99792e8*2.99792e8/(e*10e6);
  //G4double post_tar_E = fDSSSD1_energy+fDSSSD2_energy;
  //G4double beta_dsssd = sqrt(1-(1/pow(post_tar_E/rest_mass_MeV+1,2)));
  //G4double beta_tar = sqrt(1-(1/pow(ftar_energy/rest_mass_MeV+1,2)));

  G4ThreeVector ion_direction = fMCP2_pos-fMCP1_pos;
  G4double tof = fMCP2_time-fMCP1_time;
  G4double ion_veloctiy = ion_direction.mag()/tof;
  G4double beta_TOF = ion_veloctiy/lightspeed;

  //Obtain angles
  G4double target_step = (ftar_pos.getZ()-fMCP1_pos.getZ())/(fMCP2_pos.getZ()-fMCP1_pos.getZ());
  G4ThreeVector target_ion_pos = target_step*(fMCP2_pos-fMCP1_pos)+fMCP1_pos;
  G4ThreeVector ion_scatter_dir = fDSSSD2_pos-target_ion_pos;
  ion_scatter_dir.unit();
  G4ThreeVector gamma_direction = fge_pos-target_ion_pos;
  gamma_direction.unit();
  //G4ThreeVector gamma_direction = fge_pos;
  G4double angle = gamma_direction.angle(ion_scatter_dir);


  std::ifstream fin("doppler_out.dat",std::ios_base::app);
  std::string line;

  std::getline(fin,line);

  std::vector<double> properties;
  std::string dummy;


  for (int i = 0; i < line.length(); i++){
    dummy+=line.at(i);
    if (line.at(i)==' ' || i == line.length()-1 ) {
      properties.push_back(std::stod(dummy));
      dummy.clear();
    }
  }


  G4ThreeVector coulex_pos = G4ThreeVector(properties.at(0),properties.at(1),properties.at(2)); //Coulex position
  G4double coulex_angle = properties.at(3); //Coulex angle
  G4double coulex_beta = properties.at(4); //Coulex beta
  G4double coulex_energy = properties.at(5); //original coulex energy i.e excitation energy
  G4double shifted_energy = properties.at(6); //energy gamma shifted to
  G4double coulex_factor = properties.at(7); //doppler factor used to shift original gamma

  //G4cout<<"all "<<ge_energy_all<<G4endl;
  //Perform perfect reconstruction
  //G4double perfect_energy = ge_energy_all*coulex_factor;

  //G4double reso_Ge = 2.0/1000*CLHEP::MeV;

  //ge_energy_all = G4RandGauss::shoot(ge_energy_all,reso_Ge);

  //perform imperfect dopplershift
  coulex_factor = sqrt(1.0 - coulex_beta*coulex_beta)/(1-coulex_beta*cos(coulex_angle));
  coulex_energy = ge_energy_all*coulex_factor;
  G4double guild_factor = sqrt(1.0 - beta_TOF*beta_TOF)/(1-beta_TOF*cos(coulex_angle));
  G4double factor = sqrt(1.0 - coulex_beta*coulex_beta)/(1-coulex_beta*cos(angle));
  coulex_factor = sqrt(1.0 - beta_TOF*beta_TOF)/(1-beta_TOF*cos(angle));
  G4double perfect_energy = ge_energy_all*coulex_factor;
  G4double ge_e_dopp = ge_energy_all*factor;
  G4double guild = ge_energy_all*guild_factor;

  //G4cout<<" doppler correct "<<ge_e_dopp<<G4endl;
  fRunAction->AddEdep(fEdep);

  //Reconstruct with perfect variables
  G4ThreeVector coulex_incoming_direction = G4ThreeVector(properties.at(14),properties.at(15),properties.at(16)); //Direction of scattered ion
  G4ThreeVector coulex_ion_direction = G4ThreeVector(properties.at(11),properties.at(12),properties.at(13)); //Direction of scattered ion
  G4ThreeVector coulex_direction = G4ThreeVector(properties.at(8),properties.at(9),properties.at(10)); //Direction of generated gamma
  coulex_incoming_direction.unit();
  coulex_ion_direction.unit();
  coulex_direction.unit();

  G4ThreeVector beam_dir = G4ThreeVector(0,0,1);

  G4double lab_scattering_ion = beam_dir.angle(ion_scatter_dir);
  //lab_scattering_ion.unit();

  G4double scattering = coulex_ion_direction.angle(ion_direction);
  G4double scattering_lab = ion_scatter_dir.angle(ion_direction);

  if (fMCP1_pos.getX()<9999 &&fDSSSD2_pos.getX()<9999 && shifted_energy>0.5 ) {
  //G4cout<<" scattering lab "<<scattering_lab<<" ion_scatter_dir "<<ion_scatter_dir.getX()<<","<<ion_scatter_dir.getY()<<","<<ion_scatter_dir.getZ()<<" ion_direction "<<ion_direction.getX()<<","<<ion_direction.getY()<<","<<ion_direction.getZ()<<G4endl;
}
  /*
  G4cout<<"coulex_pos("<<coulex_pos.getX()<<","<<coulex_pos.getX()<<","<<coulex_pos.getX()<<")->gamma_detection("<<fge_pos.getX()<<","<<fge_pos.getY()<<","<<fge_pos.getZ()<<")"<<G4endl;
  G4cout<<"coulex_dir = ("<<coulex_direction.getX()<<","<<coulex_direction.getY()<<","<<coulex_direction.getZ()<<") & reconstruction_dir("<<gamma_direction.getX()<<","<<gamma_direction.getY()<<","<<gamma_direction.getZ()<<") angular diff = "<<coulex_direction.angle(gamma_direction)*57.29<<G4endl;
  G4cout<<G4endl;
  G4cout<<"coulex_pos("<<coulex_pos.getX()<<","<<coulex_pos.getX()<<","<<coulex_pos.getX()<<")->DSSSD2_pos("<<fDSSSD2_pos.getX()<<","<<fDSSSD2_pos.getY()<<","<<fDSSSD2_pos.getZ()<<")"<<G4endl;
  G4cout<<"coulex_ion_dir = ("<<coulex_ion_direction.getX()<<","<<coulex_ion_direction.getY()<<","<<coulex_ion_direction.getZ()<<") & reconstruction_dir_ion ("<<ion_scatter_dir.getX()<<","<<ion_scatter_dir.getY()<<","<<ion_scatter_dir.getZ()<<") angular diff = "<<coulex_ion_direction.angle(ion_scatter_dir)*57.29<<G4endl;
  G4cout<<G4endl;
  G4cout<<"reco_angle = "<<gamma_direction.angle(ion_scatter_dir)*57.29<<" real_angle = "<<coulex_direction.angle(coulex_ion_direction)*57.29<<G4endl;
  */

  //real coulex direction vs dir between ge and real position
  //G4cout<<"gamma_dir_reco("<<gamma_direction.getX()<<","<<gamma_direction.getY()<<","<<gamma_direction.getZ()<<"),ion_dir_reco("<<ion_scatter_dir.getX()<<","<<ion_scatter_dir.getY()<<","<<ion_scatter_dir.getZ()<<"),angle_reco="<<angle<<G4endl;
  //G4cout<<"target_ion_pos_reco("<<target_ion_pos.getX()<<","<<target_ion_pos.getY()<<","<<target_ion_pos.getZ()<<"),"<<"angle="<<angle<<G4endl;

  //Is safe coulex?
  int Zp = 28; int Zt = 79;
  int Ap = 64; int At = 197;

  double dist_app = 1.25*(pow(64,0.3333)+pow(197,0.333333))+5;
  int safe_coulex = 0;
  if (dist_app<0.72*(Zp*Zt)*(Ap+At)/(fDSSSD2_energy*At)*(1+1/(sin(scattering_lab/2)))) {
    safe_coulex = 1;
  }

  analysisManager->FillNtupleDColumn(0, fMCP1_pos.getX()); //MCP1_x
  analysisManager->FillNtupleDColumn(1, fMCP1_pos.getY()); //MCP1_y
  analysisManager->FillNtupleDColumn(2, fMCP1_pos.getZ()); //MCP1_z
  analysisManager->FillNtupleDColumn(3, fMCP1_time); //MCP1_timing
  analysisManager->FillNtupleDColumn(4, fMCP2_pos.getX()); //MCP2_x
  analysisManager->FillNtupleDColumn(5, fMCP2_pos.getY()); //MCP2_y
  analysisManager->FillNtupleDColumn(6, fMCP2_pos.getZ()); //MCP2_z
  analysisManager->FillNtupleDColumn(7, fMCP2_time); //MCP2_timing
  analysisManager->FillNtupleDColumn(8, tof); //MCP_TOF
  analysisManager->FillNtupleDColumn(9, fDSSSD1_energy); //DSSSD1_energy
  analysisManager->FillNtupleDColumn(10, fDSSSD1_time); //DSSSD1_timing
  analysisManager->FillNtupleDColumn(11, fDSSSD2_energy); //DSSSD2_energy
  analysisManager->FillNtupleDColumn(12, fDSSSD2_time); //DSSSD2_timing
  analysisManager->FillNtupleDColumn(13, fge_pos.getX()); //Ge_x
  analysisManager->FillNtupleDColumn(14, fge_pos.getY()); //Ge_y
  analysisManager->FillNtupleDColumn(15, fge_pos.getZ()); //Ge_z
  analysisManager->FillNtupleDColumn(16, fge_time); //Ge_timing
  analysisManager->FillNtupleDColumn(17, ge_energy_all); //Ge_energy observed
    //G4cout<<ge_energy_all*1000<<G4endl;
  if (ge_e_dopp!=0) {
    analysisManager->FillNtupleDColumn(18, ge_e_dopp); //Ge_energy doppler corrected
  }

  if (fDSSSD2_pos.getX()<9999 && fMCP1_pos.getX()<9999 && fMCP2_pos.getX()<9999 && fge_pos.getX()<9999 && ge_e_dopp>1 && shifted_energy>0.5) {
    G4cout<<"doppler "<<ge_e_dopp*1000<<" keV"<<G4endl;
  }

  analysisManager->FillNtupleDColumn(35, ftar_energy); //Gamma emmision angle

  analysisManager->FillNtupleDColumn(19, beta_TOF); //beta from TOF MCPs
  analysisManager->FillNtupleDColumn(20, angle); //angle using MCPs and DSSSDs

  analysisManager->FillNtupleDColumn(21, fmass);
  analysisManager->FillNtupleDColumn(22, fcharge);

  analysisManager->FillNtupleDColumn(23,fDSSSD2_pos.getX());
  analysisManager->FillNtupleDColumn(24,fDSSSD2_pos.getY());
  analysisManager->FillNtupleDColumn(25,fDSSSD2_pos.getZ());

  analysisManager->FillNtupleDColumn(26, coulex_energy); //excitation Energy
  //analysisManager->FillNtupleDColumn(26, fge_counter); //excitation Energy
  if (shifted_energy>0.5) {analysisManager->FillNtupleDColumn(27, shifted_energy);} //energy gamma shifted too
  analysisManager->FillNtupleDColumn(28, perfect_energy); //energy gamma shifted too
  analysisManager->FillNtupleDColumn(29, coulex_beta); //Beta from shifting step
  analysisManager->FillNtupleDColumn(30, coulex_angle); //angle from shifting step
  //analysisManager->FillNtupleDColumn(31, beta_dsssd); //Beta reoncstruced from DSSSD energy deposition

  analysisManager->FillNtupleDColumn(32, ftar_pos.getX()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(33, ftar_pos.getY()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(34, ftar_pos.getZ()); //Gamma emmision angle

  //nalysisManager->FillNtupleDColumn(36, beta_tar); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(37, angle-coulex_angle); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(38, beta_TOF-coulex_beta); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(39, factor-coulex_factor); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(40, target_ion_pos.getX()-ftar_pos.getX()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(41, target_ion_pos.getY()-ftar_pos.getY()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(42, target_ion_pos.getZ()-ftar_pos.getZ()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(43, coulex_pos.getX()-target_ion_pos.getX()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(44, coulex_pos.getY()-target_ion_pos.getY()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(45, coulex_pos.getZ()-target_ion_pos.getZ()); //Gamma emmision angle
  analysisManager->FillNtupleDColumn(46, 57.29*(coulex_ion_direction.angle(ion_scatter_dir))); //real ion direction vs dir beteen DSSSD and target
  analysisManager->FillNtupleDColumn(47, 57.29*(coulex_direction.angle(gamma_direction))); //real gamma direction vs dir between ge and target
  analysisManager->FillNtupleDColumn(48, 57.29*(coulex_angle-coulex_direction.angle(coulex_ion_direction))); //reco coulex angle vs read out
  analysisManager->FillNtupleDColumn(49, 57.29*scattering); //reco coulex angle vs read out
  analysisManager->FillNtupleDColumn(50, 57.29*scattering_lab);

  G4double one = properties.at(14);
  G4double two = properties.at(15);
  G4double three = properties.at(16);

  analysisManager->FillNtupleDColumn(51, one);
  analysisManager->FillNtupleDColumn(52, two);
  analysisManager->FillNtupleDColumn(53, three);

  analysisManager->FillNtupleDColumn(54,57.29*lab_scattering_ion);
  //analysisManager->FillNtupleDColumn(54,57.29*ftar_energy);
  analysisManager->FillNtupleDColumn(58,safe_coulex);
  analysisManager->FillNtupleDColumn(56,guild);

  analysisManager->FillNtupleDColumn(59,fdeg_pos.getX());
  analysisManager->FillNtupleDColumn(60,fdeg_pos.getY());
  analysisManager->FillNtupleDColumn(61,fdeg_pos.getZ());
  analysisManager->FillNtupleDColumn(62,fdeg_energy);
  analysisManager->FillNtupleDColumn(63,fdeg_mass);
  analysisManager->FillNtupleDColumn(64,fdeg_charge);

  analysisManager->FillNtupleDColumn(65,fz1_pos.getX());
  analysisManager->FillNtupleDColumn(66,fz1_pos.getY());
  analysisManager->FillNtupleDColumn(67,fz1_pos.getZ());
  analysisManager->FillNtupleDColumn(68,fz1_energy);
  analysisManager->FillNtupleDColumn(69,fz1_mass);
  analysisManager->FillNtupleDColumn(70,fz1_charge);

  analysisManager->FillNtupleDColumn(71,fz2_pos.getX());
  analysisManager->FillNtupleDColumn(72,fz2_pos.getY());
  analysisManager->FillNtupleDColumn(73,fz2_pos.getZ());
  analysisManager->FillNtupleDColumn(74,fz2_energy);
  analysisManager->FillNtupleDColumn(75,fz2_mass);
  analysisManager->FillNtupleDColumn(76,fz2_charge);

  analysisManager->FillNtupleDColumn(77,fz3_pos.getX());
  analysisManager->FillNtupleDColumn(78,fz3_pos.getY());
  analysisManager->FillNtupleDColumn(79,fz3_pos.getZ());
  analysisManager->FillNtupleDColumn(80,fz3_energy);
  analysisManager->FillNtupleDColumn(81,fz3_mass);
  analysisManager->FillNtupleDColumn(82,fz3_charge);

  analysisManager->FillNtupleDColumn(83,fz4_pos.getX());
  analysisManager->FillNtupleDColumn(84,fz4_pos.getY());
  analysisManager->FillNtupleDColumn(85,fz4_pos.getZ());
  analysisManager->FillNtupleDColumn(86,fz4_energy);
  analysisManager->FillNtupleDColumn(87,fz4_mass);
  analysisManager->FillNtupleDColumn(88,fz4_charge);

  analysisManager->FillNtupleDColumn(89,fz5_pos.getX());
  analysisManager->FillNtupleDColumn(90,fz5_pos.getY());
  analysisManager->FillNtupleDColumn(91,fz5_pos.getZ());
  analysisManager->FillNtupleDColumn(92,fz5_energy);
  analysisManager->FillNtupleDColumn(93,fz5_mass);
  analysisManager->FillNtupleDColumn(94,fz5_charge);

  analysisManager->AddNtupleRow();

  //}




  std::ofstream doppler_out ("doppler_out.dat");
  doppler_out<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0;



  //G4int eventID = event->GetEventID();

  //G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
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
  //


}

bool G4EventAction::ran_dead(int m, int n) {

}
