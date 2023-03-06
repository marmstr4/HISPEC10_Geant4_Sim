#ifndef G4EventAction_h
#define G4EventAction_h 1
#include <vector>

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
class G4RunAction;

/// Event action class
///

class G4EventAction : public G4UserEventAction
{
  public:
    G4EventAction(G4RunAction* runAction);
    virtual ~G4EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) { fEdep += edep; }
    ///void Addfafter_deg(G4double energy);
    void AddfMCP1_pos(G4ThreeVector position);
    void AddfMCP1_time(G4double time);
    void AddfMCP2_pos(G4ThreeVector position);
    void AddfMCP2_time(G4double time);
    void AddfDSSSD1_time(G4double time);
    void AddfDSSSD1_energy(G4double energy);
    void AddfDSSSD2_time(G4double time);
    void AddfDSSSD2_energy(G4double energy);
    void AddfDSSSD2_pos(G4ThreeVector position);
    void Addfdeg_energy(G4double energy);
    void Addfdeg_mass(G4double mass);
    void Addfdeg_charge(G4double charge);
    void Addfdeg_pos(G4ThreeVector position);

    void Addfz1_energy(G4double energy);
    void Addfz1_mass(G4double mass);
    void Addfz1_charge(G4double charge);
    void Addfz1_pos(G4ThreeVector position);

    void Addfz2_energy(G4double energy);
    void Addfz2_mass(G4double mass);
    void Addfz2_charge(G4double charge);
    void Addfz2_pos(G4ThreeVector position);

    void Addfz3_energy(G4double energy);
    void Addfz3_mass(G4double mass);
    void Addfz3_charge(G4double charge);
    void Addfz3_pos(G4ThreeVector position);

    void Addfz4_energy(G4double energy);
    void Addfz4_mass(G4double mass);
    void Addfz4_charge(G4double charge);
    void Addfz4_pos(G4ThreeVector position);

    void Addfz5_energy(G4double energy);
    void Addfz5_mass(G4double mass);
    void Addfz5_charge(G4double charge);
    void Addfz5_pos(G4ThreeVector position);

    //void Addfscint_time(G4double time);
    //void Addfscint_energy(G4double energy);
    bool ran_dead(int i, int j);
    void Addfge_time(G4double time);
    void Addfge_energy(G4double energy);
    void Addfge_pos(G4ThreeVector position);
    void Addfmass(G4double mass);
    void Addfcharge(G4double charge);
    void Addftar_pos(G4ThreeVector position);
    void Addftar_energy(G4double energy);
    G4double fge_energy[12][3];
    G4ThreeVector fge_array_pos[12][3];
    G4ThreeVector fge_array_pos2[12][3];
    G4int id_gam = 0;
    G4int id_elc = 0;
    G4double dead_rate[12][3];
    //std::vector<std::string> gamma_sample;


  private:
    G4RunAction* fRunAction;
    G4double     fEdep;
    G4ThreeVector fMCP1_pos;
    G4double fMCP1_time;
    G4ThreeVector fMCP2_pos;
    G4double fMCP2_time;
    G4double fDSSSD1_time;
    G4double fDSSSD1_energy;
    G4double fDSSSD2_time;
    G4double fDSSSD2_energy;
    G4ThreeVector fDSSSD2_pos;
    //G4double fscint_time;
    //G4double fscint_energy;
    G4double fge_time;
    G4double fge_counter;

  //  G4double fafter_deg;
    G4ThreeVector fge_pos;
    G4double fmass;
    G4double fcharge;
    G4ThreeVector ftar_pos;
    G4double ftar_energy;
    G4double fdeg_mass;
    G4double fdeg_charge;
    G4ThreeVector fdeg_pos;
    G4double fdeg_energy;

    G4double fz1_energy;
    G4double fz1_mass;
    G4double fz1_charge;
    G4ThreeVector fz1_pos;
    G4double fz2_energy;
    G4double fz2_mass;
    G4double fz2_charge;
    G4ThreeVector fz2_pos;
    G4double fz3_energy;
    G4double fz3_mass;
    G4double fz3_charge;
    G4ThreeVector fz3_pos;
    G4double fz4_energy;
    G4double fz4_mass;
    G4double fz4_charge;
    G4ThreeVector fz4_pos;
    G4double fz5_energy;
    G4double fz5_mass;
    G4double fz5_charge;
    G4ThreeVector fz5_pos;



};

inline void G4EventAction::AddfMCP1_pos(G4ThreeVector position) {
  fMCP1_pos = position;
}

inline void G4EventAction::AddfMCP1_time(G4double time) {
  fMCP1_time = time;
}

inline void G4EventAction::AddfMCP2_pos(G4ThreeVector position) {
  fMCP2_pos = position;
}

inline void G4EventAction::AddfMCP2_time(G4double time) {
  fMCP2_time = time;
}

inline void G4EventAction::AddfDSSSD1_energy(G4double energy) {
  fDSSSD1_energy += energy;
}

inline void G4EventAction::AddfDSSSD1_time(G4double time) {
  fDSSSD1_time = time;
}

inline void G4EventAction::AddfDSSSD2_energy(G4double energy) {
  fDSSSD2_energy += energy;
}

inline void G4EventAction::AddfDSSSD2_pos(G4ThreeVector position) {
  fDSSSD2_pos = position;
}

inline void G4EventAction::AddfDSSSD2_time(G4double time) {
  fDSSSD2_time = time;
}

//inline void G4EventAction::Addfscint_time(G4double time) {
  //fscint_time = time;
//}

inline void G4EventAction::Addfdeg_energy(G4double energy) {
  fdeg_energy = energy;
}

inline void G4EventAction::Addfdeg_mass(G4double mass) {
  fdeg_mass = mass;
}

inline void G4EventAction::Addfdeg_charge(G4double charge) {
  fdeg_charge = charge;
}

inline void G4EventAction::Addfdeg_pos(G4ThreeVector position) {
  fdeg_pos = position;
}

inline void G4EventAction::Addfz1_energy(G4double energy) {
  fz1_energy = energy;
}

inline void G4EventAction::Addfz1_mass(G4double mass) {
  fz1_mass = mass;
}

inline void G4EventAction::Addfz1_charge(G4double charge) {
  fz1_charge = charge;
}

inline void G4EventAction::Addfz1_pos(G4ThreeVector position) {
  fz1_pos = position;
}

inline void G4EventAction::Addfz2_energy(G4double energy) {
  fz2_energy = energy;
}

inline void G4EventAction::Addfz2_mass(G4double mass) {
  fz2_mass = mass;
}

inline void G4EventAction::Addfz2_charge(G4double charge) {
  fz2_charge = charge;
}

inline void G4EventAction::Addfz2_pos(G4ThreeVector position) {
  fz2_pos = position;
}

inline void G4EventAction::Addfz3_energy(G4double energy) {
  fz3_energy = energy;
}

inline void G4EventAction::Addfz3_mass(G4double mass) {
  fz3_mass = mass;
}

inline void G4EventAction::Addfz3_charge(G4double charge) {
  fz3_charge = charge;
}

inline void G4EventAction::Addfz3_pos(G4ThreeVector position) {
  fz3_pos = position;
}

inline void G4EventAction::Addfz4_energy(G4double energy) {
  fz4_energy = energy;
}

inline void G4EventAction::Addfz4_mass(G4double mass) {
  fz4_mass = mass;
}

inline void G4EventAction::Addfz4_charge(G4double charge) {
  fz4_charge = charge;
}

inline void G4EventAction::Addfz4_pos(G4ThreeVector position) {
  fz4_pos = position;
}

inline void G4EventAction::Addfz5_energy(G4double energy) {
  fz5_energy = energy;
}

inline void G4EventAction::Addfz5_mass(G4double mass) {
  fz5_mass = mass;
}

inline void G4EventAction::Addfz5_charge(G4double charge) {
  fz5_charge = charge;
}

inline void G4EventAction::Addfz5_pos(G4ThreeVector position) {
  fz5_pos = position;
}



inline void G4EventAction::Addfge_energy(G4double energy) {
  fge_counter = energy;
}

inline void G4EventAction::Addfge_time(G4double time) {
  fge_time = time;
}

inline void G4EventAction::Addfge_pos(G4ThreeVector position) {
  fge_pos = position;
}

inline void G4EventAction::Addfcharge(G4double charge) {
  fcharge = charge;
}

inline void G4EventAction::Addfmass(G4double mass) {
  fmass = mass;
}

inline void G4EventAction::Addftar_pos(G4ThreeVector position) {
  ftar_pos = position;
}

inline void G4EventAction::Addftar_energy(G4double energy) {
  ftar_energy = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
