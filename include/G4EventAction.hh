#ifndef G4EventAction_h
#define G4EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

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
    void Addfafter_deg(G4double energy);
    void AddfMCP1_x(G4double x);
    void AddfMCP1_y(G4double y);
    void AddfMCP1_z(G4double z);
    void AddfMCP1_time(G4double time);
    void AddfMCP2_x(G4double x);
    void AddfMCP2_y(G4double y);
    void AddfMCP2_z(G4double z);
    void AddfMCP2_time(G4double time);
    void AddfDSSSD1_time(G4double time);
    void AddfDSSSD1_energy(G4double energy);
    void AddfDSSSD2_time(G4double time);
    void AddfDSSSD2_energy(G4double energy);
    void Addfscint_time(G4double time);
    void Addfscint_energy(G4double energy);
    void Addfge_time(G4double time);
    void Addfge_energy(G4double energy);
    void Addfge_x(G4double x);
    void Addfge_y(G4double y);
    void Addfge_z(G4double y);


  private:
    G4RunAction* fRunAction;
    G4double     fEdep;
    G4double fMCP1_x;
    G4double fMCP1_y;
    G4double fMCP1_z;
    G4double fMCP1_time;
    G4double fMCP2_x;
    G4double fMCP2_y;
    G4double fMCP2_z;
    G4double fMCP2_time;
    G4double fDSSSD1_time;
    G4double fDSSSD1_energy;
    G4double fDSSSD2_time;
    G4double fDSSSD2_energy;
    G4double fscint_time;
    G4double fscint_energy;
    G4double fge_time;
    G4double fge_energy;
    G4double fafter_deg;
    G4double fge_x;
    G4double fge_y;
    G4double fge_z;

};

inline void G4EventAction::AddfMCP1_x(G4double x) {
  fMCP1_x = x;
}

inline void G4EventAction::AddfMCP1_y(G4double y) {
  fMCP1_y = y;
}

inline void G4EventAction::AddfMCP1_z(G4double z) {
  fMCP1_z = z;
}

inline void G4EventAction::AddfMCP1_time(G4double time) {
  fMCP1_time = time;
}

inline void G4EventAction::AddfMCP2_x(G4double x) {
  fMCP2_x = x;
}

inline void G4EventAction::AddfMCP2_y(G4double y) {
  fMCP2_y = y;
}

inline void G4EventAction::AddfMCP2_z(G4double z) {
  fMCP2_z = z;
}

inline void G4EventAction::AddfMCP2_time(G4double time) {
  fMCP2_time = time;
}

inline void G4EventAction::AddfDSSSD1_energy(G4double energy) {
  fDSSSD1_energy = energy;
}

inline void G4EventAction::AddfDSSSD1_time(G4double time) {
  fDSSSD1_time = time;
}

inline void G4EventAction::AddfDSSSD2_energy(G4double energy) {
  fDSSSD2_energy = energy;
}

inline void G4EventAction::AddfDSSSD2_time(G4double time) {
  fDSSSD2_time = time;
}

inline void G4EventAction::Addfscint_time(G4double time) {
  fscint_time = time;
}

inline void G4EventAction::Addfscint_energy(G4double energy) {
  fscint_energy = energy;
}

inline void G4EventAction::Addfge_energy(G4double energy) {
  fge_energy = energy;
}

inline void G4EventAction::Addfge_time(G4double time) {
  fge_time = time;
}

inline void G4EventAction::Addfge_x(G4double x) {
  fge_y = x;
}

inline void G4EventAction::Addfge_y(G4double y) {
  fge_y = y;
}

inline void G4EventAction::Addfge_z(G4double z) {
  fge_z = z;
}

inline void G4EventAction::Addfafter_deg(G4double energy) {
  fafter_deg = energy;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
