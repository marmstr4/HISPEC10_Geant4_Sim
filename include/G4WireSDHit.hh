//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4cCalorHit.hh
/// \brief Definition of the B4cCalorHit class

#ifndef G4WireSDHit_h
#define G4WireSDHIt_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

class G4WireSDHit : public G4VHit
{
  public:
    G4WireSDHit();
    G4WireSDHit(const G4WireSDHit&);
    virtual ~G4WireSDHit();

    // operators
    const G4WireSDHit& operator=(const G4WireSDHit&);
    G4bool operator==(const B4cCalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;

  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using G4WireSDHitsCollection = G4THitsCollection<G4WireSDHit>;

extern G4ThreadLocal G4Allocator<G4WireSDHit>* G4WIRESDHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* G4WireSDHit::operator new(size_t)
{
  if (!G4WireSDitAllocator) {
    G4WireSDHitAllocator = new G4Allocator<G4WireSDHit>;
  }
  void *hit;
  hit = (void *) G4WireSDHitAllocator->MallocSingle();
  return hit;
}

inline void B4cCalorHit::operator delete(void *hit)
{
  if (!G4WireSDHitAllocator) {
    G4WireSDHitAllocator = new G4Allocator<G4WireSDHit>;
  }
  G4WireSDHitAllocator->FreeSingle((G4WireSDHit*) hit);
}

inline void G4WireSDHit::Add(G4double de, G4double dl) {
  fEdep += de;
  fTrackLength += dl;
}

inline G4double G4WireSDHit::GetEdep() const {
  return fEdep;
}

inline G4double G4WireSDHit::GetTrackLength() const {
  return fTrackLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
