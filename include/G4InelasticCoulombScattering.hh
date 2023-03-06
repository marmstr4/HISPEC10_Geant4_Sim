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
//      ---------------- G4InelasticCoulombScattering header ----------------
//                 by Christian Stahl, Sep/Oct 2010
// -------------------------------------------------------------------------------
// Short description: a simple process for inelastic coulomb scattering, that is
// a Coulomb excitation reaction without handling the reaction in Geant4.
// Cross section are computed with external programm CLX and read in, the reaction
// kinematics is calculated relativistically.
// -----------------------------------------------------------------------

#ifndef G4InelasticCoulombScattering_hh
#define G4InelasticCoulombScattering_hh

// GEANT4 Headers
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTypes.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4DataInterpolation.hh"
#include "G4NucleiProperties.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4String.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <typeinfo>
#include <stdexcept>
#include "Randomize.hh"



struct isotope;
struct xsecentry;
struct xsecdata;

class BadConversion : public std::runtime_error {
 public:
   BadConversion(std::string const& s)
     : std::runtime_error(s)
     { }
};

template<typename T>
inline std::string stringify(T const& x)
{
  std::ostringstream o;
  if (!(o << x))
   throw BadConversion(std::string("stringify(")
                         + typeid(x).name() + ")");
  return o.str();
}

class G4InelasticCoulombScattering : public G4VDiscreteProcess
{
public:

  // Constructor
  G4InelasticCoulombScattering(const G4String& processName ="InelasticCoulomb");

  // Destructor
  ~G4InelasticCoulombScattering();

  G4bool IsApplicable(const G4ParticleDefinition& particle);

  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                           G4ForceCondition* condition);
  // It returns the MeanFreePath of the process for the current track :
  // (energy, material)
  // The previousStepSize and G4ForceCondition* are not used.
  // This function overloads a virtual function of the base class.
  // It is invoked by the ProcessManager of the Particle.

  void AddMaterial(G4Material* newMat, G4int pZ, G4int pA);
  // Calculates Cross-Sections and inverse mean free pathes for a couple of Projectile with Z=pZ and A=pA with the Material "newMat"
  // This routine is called for every material penetrated by the primary beam, the GetMeanFreePath - function returns interpolated values
  // from the tables created by this function

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  // It computes the final state of the process (at end of step),
  // returned as a ParticleChange object.
  // This function overloads a virtual function of the base class.
  // It is invoked by the ProcessManager of the Particle.

  G4double GetXSecClx(G4double kinE, G4int &pZ, G4int &pA, G4int &tZ, G4int &tA);
  // Calls CLX to calculate the cross section for the input parameters given in 'inputdata'
  // returns the cross section integrated over the angle theta
  // the angle dependant cross section is written into the private array 'xsecTempTable' of type xsecdata with cross section in barn/sr for scattering angles theta in CM-system

   G4double GetXSecDweiko(G4double kinE, G4int &pZ, G4int &pA, G4int &tZ, G4int &tA);

   // functions to calculate thetaMAX (scat. angle  where safe coulex condition is fulfilled)  ////////////
													 //
   G4double Ri(G4int &A);										 //
													 //
   G4double Ci(G4int &A);										 //
													 //
   G4double DFactor(G4int &pA, G4int &tA);								 //
													 //
   G4double ComToLabFaktor(G4int &pA, G4int &tA);							 //
													 //
   G4double aFactor(G4double &kinE, G4int &pA, G4int &tA, G4int &pZ, G4int &tZ);			 //
													 //
   G4double ThetaCOM(G4double &kinE, G4int &pA, G4int &tA, G4int &pZ, G4int &tZ);			 //
													 //
   G4double UsedCoulombBar(G4int &pA, G4int &tA, G4int &pZ, G4int &tZ);					 //
													 //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4LorentzVector GetEnegryMomentumConservation(){return EnMomConservation;}

  G4int GetNumberOfNeutronsInTarget(){return nOfNeutrons;}

  void SetVerboseLevel(G4int ver){VerboseLevel=ver;}

  void SetgsE(G4double en){gsE=en;}

  void SetExE(G4double en){ExE=en;}

  void SetJ0(G4double en){J0=en;}

  void SetJ1(G4double en){J1=en;}

  void Setp0(G4int p){p0=p;}

  void Setp1(G4int p){p1=p;}

  void SetLambda(G4int l){lambda=l;}

  void SetM12(G4double en){M12=en;}

  void SetM22(G4double en){M22=en;}

  void SetBinWidth(G4double en){BinWidth=en;}

  void SetMaxE(G4double en){Emax = en;}

  void SetTargetExcitation(G4double en){TargetExcitation = en;}

  void SetEnhance(G4double en){Enhance = en;}

  void SetmaxThetaCOM(G4double en){maxThetaCOM = en;}

  void SetangleBinning(G4double en){angleBinning = en;}

  void SetmaxEnergyClx(G4double en){maxEnergyClx = en;}

  void SetminEnergyDweiko(G4double en){minEnergyDweiko = en;}

  void SetSafeCoulexClx(G4bool en){SafeCoulexClx = en;}

  void SetSafeCoulexDweiko(G4bool en){SafeCoulexDweiko = en;}



private:


  // Hide assignment operator as private
  G4InelasticCoulombScattering& operator=(const G4InelasticCoulombScattering &right);

  // Copy constructor
  G4InelasticCoulombScattering(const G4InelasticCoulombScattering&);

  //instance of the Nuclear properties class
  G4NucleiProperties NucProp;

  // Working parameters
  G4LorentzVector EnMomConservation;                  // Residual of Energy/Momentum Cons.
  G4int nOfNeutrons;                                  // #of neutrons in the target nucleus

  // Variables for the reaction
  static std::vector<std::vector<isotope> > Isotopes;        // Vector with registered Isotopes
  static std::vector<G4Material*> Materials;  // Vector with registered Materials

  G4int VerboseLevel;
  G4bool WriteOutKinematics;
  G4double pAngleMaxCOM;
  G4double tAngleMaxCOM;
  G4double TwoPi;
  G4double deg2rad;
  G4double Emax;
  G4double Enhance;
  G4double maxThetaCOM;
  G4double angleBinning;
  G4double maxEnergyClx;
  G4double minEnergyDweiko;
  G4double coulombbarrier;
  G4bool SafeCoulexClx;
  G4bool SafeCoulexDweiko;
  G4Material* prevMat;
  G4int id;   //id of the Material we're in
  G4double gsE;
  G4double ExE;
  G4double J0;
  G4double J1;
  G4double K0;
  G4double K1;
  G4int p0;
  G4int p1;
  G4int lambda;
  G4double M12;
  G4double M22;
  G4double BinWidth;
  G4double TargetExcitation;
  G4double TargetExcitationCLX;	// variable for target or projectile excitation
  G4int nEbins;
  G4double* E_vals;
  G4double lowerBinDweiko;
  G4double upperBinClx;
   std::vector<xsecdata> xsecTempTable;
  std::vector<xsecdata> xsecTempTableClx;
  std::vector<xsecdata> xsecTempTableDweiko;
  //xsecdata* xsecTempTable;

};
#endif
