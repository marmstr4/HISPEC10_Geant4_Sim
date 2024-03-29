//
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
//      ---------------- G4InelasticCoulombScattering -----------------
//                 by Christian Stahl, Sep/Oct 2010
// --------------------------------------------------------------------


#include "G4InelasticCoulombScattering.hh"
#include "Randomize.hh"
#include "G4NistManager.hh"

using std::istringstream;


struct isotope
    {
        G4int Z;		//Target Z
        G4int A;		//Target A
        G4double mass;		//Target mass
        G4double density;	//Isotope density
	xsecentry* xsec_tableClx;	//Table with (interpolated) cross sections
	xsecentry* xsec_tableDweiko;	//Table with DWEIKO cross sections
	xsecentry* xsec_table;		//Table with CLX cross sections
	G4double TransConst;	//constant useful for transforming from CM to Lab system
    };

struct xsecdata{
    G4double angle;	//CM Scattering angle of the Projectile
    G4double xsec;	//make this an array if more than 1 excited states are considered
};

struct xsecentry{
    std::vector<xsecdata> xsec;
    G4double xsec_sum;
    G4double partial_inv_mfp;
};

// Initialization of static vectors
std::vector<std::vector<isotope> > G4InelasticCoulombScattering::Isotopes;     //Vector containing Z,A,density, cross section tables and summed cross section for every isotope in the known materials
std::vector<G4Material*> G4InelasticCoulombScattering::Materials;  // vector of all known materials

// Constructor
G4InelasticCoulombScattering::G4InelasticCoulombScattering(const G4String& processName):
  G4VDiscreteProcess(processName, fHadronic)
{
  //initilializing some standart values, can be changed via macro
  TwoPi=8*atan2(1,1);
  deg2rad=TwoPi/360;
  VerboseLevel = 0;
  prevMat=0;
  K0=0;
  K1=0;
  id=-1;
  //Emax=1000*MeV;
  Enhance=1;
  //Emax=6*64*MeV;
  Emax=700*MeV;
  BinWidth=20*MeV;
  coulombbarrier=290*MeV;
  //minEnergyDweiko=1000*coulombbarrier*MeV-BinWidth*5;
  minEnergyDweiko=coulombbarrier*MeV-50;
  //maxEnergyClx=790*MeV;
  maxEnergyClx=390*MeV;
  ExE=1345*keV;
  //ExE=605.9*keV;
  angleBinning=1;
  SafeCoulexDweiko=true;
  SafeCoulexClx=true;



  WriteOutKinematics = false;
  if(WriteOutKinematics){
	G4cout <<"Syntax in kinematic file 'kin.dat' (P-> Projectile, R->Recoil):"<<G4endl;
	G4cout <<"Ptheta(lab)[°] <> Rtheta(lab)[°] <> Ptheta(COM)[°] <> Eout-Ein[MeV] <> PkinE[MeV] <> RkinE[MeV]"<<G4endl;
  }

  E_vals=NULL;

}

// Destructor
G4InelasticCoulombScattering::~G4InelasticCoulombScattering() {

}


G4bool G4InelasticCoulombScattering::IsApplicable(const G4ParticleDefinition& particle)
{

  	G4int Z=static_cast<G4int>(particle.GetPDGCharge());
  	G4int A=particle.GetBaryonNumber();
  	if      (particle == *( G4Deuteron::Deuteron()     )) return true;
  	else if (particle == *( G4Alpha::Alpha()           )) return true;
  	else if (particle == *( G4Triton::Triton()         )) return true;
  	else if (particle == *( G4He3::He3()               )) return true;
  	else if (particle == *( G4GenericIon::GenericIon() )) return true;
  	else if (Z > 0 && A > 0)                              return true;
    /*
    G4int Z=static_cast<G4int>(particle.GetPDGCharge());
  	G4int A=particle.GetBaryonNumber();
    G4cout<<Z<<" "<<A<<G4endl;
    if (Z==28 && A==64) {G4cout<<"true"<<G4endl;return true;}
    */
    return false;
}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4InelasticCoulombScattering::GetMeanFreePath(const G4Track& aTrack, G4double,
                                           G4ForceCondition* Fc)
{

  //array with energies in the cross section tables; calculation of the bins for the energie gates of clx and dweiko
  if(!E_vals){
      nEbins=int(((Emax+BinWidth-ExE)/MeV)/BinWidth)+1;
      E_vals = new G4double[nEbins];
      lowerBinDweiko=-1;
      upperBinClx=-1;
      for(int i=0; i<nEbins; i++){
        E_vals[i]=ExE+i*BinWidth;
	if(lowerBinDweiko<0 && E_vals[i]>minEnergyDweiko)	lowerBinDweiko=i;
	if(upperBinClx<0 && E_vals[i]>maxEnergyClx)	upperBinClx=i;
      }
      if(minEnergyDweiko > Emax) lowerBinDweiko = nEbins;
      if(maxEnergyClx > Emax) upperBinClx = nEbins;
      G4cout << "Using DWEIKO for cross sections calculation at energies above " <<(ExE+BinWidth*lowerBinDweiko)/MeV << " MeV"<<G4endl;
      if(SafeCoulexDweiko) G4cout << "      safe coulex condition active for DWEIKO cross sections" << G4endl;
      else G4cout << "      safe coulex condition -NOT- active for DWEIKO cross sections" << G4endl;
      G4cout << "      DWEIKO cross sections cut at "<<maxThetaCOM << "deg" << G4endl;
      G4cout << "Using CLX for cross sections calculation at energies below " <<(ExE+BinWidth*upperBinClx)/MeV << " MeV"<<G4endl;
      if(SafeCoulexClx) G4cout << "      safe coulex condition active for CLX cross sections" << G4endl;
      else G4cout << "      safe coulex condition -NOT- active for CLX cross sections" << G4endl;

  }

/// Allow this process for Primaries only
 if(aTrack.GetParentID()!=0) return DBL_MAX;


  // static const G4double mProt = G4Proton::Proton()->GetPDGMass()/MeV;// proton Mass in MeV
 G4int idE;		//index of first Energy bin with Energy higher than the projectile kinetic energy
 G4double inv_MFP = 0; 	//inverse mean free path -> will contain sums of partial inverse mfp's for the present isotopes
 G4double mfp=0;	// mean free path
 G4double inv_MFP_tmp=0;

  *Fc = NotForced;   //Force Condition NotForced means: PostStepDoIt of this process is not forced unless this process limits the step, i.e. this process was chosen to happen next

/// ///////////////////////////////// ///
/// collect information on projectile
/// ///////////////////////////////// ///

  const G4DynamicParticle* projectile = aTrack.GetDynamicParticle();
  G4ParticleDefinition* projectileDefinition=projectile->GetDefinition();
  //Can process be applied to this particle?
  if( !IsApplicable(*projectileDefinition)) G4cout<<"-Warning-G4InelasticCoulombScattering::GetMeanFreePath: notImplementedParticle"<<G4endl;
	//Projectile mass and charge number, kinetic energy
  	G4int pZ=static_cast<G4int>(projectileDefinition->GetPDGCharge());
  	G4int pA=projectileDefinition->GetBaryonNumber();
    //G4int pZ = 28;
    //G4int pA = 64;
	G4double pKinE= projectile->GetKineticEnergy();

/// No coulex if projectile kin. Energy is lower than excitation energy
  if(pKinE<=ExE) return DBL_MAX;

  if(VerboseLevel>1) G4cout << "COULEX::GMFP: Projectile: Z="<<pZ<<" A="<<pA<<" ("<<projectileDefinition->GetParticleName()<<"), E="<<pKinE/MeV<<"MeV"<<G4endl;


/// ////////////////////////////////////////////
/// Block collecting information on target and calculate Coulex cross sections if necessary
/// ////////////////////////////////////////////
// check wether material is in the list in vector 'Materials'. If not, collect data on it and save list of Isotopes present in this material with their Z,A and nucleus density in the vector Isotopes[i][j] where i is the id of the material (to be found from vector Materials[i], numbering is same for both, and j refers to the different isotopes. Get number of Isotopes in Material i with Isotopes[i].size()

  G4int nE=0;

  //if (aTrack.GetMaterial()->GetName()!="Ge_Au") return DBL_MAX;


  if(prevMat==0 || prevMat != aTrack.GetMaterial()){
    prevMat = aTrack.GetMaterial();
    nE=prevMat->GetNumberOfElements();
    id=-1;
    for(int i=0; i<G4int(Materials.size()); i++){
	if(Materials[i]==prevMat){
		id=size_t(i);
		if(VerboseLevel>1) G4cout<<"COULEX::GMFP: Material " << prevMat->GetName() <<" is known (id " <<id<<") and consists of " <<nE<<" Element(s)"<<G4endl;
	}
    }
    if(id==-1){  //current Material not yet in the list

        G4cout << "COULEX::GMFP: Calculating cross sections for projectile " << projectileDefinition->GetParticleName() << G4endl;
	AddMaterial(prevMat, pZ, pA);
	if(VerboseLevel>1) G4cout<<"COULEX::GMFP: Material " << prevMat->GetName() <<" is unknown and consists of " <<nE<<" Element(s)"<<G4endl;
    }

  }else{
	 nE=prevMat->GetNumberOfElements();
	 if(VerboseLevel>1) G4cout<<"COULEX::GMFP: Still in Material " << prevMat->GetName() <<" (id " <<id<<")"<<G4endl;
  }

  if(!nE) nE=prevMat->GetNumberOfElements();

/// //////////////////////////////////////////////////////////// ///
/// Calculate mean free path for current energy by interpolation ///
/// //////////////////////////////////////////////////////////// ///

	idE =0;
	if(pKinE>Emax){
	  idE=nEbins-1;
	  //G4cout << "!!! exception: Encountered kinetic energy exceeding cross section table. Table end at " << Emax/MeV << "Mev, particle has " << pKinE/MeV << "MeV" << G4endl;
	}
	else while(pKinE>E_vals[idE]) idE++;	//store index of first Energy bin with Energy higher than pKinE

	//linear interpolation of inverse mean free path for each isotope in this material
	//with f(E)->inverse MFP at Energy E, E0 energy of the sample with first energy smaller than pKinE (index idE-1), E1 energy of the sample with first energy bigger than pKinE (index idE):
	//f(E)=f(E0)+(f(E1)-f(E0))/(E1-E0)*(E-E0)

	if(VerboseLevel>1)G4cout <<"COULEX::GMFP: Interpolating for pKinE="<< std::fixed << std::setprecision(3)<<pKinE/MeV<<"MeV: idE="<<idE<<", E[idE-1]="<<  E_vals[idE-1]/MeV<<"MeV, E[idE]="<< E_vals[idE]/MeV<<G4endl;
	for(int i=0; i<G4int((Isotopes[id]).size()); i++){
//	  G4cout<<"i:	"<<i<<"	id:	"<<id<<G4endl;
		inv_MFP_tmp=Isotopes[id][i].xsec_table[idE-1].partial_inv_mfp + (Isotopes[id][i].xsec_table[idE].partial_inv_mfp-Isotopes[id][i].xsec_table[idE-1].partial_inv_mfp)*(pKinE-E_vals[idE-1])/BinWidth;
		inv_MFP+=inv_MFP_tmp;

		if(VerboseLevel>1) G4cout << "		For Isotope Z="<<Isotopes[id][i].Z<<", A=" << Isotopes[id][i].A <<": interpolated XSec="<< std::scientific << std::setprecision(3)<< (Isotopes[id][i].xsec_table[idE-1].xsec_sum + (Isotopes[id][i].xsec_table[idE].xsec_sum-Isotopes[id][i].xsec_table[idE-1].xsec_sum)*(pKinE-E_vals[idE-1])/BinWidth)/barn<<"barn, interpolated mfp="<< std::scientific << std::setprecision(3)<<(1/inv_MFP_tmp)/micrometer<<"µm"<<G4endl;
	}
	inv_MFP*=Enhance;		//Apply Enhancement factor set via the macro
	if(VerboseLevel>1){
		G4cout << "COULEX::GMFP: MFP="<< std::scientific << std::setprecision(3) << (1/inv_MFP)/micrometer << "µm";
        	if(Enhance!=1) G4cout <<", cross section was enhanced by a factor of "<< std::fixed<<Enhance << G4endl;
	}

   if(inv_MFP>0.){
     mfp=1/inv_MFP;
   }else{
     mfp=DBL_MAX;
   }

   if(mfp>DBL_MAX) return DBL_MAX;
   else return mfp;

     //return 0.001*nanometer;  //dummy output: forcing reaction
}

/// /////////////////////////////////////////////// ///
/// AddMaterial: Comp. Xsecs for a certain material ///
/// /////////////////////////////////////////////// ///
void G4InelasticCoulombScattering::AddMaterial(G4Material* newMat, G4int pZ, G4int pA){

        Materials.push_back(newMat);
        id=Materials.size()-1;

	bool IsVacuum = false;

	if(newMat->GetName() == "Vacuum"){
	  IsVacuum = true;
	}

        G4int nE=newMat->GetNumberOfElements();

        if(VerboseLevel>0) G4cout<<"COULEX::GMFP: Material " << newMat->GetName() <<" is NOT known (new id: " <<id<<") and consists of " <<nE<<" Element(s)"<<G4endl;

        //Adding Isotope List of current material to vector Isotopes
        //Find Materials in current Target-layer
        const G4double* NOfNucPerVolume = newMat->GetVecNbOfAtomsPerVolume();
        const G4ElementVector* ElementVector = newMat->GetElementVector();
        G4IsotopeVector** IsoVec = new G4IsotopeVector*[nE]; //Array to be filled with pointers to element i's Isotope Vector
        G4double** AbundanceVec = new G4double*[nE];

//some detailed output on the material
        if(VerboseLevel>1){
          G4IsotopeVector** InfoIsoVec = new G4IsotopeVector*[nE];
          G4double** InfoAbundanceVec = new G4double*[nE];
          const G4double* InfoNOfNucPerVolume = newMat->GetVecNbOfAtomsPerVolume();
          const G4ElementVector* InfoElementVector = newMat->GetElementVector();
          for(int i=0; i<nE; i++){
                InfoIsoVec[i]=(*InfoElementVector)[i]->GetIsotopeVector();
                InfoAbundanceVec[i] = (*InfoElementVector)[i]->GetRelativeAbundanceVector();
                G4cout << "COULEX::GMFP: Element "<<i+1<<" is " << (*InfoElementVector)[i]->GetName() << " having "<< std::scientific << std::setprecision(3) << InfoNOfNucPerVolume[i]*cm3 <<" Atoms/cm³ and " << G4int((*InfoIsoVec[i]).size())<< " stable isotopes" <<G4endl;
                if(VerboseLevel>2){
                        for(int j=0; j<G4int((*InfoIsoVec[i]).size()); j++){
                                G4cout << "             " << j+1 <<": " << (*InfoIsoVec[i])[j]->GetName() << " with Z=" <<(*InfoIsoVec[i])[j]->GetZ()<<", A="<< (*InfoIsoVec[i])[j]->GetN()<<", relative Abundance "<< std::scientific << std::setprecision(3) <<InfoAbundanceVec[i][j]<< G4endl;
                        }
                }
          }
          delete []InfoIsoVec;
          delete []InfoAbundanceVec;
        }



///Create Coulex Cross section tables
        std::vector<isotope> ThisIsoVector;     //(temporary) vector of type isotope (struct defined at top) for isotopes in this Material
        isotope ThisIso;
        G4int nIsos=0;

        //count number of isotopes in this material
        for(int i=0; i<nE; i++){
                nIsos+=(*ElementVector)[i]->GetNumberOfIsotopes();
        }

        //array with angle bins for creation of the cross section histogram
        angleBinning =1;
	G4int AngleBin=(180+angleBinning)/angleBinning;

	G4double* T_vals;
	T_vals = new G4double[AngleBin];
	for(unsigned int i=0; i*angleBinning<=180; i++){
	  T_vals[i]=i*angleBinning*degree;
	}

        for(int i=0; i<nE; i++){
                IsoVec[i]=(*ElementVector)[i]->GetIsotopeVector();
                AbundanceVec[i] = (*ElementVector)[i]->GetRelativeAbundanceVector();
                for(int j=0; j<G4int((*IsoVec[i]).size()); j++){
                        ThisIso.Z=(*IsoVec[i])[j]->GetZ();
                        ThisIso.A=(*IsoVec[i])[j]->GetN();
			ThisIso.mass=NucProp.GetNuclearMass(ThisIso.A, ThisIso.Z)/MeV;
                        ThisIso.density=NOfNucPerVolume[i]*AbundanceVec[i][j];

                        if(IsVacuum) G4cout << "COULEX::GMFP: Vacuum, no cross section calculation." << G4endl;
			 else G4cout << "COULEX::GMFP: Calculating XSec Table on target Z=" << ThisIso.Z << " A=" << ThisIso.A<<" from "<< std::fixed << std::setprecision(2)<< ExE/MeV<<" to "<< (Emax+ExE)/MeV<<"MeV in steps of "<< BinWidth/MeV<<"MeV: " << G4endl;
                        if(VerboseLevel>3) G4cout  << G4endl;

                        xsecentry* ThisXsecTable = new xsecentry[nEbins];
			std::vector<std::vector<G4double> > diffXsec;
			std::vector<G4double> diffXsec_E;

                        xsecentry* ThisXsecTableClx = new xsecentry[nEbins];
			std::vector<std::vector<G4double> > diffXsecClx;
			std::vector<G4double> diffXsec_EClx;

			xsecentry* ThisXsecTableDweiko = new xsecentry[nEbins];
			std::vector<std::vector<G4double> > diffXsecDweiko;
			std::vector<G4double> diffXsec_EDweiko;


                        for(G4int k=0; k<nEbins; k++){

			  if(IsVacuum){
			    //diff. xsec table
			    xsecTempTableClx.clear();
			    xsecdata tempdatClx;
			    tempdatClx.xsec=0*barn;
			    for( int i=0; i*angleBinning<=180; i++){
				    tempdatClx.angle=i*angleBinning*degree;
				    xsecTempTableClx.push_back(tempdatClx);
			    }
			    ThisXsecTableClx[k].xsec=xsecTempTableClx;
			    ThisXsecTableClx[k].partial_inv_mfp = 0;

			    //CLX cross sections
			    ThisXsecTableClx[k].xsec_sum 	= 0.;  // write xsec table to private variable xsecTempTable, function returns integrated cross section
			    ThisXsecTableClx[k].xsec		= xsecTempTableClx;
			    ThisXsecTableClx[k].partial_inv_mfp = 0.;

			    //DWEIKO cross sections
			    ThisXsecTableDweiko[k].xsec_sum 	= 0.;  // write xsec table to private variable xsecTempTable, function returns integrated cross section
			    ThisXsecTableDweiko[k].xsec		= xsecTempTableClx;
			    ThisXsecTableDweiko[k].partial_inv_mfp = 0.;

			    //combined cross sections
			    ThisXsecTable[k].xsec_sum 	= 0.;  // write xsec table to private variable xsecTempTable, function returns integrated cross section
			    ThisXsecTable[k].xsec		= xsecTempTableClx;
			    ThisXsecTable[k].partial_inv_mfp = 0.;

			    //create vector with differential cross section for creating the histogram
			    for(unsigned int theta=0; theta<xsecTempTableClx.size(); theta++){
			      diffXsec_E.push_back(xsecTempTableClx[theta].xsec);
			    }
			    diffXsec.push_back(diffXsec_E);
			    diffXsecClx.push_back(diffXsec_E);
			    diffXsecDweiko.push_back(diffXsec_E);
			    diffXsec_E.clear();

			    continue;
			  }

			  //There are 3 tables: "clx-only"; "dweiko-only" and one with dweiko and clx ("both")
			  //G4cout << "  Coulex XSec for E=" << E_vals[k]/MeV << " MeV:";

			  if(k<=upperBinClx){
			    //fill table "clx-only"  with clx calculation only

			    ThisXsecTableClx[k].xsec_sum 	= 100*GetXSecClx(E_vals[k], pZ, pA, ThisIso.Z, ThisIso.A);  // write xsec table to private variable xsecTempTable, function returns integrated cross section
			    ThisXsecTableClx[k].xsec		= xsecTempTableClx;
          //Limit CLX to below

			    ThisXsecTableClx[k].partial_inv_mfp = ThisXsecTableClx[k].xsec_sum*ThisIso.density;
			    //G4cout << "   CLX: " << ThisXsecTableClx[k].xsec_sum/barn*1000 << " mbarn ";
          G4cout << "   CLX: " <<E_vals[k]/MeV/64<<","<<ThisXsecTableClx[k].xsec_sum/barn*1000 << " mbarn "<<G4endl;;

			    if(k<=lowerBinDweiko){
			       //fill table "both" with clx calculation
			      ThisXsecTable[k].xsec_sum 	= ThisXsecTableClx[k].xsec_sum;
			      xsecTempTable=xsecTempTableClx;
			      ThisXsecTable[k].xsec		= xsecTempTable;
			      ThisXsecTable[k].partial_inv_mfp  = ThisXsecTableClx[k].partial_inv_mfp;
			    }
			  }else{
			    //fill fill table "clx-only"  with 0
				ThisXsecTableClx[k].xsec_sum = 0;
				xsecTempTableClx.clear();
				xsecdata tempdatClx;
				tempdatClx.xsec=0*barn;

				for( int i=0; i*angleBinning<=180; i++){
				    tempdatClx.angle=i*angleBinning*degree;
				    xsecTempTableClx.push_back(tempdatClx);
				}
				ThisXsecTableClx[k].xsec=xsecTempTableClx;
				ThisXsecTableClx[k].partial_inv_mfp = 0;

			  }
			  if(k>=lowerBinDweiko){
			    //fill table "dweiko-only" with dweiko calculation only
			    ThisXsecTableDweiko[k].xsec_sum = 100000*GetXSecDweiko(E_vals[k], pZ, pA, ThisIso.Z, ThisIso.A);  // write xsec table to private variable xsecTempTable, function returns integrated cross section
			    ThisXsecTableDweiko[k].xsec		   = xsecTempTableDweiko;
			    ThisXsecTableDweiko[k].partial_inv_mfp = ThisXsecTableDweiko[k].xsec_sum*ThisIso.density;
			    G4cout << "   DWEIKO: " <<E_vals[k]/MeV/64<<","<< ThisXsecTableDweiko[k].xsec_sum/barn*1000<< " mbarn "<<G4endl;
         //G4cout << "   DWEIKO: " <<E_vals[k]/MeV<<","<<ThisXsecTableClx[k].xsec_sum/barn*1000 << " mbarn "<<G4endl;;
        //G4cout << "   DWEIKO: " <<E_vals[k]/MeV/74<<","<<ThisXsecTableDweiko[k].xsec_sum/barn*1000 << " mbarn "<<G4endl;;
			    if(k>=upperBinClx){



			      //fill table "both" with dweiko calculation only
			      ThisXsecTable[k].xsec_sum        = ThisXsecTableDweiko[k].xsec_sum;
			      xsecTempTable=xsecTempTableDweiko;
			      ThisXsecTable[k].xsec	       = xsecTempTable;
			      ThisXsecTable[k].partial_inv_mfp = ThisXsecTableDweiko[k].partial_inv_mfp;
			    }
			  }else{
			     //fill fill table "dweiko-only"  with 0
				ThisXsecTableDweiko[k].xsec_sum = 0;
				xsecTempTableDweiko.clear();
				xsecdata tempdatDweiko;
				tempdatDweiko.xsec=0*barn;

				for( int i=0; i*angleBinning<=180; i++){
				    tempdatDweiko.angle=i*angleBinning*degree;
				    xsecTempTableDweiko.push_back(tempdatDweiko);
				}
				ThisXsecTableDweiko[k].xsec=xsecTempTableDweiko;
				ThisXsecTableDweiko[k].partial_inv_mfp = 0;
			  }

			  if((k>lowerBinDweiko)&&(k<upperBinClx)){
			    // now the table "both" gets filled with a linear extrapolation of dweiko an clx


			    G4double factorClx=((upperBinClx-lowerBinDweiko)-(k-lowerBinDweiko))/(upperBinClx-lowerBinDweiko);
			    G4double factorDweiko=1-factorClx;

			    ThisXsecTable[k].xsec_sum = ThisXsecTableClx[k].xsec_sum*factorClx+ThisXsecTableDweiko[k].xsec_sum*factorDweiko;
			    xsecTempTable.clear();
			    xsecdata tempdat;

			    for( int i=0; i*angleBinning<=180; i++){
				tempdat.angle=xsecTempTableClx[i].angle;
				tempdat.xsec=xsecTempTableClx[i].xsec*factorClx+xsecTempTableDweiko[i].xsec*factorDweiko;
				xsecTempTable.push_back(tempdat);
			    }
			    ThisXsecTable[k].xsec=xsecTempTable;
			    ThisXsecTable[k].partial_inv_mfp = ThisXsecTableClx[k].partial_inv_mfp*factorClx+ThisXsecTableDweiko[k].partial_inv_mfp*factorDweiko;

			     }

			  //G4cout << "   COMBINED: "<<E_vals[k]/MeV/64<<" "<<ThisXsecTable[k].xsec_sum/barn*1000 << " mbarn " << G4endl;

			   //create vector with differential cross section for creating the histogram
			  for(unsigned int theta=0; theta<xsecTempTable.size(); theta++){
			    diffXsec_E.push_back(xsecTempTable[theta].xsec);
  			  }
			  diffXsec.push_back(diffXsec_E);
			  diffXsec_E.clear();

			  for(unsigned int theta=0; theta<xsecTempTableClx.size(); theta++){
			    diffXsec_EClx.push_back(xsecTempTableClx[theta].xsec);
  			  }
			  diffXsecClx.push_back(diffXsec_EClx);
			  diffXsec_EClx.clear();

			  for(unsigned int theta=0; theta<xsecTempTableDweiko.size(); theta++){
			    diffXsec_EDweiko.push_back(xsecTempTableDweiko[theta].xsec);
  			  }
			  diffXsecDweiko.push_back(diffXsec_EDweiko);
			  diffXsec_EDweiko.clear();

                          if(VerboseLevel>3)	    G4cout << "        For E="<< std::scientific << std::setprecision(3)<<E_vals[k]/MeV<<"MeV: Summed XSec="<< std::scientific << std::setprecision(3) << ThisXsecTableClx[k].xsec_sum/barn << "barn -> MFP="<< std::scientific << std::setprecision(3)<< (1/ThisXsecTableClx[k].partial_inv_mfp)/micrometer <<"µm"<< G4endl;

                          flush(G4cout);
                        }

			ThisIso.xsec_table=ThisXsecTable;
                        ThisIso.xsec_tableClx=ThisXsecTableClx;
                        ThisIso.xsec_tableDweiko=ThisXsecTableDweiko;

                        ThisIsoVector.push_back(ThisIso);


                        G4cout << "	done." << G4endl;

		}
        }
        Isotopes.push_back(ThisIsoVector);  //write Isotope list for current material into vector
        ThisIsoVector.clear();
}

G4VParticleChange* G4InelasticCoulombScattering::PostStepDoIt(const G4Track& track, const G4Step& step){
//variables needed here
    G4double mfp_sum=0;
    G4int idE=0;		//index of first Energy bin with Energy higher than pKinE
    G4int IsoID=0;		//id-# of picked isotope in target for the reaction
    std::vector<G4double> inv_mfp_vec;	//in entry #i: sum of interpolated inverse MFP's up to isotopes #i in this material
    G4double max_xsec=0;	//maximum differential cross section of an angle bin at this projectile energy
    G4double rnd2;
    G4int rnd1;
    G4double pAngleTheta;	//projectile polar angle in CM System
    G4double pAnglePsi;		//projectile polar angle in Lab System
    G4double tAnglePsi;		//target (recoil) polar angle in Lab System
    G4double pAnglePhi;		//projectile azimuthal in Lab System (arbitary if no polarisation)
    G4double tAnglePhi;  	//recoil azimuthal in Lab System (arbitary if no polarisation) = pAnglePhi+180°
    G4double pKinE;		//projectile kinetic Energy before scattering off target nucleus in Lab System
    G4double pKinEScat;		//projectile kinetic Energy after scattering off target nucleus in Lab System
    G4double tKinEScat;		//target (recoil) kinetic energy in Lab system
    G4double velocity = track.GetVelocity();
/// get it from the messenger!
    G4bool WriteOutKinematics=false;


/// Information on projectile
    const G4DynamicParticle* projectile = track.GetDynamicParticle();	//pointer to Projectile

    G4ParticleDefinition* projectileDefinition=projectile->GetDefinition();

    pKinE= projectile->GetKineticEnergy(); 				//projectile kinetic Energy

    G4ParticleMomentum pDir = projectile->GetMomentumDirection();	//projectile direction: It is a unit three-vector

    //G4double Momentum = projectile->GetTotalMomentum(); 		//Momentum
    G4int pZ=static_cast<G4int>(projectileDefinition->GetPDGCharge());  //projectile Z

    G4int pA=projectileDefinition->GetBaryonNumber();			//projectile A



    //if (pZ == 28 && pA == 64) {

/// Information on Target
    //Target Material cannot have changed since MFP was calculated last time (it's recalculated at material boundaries)
    // -> we are in the material prevMat with index 'id' in the Materials and Isotope vector
    // Sample on what Isotope in this Material to scatter off weighted by the MFP's for the distinct isotopes


    if(pKinE>Emax) idE=nEbins-1;
    else while(pKinE>E_vals[idE]) idE++;	//store index of first Energy bin with Energy higher than pKinE

    //linear interpolation of inverse mean free path for each isotope in this material
    //with f(E)->inverse MFP at Energy E, E0 energy of the sample with first energy smaller than pKinE (index idE-1), E1 energy of the sample with first energy bigger than pKinE (index idE):
    //f(E)=f(E0)+(f(E1)-f(E0))/(E1-E0)*(E-E0)

     for(int i=0; i<G4int((Isotopes[id]).size()); i++){
	mfp_sum+=Isotopes[id][i].xsec_table[idE-1].partial_inv_mfp + (Isotopes[id][i].xsec_table[idE].partial_inv_mfp-Isotopes[id][i].xsec_table[idE-1].partial_inv_mfp)*(pKinE-E_vals[idE-1])/BinWidth;
	inv_mfp_vec.push_back(mfp_sum);
	if(VerboseLevel>3) G4cout <<"	      invMFPsum up to of isotope #"<<i<<" ("<<Isotopes[id][i].Z<<"|"<<Isotopes[id][i].A<<"): "<< std::scientific << std::setprecision(3)<< (inv_mfp_vec[i])<<G4endl;
     }

/// pick isotope
    G4double IsoRand = mfp_sum*G4UniformRand();		//random double between 0 and mfp_sum
    if((Isotopes[id]).size()) while(IsoRand>=inv_mfp_vec[IsoID]) IsoID++;
    inv_mfp_vec.clear();

    if(VerboseLevel>3) G4cout<<"	      RND was " << IsoRand << " ->";
    if(VerboseLevel>2) G4cout<<"	      Picked target nucleus with Z="<<Isotopes[id][IsoID].Z<<" and A="<<Isotopes[id][IsoID].A << G4endl;

// array of interpolated cross section per angle
    xsecdata xsec_interpol[181];
    for(int i=0; i<181; i++){
	xsec_interpol[i].angle=Isotopes[id][IsoID].xsec_table[idE].xsec[i].angle;

	//interpolate cross section:
	xsec_interpol[i].xsec=Isotopes[id][IsoID].xsec_table[idE-1].xsec[i].xsec + (Isotopes[id][IsoID].xsec_table[idE].xsec[i].xsec-Isotopes[id][IsoID].xsec_table[idE-1].xsec[i].xsec)*(pKinE-E_vals[idE-1])/BinWidth;

	//find maximum differential cross section among all angles theta
	if(xsec_interpol[i].xsec>max_xsec)max_xsec=xsec_interpol[i].xsec;

	if(VerboseLevel>3)G4cout << "angle: " << std::scientific << std::setprecision(3)<< xsec_interpol[i].angle/degree << "°, interpolated: " << std::scientific << std::setprecision(3)<< Isotopes[id][IsoID].xsec_table[idE-1].xsec[i].xsec/millibarn<<"mb + ("<< std::scientific << std::setprecision(3)<<Isotopes[id][IsoID].xsec_table[idE].xsec[i].xsec/millibarn<<"mb - "<<Isotopes[id][IsoID].xsec_table[idE-1].xsec[i].xsec/millibarn<<"mb) * "<< std::scientific << std::setprecision(3)<<(pKinE-E_vals[idE-1])/BinWidth<<"="<< std::scientific << std::setprecision(3)<<xsec_interpol[i].xsec/millibarn<<"mb"<<G4endl;
    }

/// pick projectile scattering angle (CMS)
// notice: target angle in CMS is necessarily projectile angle - Pi
// lab angles are with respect to the direction of motion of the projectile before scattering

    //pick angle bin according to differential cross section distribution
    do{
	rnd2=max_xsec*G4UniformRand();	//random double between 0 and the maximum cross section per angle for this energy
	rnd1=int(181*G4UniformRand());	//random angle between 0 and 180
    }while(rnd2>xsec_interpol[rnd1].xsec);

    //find particle polar angle in CM system
    pAngleTheta=(xsec_interpol[rnd1].angle/degree-0.5+G4UniformRand())*degree;	//random angle within the bin just chosen
    if(pAngleTheta/degree<0) pAngleTheta=0;
    if(pAngleTheta/degree>180) pAngleTheta=180*degree;

    //if(VerboseLevel>2)
    // G4cout << "	      CM projectile polar scattering angle is "<< std::scientific << std::setprecision(3)<<pAngleTheta/degree<<"°"<<G4endl;

    pAnglePhi = 360*G4UniformRand()*degree; 	//projectile azimuthal angle phi is arbitary (no polarisation involved)
    if(pAnglePhi/degree > 180)tAnglePhi = pAnglePhi-180*degree;		//recoil will get azimuthal angle pAnglePhi+180° to conserve momentum
    else tAnglePhi = pAnglePhi+180*degree;

/// Calculate kinematics in Lab system
    //formulas taken from A.M.Baldin, W.L. Goldanksij, L.L. Reosenthal: "Kinematik der Kernreaktionen", Akademie-Verlag Berling, 1963

    G4double V2, Ei, pi2, pi, Eg, Eg2, EgCOM2, A1, A2, rho1, rho2, mi, mi2, mii, mii2, m1, m1_2, t1cos, t2cos;
    G4double pEg, tEg;

    mi = NucProp.GetNuclearMass(pA, pZ)/MeV;	//projectile mass before coulex
    mi2 = mi*mi;
    mii = Isotopes[id][IsoID].mass/MeV;	//target nucleus mass
    mii2=mii*mii;
    m1=mi+ExE/MeV;			//projectile mass after coulex (gains 'weight' by excitation)
    m1_2=m1*m1;
    Ei = pKinE/MeV+mi;			//projectile initial total energy
    pi2 = Ei*Ei - mi2;			//projectile initial momentum
    pi = std::sqrt(pi2);
    Eg=Ei+mii;				//total energy in lab system
    Eg2 = Eg*Eg;
    V2=pi2/Eg2;				//a factor
    EgCOM2 = Eg2*(1-V2);		//total energy in COM system
    A1=Eg2-pi2+m1_2-mii2;		//a factor
    A2=Eg2-pi2+mii2-m1_2;		//a factor

    rho1=std::sqrt(Ei*Ei-mi2)*A1/(Eg*std::sqrt(A1*A1-4*EgCOM2*m1_2));
    rho2=std::sqrt(Ei*Ei-mi2)*A2/(Eg*std::sqrt(A2*A2-4*EgCOM2*mii2));


    pAnglePsi  = std::atan(std::sin(pAngleTheta/rad)*std::sqrt(1-V2) / (std::cos(pAngleTheta/rad)+rho1))*rad;
    if(pAnglePsi<0)pAnglePsi+=180*degree;	//force positive angles: atan gives angles between -90° and 90°, I want between 0° and 180°
    t1cos = std::cos(pAnglePsi/rad);

    tAnglePsi  = std::atan(std::sin(pAngleTheta/rad)*std::sqrt(1-V2) / ((-1.0)*std::cos(pAngleTheta/rad)+rho2))*rad;
    t2cos = std::cos(tAnglePsi/rad);

    //if(VerboseLevel>1)
    //G4cout << "COULEX::PSDI: Projectile lab polar angle: "<< std::scientific << std::setprecision(3)<<pAnglePsi/degree<<"°, recoil lab polar angle: "<<std::scientific << std::setprecision(3)<<tAnglePsi/degree<<"°"<<G4endl<<"	      Projectile lab azimuthal angle: "<< std::scientific << std::setprecision(3)<<pAnglePhi/degree<<"°, recoil lab azimuthal angle: "<< std::scientific << std::setprecision(3)<<tAnglePhi/degree<<"°"<<G4endl;


    //projectile and recoil lab total energies after scattering

    //Take care of branching in inverse kinematics: occurs if rho1,2>1.
    //Branching means there are 2 COM angles for 1 LAB angle. The branches
    //merge at the point where the lab angle as a function of the COM angle is maximum.
    //If branching occurs, for angles below and above the critical angle (where the branches merge)
    //different signs have to be used in the energy formula.

    G4double pAngleMaxCOM = acos(-1./rho1)*rad;
  // 	  G4double pAngleMaxLAB = atan(sin(pAngleMaxCOM/rad)*sqrt(1-V2) / (cos(pAngleMaxCOM/rad)+rho1))*rad;

    G4double tAngleMaxCOM = acos(1./rho2)*rad;
  // 	  G4double tAngleMaxLAB = atan(sin(tAngleMaxCOM/rad)*sqrt(1-V2) / (((-1.)*cos(tAngleMaxCOM/rad))+rho2))*rad;
  // 	  G4double pAngleCrit = atan(sin(tAngleMaxCOM/rad)*sqrt(1-V2) / (cos(tAngleMaxCOM/rad)+rho1))*rad;
  // 	  G4double tAngleCrit = atan(sin(pAngleMaxCOM/rad)*sqrt(1-V2) / (((-1.)*cos(pAngleMaxCOM/rad))+rho2))*rad;


    if(rho1>1 && pAngleTheta>pAngleMaxCOM){		//branching for beam-like particle?
	    pEg = A1*Eg - pi*t1cos*std::sqrt(A1*A1 -4*m1_2*(Eg2 - pi2*t1cos*t1cos));
    }else{
	    pEg = A1*Eg + pi*t1cos*std::sqrt(A1*A1 -4*m1_2*(Eg2 - pi2*t1cos*t1cos));
    }
    if(rho2>1 && pAngleTheta<tAngleMaxCOM){		//branching for target-like particle?
	    tEg = A2*Eg - pi*t2cos*std::sqrt(A2*A2 -4*mii2*(Eg2 - pi2*t2cos*t2cos));
    }else{
	    tEg = A2*Eg + pi*t2cos*std::sqrt(A2*A2 -4*mii2*(Eg2 - pi2*t2cos*t2cos));
    }

    pEg /= 2*(Eg2 - pi2*t1cos*t1cos);
    tEg /= 2*(Eg2 - pi2*t2cos*t2cos);

    pKinEScat=(pEg-m1)*MeV;
    tKinEScat=(tEg-mii)*MeV;

    //if(VerboseLevel>1)
    //if(VerboseLevel>0)
    //G4cout << "COULEX::PSDI: coulex at "<< std::fixed << std::setprecision(3) << pKinE/MeV << "MeV, recoil kinetic Energy: " << tKinEScat/MeV <<"MeV, scattered projectile kinetic Energy: "<< std::scientific << std::setprecision(3)<<pKinEScat/MeV<<"MeV"<<G4endl;

     // a file showing the kinematics of the reaction
     if(WriteOutKinematics){
	std::ofstream of("kine.dat", std::ios::app); // output file
	of << pAnglePsi/degree << "	"<< tAnglePsi/degree << "	"<< pAngleTheta/degree<<"	"<<pKinEScat/MeV+tKinEScat/MeV + ExE/MeV - pKinE/MeV<<"	"<<pKinEScat/MeV<<"	"<<tKinEScat/MeV<<std::endl;
	of. close();
     }


/// Create recoiling target ion
    aParticleChange.Initialize(track);
    G4TouchableHandle trTouchable = track.GetTouchableHandle();

    /// Update projectile properties
    G4double sinpPsi = std::sin(pAnglePsi/rad);
    G4ThreeVector ProjectileNewDirection(std::cos(pAnglePhi/rad)*sinpPsi,std::sin(pAnglePhi/rad)*sinpPsi,std::cos(pAnglePsi/rad));
    //ProjectileNewDirection.rotateUz(pDir);
    //ProjectileNewDirection = G4ThreeVector(0,0,1);
    //ProjectileNewDirection = ProjectileNewDirection.unit();

    G4ParticleTable* particleTable2;
    particleTable2 = G4ParticleTable::GetParticleTable();

    G4ParticleDefinition* gamma = particleTable2->FindParticle("gamma");
    G4double theta = CLHEP::twopi*G4UniformRand(), phi = acos(2*G4UniformRand()-1);
    G4double ux = std::sin(phi)*std::cos(theta);
    G4double uy = std::sin(phi)*std::sin(theta);
    G4double uz = std::cos(phi);

    G4ThreeVector directions = G4ThreeVector(ux,uy,uz);
    //G4ThreeVector directions = G4ThreeVector(0,1,0);
    G4NistManager* nist = G4NistManager::Instance();

    if (pZ==28 && pA==64 && (track.GetMaterial() == nist->FindOrBuildMaterial("G4_Au"))) {
    //if (pZ==30 && pA==74) {

    G4cout << "COULEX::PSDI: Projectile Z="<<pZ<<", A="<<pA<<", Ekin="<< std::scientific << std::setprecision(3)<<pKinE/MeV<<"MeV"<<G4endl;

      G4NuclearLevelData* levelData = G4NuclearLevelData::GetInstance();
      const G4LevelManager* levelManager = levelData->GetLevelManager(pZ, pA);
      //G4int Nentries = levelManager->NumberOfTransitions()+1;

      directions = directions.unit();


      //G4double time_taken = (log(G4UniformRand())/log(0.5))*levelManager->LifeTime(1);
      G4double time = track.GetGlobalTime();
      //aParticleChange.ProposeGlobalTime(time+100);


      G4double lightspeed = 299.792; //mm/ns;
      G4double angle = directions.angle(ProjectileNewDirection);
      G4double beta = velocity/lightspeed;
      G4double factor = sqrt(1.0 - beta*beta)/(1-beta*cos(angle));
      G4double energy = levelManager->LevelEnergy(1)/MeV;
      //G4double energy = 0.5*MeV;
      G4double energy_gam = energy;
      energy = energy/factor;

      //G4cout<<beta<<G4endl;

      //G4ThreeVector gamma_position(position[0]+velocity*ProjectileNewDirection[0]*time_taken,position[1]+velocity*ProjectileNewDirection[1]*time_taken,position[2]+velocity*ProjectileNewDirection[2]*time_taken);
      G4ThreeVector gamma_position = track.GetPosition();

      std::ofstream doppler_out ("doppler_out.dat");
      doppler_out<<gamma_position.getX()<<" "<<gamma_position.getY()<<" "<<gamma_position.getZ()<<" "<<angle<<" "<<beta<<" "<<energy_gam<<" "<<energy<<" "<<factor<<" "<<directions.getX()<<" "<<directions.getY()<<" "<<directions.getZ()<<" "<<ProjectileNewDirection.getX()<<" "<<ProjectileNewDirection.getY()<<" "<<ProjectileNewDirection.getZ()<<" "<<pKinEScat/MeV<<" "<<pAnglePsi/degree<<" "<<pAngleTheta/degree;

  	  G4DynamicParticle* gamma_part = new G4DynamicParticle(gamma, directions, energy);

      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(gamma_part, gamma_position);


    G4double time_taken = 0.00333;
    G4ThreeVector ion_position = track.GetPosition();
    G4ThreeVector ion_step(ion_position[0]+velocity*ProjectileNewDirection[0]*time_taken,ion_position[1]+velocity*ProjectileNewDirection[1]*time_taken,ion_position[2]+velocity*ProjectileNewDirection[2]*time_taken);




  }

  aParticleChange.ProposeEnergy(pKinEScat);
  aParticleChange.ProposeLocalEnergyDeposit(pKinE-pKinEScat);
  //aParticleChange.ProposeLocalTime(100);
  aParticleChange.ProposeMomentumDirection(ProjectileNewDirection);

  return G4VDiscreteProcess::PostStepDoIt(track,step);

}





//Functions for calculation of maximum safe coulex angle http://web-docs.gsi.de/~wolle/EB_at_GSI/FRS-WORKING/PHYSICS/GRAZING/grazing.html

 G4double G4InelasticCoulombScattering::Ri(G4int &A){
  return  1.28*pow(A,(1.0/3.0))-0.76+0.8*(1/pow(A,(1.0/3.0)));
 }

G4double G4InelasticCoulombScattering::Ci(G4int &A){
  return Ri(A)*(1-(1/pow(Ri(A),2)));
 }

G4double G4InelasticCoulombScattering::DFactor(G4int &pA, G4int &tA){
  return  Ci(pA)+Ci(tA);
 }
G4double G4InelasticCoulombScattering::ComToLabFaktor(G4int &pA, G4int &tA){
   return (pA+tA)/tA;
 }
G4double G4InelasticCoulombScattering::aFactor(G4double &kinE, G4int &pA, G4int &tA, G4int &pZ, G4int &tZ){
 return ((1.44*pZ*tZ*(931.5+(kinE/MeV)/(pA*ComToLabFaktor(pA, tA))))/(pA*(1863*(kinE/MeV)/(pA*ComToLabFaktor(pA, tA))+pow(((kinE/MeV)/(pA*ComToLabFaktor(pA, tA))),2.0))));
}

G4double G4InelasticCoulombScattering::ThetaCOM(G4double &kinE, G4int &pA, G4int &tA, G4int &pZ, G4int &tZ){
  return  2*asin(aFactor(kinE, pA, tA, pZ, tZ)/(DFactor(pA, tA)-aFactor(kinE, pA, tA, pZ, tZ)))*rad;
 }

G4double G4InelasticCoulombScattering::UsedCoulombBar(G4int &pA, G4int &tA, G4int &pZ, G4int &tZ){
 return ((1.44*pZ*tZ)/DFactor(pA, tA))*((tA+pA)/tA);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


G4double G4InelasticCoulombScattering::GetXSecClx(G4double kinE, G4int &pZ, G4int &pA, G4int &tZ, G4int &tA){
 	xsecTempTableClx.clear();

	G4int ThetaCOMInt;

	TargetExcitationCLX=1+TargetExcitation;	// variable for target or projectile excitation

	//stringstream for conversion of strings to double/integer
	istringstream isst;

	std::ofstream of("clx.dat", std::ios::out); // output file
	of <<  std::endl;		//empty header to input file
	of << "00000000" << std::endl;	//output options to input filename
	of << "2	4	10" << std::endl;	//number of states to be considered, a useless number(max multipol?), max. M quantum number, default is 10 (must be greater then ground state spin
	of << "0	0	0" << std::endl;	//accuracy (use standart e-8), max Param XI (use standert 6), E1 polarization (use standart 0.005)
	//infos on projectile
	of << pZ << "	" << pA << std::endl;
	//infos on target
	of << tZ << "	" << tA<< std::endl;
	of << 1 << std::endl; //1 for projectile excitation, 2 for target excitation
	of << kinE/MeV<< std::endl;

	// The angle for safe coulex is calcutated by "UsedCoulombBar".
	// Beyond the coulomb barrier cross sections calculated by CLX are used if the safty condition is fulfilled.
	// This can be controlled by the SafeAngleDweiko-botton.

    //G4cout<<UsedCoulombBar(pA, tA, pZ, tZ)<<" "<<pA<<" "<<pZ<<" "<<tA<<" "<<tZ<<G4endl;

	if((kinE/MeV>UsedCoulombBar(pA, tA, pZ, tZ))&&(SafeCoulexClx==true)){
	   G4double ThetaRound= ThetaCOM(kinE, pA, tA, pZ, tZ)/degree+1.;
	   ThetaCOMInt=(G4int)ThetaRound;
	}else{
	 ThetaCOMInt=180;
	}


/*
130Xe
of << "1	"<<179<<"	"<< 1 << std::endl; //Theta_min, Theata_max, Delta_Theta. If you make it dynamic: change dimension of xsec-array!
// infos on nuclear states
//gs
of << "1	" <<0<< "	"<<0<<"	"<<1<<"	"<<0.0 << std::endl;
//first excited state
of << "2	"<<2.0<<"	"<<0.5361<<"	"<<1<<"	"<<0.0 << std::endl;
of << "3	"<<4.0<<"	"<<0.699<<"	"<<1<<"	"<<0.0 << std::endl;
of << "4	"<<6.0<<"	"<<0.739<<"	"<<1<<"	"<<0.0 << std::endl;
*/

// 	G4cout << "ThetaCOMInt=" << ThetaCOMInt << G4endl;

of << "1	"<<180<<"	"<< 1 << std::endl; //Theta_min, Theata_max, Delta_Theta. If you make it dynamic: change dimension of xsec-array!
// infos on nuclear states
//gs
of << "1	" <<0<< "	"<<0<<"	"<<1<<"	"<<0.0 << std::endl;
//first excited state
of << "2	"<<2.0<<"	"<<1.345<<"	"<<1<<"	"<<0.0 << std::endl;
//of << "2	"<<2.0<<"	"<<0.6059<<"	"<<1<<"	"<<0.0 << std::endl;

//of << "3	" <<4.0<< "	"<<3.308<<"	"<<1<<"	"<<0.0 << std::endl;
//of << "4	" <<4.0<< "	"<<10.356<<"	"<<1<<"	"<<0.0 << std::endl;
//of << "5	" <<4.0<< "	"<<11.0967<<"	"<<1<<"	"<<0.0 << std::endl;
//of << "4	"<<6.0<<"	"<<0.5437<<"	"<<1<<"	"<<0.0 << std::endl;
//of << "5	"<<8.0<<"	"<<0.9026<<"	"<<1<<"	"<<0.0 << std::endl;


/*
of << "1	"<<179<<"	"<< 1 << std::endl; //Theta_min, Theata_max, Delta_Theta. If you make it dynamic: change dimension of xsec-array!
// infos on nuclear states
//gs
of << "1	" <<0<< "	"<<0<<"	"<<1<<"	"<<0.0 << std::endl;
//first excited state
of << "2	"<<0.0<<"	"<<6.0494<<"	"<<1<<"	"<<0.0 << std::endl;
of << "3	"<<2.0<<"	"<<6.9171<<"	"<<1<<"	"<<0.0 << std::endl;
of << "4	"<<3.0<<"	"<<6.12989<<"	"<<-1<<"	"<<0.0 << std::endl;

of << "5	"<<1.0<<"	"<<7.11685<<"	"<<-1<<"	"<<0.0 << std::endl;
*/

//Matrix Elements, Assuming electric transition of multipole order lambda
//of <<"1	2	"<<2.065<<"	"<<2<< std::endl;// M(E<lambda> state1->state2)
//of <<"2	1	"<<0.611<<"	"<<3<< std::endl;
of <<"1	2	"<<0.11857<<"	"<<2<< std::endl;// M(E<lambda> state1->state2)
//of <<"1	2	"<<0.68<<"	"<<2<< std::endl;
//of <<"2	2	"<<0.6093<<"	"<<2<< std::endl;// M(E<lambda> state1->state2)
//of <<"2	3	"<<0.3529<<"	"<<2<< std::endl;// M(E<lambda> state1->state2)
//of <<"3	3	"<<0.4588<<"	"<<2<< std::endl;// M(E<lambda> state1->state2)
//of <<"2	1	"<<0.20405<<"	"<<2<< std::endl;// M(E<lambda> state1->state2)
//of <<"3	1	"<<0.566<<"	"<<2<< std::endl;
//of <<"4	2	"<<4.1688<<"	"<<2<< std::endl;
//of <<"5	2	"<<0.527<<"	"<<2<< std::endl;
//of <<"2	3	"<<3.280<<"	"<<2 << std::endl;	//M(E<lambda> state2->state2) (quadrupole moment for lambda=2)
//of <<"3	4	"<<4.298<<"	"<<2 << std::endl;
//of <<"4	5	"<<5.147<<"	"<<2 << std::endl;

	of.close();

	//run CLX
	(void)(system("rm clxout")+1);
	(void)(system("rm clxdcy.dat")+1);
	(void)(system("./clx > clxout")+1);


/// read from CLX output

	std::ifstream fi;
	fi.open("clxout"); // open clx outputile

	G4int anglecount=0;

	G4String OneLine, anglestring, xsecstring;
	G4double angleClx, xsecClx, xsec_sumClx=0.;
	G4String numbers = "0123456789";

	xsecdata tempdatClx;

	while(!fi.eof()){
	  getline(fi, OneLine);
	  if(OneLine.find("SCATTERING ANGLE ")!=std::string::npos){
	    //read scattering angle

	    anglestring =OneLine.substr(33, 41);
	    anglestring.append(" ");  // deliminator for the number

	    isst.clear();
	    isst.str(anglestring);
	    isst >> angleClx;
	    tempdatClx.angle=angleClx*degree;

	    getline(fi, OneLine);
	    getline(fi, OneLine);
	    getline(fi, OneLine);
	    getline(fi, OneLine);
	    getline(fi, OneLine);
	    //read xsec
	    ///loop from here if more than 1 excited state should be considered!
	    getline(fi, OneLine);
	    //lastnum=OneLine.find_last_of(numbers);

	    xsecstring = OneLine.substr(29, 39);

	    if(xsecstring.find("nan")!=std::string::npos) xsecClx = 0;
//	    else if(angleClx < (-0.345*kinE/MeV +210)) xsecClx = 0;
	    else{
	      	xsecstring.append(" ");  // deliminator for the number
		isst.clear();
		isst.str(xsecstring);
		isst >> xsecClx;
	    }

	    //numerical integration of the cross section over the solid angle element
	    //(excluding phi-integration: multiply by 2Pi on return)
//Here to edit cross section tempdatClx.xsec=1000*fabs(xsecClx)*barn*std::sin(tempdatClx.angle/rad)*deg2rad*angleBinning*TwoPi;


	    tempdatClx.xsec=fabs(xsecClx)*barn*std::sin(tempdatClx.angle/rad)*deg2rad*angleBinning*TwoPi;
	    xsecTempTableClx.push_back(tempdatClx);
	    xsec_sumClx+=tempdatClx.xsec;
	    ///loop until here
	    anglecount++;
	  }
	}

	tempdatClx.xsec=0*barn;
	for( int i=0; i*angleBinning<=180; i++){
	  //if((i*angleBinning) > ThetaCOMInt){
	    tempdatClx.angle=i*angleBinning*degree;
	    xsecTempTableClx.push_back(tempdatClx);

      //G4cout<<tempdatClx.xsec/barn<<" "<<tempdatClx.angle/degree<<" "<<i*angleBinning<<G4endl;

	  //}
	}

	fi.close();

	return xsec_sumClx;

}



G4double G4InelasticCoulombScattering::GetXSecDweiko(G4double kinE, G4int &pZ, G4int &pA, G4int &tZ, G4int &tA){
 	xsecTempTableDweiko.clear();

	//stringstream for conversion of strings to double/integer
	istringstream isst;
if (tA==1) {tA=2;}
	std::ofstream of("dweiko.in", std::ios::out); // output file
	of <<  std::endl;		//empty header to input file
	of << pA <<"      "<< pZ <<"    "<< 2 <<"      "<< 1 <<"      "<< (kinE/MeV)/pA  << std::endl;

	// IW=0(1) for projectile (target) excitation.
	// IOPM=1(0) for output (none) of optical model potentials.
	// IOELAS=(0)[1]{2} for (no) [center of mass] {laboratory} elastic scattering cross section.
	// IOINEL=(0)[1]{2} for (no) [center of mass] {laboratory} inelastic scattering cross section.
	// IOGAM=(0)[1]{2} for (no) output [output of statistical tensors] {output of  gamma-ray angular distributions}
	// IW=0(1)   IOPM=0(1)    IOELAS=0(1)[2]    IOINEL=0(1)[2]   IOGAM=0(1)[2]
	of <<0<<"         "<<0<<"            "<<0<<"                 "<<1<<"                "<<0<< std::endl;
	// NB=number of impact parameter points (NB <= NBMAX). currently 200
	// ACCUR=accuracy required for time integration at each impact parameter.
	// BMIN=minimum impact parameter (enter 0 for default)
	// IOB=1(0) prints (does not print) out impact parameter probabilities.
	// NB         ACCUR        BMIN[fm]     IOB=1(0)
	of <<200<<"        "<<0.001<<"        "<<0<<"           "<<0<< std::endl;
	// OMP switch:
	// IOPW=0 (no OMP)          IOPNUC=0 (no nuclear)
	// 1 (Woods-Saxon)            1 (vibrational excitations)
	// 2 (read)
	// 3 (t-rho-rho folding potential)
	// 4 (M3Y folding potential)
	// IOPW                       IOPNUC
	of <<"       "<<1<<"                          "<<0<< std::endl;

	// If IOPW=1, enter V0_ws [VI_ws] = real part [imaginary] (>0, both) of Woods-Saxon.
        // r0_ws [r0I_ws] = radius parameter (R_ws = r0 * (ap^1/3 + at^1/3).
        // d_ws [dI_ws]  = diffuseness.
	// If IOPW is not equal to 1, place a '#' sign at the beginning of this line, or delete it.
	// V0 [MeV]     r0[fm]      d[fm]      VI [MeV]     r0_I [fm]     dI [fm]
  of <<50.<<"          "<<1.067<<"       "<<0.8<<"        "<<58.<<"          "<<1.067<<"         "<<0.8<< std::endl;

	//The angle for safe coulex is calcutated by "UsedCoulombBar". Use the double angle range for Dweiko only.
	// This can be controlled by the SafeAngleDweiko-botton.
	G4int ThetaCOMInt;
	if((kinE/MeV>UsedCoulombBar(pA, tA, pZ, tZ))&&(SafeCoulexDweiko==true)){
	   G4double ThetaRound= ThetaCOM(kinE, pA, tA, pZ, tZ)/degree+1.;
	   ThetaCOMInt=(G4int) ThetaRound;
	}else{
	 ThetaCOMInt=maxThetaCOM;
	}
	// If IOELAS=1,2 or IOINEL=1,2 enter here THMAX, maximum angle (in degrees and in the
	// center of mass), and NTHETA (<= NGRID), the number of points in scatering angle.
	// If IOELAST or IOINEL are not 1, or 2, place a '//' sign at the beginning of this line
	// THMAX       NTHETA
  //ThetaCOMInt=180;


	//of <<ThetaCOMInt <<"	"<<ThetaCOMInt/angleBinning<< std::endl;
  of <<45. <<"          "<<100<< std::endl;

	// If IOINEL=1 enter the state (JINEL) for the inelastic angular distribution.
	// If IOINEL is not 1, or 2, place a '//' sign at the beginning of this line

	// JINEL
  of <<2<< std::endl;
	// NST: number of nuclear levels (<= NSTMAX).
	of <<2<<std::endl;
	// state label (I), energy (EX), and  spin (SPIN).
	// I ranges from 1 to NST.
	// I     Ex[MeV]    SPIN
  of <<1 <<"        "<<0.000<<"      "<<0<< std::endl;
	//first excited state
  of << 2<<"        "<<1.345<<"	    "<<2.0<<std::endl;
  //of << 2<<"        "<<0.6059<<"	    "<<2.0<<std::endl;
	// Reduced matrix elements for E1, E2, E3, M1 and M2 excitations:
	// <I_j||O(E/M;L)||I_i>,      j > i , for the electromagnetic transitions.
	// To stop, add a row of zeros at the end of this list

	//if((lambda==1)&&((p0*p1)==-1)){
	// i -> j   *E1[e fm]*  E2[e fm^2]  E3[e fm^3]   M1[e fm]   M2[e fm^2]
	of <<1<<"    "<<2<<"   "<<0<<"         "<<0.11857<<"      "<<0<<"            "<<0<<"          "<<0<<std::endl;
  //of <<1<<"    "<<2<<"   "<<0<<"         "<<0.68<<"      "<<0<<"            "<<0<<"          "<<0<<std::endl;
  of <<0<<"    "<<0<<"   "<<0<<"         "<<0<<"      "<<0<<"            "<<0<<"          "<<0<<std::endl;
  of <<std::endl;
	//of <<"0	"<<"0	"<<"0 "<<"	0	"<<"0	"<<"0	"<<"0	"<< std::endl;
	//}
  /*
	if((lambda==1)&&((p0*p1)==1)){
	// i -> j   E1[e fm]  E2[e fm^2]  E3[e fm^3]   *M1[e fm]*   M2[e fm^2]
	of <<"1	"<<"2	"<<"0	"<<"0	"<<"0	"<<M12<<"	0	"<< std::endl;
	of <<"2	"<<"2	"<<"0	"<<"0	"<<"0	"<<M22<<"	0	"<< std::endl;
	of <<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<< std::endl;
	}
	if((lambda==2)&&((p0*p1)==1)){
	// i -> j   E1[e fm]  *E2[e fm^2]*  E3[e fm^3]   M1[e fm]   M2[e fm^2]
	of <<"1	"<<"2	"<<"0	"<<M12<<"	0	"<<"0	"<<"0	"<< std::endl;
	of <<"2	"<<"2	"<<"0	"<<M22<<"	0	"<<"0	"<<"0	"<< std::endl;
	of <<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<< std::endl;
	}
	if((lambda==2)&&((p0*p1)==-1)){
	// i -> j   E1[e fm]  E2[e fm^2]  E3[e fm^3]   M1[e fm]   *M2[e fm^2]*
	of <<"1	"<<"2	"<<"0	"<<"0	"<<"0	"<<"0	"<<M12<< std::endl;
	of <<"2	"<<"2	"<<"0	"<<"0	"<<"0	"<<"0	"<<M22<< std::endl;
	of <<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<< std::endl;
	}
	if((lambda==3)&&((p0*p1)==-1)){
	// i -> j   E1[e fm]  E2[e fm^2]  E3[e fm^3]   M1[e fm]   *M2[e fm^2]*
	of <<"1	"<<"2	"<<"0	"<<"0	"<<M12<<"	0	"<<"0	"<< std::endl;
	of <<"2	"<<"2	"<<"0	"<<"0	"<<M22<<"	0	"<<"0	"<< std::endl;
	of <<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<<"0	"<< std::endl;


	}
*/

	of.close();

	///run DWEIKO
	//ugly construction to avoid compiler warning about unsused return value..
	(void)(system("rm dweiko_inel.out")+1);
	(void)(system("rm dweiko.out")+1);
	(void)(system("./dweiko")+1);


	/// read from DWEIKO output

	std::ifstream fi;
	fi.open("dweiko_inel.out"); // open dweiko outputile

	G4int anglecount=0;

	G4String OneLine, anglestring, xsecstring;
	G4double angleDweiko, xsecDweiko;
	G4String numbers = "0123456789";

	//Dweiko doesn't calculate at theta=0, so we have to push back an entry for angle 0 with 0 cross-sections
	xsecdata tempdatDweiko;
	tempdatDweiko.xsec=0*barn;
	tempdatDweiko.angle=0*degree;
	xsecTempTableDweiko.push_back(tempdatDweiko);

	//process the DWEIKO output file
	while(!fi.eof()){

 	  if(anglecount==0){
 	     getline(fi, OneLine);
 	     getline(fi, OneLine);
 	  }
	  else{
	    getline(fi, OneLine);
	  }

	  if(OneLine.size() < 11) continue;

	  anglestring =OneLine.substr(1, 10);
	  xsecstring = OneLine.substr(11, OneLine.size());

	  if(anglestring.find("nan")!=std::string::npos) angleDweiko = 0;
	  else{
	      	anglestring.append(" ");  // deliminator for the number
		isst.str(anglestring);
		isst >> angleDweiko;
	  }
	  if(xsecstring.find("nan")!=std::string::npos) xsecDweiko = 0;
	  else{
	      	xsecstring.append(" ");  // deliminator for the number
		isst.str(xsecstring);
		isst >> xsecDweiko;
	  }

	  //numerical integration of the cross section over the solid angle
	  //(excluse phi-integration: multiply by 2Pi on return)
 	  tempdatDweiko.angle=angleDweiko*degree;
 	  tempdatDweiko.xsec=xsecDweiko*millibarn*std::sin(tempdatDweiko.angle/rad)*deg2rad*angleBinning*TwoPi;
 	  xsecTempTableDweiko.push_back(tempdatDweiko);
 	  ///loop until here
 	  anglecount++;
	}

	tempdatDweiko.xsec=0*barn;
	for( int i=0; i*angleBinning<=180; i++){
	  if((i*angleBinning) > ThetaCOMInt){
	    tempdatDweiko.angle=i*angleBinning*degree;
	    xsecTempTableDweiko.push_back(tempdatDweiko);
	  }
	}

	fi.close();

	fi.open("dweiko.out"); // open other dweiko outputile
	G4String LastLine, ThisLine;
	if(fi.is_open()) {
	  while(!fi.eof()){
	      LastLine = ThisLine;
 	     getline(fi, ThisLine);

	  }
	}

	fi.close();

        G4String TotalXSec;
        TotalXSec = LastLine.substr(17, LastLine.size());

	isst.str(TotalXSec);
	G4double XSec_Sum;
	isst >> XSec_Sum;
	XSec_Sum*=millibarn;

	return XSec_Sum;
}
