/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/
//      ------------ BDSSynchrotronRadiation physics process --------
//                     by Grahame Blair, 18 October 2001

#include "BDSGlobalConstants.hh" // must be first in include list
#include "BDSSynchrotronRadiation.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "BDSAcceleratorComponent.hh"


extern G4int event_number;

typedef list<BDSAcceleratorComponent*>  BDSBeamline;
extern BDSBeamline theBeamline;


BDSSynchrotronRadiation::BDSSynchrotronRadiation(const G4String& processName)
  : G4VDiscreteProcess(processName)
     // initialization
{
  nExpConst=5*fine_structure_const/(2*sqrt(3))/electron_mass_c2;
  CritEngFac=3./2.*hbarc/pow(electron_mass_c2,3);

} 
 

G4VParticleChange* 
BDSSynchrotronRadiation::PostStepDoIt(const G4Track& trackData,
				     const G4Step& stepData)
{
  aParticleChange.Initialize(trackData);

  
  // tmp
  //    return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);


  G4double eEnergy=trackData.GetTotalEnergy();

  G4double R=BDSLocalRadiusOfCurvature;

  G4double NewKinEnergy = trackData.GetKineticEnergy();
 
  G4double GamEnergy=0;


  if(abs(R)>0)
    GamEnergy=SynGenC(BDSGlobals->GetSynchLowX())*
      CritEngFac*pow(eEnergy,3)/abs(R);

  NewKinEnergy -= GamEnergy;

  if(GamEnergy>0)
    {
      if((BDSGlobals->GetSynchTrackPhotons())&&
	 (GamEnergy>BDSGlobals->GetSynchLowGamE()) )
	{
	  G4DynamicParticle* aGamma= 
	    new G4DynamicParticle (G4Gamma::Gamma(), 
				   trackData.GetMomentumDirection(),
				   GamEnergy);

	  aParticleChange.SetNumberOfSecondaries(1);
	  aParticleChange.AddSecondary(aGamma); 
	}

      /* remove logging for now:
      BDSBeamline::const_iterator iBeam;
      G4double zpos=trackData.GetPosition().z();

      for(iBeam=theBeamline.begin();
	  iBeam!=theBeamline.end() && zpos>=(*iBeam)->GetZUpper(); 
	  iBeam++){}

      if(*iBeam)(*iBeam)->AddSynchEnergyLoss(GamEnergy);
 
      */


     /*
      G4cout.precision(8);
      G4cout<<" eGam="<<GamEnergy/GeV<<G4endl;
      G4cout<<" trackdata:"<<trackData.GetVolume()->GetName()<<G4endl;
      G4cout<<" beamline:"<<(*iBeam)->GetName()<<" z="<<zpos/m<<G4endl;
      */

      if (NewKinEnergy > 0.)
	{
	  //
	  // Update the incident particle 
	    aParticleChange.SetEnergyChange(NewKinEnergy);
	} 
      else
	{ 
	  aParticleChange.SetEnergyChange( 0. );
	  aParticleChange.SetLocalEnergyDeposit (0.);
	  G4double charge= trackData.GetDynamicParticle()->GetCharge();
	  if (charge<0.) aParticleChange.SetStatusChange(fStopAndKill);
	  else       aParticleChange.SetStatusChange(fStopButAlive);
	}
    }
  return G4VDiscreteProcess::PostStepDoIt(trackData,stepData);
}


BDSSynchrotronRadiation::~BDSSynchrotronRadiation(){
}

//--------------------------------------------------------------------
// routines from H. Burkhardt (Lep Note 632)
G4double BDSSynchrotronRadiation::SynGenC(G4double xmin)
// C++ version
// synchrotron radiation photon spectrum generator
// returns photon energy in units of the critical energy
// xmin is the lower limit for the photon energy to be generated
//          see H.Burkhardt, LEP Note 632
{ static bool DoInit=true;
  static G4double a1,a2,c1,xlow,ratio,LastXmin=-1.;
  
  if(DoInit | xmin!=LastXmin)
  { DoInit=false;
     LastXmin=xmin;
     if(xmin>1.) xlow=xmin; else xlow=1.;
     // initialize constants used in the approximate expressions
     // for SYNRAD   (integral over the modified Bessel function K5/3)
    a1=SynRadC(1.e-38)/pow(1.e-38,-2./3.); // = 2**2/3 GAMMA(2/3)
    a2=SynRadC(xlow)/exp(-xlow);
    c1=pow(xmin,1./3.);
     // calculate the integrals of the approximate expressions
    if(xmin<1.) // low and high approx needed
    { G4double sum1ap=3.*a1*(1.-pow(xmin,1./3.)); //     integral xmin --> 1
      G4double sum2ap=a2*exp(-1.);                    //     integral 1 --> infin
      ratio=sum1ap/(sum1ap+sum2ap);
    }
    else // only high approx needed
    {  ratio=0.;     //     generate only high energies using approx. 2
       a2=SynRadC(xlow)/exp(-xlow);
    }
  }
  // Init done, now generate
  G4double appr,exact,result;
  do
  { if(G4UniformRand()<ratio)           // use low energy approximation
     { result=c1+(1.-c1)*G4UniformRand();
       result*=result*result;  // take to 3rd power;
       exact=SynRadC(result);
       appr=a1*pow(result,-2./3.);
     }
     else                      // use high energy approximation
     { result=xlow-log(G4UniformRand());
       exact=SynRadC(result);
       appr=a2*exp(-result);
     }
  }
  while(exact<appr*G4UniformRand());    // reject in proportion of approx
  return result;               // result now exact spectrum with unity weight
}

G4double BDSSynchrotronRadiation::SynRadC(G4double x)
/*   x :    energy normalized to the critical energy
     returns function value SynRadC   photon spectrum dn/dx
     (integral of modified 1/3 order Bessel function)
     principal: Chebyshev series see H.H.Umstaetter CERN/PS/SM/81-13 10-3-1981
     see also my LEP Note 632 of 12/1990
     converted to C++, H.Burkhardt 21-4-1996    */
{ G4double synrad=0.;
  if(x>0. & x<800.) // otherwise result synrad remains 0
  { if(x<6.)
     { G4double a,b,z;
       z=x*x/16.-2.;
       b=          .00000000000000000012;
       a=z*b  +    .00000000000000000460;
       b=z*a-b+    .00000000000000031738;
       a=z*b-a+    .00000000000002004426;
       b=z*a-b+    .00000000000111455474;
       a=z*b-a+    .00000000005407460944;
       b=z*a-b+    .00000000226722011790;
       a=z*b-a+    .00000008125130371644;
       b=z*a-b+    .00000245751373955212;
       a=z*b-a+    .00006181256113829740;
       b=z*a-b+    .00127066381953661690;
       a=z*b-a+    .02091216799114667278;
       b=z*a-b+    .26880346058164526514;
       a=z*b-a+   2.61902183794862213818;
       b=z*a-b+  18.65250896865416256398;
       a=z*b-a+  92.95232665922707542088;
       b=z*a-b+ 308.15919413131586030542;
       a=z*b-a+ 644.86979658236221700714;
       G4double p;
       p=.5*z*a-b+  414.56543648832546975110;
       a=          .00000000000000000004;
       b=z*a+      .00000000000000000289;
       a=z*b-a+    .00000000000000019786;
       b=z*a-b+    .00000000000001196168;
       a=z*b-a+    .00000000000063427729;
       b=z*a-b+    .00000000002923635681;
       a=z*b-a+    .00000000115951672806;
       b=z*a-b+    .00000003910314748244;
       a=z*b-a+    .00000110599584794379;
       b=z*a-b+    .00002581451439721298;
       a=z*b-a+    .00048768692916240683;
       b=z*a-b+    .00728456195503504923;
       a=z*b-a+    .08357935463720537773;
       b=z*a-b+    .71031361199218887514;
       a=z*b-a+   4.26780261265492264837;
       b=z*a-b+  17.05540785795221885751;
       a=z*b-a+  41.83903486779678800040;
       G4double q;
       q=.5*z*a-b+28.41787374362784178164;
       G4double y;
       G4double const twothird=2./3.;
       y=pow(x,twothird);
       synrad=(p/y-q*y-1.)*1.81379936423421784215530788143;
     }
     else // 6 < x < 174
     { G4double a,b,z;
       z=20./x-2.;
       a=      .00000000000000000001;
       b=z*a  -.00000000000000000002;
       a=z*b-a+.00000000000000000006;
       b=z*a-b-.00000000000000000020;
       a=z*b-a+.00000000000000000066;
       b=z*a-b-.00000000000000000216;
       a=z*b-a+.00000000000000000721;
       b=z*a-b-.00000000000000002443;
       a=z*b-a+.00000000000000008441;
       b=z*a-b-.00000000000000029752;
       a=z*b-a+.00000000000000107116;
       b=z*a-b-.00000000000000394564;
       a=z*b-a+.00000000000001489474;
       b=z*a-b-.00000000000005773537;
       a=z*b-a+.00000000000023030657;
       b=z*a-b-.00000000000094784973;
       a=z*b-a+.00000000000403683207;
       b=z*a-b-.00000000001785432348;
       a=z*b-a+.00000000008235329314;
       b=z*a-b-.00000000039817923621;
       a=z*b-a+.00000000203088939238;
       b=z*a-b-.00000001101482369622;
       a=z*b-a+.00000006418902302372;
       b=z*a-b-.00000040756144386809;
       a=z*b-a+.00000287536465397527;
       b=z*a-b-.00002321251614543524;
       a=z*b-a+.00022505317277986004;
       b=z*a-b-.00287636803664026799;
       a=z*b-a+.06239591359332750793;
       G4double p;
       p=.5*z*a-b    +1.06552390798340693166;
       G4double const pihalf=pi/2.;
       synrad=p*sqrt(pihalf/x)/exp(x);
     }
  }
  return synrad;
}


