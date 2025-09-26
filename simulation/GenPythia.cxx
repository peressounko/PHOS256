// Generates Pythia event
//

#include <math.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TVirtualMC.h>

#include "GenPythia.h"

GenPythia::GenPythia(const char* xmlDir)
{
  // Constructor with an xmlDir (eg "../xmldoc"
  fgXmldocPath = xmlDir;
  fPythia = new Pythia8::Pythia(xmlDir);
}

void GenPythia::SetPythiaSeed(UInt_t seed)
{
  //
  // set seed in PYTHIA8
  // NB. 900000000 is the maximum seed (0 is not allowed)
  //
  fPythia->readString("Random:setSeed = on");
  fPythia->readString(Form("Random:seed = %d", (seed % 900000000) + 1));

  if (fPythiaPartonLevel) {
    fPythiaPartonLevel->readString("Random:setSeed = on");
    fPythiaPartonLevel->readString(Form("Random:seed = %d", (seed % 900000000) + 1));
  }
}

// Pythia8::Pythia* AliTPythia8::GetPythiaPartonLevel()
// {
//   // Return the parton level object, used when events are superpositioned
//   // Basic initialization is done (can be changed by the user by retrieving the instance
//   // NOTE that the standard call AliTPythia8::Instnce() gives a different object, so check where the initialization should go

//   if (!fPythiaPartonLevel) {
//     if (fgXmldocPath != 0)
//       fPythiaPartonLevel    = new Pythia8::Pythia(fgXmldocPath);
//     else
//       fPythiaPartonLevel    = new Pythia8::Pythia();

//     // special settings to avoid color reconnection and hadronization
//     fPythiaPartonLevel->readString("HadronLevel:all = off");
//     fPythiaPartonLevel->readString("ColourReconnection:reconnect = off");
//   }

//   return fPythiaPartonLevel;
// }

//___________________________________________________________________________
bool GenPythia::Initialize(int idAin, int idBin, double eLab)
{
  // Initialization
  // AddParticlesToPdgDataBase();
  // UpdateParticleProperties();
  if (idAin < 100000000 && idBin < 100000000) { // pp
    fPythia->readString("Tune:pp = 5 ! 4C tune");
    fPythia->readString("Beams:frameType = 2");
    fPythia->readString(Form("Beams:eA = %13.4f", eLab));
    fPythia->readString("Beams:eB = 0.93827208944");
    fPythia->readString(Form("Beams:idA = %10d", idAin));
    fPythia->readString(Form("Beams:idB = %10d", idBin));
  } else {
    // Setup the beams.
    fPythia->readString("Tune:pp = 5 ! 4C tune");
    fPythia->readString(Form("Beams:idA = %10d", idAin));
    fPythia->readString(Form("Beams:idB = %10d", idBin));
    fPythia->readString("Beams:frameType = 2");
    fPythia->readString(Form("Beams:eA = %13.4f", eLab));
    // get mass of the target
    // int A = (idBin % 10000) / 10;
     // If the particle energy is smaller than its mass it is assumed to be at rest. 
    // fPythia->readString(Form("Beams:eB =  0.")); //, A * 0.93827208944));
    // // Initialize the Angantyr model to fit the total and semi-includive
    // // cross sections in Pythia within some tolerance.
    // fPythia->readString(
    //   "HeavyIon:SigFitErr = "
    //   "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
    // // These parameters are typicall suitable for sqrt(S_NN)=5TeV
    // fPythia->readString("HeavyIon:SigFitDefPar = 2.15,17.24,0.33");
    // // A simple genetic algorithm is run for 20 generations to fit the
    // // parameters.
    // fPythia->readString("HeavyIon:SigFitNGen = 20");
  }

  fPythia->readString("Random:setSeed = on");
  fPythia->readString("HadronLevel:Hadronize = on");

  fPythia->readString("HardQCD:all = on");
  fPythia->readString("SoftQCD:all = on");

  // if (fPythiaPartonLevel)
  //   if (fPythiaPartonLevel->init() == kFALSE)
  //     return false;

  return fPythia->init();
}

//_____________________________________________________________________________
void GenPythia::Generate()
{
  while (!fPythia->next())
    ; // Generate while produce correct event

  // Polarization
  double polx = 0.;
  double poly = 0.;
  double polz = 0.;

  // Option: to be tracked
  int toBeDone = 1;
  int ntr;
  std::vector<int> labels(fPythia->event.size());

  if (fPythia->event.size() == 0) {
    return;
  }

  fNofPrimaries = 0;
  for (int i = 0; i < fPythia->event.size(); i++) {
    if (fPythia->event[i].id() == 90) {
      labels.push_back(-1); // no such primary
      continue;
    }
    toBeDone = 1;
    if (fStore == kFinalOnly) {
      if (!fPythia->event[i].isFinal()) {
        labels.push_back(-1); // no such primary
        continue;
      }
    } else {
      if (!fPythia->event[i].isFinal()) {
        toBeDone = 0;
      }
    }
    if (abs(fPythia->event[i].id()) > 10000) { // Do not put fragments [100ZZZAAA0]
      labels.push_back(-1); // no such primary
      continue;
    }
    // Check acceptance
    const double radToDeg = 180. / M_PI;
    double phi = radToDeg * std::atan2(fPythia->event[i].py(), fPythia->event[i].px());
    double pt = std::sqrt(fPythia->event[i].px() * fPythia->event[i].px() + fPythia->event[i].py() * fPythia->event[i].py());
    double theta = radToDeg * std::atan2(pt, fPythia->event[i].pz());
    if (phi < fPhiMin || phi > fPhiMax || theta < fThetaMin || theta > fThetaMax) {
      labels.push_back(-1); // no such primary
      continue;
    }

    // Add particle to stack
    /// \param toBeDone  1 if particles should go to tracking, 0 otherwise
    /// \param parent    number of the parent track, -1 if track is primary
    /// \param pdg       PDG encoding
    /// \param px        particle momentum - x component [GeV/c]
    /// \param py        particle momentum - y component [GeV/c]
    /// \param pz        particle momentum - z component [GeV/c]
    /// \param e         total energy [GeV]
    /// \param vx        position - x component [cm]
    /// \param vy        position - y component  [cm]
    /// \param vz        position - z component  [cm]
    /// \param tof       time of flight [s]
    /// \param polx      polarization - x component
    /// \param poly      polarization - y component
    /// \param polz      polarization - z component
    /// \param mech      creator process VMC code
    /// \param ntr       track number (is filled by the stack)
    /// \param weight    particle weight
    /// \param is        generation status code
    fStack->PushTrack(toBeDone,
                      labels[fPythia->event[i].mother1()],
                      fPythia->event[i].id(),
                      fPythia->event[i].px(),    // [GeV/c]
                      fPythia->event[i].py(),    // [GeV/c]
                      fPythia->event[i].pz(),    // [GeV/c]
                      fPythia->event[i].e(),     // [GeV]
                      fPythia->event[i].xProd(), // [mm]
                      fPythia->event[i].yProd(), // [mm]
                      fPythia->event[i].zProd(), // [mm]
                      fPythia->event[i].tProd(), // [mm/c]
                      polx,
                      poly,
                      polz,
                      kPPrimary,
                      ntr,
                      1.,
                      0);
    labels.push_back(ntr);
    fNofPrimaries++;
  }
}
