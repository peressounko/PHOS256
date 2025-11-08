// Generates N particles
//

#include <math.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TVirtualMC.h>
// #include <TVirtualMCApplication.h>
// #include <TVirtualMCStack.h>

#include "GenBox.h"
#include "Stack.h"

//_____________________________________________________________________________
GenBox::GenBox(TVirtualMCStack* stack)
  : TObject(),
    fStack(stack),
    fNofPrimaries(1),
    fIpart(22),
    fPtMin(0.1),
    fPtMax(1.),
    fThetaMin(10.),
    fThetaMax(20.),
    fPhiMin(80.),
    fPhiMax(110.)
{
}

//_____________________________________________________________________________
void GenBox::Generate()
{

  double polar[3] = {0, 0, 0};

  double origin[3] = {0., 0., 0.};
  double time = 0.;
  double p[3];

  // Polarization
  double polx = 0.;
  double poly = 0.;
  double polz = 0.;

  // Option: to be tracked
  int toBeDone = 1;

  int ntr;
  double m = TDatabasePDG::Instance()->GetParticle(fIpart)->Mass();

  const double degToRad = M_PI / 180.;

  for (int i = 0; i < fNofPrimaries; i++) {
    double random[3];
    gRandom->RndmArray(3, random);

    // usniform solid angle d(cos Theta) d(phi)
    double pmom = fPtMin + random[0] * (fPtMax - fPtMin);
    double cosMax = std::cos(degToRad * fThetaMin);
    double cosMin = std::cos(degToRad * fThetaMax);
    double cosTheta = cosMin + random[1] * (cosMax - cosMin);
    double phi = degToRad * (fPhiMin + random[2] * (fPhiMax - fPhiMin));
    double pt = pmom * std::sqrt(1. - cosTheta * cosTheta);
    double e = std::sqrt(pmom * pmom + m * m);

    p[0] = pt * std::cos(phi);
    p[1] = pt * std::sin(phi);
    p[2] = pmom * cosTheta;

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
    fStack->PushTrack(toBeDone, -1, fIpart, p[0], p[1], p[2], e, origin[0], origin[1], origin[2], time, polx,
                      poly, polz, kPPrimary, ntr, 1., 0);
    static_cast<Stack*>(TVirtualMC::GetMC()->GetStack())->StoreTrack(i);
  }
}
