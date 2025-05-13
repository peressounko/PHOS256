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

    double pt = fPtMin + random[0] * (fPtMax - fPtMin);
    double theta = degToRad * (fThetaMin + random[1] * (fThetaMax - fThetaMin));
    double phi = degToRad * (fPhiMin + random[2] * (fPhiMax - fPhiMin));
    double pmom = pt / TMath::Sin(theta);
    double e = TMath::Sqrt(pmom * pmom + m * m);

    p[0] = pt * TMath::Cos(phi);
    p[1] = pt * TMath::Sin(phi);
    p[2] = pmom * TMath::Cos(theta);

    // Add particle to stack
    fStack->PushTrack(toBeDone, -1, fIpart, p[0], p[1], p[1], e, origin[0], origin[1], origin[0], time, polx,
                      poly, polz, kPPrimary, ntr, 1., 0);
  }
}
