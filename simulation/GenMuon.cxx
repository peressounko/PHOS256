// Generates one muon per event
// with realistic parameterisation of
// spectrum and angular distribution
// https://arxiv.org/pdf/1606.06907
// simplified for small energies

#include <math.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TVirtualMC.h>
// #include <TVirtualMCApplication.h>
#include <TVirtualMCStack.h>

#include "GenMuon.h"
#include "Stack.h"

//_____________________________________________________________________________
GenMuon::GenMuon(TVirtualMCStack* stack)
  : GenBox(stack)
{
}

//_____________________________________________________________________________
void GenMuon::Generate()
{
  // Generate one muon passing trigger plates
  //  we assume z axis as zenith axis
  //  muon spectrum and angular distributions are taken from
  //  https://arxiv.org/pdf/1606.06907
  //  we first generate muon passing first trigger plate in random place
  //  then require passing second plate if SetSecondTrigger() is called

  // Init
  if (fIntMax == 0.) {
    fIntMax = std::pow(kEmin + kE0, -knEn);
    fIpart = 13; // muon
  }

  double polar[3] = {0, 0, 0};

  double origin[3] = {0};
  double time = 0.;
  double p[3] = {0};
  double en = 0.;

  // Polarization
  double polx = 0.;
  double poly = 0.;
  double polz = 0.;

  // Option: to be tracked
  int toBeDone = 1;

  int ntr;
  const double m = 0.1056583755;

  const double degToRad = M_PI / 180.;

  bool passedTrig = false;
  while (!passedTrig) {
    double random[3];
    gRandom->RndmArray(3, random);

    // Generate energy
    //  (𝐸0 + 𝐸)−(𝑛+1)cos(theta)^n/(1 + 𝐸/854) : 854 GeV! therefore neglect factor (1+E/854)
    en = kEmin - kE0 + std::pow(fIntMax * (1. - random[0]), -1. / knEn);
    double pmom = 0.;
    if (en > m) {
      pmom = std::sqrt(en * en - m * m);
    } else {
      continue;
    }

    // Zenith angle (angle between particle and vertical)
    //[D. Pagano et al, Nucl.Instrum.Meth. A, 1014 (2021), 165732
    // 𝐼 = 𝐼0 𝑐𝑜𝑠𝑛(𝑝), 𝑛(𝑝) = 2.856 - 0.655ln(𝑝 𝑝0 ),
    // 𝑝0 = 1.0 GeV/c with p > 0.040 GeV/c
    const double p0 = 1.;

    double ntheta = 2.856 - 0.655 * std::log(pmom / p0);
    double theta = std::acos(std::pow(random[1], 1. / (ntheta + 1.)));
    double phi = random[2] * 2 * M_PI;

    p[0] = pmom * std::sin(theta) * std::cos(phi);
    p[1] = pmom * std::sin(theta) * std::sin(phi);
    p[2] = pmom * std::cos(theta);

    // test triggers
    // generate position in top trigger
    //  and test if it will pass through bottom
    gRandom->RndmArray(2, random);
    origin[0] = fTr1X0 + 2. * (random[0] - 0.5) * fTr1Dx;
    origin[1] = fTr1Y0 + 2. * (random[1] - 0.5) * fTr1Dy;
    origin[2] = 0;

    if (fUseSecondTr) {

      double scale = 0.;
      if (p[2] > 0) {
        scale = fDistTrigZ / p[2];
      }
      double xNew = origin[0] + p[0] * scale;
      double yNew = origin[1] + p[1] * scale;
      if (std::abs(xNew - fTr2X0) < fTr2Dx && std::abs(yNew - fTr2Y0) < fTr2Dy) {
        passedTrig = true;
      } else {
        passedTrig = false;
      }

    } else {
      passedTrig = true;
    }
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
  fStack->PushTrack(toBeDone, -1, fIpart, p[0], p[1], p[2], en, origin[0], origin[1], origin[2], time, polx,
                    poly, polz, kPPrimary, ntr, 1., 0);
  static_cast<Stack*>(TVirtualMC::GetMC()->GetStack())->StoreTrack(0); // single track
}
