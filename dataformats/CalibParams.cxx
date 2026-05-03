#include "Geometry.h"
#include <TH2.h>

#include <iostream>
#include "CalibParams.h"

CalibParams::CalibParams(int /*dummy*/)
{
  // produce reasonable objest for test purposes
  mGainCalib.fill(0.015);
  mHGLGRatio.fill(16.);
  mHGTimeCalib.fill(0.);
  mLGTimeCalib.fill(0.);
}

bool CalibParams::SetGain(TH2* h)
{
  const int MAXX = 16,
            MAXZ = 16;
  if (!h) {
    std::cout << " Error: no input histogam" << std::endl;
    return false;
  }

  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    std::cout << " Error: Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ << std::endl;
    return false;
  }

  int relid[2] = {1, 1};
  int absId;
  for (int ix = 1; ix <= MAXX; ix++) {
    relid[0] = ix;
    for (int iz = 1; iz <= MAXZ; iz++) {
      relid[1] = iz;
      Geometry::RelToAbsId(relid, absId);
      mGainCalib[absId - OFFSET] = h->GetBinContent(ix, iz);
    }
  }
  return true;
}

bool CalibParams::SetHGLGRatio(TH2* h)
{
  const int MAXX = 64,
            MAXZ = 56;
  if (!h) {
    std::cout << " Error: no input histogam" << std::endl;
    return false;
  }

  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    std::cout << " Error: Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ << std::endl;
    return false;
  }

  int relid[2] = {1, 1};
  int absId;
  for (int ix = 1; ix <= MAXX; ix++) {
    relid[0] = ix;
    for (int iz = 1; iz <= MAXZ; iz++) {
      relid[1] = iz;
      Geometry::RelToAbsId(relid, absId);
      mHGLGRatio[absId - OFFSET] = h->GetBinContent(ix, iz);
    }
  }
  return true;
}

bool CalibParams::SetHGTimeCalib(TH2* h)
{
  const int MAXX = 64,
            MAXZ = 56;
  if (!h) {
    std::cout << " Error: no input histogam" << std::endl;
    return false;
  }

  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    std::cout << " Error: Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ << std::endl;
    return false;
  }

  int relid[2] = {1, 1};
  int absId;
  for (int ix = 1; ix <= MAXX; ix++) {
    relid[0] = ix;
    for (int iz = 1; iz <= MAXZ; iz++) {
      relid[1] = iz;
      Geometry::RelToAbsId(relid, absId);
      mHGTimeCalib[absId - OFFSET] = h->GetBinContent(ix, iz);
    }
  }
  return true;
}

bool CalibParams::SetLGTimeCalib(TH2* h)
{
  const int MAXX = 64,
            MAXZ = 56;
  if (!h) {
    std::cout << " Error no input histogam" << std::endl;
    return false;
  }

  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    std::cout << " Error: Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ << std::endl;
    return false;
  }
  int relid[2] = {1, 1};
  int absId;
  for (int ix = 1; ix <= MAXX; ix++) {
    relid[0] = ix;
    for (int iz = 1; iz <= MAXZ; iz++) {
      relid[1] = iz;
      Geometry::RelToAbsId(relid, absId);
      if (absId - OFFSET < 0) { // non-existing part of a module 1
        continue;
      }
      mLGTimeCalib[absId - OFFSET] = h->GetBinContent(ix, iz);
    }
  }
  return true;
}
