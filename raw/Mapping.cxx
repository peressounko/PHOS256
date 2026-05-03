// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Mapping.cxx
/// \author Dmitri Peresunko

#include <fstream>
#include <iostream>
#include "TSystem.h"
#include "Mapping.h"
#include "Geometry.h"

//_______________________________________________________
Mapping::ErrorStatus Mapping::hwToAbsId(short hwAddr, short& absId, CaloFlag& caloFlag)
{

  if (!mInitialized) {
    SetMapping();
  }

  if (hwAddr < 0 || hwAddr >= NMaxHWAddress) {
    return kWrongHWAddress;
  }

  // transform
  absId = mAbsId[hwAddr];
  caloFlag = mCaloFlag[hwAddr];

  return kOK;
}
//_______________________________________________________
int Mapping::GetXcell(short hw)
{
  short absId;
  CaloFlag caloFlag;
  hwToAbsId(hw, absId, caloFlag);
  int relid[2] = {0};
  Geometry::AbsToRelNumbering(absId, relid);
  return relid[0];
}
//_______________________________________________________
int Mapping::GetZcell(short hw)
{
  short absId;
  CaloFlag caloFlag;
  hwToAbsId(hw, absId, caloFlag);
  int relid[2] = {0};
  Geometry::AbsToRelNumbering(absId, relid);
  return relid[1];
}
//_______________________________________________________
int Mapping::GetCaloFlag(short hw)
{
  short absId;
  CaloFlag caloFlag;
  hwToAbsId(hw, absId, caloFlag);
  return caloFlag;
}
//_______________________________________________________
Mapping::ErrorStatus Mapping::absIdTohw(short absId, short caloFlag, short& hwAddr)
{

  if (caloFlag < 0 || caloFlag > 2) {
    hwAddr = 0;
    return kWrongCaloFlag;
  }
  if (caloFlag < 2) {
    if (absId <= 0 || absId > NCHANNELS) {
      hwAddr = 0;
      return kWrongAbsId;
    }
  }

  if (!mInitialized) {
    std::cout << "Mapping not initialized" << std::endl;
    return kNotInitialized;
  }

  hwAddr = mAbsToHW[absId - 1][caloFlag];
  return kOK;
}
//_______________________________________________________
Mapping::ErrorStatus Mapping::SetMapping(std::string_view path)
{
  // Read mapping from data file
  mPath = path;

  std::string p;
  if (mPath.empty()) { // use default path
    p = gSystem->Getenv("PHOS256_ROOT");
    p += "/share/Mapping";
  } else {
    p = mPath.data();
  }

  short numberOfChannels = 0;
  short maxHWAddress = 0;
  std::string fname = p + "/Mapping.data";
  std::ifstream fIn(fname);
  if (!fIn.is_open()) {
    std::cout << "Missing mapping file " << fname << std::endl;
    abort();
    return kNotInitialized;
  }
  if (!(fIn >> numberOfChannels)) {
    std::cout << "Syntax of mapping file " << p << " is wrong: no numberOfChannels" << std::endl;
    abort();
    return kNotInitialized;
  }
  if (numberOfChannels != NHWPERDDL) {
    std::cout << "Unexpected number of channels: " << numberOfChannels << " expecting " << NHWPERDDL << " file " << p << " is wrong: no numberOfChannels" << std::endl;
    abort();
    return kNotInitialized;
  }
  if (!(fIn >> maxHWAddress)) {
    std::cout << "Syntax of mapping file " << p << " is wrong: no maxHWAddress" << std::endl;
    abort();
    return kNotInitialized;
  }
  if (maxHWAddress >= NMaxHWAddress) {
    std::cout << "Maximal HW address in file " << maxHWAddress << "larger than array size " << NMaxHWAddress << "for " << p << std::endl;
    abort();
    return kNotInitialized;
  }
  for (short ich = 0; ich < numberOfChannels; ich++) { // 512 = 2*256 channels
    int hwAddress;
    if (!(fIn >> hwAddress)) {
      std::cout << "Syntax of mapping file " << p << " is wrong: no HWadd for ch " << ich << std::endl;
      abort();
      return kNotInitialized;
    }
    if (hwAddress > maxHWAddress) {
      std::cout << "Hardware (ALTRO) adress (" << hwAddress << ") outside the range (0 -> " << maxHWAddress << ") !" << std::endl;
      return kNotInitialized;
    }
    int row, col, caloFlag;
    if (!(fIn >> row >> col >> caloFlag)) {
      std::cout << "Syntax of mapping file " << p << " is wrong:  no (raw col caloFlag)" << std::endl;
      abort();
      return kNotInitialized;
    }

    if (caloFlag < 0 || caloFlag > 2) {
      std::cout << "Wrong CaloFlag value found (" << caloFlag << "). Should be 0, 1 !" << std::endl;
      abort();
      return kNotInitialized;
    }

    // convert col,raw caloFlag to AbsId
    // module numbering:
    // absId:
    // start from 1 till 16*16=256. Numbering in each module starts at bottom left and first go in z direction:
    //  16   32     256
    //  ...  ...    ...
    //  1    17 ... 241
    //  relid[2]: (iphi[1...16], iz[1...16])
    int absId = 0;
    if (caloFlag < 2) { // readout channels
      int relid[2] = {row, col};
      Geometry::RelToAbsId(relid, absId);
    }

    mAbsId[hwAddress] = absId;
    mCaloFlag[hwAddress] = (CaloFlag)caloFlag;
    mAbsToHW[absId - 1][caloFlag] = hwAddress;
  }
  fIn.close();
  mInitialized = true;
  return kOK;
}
