/// \class Mapping
/// \brief Checks validity of hardware address (HW) and transform it to digit AbsId index
///
/// \author Dmitri Peresunko
/// \since Jan.2020
///

#ifndef MAPPING_H_
#define MAPPING_H_

#include <string>
#include <string_view>
#include "Rtypes.h" //for ClassDef

class Mapping
{
 public:
  enum ErrorStatus { kOK,
                     kWrongHWAddress,
                     kWrongAbsId,
                     kWrongCaloFlag,
                     kNotInitialized };
  static constexpr short NCHANNELS = 448;      ///< Number of channels starting from 1
  static constexpr short NHWPERDDL = 832;      ///< Number of HW addressed per DDL 13 FEC (first FEC missing)
  static constexpr short NMaxHWAddress = 1872; ///< Maximal HW address+1 (size of array)

  enum CaloFlag { kLowGain,
                  kHighGain };

  Mapping() = default;
  ~Mapping() = default;

  /// \brief convert hardware address to absId and caloFlag
  ErrorStatus hwToAbsId(short hw, short& absId, CaloFlag& caloFlag);
  /// \brief convert absId and caloflag to hardware address and ddl
  ErrorStatus absIdTohw(short absId, short caloFlag, short& hwAddr);

  int GetXcell(short hw);
  int GetZcell(short hw);
  int GetCaloFlag(short hw);

  ErrorStatus SetMapping(std::string_view path = "");

 protected:
  /// \brief Construct vector for conversion only if necessary
  ErrorStatus constructAbsToHWMatrix();

 private:
  std::string mPath = "";                         ///< path to mapping files
  bool mInitialized = false;                      ///< If conversion tables created
  bool mInvInitialized = false;                   ///< If inverse conversion tables created
  short mAbsId[NMaxHWAddress] = {0};              ///< Conversion table (ddl,branch,fec,chip,channel) to absId
  CaloFlag mCaloFlag[NMaxHWAddress] = {kLowGain}; ///< Conversion table (ddl,branch,fec,chip,channel) to absId
  short mAbsToHW[NCHANNELS][3] = {0};             ///< Conversion table [AbsId][caloFlag] to hw address

  ClassDefNV(Mapping, 1)
}; // End of Mapping

#endif
