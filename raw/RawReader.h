#ifndef RAWREADER_H
#define RAWREADER_H
///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a date file or event.
///
///////////////////////////////////////////////////////////////////////////////

#include "Mapping.h"
#include <string>
#include <string_view>
#include "Rtypes.h" //for ClassDef

#include "event.h"

class TH2F;

class RawDataHeader;
class RawDataHeaderV3;

class RawReader
{
 public:
  RawReader(std::basic_string_view<char> fileName);
  virtual ~RawReader();

  // Operations with events etc
  bool Reset();
  bool NextEvent();   // read next event
  bool NextChannel(); // read next channel

  // Global getters
  uint GetEventType() const; // event type
  uint GetRunNumber() const;
  uint GetTimestamp() const;
  uint GetOrbitID() const;
  ushort GetBCID() const;
  ULong64_t GetEventIdAsLong() const;
  const uint* GetTriggerPattern() const;
  ULong64_t GetClassMask() const;
  const uint* GetEventId() const;
  int GetDataSize() const;
  uint GetPeriod() const;

  TH2F* GetErrorSummary() const { return fErrorSummary; };

  bool NextDDL();
  bool NextBunch();

  void EventSelectionType(int t = 7) { fSelectEventType = t; }

  const UChar_t* GetChannelPayload() const;
  UChar_t GetNPretriggerSamples() const { return (fAltroCFG2 >> 20) & 0xF; }
  UShort_t GetNSamplesPerCh() const { return (fAltroCFG2 >> 10) & 0x3FF; }
  UInt_t GetAltroCFG1() const { return fAltroCFG1; }
  const UShort_t* GetSignals() const { return fBunchDataPointer; } // Provide access to altro data itself
  UInt_t GetStartTimeBin() const { return fStartTimeBin; }         // Provide the index if the first time-bin in current bunch
  Int_t GetBunchLength() const { return fBunchLength; }            // Provide the current bunch length

  // Calo getters
  int GetHWAddress() const { return fHWAddress; }
  int GetModule() const { return fModule; }
  int GetCellX() const { return fXcell; }
  int GetCellZ() const { return fZcell; }

  enum CaloFlag { kLowGain = 0,
                  kHighGain = 1,
                  kTRUData = 2,
                  kLEDMonData = 3 };
  bool IsLowGain() const { return (fCaloFlag == kLowGain); }
  bool IsHighGain() const { return (fCaloFlag == kHighGain); }
  // bool  IsTRUData()      const {return (fCaloFlag == kTRUData)   ;}
  // bool  IsLEDMonData()   const {return (fCaloFlag == kLEDMonData);}

  int GetCaloFlag() const { return fCaloFlag; }

  void SetMapping(std::string fn = "./Mapping.data") { fMapping.SetMapping(fn); }

 protected:
  void ApplyMapping();
  bool ReadHeader();
  bool ReadNextData();
  bool ReadNextInt(uint& data);
  bool ReadNextShort(ushort& data);
  bool ReadNextChar(UChar_t& data);
  bool ReadNext(UChar_t* data, int size);
  int CheckData() const;
  uint SwapWord(uint x) const;
  ushort SwapShort(ushort x) const;
  uint Get32bitWord(int index) const;
  bool ReadRCUTrailer(int rcuVer);

  int GetEquipmentSize() const;
  int GetEquipmentType() const;
  int GetEquipmentId() const;
  const unsigned int* GetEquipmentAttributes() const;
  int GetEquipmentElementSize() const;
  int GetEquipmentHeaderSize() const;
  bool IsEventSelected() const;

  static constexpr int kMaxNTimeBins = 1024;
  enum { kErrMagic = 1,
         kErrNoDataHeader = 2,
         kErrSize = 4,
         kErrOutOfBounds = 8,
         kErrVersion = 16 };

  FILE* fFile = nullptr;                       // DATE file
  eventHeaderStruct* fEvent = nullptr;         // raw data super event
  eventHeaderStruct* fSubEvent = nullptr;      // raw data sub event
  equipmentHeaderStruct* fEquipment = nullptr; // raw data equipment header

  UChar_t* fPosition = nullptr; // current position in the raw data
  UChar_t* fEnd = nullptr;      // end position of the current data block

  int fErrorCode = 0;               // code of last error
  RawDataHeader* fHeader = nullptr; // current data header
  int fCount = 0;                   // counter of bytes to be read for current DDL

  // int            fSelectEquipmentType;  // type of selected equipment (<0 = no selection)
  // int            fSelectMinEquipmentId; // minimal index of selected equipment (<0 = no selection)
  // int            fSelectMaxEquipmentId; // maximal index of selected equipment (<0 = no selection)
  // bool           fSkipInvalid;          // skip invalid data
  int fSelectEventType = 7;           // type of selected events (<0 = no selection)
  ULong64_t fSelectTriggerMask = 0;   // trigger mask for selecting events (0 = no selection)
  ULong64_t fSelectTriggerMask50 = 0; // trigger maskNext50 for selecting events (0 = no selection)

  int fEventNumber = -1; // current event number

  RawDataHeader* fHeaderSwapped = nullptr; // temporary buffer for swapping header on PowerPC
  UChar_t* fData;                          // raw data

  int fDDLNumber;       // index of current DDL number
  int fPrevDDLNumber;   // index of previous DDL number
  int fRCUId;           // current RCU identifier
  int fPrevRCUId;       // previous RCU identifier
  short fHWAddress;     // current hardware address
  short fPrevHWAddress; // previous hardware address
  int fTime;            // index of current time bin
  int fPrevTime;        // index of previous time bin
  int fSignal;          // signal in ADC counts
  int fTimeBunch;       // total length of the current time bunch

  int fStreamPosition;     // current (10 bit) position in fData
  int fStreamCount;        // counter of words to be read for current trailer
  int fChannelPayloadSize; //
  int fBunchLength;        // remaining number of signal bins in the current bunch

  UChar_t* fRCUTrailerData; // pointer to RCU trailer data
  int fRCUTrailerSize;      // size of RCU trailer data in bytes

  // RCU trailer contents
  uint fFECERRA;       // contains errors related to ALTROBUS transactions
  uint fFECERRB;       // contains errors related to ALTROBUS transactions
  ushort fERRREG2;     // contains errors related to ALTROBUS transactions or trailer of ALTRO channel block
  ushort fERRREG3;     // contains number of altro channels skipped due to an address mismatch
  ushort fERRREG4;     // contains number of altro channels skipped due to a block length mismatch
  ushort fActiveFECsA; // bit pattern of active FECs in branch A
  ushort fActiveFECsB; // bit pattern of active FECs in branch B
  uint fAltroCFG1;     // ALTROCFG1 register
  uint fAltroCFG2;     // ALTROCFG2 and ALTROIF registers

  int fChannelStartPos; // start index of the current channel
  // int            fPosition;     // current position (32-bit words) in fData
  // int            fCount;        //
  int fStartTimeBin; //
  // int            fBunchLength;  //

  bool fBadChannel; //
  int fPayloadSize; //

  // int            fChannelPayloadSize; //

  UShort_t fBunchData[kMaxNTimeBins]; // cache for the decoded altro payload
  UShort_t* fBunchDataPointer;        // pointer to the current bunch samples
  int fBunchDataIndex;                // current position in the payload

  // UChar_t*         fRCUTrailerData; // pointer to RCU trailer data
  // int            fRCUTrailerSize; // size of RCU trailer data in bytes

  // // RCU trailer contents
  // Uint           fFECERRA;      // contains errors related to ALTROBUS transactions
  // Uint           fFECERRB;      // contains errors related to ALTROBUS transactions
  // UShort_t         fERRREG2;      // contains errors related to ALTROBUS transactions or trailer of ALTRO channel block
  // Uint           fERRREG3;      // contains number of altro channels skipped due to an address mismatch
  // UShort_t         fActiveFECsA;  // bit pattern of active FECs in branch A
  // UShort_t         fActiveFECsB;  // bit pattern of active FECs in branch B
  // Uint           fAltroCFG1;    // ALTROCFG1 register
  // Uint           fAltroCFG2;    // ALTROCFG2 and ALTROIF registers

  // AliAltroRawStream* fOldStream;  // streamer for old altro format

  bool fCheckAltroPayload; // check altro payload correctness or not?
  int fFormatVersion;

  int fModule = 0;   // index of current module
  int fXcell = 0;    // index of current row
  int fZcell = 0;    // index of current column
  int fCaloFlag = 0; // low (0) or (1) high gain; see enum EAliCaloFlag above
  Mapping fMapping;
  std::string fMappingPath = "";

  TH2F* fErrorSummary = nullptr;

  ClassDef(RawReader, 0) // class for reading raw digits from a root file
};

#endif
