#include "event.h"
#include "Mapping.h"
#include "RawReader.h"
#include "RawDataHeader.h"

#include <iostream>

ClassImp(RawReader)
  //_______________________________________________________
  RawReader::RawReader(std::basic_string_view<char> fileName) : fFile(nullptr),
                                                                fEvent(nullptr),
                                                                fSubEvent(nullptr),
                                                                fEquipment(nullptr),
                                                                fPosition(nullptr),
                                                                fEnd(nullptr),
                                                                fModule(-1),
                                                                fXcell(-1),
                                                                fZcell(-1),
                                                                fCaloFlag(0),
                                                                fStreamPosition(0)
{
  // create an object to read digits from the given date event

  fFile = fopen(fileName.data(), "rb");
  if (!fFile) {
    std::cout << "ERROR: RawReader: could not open file " << fileName << std::endl;
    return;
  }

  eventHeaderStruct header;
  unsigned int headerSize = sizeof(eventHeaderStruct);
  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    UChar_t* buffer = new UChar_t[header.eventSize];
    fseek(fFile, -(long)headerSize, SEEK_CUR);
    if (fread(buffer, 1, header.eventSize, fFile) != header.eventSize)
      break;
    fEvent = reinterpret_cast<eventHeaderStruct*>(buffer);
    break;
  }

  // create an object to read PHOS/EMCAL raw digits
  // SelectRawData(calo); <---------------
}
//_______________________________________________________
RawReader::~RawReader()
{
  // destructor

  if (fEvent)
    delete[] fEvent;
  if (fFile) {
    fclose(fFile);
  }
}
//_______________________________________________________
uint RawReader::GetEventType() const
{
  // get the type from the event header

  if (!fEvent)
    return 0;
  return fEvent->eventType;
}
//_______________________________________________________
uint RawReader::GetRunNumber() const
{
  // get the run number from the event header

  if (!fEvent)
    return 0;
  return fEvent->eventRunNb;
}
//_______________________________________________________
const uint* RawReader::GetEventId() const
{
  // get the event id from the event header

  if (!fEvent)
    return nullptr;
  return fEvent->eventId;
}
//_______________________________________________________
const uint* RawReader::GetTriggerPattern() const
{
  // get the trigger pattern from the event header

  if (!fEvent)
    return nullptr;
  return fEvent->eventTriggerPattern;
}
// //_______________________________________________________
// const uint* RawReader::GetDetectorPattern() const
// {
// // get the detector pattern from the event header

//   if (!fEvent) return nullptr;
//   return fEvent->eventDetectorPattern;
// }
//_______________________________________________________
// const uint* RawReader::GetAttributes() const
// {
// // get the type attributes from the event header

//   if (!fEvent) return nullptr;
//   return fEvent->eventTypeAttribute;
// }
// //_______________________________________________________
// const unsigned int* RawReader::GetSubEventAttributes() const
// {
// // get the type attributes from the sub event header

//   if (!fSubEvent) return nullptr;
//   return fSubEvent->eventTypeAttribute;
// }
// //_______________________________________________________
// unsigned int RawReader::GetLDCId() const
// {
// // get the LDC Id from the event header

//   if (!fSubEvent) return 0;
//   return fSubEvent->eventLdcId;
// }
// //_______________________________________________________
// unsigned int RawReader::GetGDCId() const
// {
// // get the GDC Id from the event header

//   if (!fEvent) return 0;
//   return fEvent->eventGdcId;
// }
//_______________________________________________________
uint RawReader::GetTimestamp() const
{
  // get the timestamp from the event header

  if (!fEvent)
    return 0;
  return fEvent->eventTimestamp;
}
//_______________________________________________________
int RawReader::GetEquipmentSize() const
{
  // get the size of the equipment (including the header)

  if (!fEquipment)
    return 0;
  if (fSubEvent->eventVersion <= 0x00030001) {
    return fEquipment->equipmentSize + sizeof(equipmentHeaderStruct);
  } else {
    return fEquipment->equipmentSize;
  }
}
//_______________________________________________________
int RawReader::GetEquipmentType() const
{
  // get the type from the equipment header

  if (!fEquipment)
    return -1;
  return fEquipment->equipmentType;
}
//_______________________________________________________
int RawReader::GetEquipmentId() const
{
  // get the ID from the equipment header

  if (!fEquipment)
    return -1;
  return fEquipment->equipmentId;
}
//_______________________________________________________
const unsigned int* RawReader::GetEquipmentAttributes() const
{
  // get the attributes from the equipment header

  if (!fEquipment)
    return nullptr;
  return fEquipment->equipmentTypeAttribute;
}
//_______________________________________________________
int RawReader::GetEquipmentElementSize() const
{
  // get the basic element size from the equipment header

  if (!fEquipment)
    return 0;
  return fEquipment->equipmentBasicElementSize;
}
//_______________________________________________________
int RawReader::GetEquipmentHeaderSize() const
{
  // Get the equipment header size
  // 28 bytes by default
  return sizeof(equipmentHeaderStruct);
}
//_______________________________________________________
bool RawReader::ReadHeader()
{
  // read a data header at the current position
  // returns false if the data header could not be read

  fErrorCode = 0;

  fHeader = nullptr;
  if (!fEvent)
    return false;
  // check whether there are sub events
  if (fEvent->eventSize <= fEvent->eventHeadSize)
    return false;

  do {
    // skip payload (if event was not selected)

    if (fCount > 0)
      fPosition += fCount;

    // get the first or the next equipment if at the end of an equipment
    if (!fEquipment || (fPosition >= fEnd)) {
      fEquipment = nullptr;

      // get the first or the next sub event if at the end of a sub event
      if (!fSubEvent || (fPosition >= ((UChar_t*)fSubEvent) + fSubEvent->eventSize)) {

        // check for end of event data
        if (fPosition >= ((UChar_t*)fEvent) + fEvent->eventSize) {
          return false;
        }
        if (!TEST_SYSTEM_ATTRIBUTE(fEvent->eventTypeAttribute, ATTR_SUPER_EVENT)) {
          fSubEvent = fEvent; // no super event
        } else if (fSubEvent) {
          fSubEvent = (eventHeaderStruct*)(((UChar_t*)fSubEvent) +
                                           fSubEvent->eventSize);
        } else {
          fSubEvent = (eventHeaderStruct*)(((UChar_t*)fEvent) +
                                           fEvent->eventHeadSize);
        }

        // check the magic word of the sub event
        if (fSubEvent->eventMagic != EVENT_MAGIC_NUMBER) {
          std::cout << "ERROR: ReadHeader: wrong magic number in sub event! run:" << fSubEvent->eventRunNb << "  event: " << fSubEvent->eventId[0] << " " << fSubEvent->eventId[1] << std::endl;
          fErrorCode = kErrMagic;
          return false;
        }

        // continue if no data in the subevent
        if (fSubEvent->eventSize == fSubEvent->eventHeadSize) {
          fPosition = fEnd = ((UChar_t*)fSubEvent) + fSubEvent->eventSize;
          fCount = 0;
          continue;
        }

        fEquipment = (equipmentHeaderStruct*)(((UChar_t*)fSubEvent) + fSubEvent->eventHeadSize);

      } else {
        fEquipment = (equipmentHeaderStruct*)fEnd;
      }

      fCount = 0;
      fPosition = ((UChar_t*)fEquipment) + sizeof(equipmentHeaderStruct);
      if (fSubEvent->eventVersion <= 0x00030001) {
        fEnd = fPosition + fEquipment->equipmentSize;
      } else {
        fEnd = ((UChar_t*)fEquipment) + fEquipment->equipmentSize;
      }
    }

    // continue with the next sub event if no data left in the payload
    if (fPosition >= fEnd)
      continue;

    // check that there are enough bytes left for the data header
    if (fPosition + sizeof(RawDataHeader) > fEnd) {
      std::cout << "ERROR: ReadHeader: could not read data header data!" << std::endl;
      fCount = 0;
      fPosition = fEnd;
      fErrorCode = kErrNoDataHeader;
      continue;
    }

    // "read" the data header
    fHeader = reinterpret_cast<RawDataHeader*>(fPosition);
    // Now check the version of the header
    int version = 0;
    if (fHeader) {
      version = static_cast<int>(fHeader->GetVersion());
    }
    if (version != 3) {
      std::cout << "ERROR: ReadHeader: wrong version=" << version << ", expect 3" << std::endl;
      fCount = 0;
      fPosition = fEnd;
      fErrorCode = kErrVersion;
      continue;
    }

    if ((fPosition + fHeader->fSize) != fEnd) {
      if ((fHeader->fSize != 0xFFFFFFFF) && (fEquipment->equipmentId != 4352)) {
        std::cout << "WARNING:  ReadHeader: raw data size found in the header is wrong (" << fHeader->fSize << "!=" << fEnd - fPosition << ")! Using the equipment size instead !" << std::endl;
      }
      fHeader->fSize = fEnd - fPosition;
    }
    fPosition += sizeof(RawDataHeader);

    if (fHeader && (fHeader->fSize != 0xFFFFFFFF)) {
      fCount = fHeader->fSize - sizeof(RawDataHeader);

      // check consistency of data size in the header and in the sub event
      if (fPosition + fCount > fEnd) {
        std::cout << "ERROR: ReadHeader: size in data header exceeds event size!" << std::endl;
        fCount = 0;
        fPosition = fEnd;
        fErrorCode = kErrSize;
        continue;
      }

    } else {
      fCount = fEnd - fPosition;
    }

  } while (!fEquipment);

  return true;
}
//_______________________________________________________
bool RawReader::ReadNextData()
{
  // reads the next payload at the current position
  // returns false if the data could not be read

  fErrorCode = 0;
  while (fCount == 0) {
    if (!ReadHeader())
      return false;
  }
  fData = fPosition;
  fPosition += fCount;
  fCount = 0;
  return true;
}
//_______________________________________________________
bool RawReader::ReadNext(UChar_t* data, int size)
{
  // reads the next block of data at the current position
  // returns false if the data could not be read

  fErrorCode = 0;
  if (fPosition + size > fEnd) {
    std::cout << "ERROR: ReadNext: could not read data!" << std::endl;
    fErrorCode = kErrOutOfBounds;
    return false;
  }
  memcpy(data, fPosition, size);
  fPosition += size;
  fCount -= size;
  return true;
}
//_______________________________________________________
bool RawReader::Reset()
{
  // reset the current position to the beginning of the event

  fSubEvent = nullptr;
  fEquipment = nullptr;
  fCount = 0;
  fPosition = nullptr;
  fEnd = nullptr;
  fHeader = nullptr;

  // reset PHOS raw stream params
  // Complete reset of raw stream params

  fHWAddress = -1;
  fChannelStartPos = -1;
  fStreamCount = -1;
  fBunchLength = fStartTimeBin = -1;
  fBadChannel = false;
  fPayloadSize = -1;
  fChannelPayloadSize = -1;
  fBunchDataPointer = NULL;
  fBunchDataIndex = -1;

  fRCUTrailerData = NULL;
  fRCUTrailerSize = 0;

  fFECERRA = fFECERRB = fERRREG2 = fERRREG3 = fActiveFECsA = fActiveFECsB = fAltroCFG1 = fAltroCFG2 = 0;
  fModule = fXcell = fZcell = 0;
  fCaloFlag = 0;

  return true;
}

bool RawReader::NextEvent()
{
  // go to the next event in the date file

  if (!fFile) {
    return false;
  }

  Reset();
  eventHeaderStruct header;
  unsigned int headerSize = sizeof(eventHeaderStruct);
  if (fEvent)
    delete[] fEvent;
  fEvent = &header;

  while (fread(&header, 1, headerSize, fFile) == headerSize) {
    if (!IsEventSelected()) {
      fseek(fFile, header.eventSize - headerSize, SEEK_CUR);
      continue;
    }
    UChar_t* buffer = new UChar_t[header.eventSize];
    fseek(fFile, -(long)headerSize, SEEK_CUR);
    if (fread(buffer, 1, header.eventSize, fFile) != header.eventSize) {
      std::cout << "ERROR: NextEvent: could not read event from file" << std::endl;
      delete[] buffer;
      break;
    }
    fEvent = (eventHeaderStruct*)buffer;
    fEventNumber++;

    // Read the only DDL (from AliAltroRawStreamV3::NextDDL())
    //  Read the next DDL payload (CDH + RCU trailer)
    //  Updates the information which is coming from these
    //  two sources
    fFormatVersion = 0;
    fStreamPosition = 0;
    // Get next DDL payload
    // return wtih false in case no more data payloads
    // are found
    do {
      if (!ReadNextData())
        return false;
    } while (GetDataSize() == 0);

    // fDDLNumber = GetDDLID();
    fChannelPayloadSize = -1;
    fChannelStartPos = -1;

    fFormatVersion = static_cast<int>(fHeader->GetAttributes() & 0xF); // altro format version

    if (fFormatVersion < 2) {
      // old altro format data
      std::cout << "ERROR: NextEvent: wrong format version " << fFormatVersion << " expect 3" << std::endl;
      return false;
    }
    return ReadRCUTrailer(fFormatVersion);
  }

  fEvent = nullptr;

  return false;
}

int RawReader::CheckData() const
{
  // check the consistency of the data

  if (!fEvent)
    return 0;
  // check whether there are sub events
  if (fEvent->eventSize <= fEvent->eventHeadSize)
    return 0;

  eventHeaderStruct* subEvent = nullptr;
  UChar_t* position = 0;
  UChar_t* end = 0;
  int result = 0;

  while (true) {
    // get the first or the next sub event if at the end of a sub event
    if (!subEvent || (position >= end)) {

      // check for end of event data
      if (position >= ((UChar_t*)fEvent) + fEvent->eventSize) {
        return result;
      }
      if (!TEST_SYSTEM_ATTRIBUTE(fEvent->eventTypeAttribute,
                                 ATTR_SUPER_EVENT)) {
        subEvent = fEvent; // no super event
      } else if (subEvent) {
        subEvent = reinterpret_cast<eventHeaderStruct*>(((UChar_t*)subEvent) + subEvent->eventSize);
      } else {
        subEvent = reinterpret_cast<eventHeaderStruct*>(((UChar_t*)fEvent) + fEvent->eventHeadSize);
      }

      // check the magic word of the sub event
      if (subEvent->eventMagic != EVENT_MAGIC_NUMBER) {
        result |= kErrMagic;
        return result;
      }

      position = ((UChar_t*)subEvent) + subEvent->eventHeadSize + sizeof(equipmentHeaderStruct);
      end = ((UChar_t*)subEvent) + subEvent->eventSize;
    }

    // continue with the next sub event if no data left in the payload
    if (position >= end)
      continue;

    // check that there are enough bytes left for the data header
    if (position + sizeof(RawDataHeader) > end) {
      result |= kErrNoDataHeader;
      position = end;
      continue;
    }

    // Here we have to check if we have header v2 or v3
    // check consistency of data size in the data header and in the sub event
    RawDataHeader* header = reinterpret_cast<RawDataHeader*>(position);
    UChar_t version = header->GetVersion();
    if (version == 3) {
      if ((position + header->fSize) != end) {
        if (header->fSize != 0xFFFFFFFF) {
          std::cout << "WARNING: CheckData raw data size found in the header V3 is wrong (" << header->fSize << "!=" << end - position << ")! Using the equipment size instead !" << std::endl;
        }
        header->fSize = end - position;
        result |= kErrSize;
      }
    }
    position = end;
  }
  return 0;
}

// Frow CaloStreamV3

///////////////////////////////////////////////////////////////////////////////
//
// This class provides access to PHOS/EMCAL digits in raw data.
//
// It loops over all PHOS/EMCAL digits in the raw data given by the AliRawReader.
// The Next method goes to the next digit. If there are no digits left
// it returns false.
// Several getters provide information about the current digit.
// usage:
//    AliRawReader *reader = AliRawReader::Create(fileName);
//    AliCaloRawStreamV3 *stream = new AliCaloRawStreamV3(reader,calo);
//    while (reader->NextEvent())
//      while (stream->NextDDL())
//        while (stream->NextChannel()) ...
///
/// Yuri Kharlov. 23 June 2009
///////////////////////////////////////////////////////////////////////////////

// //_____________________________________________________________________________
// RawReader::~RawReader()
// {
// // destructor

//   if (!fExternalMapping)
//     for(int i = 0; i < fNModules*fNRCU*fNSides; i++)
//       delete fMapping[i];
// }

//_____________________________________________________________________________
bool RawReader::NextChannel()
{

  // Read the next Altro channel from the
  // raw-data stream
  // Updates the channel hardware address member and
  // channel data size. Sets the error flag in case
  // RCU signals readout error in this channel
  // Read next PHOS signal
  // Apply the PHOS altro mapping to get
  // the module,row and column indeces

  // if (fOldStream) {
  //   bool status = fOldStream->NextChannel();
  //   if (status) {
  //     fHWAddress = fOldStream->GetHWAddress();
  //     fChannelPayloadSize = fOldStream->GetChannelPayloadSize();
  //   }
  //   ApplyMapping();
  //   return status;
  // }

  int channelStartPos = fStreamPosition;
  fChannelStartPos = -1;
  fStreamCount = -1;
  fBadChannel = false;
  fBunchDataIndex = 0;
  fBunchLength = -1;

  uint word = 0;
  do {
    word = Get32bitWord(fStreamPosition++);
    if (fStreamPosition > fPayloadSize)
      return false;
  } while ((word >> 30) != 1);

  // check for readout errors
  fBadChannel = (word >> 29) & 0x1;

  // extract channel payload and hw address
  fStreamCount = (word >> 16) & 0x3FF;
  fChannelPayloadSize = fStreamCount;
  fHWAddress = word & 0xFFF;

  // Now unpack the altro data
  // Revert the order of the samples
  // inside the bunch so that the
  // first time is first in the samples
  // array
  int isample = 0;
  int nwords = (fStreamCount + 2) / 3;
  for (int iword = 0; iword < nwords; iword++) {
    word = Get32bitWord(fStreamPosition++);
    if ((word >> 30) != 0) {
      // Unexpected end of altro channel payload
      std::cout << "WARNING: NextChannel: Unexpected end of payload in altro channel payload! DDL=" << fDDLNumber << ", Address=" << fHWAddress << ", word=" << word << std::endl;
      fStreamCount = -1;
      fStreamPosition--;
      return false;
    }
    fBunchData[isample++] = (word >> 20) & 0x3FF;
    fBunchData[isample++] = (word >> 10) & 0x3FF;
    fBunchData[isample++] = word & 0x3FF;
  }

  fChannelStartPos = channelStartPos;
  ApplyMapping();
  return true;
}

//_____________________________________________________________________________
void RawReader::ApplyMapping()
{

  fXcell = fMapping.GetXcell(fHWAddress);
  fZcell = fMapping.GetZcell(fHWAddress);
  fCaloFlag = fMapping.GetCaloFlag(fHWAddress);
}

//_____________________________________________________________________________
bool RawReader::NextDDL()
{
  // // Read the next DDL payload (CDH + RCU trailer)
  // // Updates the information which is coming from these
  // // two sources
  //   fFormatVersion = 0;
  //   fStreamPosition = 0;
  //   // Get next DDL payload
  //   // return wtih false in case no more data payloads
  //   // are found
  //   do {
  //     if (!ReadNextData(fData)) return false;
  //   } while (GetDataSize() == 0);

  //   fDDLNumber = GetDDLID();
  //   fChannelPayloadSize = -1;
  //   fChannelStartPos = -1;

  //   fFormatVersion = (GetBlockAttributes() & 0xF);

  //   return ReadRCUTrailer(fFormatVersion);
  return true;
}

//_____________________________________________________________________________
bool RawReader::NextBunch()
{
  // Read next altro bunch from the
  // raw-data stream.
  // Updates the start/end time-bins
  // and the array with altro samples
  // if (fOldStream) {
  //   bool status = fOldStream->NextBunch(fBunchData,fBunchLength,fStartTimeBin);
  //   if (status) fBunchDataPointer = &fBunchData[0];
  //   else fBunchDataPointer = NULL;
  //   return status;
  // }

  int prevTimeBin = (fBunchLength > 0) ? fStartTimeBin - fBunchLength + 1 : 1024;
  fBunchLength = fStartTimeBin = -1;
  fBunchDataPointer = nullptr;

  if ((fBunchDataIndex >= fStreamCount) || fBadChannel) {
    return false;
  }

  fBunchLength = fBunchData[fBunchDataIndex];
  if (fBunchLength <= 2) {
    // Invalid bunch size
    // AliDebug(1,Form("Too short bunch length (%d) @ %d in Address=0x%x (DDL=%03d)!",
    //    fBunchLength,fBunchDataIndex,fHWAddress,fDDLNumber));
    // AddMinorErrorLog(kAltroBunchHeadErr,Form("hw=0x%x",fHWAddress));
    // if (AliDebugLevel() > 0) HexDumpChannel();
    fStreamCount = fBunchLength = -1;
    return false;
  }
  if ((fBunchDataIndex + fBunchLength) > fStreamCount) {
    // Too long bunch detected
    // AliDebug(1,Form("Too long bunch detected @ %d in Address=0x%x (DDL=%03d) ! Expected <= %d 10-bit words, found %d !", fBunchDataIndex,
    //     fHWAddress,fDDLNumber,fStreamCount-fBunchDataIndex,fBunchLength));
    // AddMinorErrorLog(kAltroBunchHeadErr,Form("hw=0x%x",fHWAddress));
    fStreamCount = fBunchLength = -1;
    return false;
  }
  fBunchDataIndex++;
  fBunchLength -= 2;

  fStartTimeBin = fBunchData[fBunchDataIndex++];
  if (fCheckAltroPayload) {
    if ((fStartTimeBin - fBunchLength + 1) < 0) {
      // AliWarning(Form("Invalid start time-bin @ %d in Address=0x%x (DDL=%03d)! (%d-%d+1) < 0", fBunchDataIndex-1,
      //       fHWAddress,fDDLNumber,fStartTimeBin,fBunchLength));

      // AddMinorErrorLog(kAltroPayloadErr,Form("hw=0x%x",fHWAddress));
      // if (AliDebugLevel() > 0) HexDumpChannel();
      fStreamCount = fBunchLength = -1;
      return false;
    }
    if (fStartTimeBin >= prevTimeBin) {
      // AliWarning(Form("Invalid start time-bin @ %d in Address=0x%x (DDL=%03d)! (%d>=%d)", fBunchDataIndex-1,
      //       fHWAddress,fDDLNumber,fStartTimeBin,prevTimeBin));
      // AddMinorErrorLog(kAltroPayloadErr,Form("hw=0x%x",fHWAddress));
      fStreamCount = fBunchLength = -1;
      return false;
    }
  }

  fBunchDataPointer = &fBunchData[fBunchDataIndex];

  fBunchDataIndex += fBunchLength;

  return true;
}

//_____________________________________________________________________________
const UChar_t* RawReader::GetChannelPayload() const
{
  // returns raw channel data, length 4+(fChannelPayloadSize+2)/3*4
  if (fChannelStartPos < 0 || fChannelPayloadSize < 0)
    return nullptr;
  int channelSize = 1 + (fChannelPayloadSize + 2) / 3; // nof 32bit words
  if (fStreamPosition < fChannelStartPos + channelSize)
    return nullptr;
  return fData + (fChannelStartPos << 2);
}

//_____________________________________________________________________________
uint RawReader::Get32bitWord(int index) const
{
  // This method returns the 32 bit word at a given
  // position inside the raw data payload.
  // The 'index' points to the beginning of the word.
  // The method is supposed to be endian (platform)
  // independent.

  index = (index << 2);
  uint word = 0;
  word |= fData[index++];
  word |= fData[index++] << 8;
  word |= fData[index++] << 16;
  word |= fData[index++] << 24;

  return word;
}

///_____________________________________________________________________________
bool RawReader::ReadRCUTrailer(int rcuVer)
{
  // Read the RCU trailer according
  // to the RCU formware version
  // specified in CDH
  // Cross-check with version found in the
  // trailer

  fRCUTrailerData = nullptr;
  fRCUTrailerSize = 0;
  fPayloadSize = -1;

  int index = GetDataSize() / 4;

  // First read 32-bit word with the
  // trailer size (7 bits), RCU ID (9 bits) and
  // RCU firmware version (8 bits?)
  // The two major bit should be 11 (identifies
  // the end of the trailer)
  uint word = Get32bitWord(--index);

  if ((word >> 30) != 3) {
    // fRawReader->AddFatalErrorLog(kRCUTrailerErr,"");
    std::cout << "ERROR: ReadRCUTrailer Last RCU trailer word not found!" << std::endl;
    return false;
  }

  UChar_t ver = (word >> 16) & 0xFF;
  if (ver != rcuVer) {
    // Wrong RCU firmware version detected
    // fRawReader->AddMajorErrorLog(kRCUVerErr,Form("%d!=%d",ver,rcuVer));
  }

  fRCUId = (int)((word >> 7) & 0x1FF);
  int trailerSize = (word & 0x7F);

  // Now read the beginning of the trailer
  // where the payload size is written
  if (trailerSize < 2) {
    // fRawReader->AddMajorErrorLog(kRCUTrailerErr,Form("tr=%d bytes",trailerSize*4));
    // AliWarning(Form("Invalid trailer size found (%d bytes) !",
    //     trailerSize*4));
    return false;
  }

  trailerSize -= 2;
  fRCUTrailerSize = trailerSize * 4;

  for (; trailerSize > 0; trailerSize--) {
    word = Get32bitWord(--index);
    if ((word >> 30) != 2) {
      // fRawReader->AddMinorErrorLog(kRCUTrailerErr,"missing 10");
      continue;
    }
    int parCode = (word >> 26) & 0xF;
    int parData = word & 0x3FFFFFF;
    switch (parCode) {
      case 1:
        // ERR_REG1
        fFECERRA = ((parData >> 13) & 0x1FFF) << 7;
        fFECERRB = ((parData & 0x1FFF)) << 7;
        break;
      case 2:
        // ERR_REG2
        fERRREG2 = parData & 0x1FF;
        break;
      case 3:
        // ERR_REG3
        fERRREG3 = parData & 0x1FFFFFF;
        break;
      case 4:
        // FEC_RO_A
        fActiveFECsA = parData & 0xFFFF;
        break;
      case 5:
        // FEC_RO_B
        fActiveFECsB = parData & 0xFFFF;
        break;
      case 6:
        // RDO_CFG1
        fAltroCFG1 = parData & 0xFFFFF;
        break;
      case 7:
        // RDO_CFG2
        fAltroCFG2 = parData & 0x1FFFFFF;
        break;
      default:
        // fRawReader->AddMinorErrorLog(kRCUTrailerErr,"undef word");
        break;
    }
  }

  if (index < 1) {
    // fRawReader->AddMajorErrorLog(kRCUTrailerErr,Form("tr=%d raw=%d bytes",
    //              fRCUTrailerSize+8,
    //              fRawReader->GetDataSize()));
    // AliWarning(Form("Invalid trailer size found (%d bytes) ! The size is bigger than the raw data size (%d bytes)!",
    //     fRCUTrailerSize,
    //     fRawReader->GetDataSize()));
  }

  fRCUTrailerData = fData + index * 4;

  // Now read the payload size
  // (First word in the RCU trailer)
  fPayloadSize = Get32bitWord(--index) & 0x3FFFFFF;

  if ((GetDataSize() - fRCUTrailerSize - 8) != (fPayloadSize * 4)) {
    // AddMajorErrorLog(kRCUTrailerSizeErr,Form("h=%d tr=%d rcu=%d bytes",
    //            fRawReader->GetDataSize(),
    //            fRCUTrailerSize+8,
    //            fPayloadSize*4));
    // AliWarning(Form("Inconsistent raw data size ! Raw data size - %d bytes (from CDH), RCU trailer - %d bytes, raw data size (from RCU trailer) - %d bytes !",
    //     fRawReader->GetDataSize(),
    //     fRCUTrailerSize+8,
    //     fPayloadSize*4));
  }
  return true;
}

// //_____________________________________________________________________________
// int RawReader::GetBranch() const
// {
//   // The method provides the RCU branch index (0 or 1)
//   // for the current hardware address.
//   // In case the hardware address has not been yet
//   // initialized, the method returns -1
//   if (fHWAddress == -1) return -1;

//   return ((fHWAddress >> 11) & 0x1);
// }

// //_____________________________________________________________________________
// int RawReader::GetFEC() const
// {
//   // The method provides the front-end card index
//   // for the current hardware address.
//   // In case the hardware address has not been yet
//   // initialized, the method returns -1
//   if (fHWAddress == -1) return -1;

//   return ((fHWAddress >> 7) & 0xF);
// }

// //_____________________________________________________________________________
// int RawReader::GetAltro() const
// {
//   // The method provides the altro chip index
//   // for the current hardware address.
//   // In case the hardware address has not been yet
//   // initialized, the method returns -1
//   if (fHWAddress == -1) return -1;

//   return ((fHWAddress >> 4) & 0x7);
// }

// //_____________________________________________________________________________
// int RawReader::GetChannel() const
// {
//   // The method provides the channel index
//   // for the current hardware address.
//   // In case the hardware address has not been yet
//   // initialized, the method returns -1
//   if (fHWAddress == -1) return -1;

//   return (fHWAddress & 0xF);
// }

// //_____________________________________________________________________________
// bool RawReader::GetRCUTrailerData(UChar_t*& data) const
// {
//   // Return a pointer to the RCU trailer
//   // data. Should be called always after
//   // the RCU trailer was already processed
//   // in the GetPosition() method
//   if (!fRCUTrailerSize || !fRCUTrailerData) {
//     std::cout << "ERROR: GetRCUTrailerData: No valid RCU trailer data is found !"<< std::endl;
//     data = nullptr;
//     return false;
//   }

//   data = fRCUTrailerData;

//   return true;
// }

// //_____________________________________________________________________________
// void RawReader::PrintRCUTrailer() const
// {
//   // Prints the contents of
//   // the RCU trailer data
//   printf("RCU trailer (Format version %d):\n"
//    "==================================================\n",  GetFormatVersion());
//   printf("FECERRA:                                   0x%x\n", fFECERRA);
//   printf("FECERRB:                                   0x%x\n", fFECERRB);
//   printf("ERRREG2:                                   0x%x\n", fERRREG2);
//   printf("#channels skipped due to address mismatch: %d\n",GetNChAddrMismatch());
//   printf("#channels skipped due to bad block length: %d\n",GetNChLengthMismatch());
//   printf("Active FECs (branch A):                    0x%x\n", fActiveFECsA);
//   printf("Active FECs (branch B):                    0x%x\n", fActiveFECsB);
//   printf("Baseline corr:                             0x%x\n",GetBaselineCorr());
//   printf("Number of presamples:                      %d\n", GetNPresamples());
//   printf("Number of postsamples:                     %d\n",GetNPostsamples());
//   printf("Second baseline corr:                      %d\n",GetSecondBaselineCorr());
//   printf("GlitchFilter:                              %d\n",GetGlitchFilter());
//   printf("Number of non-ZS postsamples:              %d\n",GetNNonZSPostsamples());
//   printf("Number of non-ZS presamples:               %d\n",GetNNonZSPresamples());
//   printf("Number of ALTRO buffers:                   %d\n",GetNAltroBuffers());
//   printf("Number of pretrigger samples:              %d\n",GetNPretriggerSamples());
//   printf("Number of samples per channel:             %d\n",GetNSamplesPerCh());
//   printf("Sparse readout:                            %d\n",GetSparseRO());
//   printf("Sampling time:                             %e s\n",GetTSample());
//   printf("L1 Phase:                                  %e s\n",GetL1Phase());
//   printf("AltroCFG1:                                 0x%x\n",GetAltroCFG1());
//   printf("AltroCFG2:                                 0x%x\n",GetAltroCFG2());
//   printf("==================================================\n");
// }

// //_____________________________________________________________________________
// double RawReader::GetTSample() const
// {
//   // Returns the sampling time
//   // in seconds. In case the rcu trailer
//   // was note read, return an invalid number (0)

//   if (!fRCUTrailerData) return 0.;

//   const double kLHCTimeSample = 25.0e-9; // LHC clocks runs at 40 MHz
//   UChar_t fq = (fAltroCFG2 >> 5) & 0xF;
//   double tSample;
//   switch (fq) {
//   case 0:
//     // 20 MHz
//     tSample = 2.0*kLHCTimeSample;
//     break;
//   case 1:
//     // 10 Mhz
//     tSample = 4.0*kLHCTimeSample;
//     break;
//   case 2:
//     // 5 MHz
//     tSample = 8.0*kLHCTimeSample;
//     break;
//   default:
//     AliWarning(Form("Invalid sampling frequency value %d !",
//           fq));
//     tSample = 0.;
//     break;
//   }

//   return tSample;
// }

// //_____________________________________________________________________________
// double RawReader::GetL1Phase() const
// {
//   // Returns the L1 phase w.r.t to the
//   // LHC clock
//   if (!fRCUTrailerData) return 0.;

//   const double kLHCTimeSample = 25.0e-9; // LHC clocks runs at 40 MHz
//   double phase = ((double)(fAltroCFG2 & 0x1F))*kLHCTimeSample;

//   double tSample = GetTSample();
//   if (phase >= tSample) {
//     AliWarning(Form("Invalid L1 trigger phase (%f >= %f) !",
//         phase,tSample));
//     phase = 0.;
//   }

//   return phase;
// }

// //_____________________________________________________________________________
// UChar_t *RawReader::GetRCUPayloadInSOD() const
// {
//   // Get a pointer to the data in case
//   // of SOD events
//   if (fRawReader) {
//     if (fRawReader->GetType() == AliRawEventHeaderBase::kStartOfData) {
//       return fData;
//     }
//   }
//   return NULL;
// }

// //_____________________________________________________________________________
// int RawReader::GetRCUPayloadSizeInSOD() const
// {
//   // Get the size of the RCU data in case
//   // of SOD events
//   if (fRawReader) {
//     if (fRawReader->GetType() == AliRawEventHeaderBase::kStartOfData) {
//       return fPayloadSize;
//     }
//   }
//   return -1;
// }
//_____________________________________________________________________________
uint RawReader::GetPeriod() const
{
  const UInt_t* id = GetEventId();
  return id ? (((id)[0] >> 4) & 0x0fffffff) : 0;
}
//_____________________________________________________________________________
unsigned int RawReader::GetOrbitID() const
{
  const unsigned int* id = GetEventId();
  return id ? ((((id)[0] << 20) & 0xf00000) | (((id)[1] >> 12) & 0xfffff)) : 0;
}
//_____________________________________________________________________________
ushort RawReader::GetBCID() const
{
  const unsigned int* id = GetEventId();
  return id ? ((id)[1] & 0x00000fff) : 0;
}
//_____________________________________________________________________________
ULong64_t RawReader::GetEventIdAsLong() const
{
  return (((ULong64_t)GetPeriod() << 36) |
          ((ULong64_t)GetOrbitID() << 12) |
          (ULong64_t)GetBCID());
}
//_____________________________________________________________________________
ULong64_t RawReader::GetClassMask() const
{
  const unsigned int* pattern = GetTriggerPattern();
  return pattern ? (((ULong64_t)pattern[1] & 0x3ffff) << 32) | (pattern[0]) : 0;
}
//_____________________________________________________________________________
int RawReader::GetDataSize() const
{
  if (fHeader) {
    if (fHeader->fSize != 0xFFFFFFFF) {
      return fHeader->fSize - sizeof(RawDataHeader);
    } else {
      return GetEquipmentSize() - GetEquipmentHeaderSize() - sizeof(RawDataHeader);
    }
  } else {
    return GetEquipmentSize() - GetEquipmentHeaderSize();
  }
}
//_____________________________________________________________________________
bool RawReader::IsEventSelected() const
{
  // apply the event selection (if any)

  // First check the event type
  if (fSelectEventType >= 0) {
    if (GetEventType() != (UInt_t)fSelectEventType)
      return false;
  }

  // // Check the list of detectors
  // UInt_t clmask = GetDetectorPattern()[0];
  // if (fSelectDetectorExpr && (fSelectDetectorExpr&clmask)!=fSelectDetectorExpr) return false;
  // if (fExcludeDetectorExpr && (fExcludeDetectorExpr&clmask)) return false;

  // // Then check the trigger pattern and compared it
  // // to the required trigger mask
  // if (fSelectTriggerMask!=0 || fSelectTriggerMask50!=0) {
  //   if ( !(GetClassMask()&fSelectTriggerMask) && !(GetClassMaskNext50() & fSelectTriggerMask50)) return false;
  // }

  // if (  fIsTriggerClassLoaded && (!fSelectTriggerExpr.IsNull() || !fExcludeTriggerExpr.IsNull()  || !fExcludeClusterExpr.IsNull())) {
  //   TString expr(fSelectTriggerExpr);
  //   ULong64_t mask   = GetClassMask();
  //   ULong64_t maskNext50 = GetClassMaskNext50();
  //   //    if (mask) {
  //   for(Int_t itrigger = 0; itrigger < 50; itrigger++){
  //     if (mask & (1ull << itrigger)){
  // if(fVeto[itrigger]) return false;
  // expr.ReplaceAll(Form("[%d]",itrigger),"1");
  //     }else{
  // expr.ReplaceAll(Form("[%d]",itrigger),"0");
  //     }
  //   }
  //     //    if (maskNext50) {
  //   for(Int_t itrigger = 0; itrigger < 50; itrigger++){
  //     if (maskNext50 & (1ull << itrigger)){
  // if(fVeto[itrigger+50]) return false;
  // expr.ReplaceAll(Form("[%d]",itrigger+50),"1");
  //     }else{
  // expr.ReplaceAll(Form("[%d]",itrigger+50),"0");
  //     }
  //   }

  //   //
  //   // Possibility to introduce downscaling
  //   if(!fSelectTriggerExpr.IsNull()){
  //     // the gROOT->ProcessLineFast just evaluates C line, effectively a number here. If it gets > 32 0 in the end
  //     // or spaces, the evaluation will be wrong!
  //     expr.ReplaceAll(" ","");
  //     TPRegexp("0+").Substitute(expr,"0","g");
  //     TPRegexp("(%\\s*\\d+)").Substitute(expr,Form("&& !(%d$1)",GetEventIndex()),"g"); // RS: what is this line ?
  //     Int_t error;
  //     Bool_t result = gROOT->ProcessLineFast(expr.Data(),&error);
  //     if ( error == TInterpreter::kNoError)
  // return result;
  //     else
  // return false;
  //   }
  // }

  return true;
}
