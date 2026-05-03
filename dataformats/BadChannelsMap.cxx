#include <TH2.h>
#include <iostream>
#include "Geometry.h"
#include "BadChannelsMap.h"

BadChannelsMap::BadChannelsMap(int /*dummy*/)
{

  // Mark few channels as bad for test peurposes
  for (int i = 0; i < 16; i++) {
    int channelID = 1 + i * 17;
    mBadCells.set(channelID - OFFSET);
  }
}

void BadChannelsMap::GetHistogramRepresentation(TH2* h) const
{
  const int MAXX = 16,
            MAXZ = 16;
  if (!h) {
    std::cout << "Error:provide histogram to be filled" << std::endl;
  }
  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    std::cout << "Error:Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ << std::endl;
    return;
  }

  h->Reset();
  int relid[2] = {1, 1};
  int absId;
  int xmin = 1;
  for (int ix = 1; ix <= MAXX; ix++) {
    relid[0] = ix;
    for (int iz = 1; iz <= MAXZ; iz++) {
      relid[1] = iz;
      Geometry::RelToAbsId(relid, absId);
      if (!IsChannelGood(absId)) {
        h->SetBinContent(ix, iz, 1);
      }
    }
  }
}

void BadChannelsMap::PrintStream(std::ostream& stream) const
{
  // first sort bad channel IDs
  stream << "Number of bad cells:  " << mBadCells.count() << "\n";
  for (std::size_t cellID = 0; cellID < mBadCells.size(); cellID++) {
    if (mBadCells.test(cellID)) {
      stream << cellID + OFFSET << "\n";
    }
  }
}

std::ostream& operator<<(std::ostream& stream, const BadChannelsMap& bcm)
{
  bcm.PrintStream(stream);
  return stream;
}
