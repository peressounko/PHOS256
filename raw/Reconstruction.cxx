#include <algorithm>
#include "Digit.h"
#include "Reconstruction.h"

//_______________________________________
void Reconstruction::Run(std::string option)
{

  if (option.find("c") != std::string::npos || option.find("C") != std::string::npos) {
    fMakeClu = true;
  }

  Init();
  // Prepare RawReader
  fRawReader->Reset();

  int iev = 0, ievGood = 0;

  // Loop over events
  while (fRawReader->NextEvent()) {
    ++iev;
    // hEvCounter->Fill(fRawReader->GetEventType());
    if (fRawReader->GetEventType() != fEventSelection) {
      continue;
    }
    ++ievGood;
    Reset();

    MakeDigits();

    if (fMakeClu) {
      fClusterizer->ProcessEvent();
    }

    fTree->Fill();
    if (iev % 1000 == 0) {
      printf("Processed %d events and %d passed event selection\n", iev, ievGood);
    }
  }

  printf("Processed %d events and %d passed event selection\n", iev, ievGood);

  // WriteOutput
  fOutFile->cd();
  fTree->Write();
  fOutFile->Close();
}
//_______________________________________
bool Reconstruction::MakeDigits()
{

  fDigits->Clear();
  // Construct temporary list of HG and LG digits
  fCells.clear();

  // Get ZS parameters
  int nPreSamples = fRawReader->GetNPretriggerSamples();
  uint value = fRawReader->GetAltroCFG1();
  bool zeroSuppressionEnabled = (value >> 15) & 0x1;
  fFitter->SetPedSubtract(!zeroSuppressionEnabled);
  // fFitter.SetNPresamples(fRawReader->GetNPresamples());
  fFitter->SetNPresamples(nPreSamples);
  int offset = 0;
  int threshold = 0;
  if (zeroSuppressionEnabled) {
    offset = (value >> 10) & 0xf;
    threshold = value & 0x3ff;
  }
  while (fRawReader->NextChannel()) {

    // Get geometry of the channel
    int cellX = fRawReader->GetCellX(); // 1...16
    int cellZ = fRawReader->GetCellZ(); // 1...16
    int caloFlag = fRawReader->GetCaloFlag();
    if (cellX < 1 || cellX > 16 || cellZ < 1 || cellZ > 16 || (caloFlag != 0 && caloFlag != 1)) {
      continue;
    }

    int nBunches = 0;
    float maxAmp = 0;
    float t = 0; // time
    double meanAmp = 0;
    double rmsAmp = 0;
    int nsam = 0;
    int nMax = 0;
    short status = Fitter::kNotEvaluated;
    while (fRawReader->NextBunch()) {
      nBunches++;
      const unsigned short* sig = fRawReader->GetSignals();
      int sigStart = fRawReader->GetStartTimeBin();
      int sigLength = fRawReader->GetBunchLength();

      status = fFitter->Evaluate(sig, sigLength);

      for (int i = 0; i < sigLength; i++) {
        if (sig[i] > maxAmp) {
          maxAmp = sig[i];
          nMax = 0;
        }
      }
      maxAmp = fFitter->GetAmp();
      // maxAmp -= offset;
      t = fFitter->GetTime();
    }
    int absId = 0;
    int relid[2] = {cellX, cellZ};
    Geometry::RelToAbsId(relid, absId);
    fCells.emplace_back(absId, (caloFlag == 1), maxAmp, t, status);
  } // end of NextChannel()

  // Combine HG and LG, purge
  // Calibrate and remove (put amp to zero) bad channels
  std::sort(fCells.begin(), fCells.end()); // sort according to absId, then first HG then LG
  int iDigit = 0;
  auto cell = fCells.begin();
  while (cell != fCells.end()) {
    auto nextCell = cell;
    ++nextCell;

    float amp = cell->amp;
    float time = cell->time;
    bool isHG = cell->isHG;

    if (nextCell != fCells.end() && cell->absId == nextCell->absId) { // combine HG and LG
      if (cell->status != Fitter::kOK) {                              // use LG
        amp = nextCell->amp * 16;                                     // TODO: true calibration!!!!!!!!!!!!!!!!!!!!
        time = nextCell->time;
        isHG = nextCell->isHG;
      }
      cell = nextCell; // pass both cells
    }
    float energy = amp; // Calibrate(cell->absId,isHG, amp);
    // time   = CalibrateT(cell->absId,isHG,time);

    if (fBadMap->IsChannelGood(cell->absId)) {
      new ((*fDigits)[iDigit++]) Digit(cell->absId, energy, time, -1);
    }
    ++cell;
  }

  return true;
}
//_______________________________________
void Reconstruction::Init()
{

  if (fRawReader)
    return;

  fRawReader = new RawReader(fRawfilename);

  fOutFile = TFile::Open("PHOSReco.root", "recreate");
  fTree = new TTree("PHOS256", "Reconstruction tree");

  fDigits = new TClonesArray("Digit", 256);
  fTree->Branch("Digits", "TClonesArray", fDigits, 32000, 99);

  if (fMakeClu) {
    fClusterizer = new Clusterizer();
    fClusterizer->Init();
    fClusters = new TObjArray();
    fClusterizer->SetDigits(fDigits);
    fClusterizer->SetClusters(fClusters);
    fTree->Branch("Clusters", "TObjArray", &fClusters, 32000, 2);
  }

  if (!fBadMap) {
    fBadMap = new BadChannelsMap(); // test default
  }

  fCells.reserve(2 * NCELLS);
  fDigits->Expand(NCELLS);

  fFitter = new Fitter();
}

//_______________________________________
void Reconstruction::Reset()
{

  fCells.clear();
  fCells.reserve(2 * NCELLS);
  fDigits->Clear();
  fDigits->Expand(NCELLS);

  if (fMakeClu) {
    fClusters->Delete();
  }
}
//_______________________________________
float Reconstruction::Calibrate(int absId, bool isHG, float amp)
{
  return amp * 0.015;
}
//_______________________________________
float Reconstruction::CalibrateT(int absId, bool isHG, float time)
{
  return 100.e-9 * time;
}
