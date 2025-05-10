///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running generation, simulation and digitization                 //
//                                                                           //
// Hits and digits are created by typing:                                    //
//                                                                           //
//   Simulation sim;                                                         //
//   sim.Run();                                                              //
//                                                                           //
// The Run method returns kTRUE in case of successful execution.             //
// The number of events can be given as argument to the Run method or it     //
// can be set by                                                             //
//                                                                           //
//   sim.SetNumberOfEvents(n);                                               //
//                                                                           //
// The configuration file Config.C                                           //
//                                                                           //
// The generation of particles and the simulation of detector hits can be    //
// switched on or off by                                                     //
//                                                                           //
//   sim.SetRunGeneration(kTRUE);   // generation of primary particles       //
//   sim.SetRunSimulation(kFALSE);  // but no tracking                       //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TVirtualMCApplication.h>
#include <TDatime.h>
#include <TInterpreter.h>

#include "AliAlignObj.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliCDBStorage.h"
#include "AliCTPRawData.h"
#include "AliCentralTrigger.h"
#include "AliCentralTrigger.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliDigitizer.h"
#include "AliESDEvent.h"
#include "AliFileUtilities.h"
#include "AliGRPObject.h"
#include "AliGenEventHeader.h"
#include "AliGenerator.h"
#include "AliGeomManager.h"
#include "AliHLTSimulation.h"
#include "AliHeader.h"
#include "AliLego.h"
#include "AliLegoGenerator.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliModule.h"
#include "AliPDG.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"
#include "AliRun.h"
#include "AliDigitizationInput.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "Simulation.h"
#include "AliSysInfo.h"
#include "AliVertexGenFile.h"
#include "AliLumiTools.h"
#include <TGraph.h>
#include <fstream>

using std::ofstream;
ClassImp(Simulation)

  //_____________________________________________________________________________
  Simulation::Simulation() :

                             fRunGeneratorOnly(kFALSE),
                             fRunGeneration(kTRUE),
                             fRunSimulation(kTRUE),
                             fLoadAlignFromCDB(kTRUE),
                             fLoadAlObjsListOfDets("ALL"),
                             fMakeSDigits("ALL"),
                             fMakeDigits("ALL"),
                             fTriggerConfig(""),
                             fMakeDigitsFromHits(""),
                             fWriteRawData(""),
                             fRawDataFileName(""),
                             fDeleteIntermediateFiles(kFALSE),
                             fWriteSelRawData(kFALSE),
                             fStopOnError(kFALSE),
                             fUseMonitoring(kFALSE),
                             fNEvents(1),
                             fConfigFileName(configFileName),
                             fGAliceFileName("galice.root"),
                             fEventsPerFile(),
                             fBkgrdFileNames(NULL),
                             fAlignObjArray(NULL),
                             fUseBkgrdVertex(kTRUE),
                             fRegionOfInterest(kFALSE),
                             fCDBUri(""),
                             fQARefUri(""),
                             fSpecCDBUri(),
                             fRun(-1),
                             fSeed(0),
                             fInitCDBCalled(kFALSE),
                             fInitRunNumberCalled(kFALSE),
                             fSetRunNumberFromDataCalled(kFALSE),
                             fEmbeddingFlag(kFALSE),
                             fLego(NULL),
                             fKey(0),
                             fUseVertexFromCDB(0),
                             fUseMagFieldFromGRP(0),
                             fGRPWriteLocation(Form("local://%s", gSystem->pwd())),
                             fUseDetectorsFromGRP(kTRUE),
                             fUseTimeStampFromCDB(0),
                             fTimeStart(0),
                             fTimeEnd(0),
                             fLumiDecayH(-1.), // by default, use lumi from CTP
                             fOrderedTimeStamps(),
                             fQADetectors("ALL"),
                             fQATasks("ALL"),
                             fRunQA(kTRUE),
                             fEventSpecie(AliRecoParam::kDefault),
                             fWriteQAExpertData(kTRUE),
                             fGeometryFile(),
                             fRunHLT(fgkRunHLTAuto),
                             fpHLT(NULL),
                             fWriteGRPEntry(kTRUE)
{
  // create simulation object with default parameters
  fgInstance = this;
}

//_____________________________________________________________________________
void Simulation::SetNumberOfEvents(int nEvents)
{
  // set the number of events for one run
  fNEvents = nEvents;
}

//_____________________________________________________________________________
void Simulation::SetSeed(int seed)
{
  // sets seed number
  // Activate it later within the Run() method

  fSeed = seed;
}

//_____________________________________________________________________________
bool Simulation::Run(int nEvents)
{
  // run the generation, simulation and digitization
  if (nEvents > 0)
    fNEvents = nEvents;

  InitDB();
  InitGeometry();

  gRandom->SetSeed(fSeed);

  // generation and simulation -> hits
  if (!RunSimulation()) {
    if (fStopOnError)
      return false;
  }

  // hits -> Digits
  if (!RunDigitization()) {
    return kFALSE;
  }

  return true;
}

//_____________________________________________________________________________
bool Simulation::RunSimulation(int nEvents)
{
  // run the generation and simulation
  if (nEvents <= 0)
    nEvents = fNEvents;

  // Init magnetic Field

  // Execute Config.C
  TInterpreter::EErrorCode interpreterError = TInterpreter::kNoError;
  gROOT->LoadMacro(fConfigFileName.data());
  Long_t interpreterResult = gInterpreter->ProcessLine(gAlice->GetConfigFunction(), &interpreterError);
  if (interpreterResult != 0 || interpreterError != TInterpreter::kNoError) {
    std::cout << "execution of config file " << fConfigFileName << " failed with error " interpreterError << std::endl;
    return false;
  }

  AliRunLoader::Instance()->CdGAFile();

  AliPDG::AddParticlesToPdgDataBase();

  TVirtualMC::GetMC()->SetMagField(TGeoGlobalMagField::Instance()->GetField());
  gAlice->GetMCApp()->Init();

  // Must be here because some MCs (G4) adds detectors here and not in Config.C
  gAlice->InitLoaders();
  AliRunLoader::Instance()->MakeTree("E");
  AliRunLoader::Instance()->LoadKinematics("RECREATE");
  AliRunLoader::Instance()->LoadTrackRefs("RECREATE");
  AliRunLoader::Instance()->LoadHits("all", "RECREATE");

  // Save stuff at the beginning of the file to avoid file corruption
  AliRunLoader::Instance()->CdGAFile();
  gAlice->Write();
  gAlice->SetEventNrInRun(-1); // important - we start Begin event from increasing current number in run
  AliSysInfo::AddStamp("RunSimulation_InitLoaders");

  AliRunLoader* runLoader = AliRunLoader::Instance();
  SetGAliceFile(runLoader->GetFileName());

  if (!gAlice->GetMCApp()->Generator()) {
    std::cout << "gAlice has no generator object. " << std::endl;
    return false;
  }

  if (!fRunSimulation) {
    gAlice->GetMCApp()->Generator()->SetTrackingFlag(0);
  }

  // Create the Root Tree with one branch per detector
  // Hits moved to begin event -> now we are crating separate tree for each event
  TVirtualMC::GetMC()->ProcessRun(nEvents);

  // End of this run, close files
  if (nEvents > 0)
    FinishRun();

  std::cout << "Stop_ProcessRun" << std::endl;

  return true;
}

//_____________________________________________________________________________
bool Simulation::RunGeneratorOnly()
{
  // Execute Config.C
  InitCDB();
  InitRunNumber();
  if (fUseMagFieldFromGRP) {
    AliGRPManager grpM;
    grpM.ReadGRPEntry();
    grpM.SetMagField();
    AliInfo("Field is locked now. It cannot be changed in Config.C");
  }

  TInterpreter::EErrorCode interpreterError = TInterpreter::kNoError;
  gROOT->LoadMacro(fConfigFileName.Data());
  Long_t interpreterResult = gInterpreter->ProcessLine(gAlice->GetConfigFunction(), &interpreterError);
  if (interpreterResult != 0 || interpreterError != TInterpreter::kNoError) {
    AliFatal(Form("execution of config file \"%s\" failed with error %d", fConfigFileName.Data(), (int)interpreterError));
  }

  // Setup the runloader and generator, check if everything is OK
  AliRunLoader* runLoader = AliRunLoader::Instance();
  AliGenerator* generator = gAlice->GetMCApp()->Generator();
  if (!runLoader) {
    AliError(Form(
      "gAlice has no run loader object. "
      "Check your config file: %s",
      fConfigFileName.Data()));
    return kFALSE;
  }
  if (!generator) {
    AliError(Form(
      "gAlice has no generator object. "
      "Check your config file: %s",
      fConfigFileName.Data()));
    return kFALSE;
  }

  runLoader->LoadKinematics("RECREATE");
  runLoader->MakeTree("E");

  // Create stack and header
  runLoader->MakeStack();
  AliStack* stack = runLoader->Stack();
  AliHeader* header = runLoader->GetHeader();

  // Intialize generator
  generator->Init();
  generator->SetStack(stack);

  // Run main generator loop

  for (int iev = 0; iev < fNEvents; iev++) {
    // Initialize event
    header->Reset(0, iev);
    runLoader->SetEventNumber(iev);
    stack->Reset();
    runLoader->MakeTree("K");

    // Generate event
    generator->Generate();

    // Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());
    stack->FinishEvent();
    header->SetStack(stack);
    runLoader->TreeE()->Fill();
    runLoader->WriteKinematics("OVERWRITE");
  }

  // Finalize
  generator->FinishRun();
  // Write file
  runLoader->WriteHeader("OVERWRITE");
  generator->Write();
  runLoader->Write();

  return kTRUE;
}

//_____________________________________________________________________________
bool Simulation::RunDigitization(const char* detectors,
                                 const char* excludeDetectors)
{
  // run the digitization and produce digits from sdigits
  AliCodeTimerAuto("", 0)

    // initialize CDB storage, run number, set CDB lock
    InitCDB();
  if (!SetRunNumberFromData())
    if (fStopOnError)
      return kFALSE;
  SetCDBLock();

  delete AliRunLoader::Instance();
  delete gAlice;
  gAlice = NULL;

  int nStreams = 1;
  if (fBkgrdFileNames)
    nStreams = fBkgrdFileNames->GetEntriesFast() + 1;
  int signalPerBkgrd = GetNSignalPerBkgrd();
  AliDigitizationInput digInp(nStreams, signalPerBkgrd);
  // digInp.SetEmbeddingFlag(fEmbeddingFlag);
  digInp.SetRegionOfInterest(fRegionOfInterest);
  digInp.SetInputStream(0, fGAliceFileName.Data());
  for (int iStream = 1; iStream < nStreams; iStream++) {
    const char* fileName = ((TObjString*)(fBkgrdFileNames->At(iStream - 1)))->GetName();
    digInp.SetInputStream(iStream, fileName);
  }
  TObjArray detArr;
  detArr.SetOwner(kTRUE);
  TString detStr = detectors;
  TString detExcl = excludeDetectors;
  if (!static_cast<AliStream*>(digInp.GetInputStream(0))->ImportgAlice()) {
    AliError("Error occured while getting gAlice from Input 0");
    return kFALSE;
  }
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader(digInp.GetInputStream(0)->GetFolderName());
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  //
  if (fUseDetectorsFromGRP) {
    AliInfo("Will run only for detectors seen in the GRP");
    DeactivateDetectorsAbsentInGRP(detArray);
  }
  //
  for (int iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*)detArray->At(iDet);
    if (!det || !det->IsActive())
      continue;
    if (!IsSelected(det->GetName(), detStr) || IsSelected(det->GetName(), detExcl))
      continue;
    AliDigitizer* digitizer = det->CreateDigitizer(&digInp);
    if (!digitizer || !digitizer->Init()) {
      AliError(Form("no digitizer for %s", det->GetName()));
      if (fStopOnError)
        return kFALSE;
      else
        continue;
    }
    detArr.AddLast(digitizer);
    AliInfo(Form("Created digitizer from SDigits -> Digits for %s", det->GetName()));
  }
  //
  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", detStr.Data()));
    if (fStopOnError)
      return kFALSE;
  }
  //
  int ndigs = detArr.GetEntriesFast();
  int eventsCreated = 0;
  AliRunLoader* outRl = digInp.GetOutRunLoader();
  while ((eventsCreated++ < fNEvents) || (fNEvents < 0)) {
    if (!digInp.ConnectInputTrees())
      break;
    digInp.InitEvent(); // this must be after call of Connect Input tress.
    if (outRl)
      outRl->SetEventNumber(eventsCreated - 1);
    static_cast<AliStream*>(digInp.GetInputStream(0))->ImportgAlice(); // use gAlice of the first input stream
    for (int id = 0; id < ndigs; id++) {
      ((AliDigitizer*)detArr[id])->Digitize("");
      AliSysInfo::AddStamp(Form("Digit_%s_%d", detArr[id]->GetName(), eventsCreated), 0, 2, eventsCreated);
    }
    digInp.FinishEvent();
  };
  digInp.FinishGlobal();
  //
  return kTRUE;
}

//_____________________________________________________________________________
AliRunLoader* Simulation::LoadRun(const char* mode) const
{
  // delete existing run loaders, open a new one and load gAlice

  delete AliRunLoader::Instance();
  AliRunLoader* runLoader =
    AliRunLoader::Open(fGAliceFileName.Data(),
                       AliConfig::GetDefaultEventFolderName(), mode);
  if (!runLoader) {
    AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
    return NULL;
  }
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    AliError(Form("no gAlice object found in file %s",
                  fGAliceFileName.Data()));
    return NULL;
  }
  return runLoader;
}

//_____________________________________________________________________________
bool Simulation::IsSelected(TString detName, TString& detectors) const
{
  // check whether detName is contained in detectors
  // if yes, it is removed from detectors

  // check if all detectors are selected
  if ((detectors.CompareTo("ALL") == 0) ||
      detectors.BeginsWith("ALL ") ||
      detectors.EndsWith(" ALL") ||
      detectors.Contains(" ALL ")) {
    detectors = "ALL";
    return kTRUE;
  }

  // search for the given detector
  bool result = kFALSE;
  if ((detectors.CompareTo(detName) == 0) ||
      detectors.BeginsWith(detName + " ") ||
      detectors.EndsWith(" " + detName) ||
      detectors.Contains(" " + detName + " ")) {
    detectors.ReplaceAll(detName, "");
    result = kTRUE;
  }

  // clean up the detectors string
  while (detectors.Contains("  "))
    detectors.ReplaceAll("  ", " ");
  while (detectors.BeginsWith(" "))
    detectors.Remove(0, 1);
  while (detectors.EndsWith(" "))
    detectors.Remove(detectors.Length() - 1, 1);

  return result;
}

//_____________________________________________________________________________
int Simulation::ConvertRaw2SDigits(const char* rawDirectory, const char* esdFileName, int N, int nSkip)
{
  //
  // Steering routine  to convert raw data in directory rawDirectory/ to fake SDigits.
  // These can be used for embedding of MC tracks into RAW data using the standard
  // merging procedure.
  //
  // If an ESD file is given the reconstructed vertex is taken from it and stored in the event header.
  //
  if (!gAlice) {
    AliError("no gAlice object. Restart aliroot and try again.");
    return kFALSE;
  }
  if (gAlice->Modules()->GetEntries() > 0) {
    AliError("gAlice was already run. Restart aliroot and try again.");
    return kFALSE;
  }

  AliInfo(Form("initializing gAlice with config file %s", fConfigFileName.Data()));

  gAlice->Announce();

  gROOT->LoadMacro(fConfigFileName.Data());
  gInterpreter->ProcessLine(gAlice->GetConfigFunction());

  if (AliCDBManager::Instance()->GetRun() >= 0) {
    SetRunNumber(AliCDBManager::Instance()->GetRun());
  } else {
    AliWarning("Run number not initialized!!");
  }

  AliRunLoader::Instance()->CdGAFile();

  AliPDG::AddParticlesToPdgDataBase();

  TVirtualMC::GetMC()->SetMagField(TGeoGlobalMagField::Instance()->GetField());

  gAlice->GetMCApp()->Init();

  // Must be here because some MCs (G4) adds detectors here and not in Config.C
  gAlice->InitLoaders();
  AliRunLoader::Instance()->MakeTree("E");
  AliRunLoader::Instance()->LoadKinematics("RECREATE");
  AliRunLoader::Instance()->LoadTrackRefs("RECREATE");
  AliRunLoader::Instance()->LoadHits("all", "RECREATE");

  //
  // Save stuff at the beginning of the file to avoid file corruption
  AliRunLoader::Instance()->CdGAFile();
  gAlice->Write();
  //
  //  Initialize CDB
  InitCDB();
  // AliCDBManager* man = AliCDBManager::Instance();
  // man->SetRun(0); // Should this come from rawdata header ?

  int iDet;
  //
  // Get the runloader
  AliRunLoader* runLoader = AliRunLoader::Instance();
  //
  // Open esd file if available
  TFile* esdFile = 0;
  TTree* treeESD = 0;
  AliESDEvent* esd = 0;
  if (esdFileName && (strlen(esdFileName) > 0)) {
    esdFile = TFile::Open(esdFileName);
    if (esdFile) {
      esd = new AliESDEvent();
      esdFile->GetObject("esdTree", treeESD);
      if (treeESD) {
        esd->ReadFromTree(treeESD);
        if (nSkip > 0) {
          AliInfo(Form("Asking to skip first %d ESDs events", nSkip));
        } else {
          nSkip = 0;
        }
      }
    }
  }

  //
  // Create the RawReader
  TString fileName(rawDirectory);
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  if (!rawReader)
    return (kFALSE);

  //     if (!fEquipIdMap.IsNull() && fRawReader)
  //       fRawReader->LoadEquipmentIdsMap(fEquipIdMap);
  //
  // Get list of detectors
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  if (fUseDetectorsFromGRP) {
    AliInfo("Will run only for detectors seen in the GRP");
    DeactivateDetectorsAbsentInGRP(detArray);
  }
  //
  // Get Header
  AliHeader* header = runLoader->GetHeader();
  // Event loop
  int nev = 0;
  int prevEsdID = nSkip - 1;
  while (kTRUE) {
    if (!(rawReader->NextEvent()))
      break;
    runLoader->SetEventNumber(nev);
    runLoader->GetHeader()->Reset(rawReader->GetRunNumber(),
                                  nev, nev);
    runLoader->GetEvent(nev);
    AliInfo(Form("We are at event %d", nev));
    //
    // Detector loop
    TString detStr = fMakeSDigits;
    for (iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
      AliModule* det = (AliModule*)detArray->At(iDet);
      if (!det || !det->IsActive())
        continue;
      if (IsSelected(det->GetName(), detStr)) {
        AliInfo(Form("Calling Raw2SDigits for %s", det->GetName()));
        det->Raw2SDigits(rawReader);
        rawReader->Reset();
      }
    } // detectors

    //
    //  If ESD information available obtain reconstructed vertex and store in header.
    if (treeESD) {
      int rawID = rawReader->GetEventIndex();
      ULong64_t rawGID = rawReader->GetEventIdAsLong();

      int esdID = nSkip + rawID;
      if (esdID > treeESD->GetEntriesFast())
        esdID = treeESD->GetEntriesFast();
      bool bFound = kFALSE;
      while (esdID > prevEsdID) {
        treeESD->GetEvent(esdID);
        if (ULong64_t(esd->GetHeader()->GetEventIdAsLong()) == rawGID) {
          bFound = kTRUE;
          prevEsdID = esdID;
          break; // found!
        }
        esdID--;
      }
      if (!bFound) {
        AliInfo("Failed to find event ... skipping");
        continue;
      }

      AliInfo(Form("Selected event %d correspond to event %d in raw and to %d in esd", nev, rawReader->GetEventIndex(), prevEsdID));
      const AliESDVertex* esdVertex = esd->GetPrimaryVertex();
      Double_t position[3];
      esdVertex->GetXYZ(position);
      AliGenEventHeader* mcHeader = new AliGenEventHeader("ESD");
      TArrayF mcV;
      mcV.Set(3);
      for (int i = 0; i < 3; i++)
        mcV[i] = position[i];
      mcHeader->SetPrimaryVertex(mcV);
      header->Reset(0, nev);
      header->SetGenEventHeader(mcHeader);
      AliInfo(Form("***** Saved vertex %f %f %f \n", position[0], position[1], position[2]));
    }
    //
    //      Finish the event
    runLoader->TreeE()->Fill();
    AliInfo(Form("Finished event %d", nev));
    nev++;
    if (N > 0 && nev >= N)
      break;
  } // events

  delete rawReader;
  //
  //  Finish the run
  runLoader->CdGAFile();
  runLoader->WriteHeader("OVERWRITE");
  runLoader->WriteRunLoader();

  return nev;
}

//_____________________________________________________________________________
void Simulation::FinishRun()
{
  //
  // Called at the end of the run.
  //

  if (IsLegoRun()) {
    AliDebug(1, "Finish Lego");
    AliRunLoader::Instance()->CdGAFile();
    fLego->FinishRun();
  }

  // Clean detector information
  TIter next(gAlice->Modules());
  AliModule* detector;
  while ((detector = dynamic_cast<AliModule*>(next()))) {
    AliDebug(2, Form("%s->FinishRun()", detector->GetName()));
    detector->FinishRun();
  }

  AliDebug(1, "AliRunLoader::Instance()->WriteHeader(OVERWRITE)");
  AliRunLoader::Instance()->WriteHeader("OVERWRITE");

  // Write AliRun info and all detectors parameters
  AliRunLoader::Instance()->CdGAFile();
  gAlice->Write(0, TObject::kOverwrite);                   // write AliRun
  AliRunLoader::Instance()->Write(0, TObject::kOverwrite); // write RunLoader itself
  //
  if (gAlice->GetMCApp())
    gAlice->GetMCApp()->FinishRun();
  AliRunLoader::Instance()->Synchronize();
}

//_____________________________________________________________________________
int Simulation::GetDetIndex(const char* detector)
{
  // return the detector index corresponding to detector
  int index = -1;
  for (index = 0; index < fgkNDetectors; index++) {
    if (strcmp(detector, fgkDetectorName[index]) == 0)
      break;
  }
  return index;
}
