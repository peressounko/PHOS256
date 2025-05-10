#ifndef SIMULATION_H
#define SIMULATION_H

//
// class for running generation, simulation and digitization
// Hits and digits are created by typing:
//   Simulation sim;
//   sim.Run();

#include <TNamed.h>

class Simulation
{
 public:
  Simulation() = default;
  ~Simulation() = default;

  static Simulation* Instance()
  {
    if (!fgInstance)
      fgInstance = new Simulation();
    return fgInstance;
  }

  void SetNumberOfEvents(int nEvents);

  void SetRunNumber(int run);
  void SetSeed(int seed);

  void SetCalibration(std::string filename);

  bool Run(int nEvents = 0);

  bool RunSimulation(int nEvents = 0);
  bool RunDigitization();

  // Sets the name of the file from which the geometry is loaded
  void SetGeometryFile(std::string filename) { fGeometryFile = filename; }
  bool IsGeometryFromFile() const { return !fGeometryFile.size() == 0; }
  void FinishRun();

 private:
  void InitDB();
  void InitRunNumber();

  static Simulation* fgInstance = nullptr; // Static pointer to object

  bool fRunGeneration;           // generate prim. particles or not
  bool fRunSimulation;           // simulate detectors (hits) or not
  TString fLoadAlObjsListOfDets; // Load alignment data from CDB for these detectors
  TString fMakeDigits;           // create digits for these detectors

  int fNEvents;            // number of events
  TString fConfigFileName; // name of the config file

  int fRun;                  //! Run number, will be passed to CDB and gAlice!!
  int fSeed;                 //! Seed for random number generator
  bool fInitCDBCalled;       //! flag to check if CDB storages are already initialized
  bool fInitRunNumberCalled; //! flag to check if run number is already initialized

  bool fUseDetectorsFromGRP; // do not simulate detectors absent in the GRP

  TString fGeometryFile; // Geometry file

  ClassDefNV(Simulation, 1) // class for running generation, simulation and digitization
};

#endif
