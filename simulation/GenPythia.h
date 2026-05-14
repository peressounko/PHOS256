#ifndef GENPYTHIA_H
#define GENPYTHIA_H

#include "GenBox.h"

#include "Pythia8/Pythia.h"

class GenPythia : public GenBox
{
 public:
  GenPythia() = default;
  GenPythia(const char* xmlDir);
  ~GenPythia() = default;

  bool Initialize(int idAin, int idBin, double eLab);

  // methods
  virtual void Generate();

  // set methods
  void SetPythiaSeed(UInt_t seed) { fSeed = seed; }
  void SetStore(int kind) { fStore = kind; }

 private:
  enum store { kAll = 0,
               kFinalOnly };
  // data members
  Pythia8::Pythia* fPythia = nullptr; //! The pythia8 instance
  // Pythia8::Pythia         *fPythiaPartonLevel = nullptr;     //! The pythia8 instance for the parton level object (used for superposition of events)
  std::string fgXmldocPath; //  path to xmldoc
  int fStore = kAll;        //  store all particles or final only
  UInt_t fSeed = 0;         //  Zero seed triggers system clock initialization. Valid range covers 1 to 900,000,000.

  ClassDef(GenPythia, 1) // GenPythia
};
#endif // GENPYTHIA_H
