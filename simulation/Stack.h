#ifndef STACK_H
#define STACK_H
#include <TVirtualMCStack.h>
#include <stack>

class TParticle;
class TClonesArray;

class Stack : public TVirtualMCStack
{
 public:
  Stack(int size);
  Stack() = default;
  virtual ~Stack() = default;

  // methods
  virtual void PushTrack(int toBeDone, int parent, int pdg, Double_t px,
                         Double_t py, Double_t pz, Double_t e, Double_t vx, Double_t vy, Double_t vz,
                         Double_t tof, Double_t polx, Double_t poly, Double_t polz, TMCProcess mech,
                         int& ntr, Double_t weight, int is);
  virtual TParticle* PopNextTrack(int& track);
  virtual TParticle* PopPrimaryForTracking(int i);
  virtual void Print(Option_t* option = "") const;
  void Reset();

  // set methods
  virtual void SetCurrentTrack(int track);
  void StoreTrack(int track); // mark trask to be stored
  void Purge();               // Remove all tracks not marked to be stored
  int GetNewLabel(int oldLab)
  {
    if (oldLab < 0 || oldLab >= fLabels.size()) {
      return -1;
    } else
      return fLabels[oldLab];
  }

  // get methods
  virtual int GetNtrack() const;
  virtual int GetNprimary() const;
  virtual TParticle* GetCurrentTrack() const;
  virtual int GetCurrentTrackNumber() const;
  virtual int GetCurrentParentTrackNumber() const;
  TParticle* GetParticle(int id) const;
  TClonesArray* GetParticles() { return fParticles; }

 private:
  enum { kKeep = 1,
         kToBeDone = 2,
         kTransport = 4 };
  TClonesArray* fParticles = nullptr; ///< The array of particle (persistent)
  TParticle* fCurrentTrack = nullptr; ///< The current track number
  int fNPrimary = 0;                  ///< The number of primaries
  int fCurrentTrackNumber = 0;
  std::vector<int> fLabels; ///< remapping of labels

  ClassDef(Stack, 1) // Stack
};

#endif // STACK_H
