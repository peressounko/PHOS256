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

  // get methods
  virtual int GetNtrack() const;
  virtual int GetNprimary() const;
  virtual TParticle* GetCurrentTrack() const;
  virtual int GetCurrentTrackNumber() const;
  virtual int GetCurrentParentTrackNumber() const;
  TParticle* GetParticle(int id) const;

 private:
  // data members
  std::stack<TParticle*> fStack;      //!< The stack of particles (transient)
  TClonesArray* fParticles = nullptr; ///< The array of particle (persistent)
  int fCurrentTrack = -1;             ///< The current track number
  int fNPrimary = 0;                  ///< The number of primaries

  ClassDef(Stack, 1) // Stack
};

#endif // STACK_H
