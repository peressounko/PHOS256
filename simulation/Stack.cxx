/// Geant4 example A01 adapted to Virtual Monte Carlo \n
/// The implementation of the Stack taken from the E03 example.

#include <TClonesArray.h>
#include <TParticle.h>

#include "Stack.h"

using namespace std;

/// \cond CLASSIMP
ClassImp(Stack)
  /// \endcond

  //_____________________________________________________________________________
  Stack::Stack(int size)
  : fParticles(0), fCurrentTrack(-1), fNPrimary(0)
{
  /// Standard constructor
  /// \param size  The stack size

  fParticles = new TClonesArray("TParticle", size);
}

// //_____________________________________________________________________________
// Stack::~Stack()
// {
//   /// Destructor

//   if (fParticles)
//     fParticles->Delete();
//   delete fParticles;
// }

// private methods

// public methods

//_____________________________________________________________________________
void Stack::PushTrack(int toBeDone, int parent, int pdg,
                      Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vx, Double_t vy,
                      Double_t vz, Double_t tof, Double_t polx, Double_t poly, Double_t polz,
                      TMCProcess mech, int& ntr, Double_t weight, int is)
{
  /// Create a new particle and push into stack;
  /// adds it to the particles array (fParticles) and if not done to the
  /// stack (fStack).
  /// Use TParticle::fMother[1] to store Track ID.
  /// \param toBeDone  1 if particles should go to tracking, 0 otherwise
  /// \param parent    number of the parent track, -1 if track is primary
  /// \param pdg       PDG encoding
  /// \param px        particle momentum - x component [GeV/c]
  /// \param py        particle momentum - y component [GeV/c]
  /// \param pz        particle momentum - z component [GeV/c]
  /// \param e         total energy [GeV]
  /// \param vx        position - x component [cm]
  /// \param vy        position - y component  [cm]
  /// \param vz        position - z component  [cm]
  /// \param tof       time of flight [s]
  /// \param polx      polarization - x component
  /// \param poly      polarization - y component
  /// \param polz      polarization - z component
  /// \param mech      creator process VMC code
  /// \param ntr       track number (is filled by the stack
  /// \param weight    particle weight
  /// \param is        generation status code

  const int kFirstDaughter = -1;
  const int kLastDaughter = -1;

  int trackId = GetNtrack();
  TParticle* particle = new ((*fParticles)[trackId]) TParticle(pdg, is, parent,
                                                               trackId, kFirstDaughter, kLastDaughter, px, py, pz, e, vx, vy, vz, tof);

  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(mech);
  fNPrimary++;

  // printf("toBeDone=%d, parent=%d, pdg=%d, fNPrimary=%d",toBeDone,parent,pdg, fNPrimary);
  // particle->Print() ;

  // if (parent < 0)
  // if (toBeDone)

  if (toBeDone)
    fStack.push(particle);

  ntr = GetNtrack() - 1;
}

//_____________________________________________________________________________
TParticle* Stack::PopNextTrack(int& itrack)
{
  /// Get next particle for tracking from the stack.
  /// \return       The popped particle object
  /// \param track  The index of the popped track

  itrack = -1;
  if (fStack.empty())
    return static_cast<TParticle*>(nullptr);

  TParticle* particle = fStack.top();
  fStack.pop();

  if (!particle)
    return static_cast<TParticle*>(nullptr);

  fCurrentTrack = particle->GetSecondMother();
  itrack = fCurrentTrack;

  return particle;
}

//_____________________________________________________________________________
TParticle* Stack::PopPrimaryForTracking(int i)
{
  /// Return \em i -th particle in fParticles.
  /// \return   The popped primary particle object
  /// \param i  The index of primary particle to be popped

  if (i < 0 || i >= fNPrimary){
    Error("GetPrimaryForTracking", "Index %d out of range %d", i, fNPrimary);
    return static_cast<TParticle*>(nullptr);
  }

  return (TParticle*)fParticles->At(i);
}

//_____________________________________________________________________________
void Stack::Print(Option_t* /*option*/) const
{
  /// Print info for all particles.

  cout << "Stack Info  " << endl;
  cout << "Total number of particles:   " << GetNtrack() << endl;
  cout << "Number of primary particles: " << GetNprimary() << endl;

  for (int i = 0; i < GetNtrack(); i++)
    GetParticle(i)->Print();
}

//_____________________________________________________________________________
void Stack::Reset()
{
  /// Delete contained particles, reset particles array and stack.

  fCurrentTrack = -1;
  fNPrimary = 0;
  fParticles->Clear();
}

//_____________________________________________________________________________
void Stack::SetCurrentTrack(int track)
{
  /// Set the current track number to a given value.
  /// \param  track The current track number

  fCurrentTrack = track;
}

//_____________________________________________________________________________
int Stack::GetNtrack() const
{
  /// \return  The total number of all tracks.

  return fParticles->GetEntriesFast();
}

//_____________________________________________________________________________
int Stack::GetNprimary() const
{
  /// \return  The total number of primary tracks.

  return fNPrimary;
}

//_____________________________________________________________________________
TParticle* Stack::GetCurrentTrack() const
{
  /// \return  The current track particle

  TParticle* current = GetParticle(fCurrentTrack);

  if (!current)
    Warning("GetCurrentTrack", "Current track not found in the stack");

  return current;
}

//_____________________________________________________________________________
int Stack::GetCurrentTrackNumber() const
{
  /// \return  The current track number

  return fCurrentTrack;
}

//_____________________________________________________________________________
int Stack::GetCurrentParentTrackNumber() const
{
  /// \return  The current track parent ID.

  TParticle* current = GetCurrentTrack();

  if (current)
    return current->GetFirstMother();
  else
    return -1;
}

//_____________________________________________________________________________
TParticle* Stack::GetParticle(int id) const
{
  /// \return   The \em id -th particle in fParticles
  /// \param id The index of the particle to be returned

  if (id < 0 || id >= fParticles->GetEntriesFast())
    Fatal("GetParticle", "Index %d out of range= %d", id, fParticles->GetEntriesFast());

  return (TParticle*)fParticles->At(id);
}
//_____________________________________________________________________________
void Stack::StoreTrack(int track) // mark trask to be stored
{
  // Mark track to be stored (set FirstDaughter field to 0 instead of def. -1)
  // track all all his ancestors
  if (track < 0 || track >= fParticles->GetEntriesFast()) {
    return;
  }
  TParticle* p = static_cast<TParticle*>(fParticles->At(track));
  if (p->GetFirstDaughter() == 0) { // already marked
    return;
  }
  p->SetFirstDaughter(0);
  track = p->GetMother(0);
  while (track >= 0) {
    p = static_cast<TParticle*>(fParticles->At(track));
    if (p->GetFirstDaughter() == 0) { // already marked
      return;
    }
    p->SetFirstDaughter(0);
    track = p->GetMother(0);
  }
  return;
}
//_____________________________________________________________________________
void Stack::Purge()
{
  // Remove all tracks not marked to be stored
  if (fParticles->GetEntriesFast() == 0) { // nothing to do
    return;
  }
  fLabels.clear();
  fLabels.reserve(fParticles->GetEntriesFast());
  int iToKeep = 0;
  TParticle* pNew = nullptr;
  // printf("PURGE: particls=%d \n",fParticles->GetEntriesFast()) ;
  for (int i = 0; i < fParticles->GetEntriesFast(); i++) {
    TParticle* p = static_cast<TParticle*>(fParticles->At(i));
    // printf("i=%d, E=%f, toKeep=%d \n",i,p->Energy(),p->GetFirstDaughter()) ;
    if (p->GetFirstDaughter() == 0) { // keep
      fLabels.push_back(iToKeep);     // new position
      if (i == iToKeep) {             // no need to copy
        if (p->GetMother(0) >= 0) {
          p->SetMother(0, fLabels[p->GetMother(0)]);
        }
      } else { // copy to new position
        pNew = static_cast<TParticle*>(fParticles->At(iToKeep));
        pNew = p;
        if (p->GetMother(0) >= 0) {
          pNew->SetMother(0, fLabels[p->GetMother(0)]);
        }
      }
      ++iToKeep;
    } else { // do not store
      fLabels.push_back(-1);
    }
  }
  // remove tail
  fParticles->RemoveRange(iToKeep, fParticles->GetEntriesFast() - 1);
  fParticles->Compress();
}