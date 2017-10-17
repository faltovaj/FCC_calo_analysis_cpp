#ifndef SHOWERPROFILES_H
#define SHOWERPROFILES_H

#include "BaseAnalysis.h"

#include "TObject.h"
#include "TH1F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class ShowerProfiles: public BaseAnalysis {

 public:
  ShowerProfiles(double aEnergy, double aSf);
  ~ShowerProfiles();

  virtual void processEvent(podio::EventStore& store, int aEventId, bool verbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;

  /// Hit energy
  TH1F* h_hitEnergy;
  /// Calibrated hit energy (SF)
  TH1F* h_cellEnergy;
  /// Radial profile, direction of the shower from hits in the first layer
  TH1F* h_radialProfile_x0;
  /// Long. profile, direction of the shower from hits in the first layer
  TH1F* h_longProfile_x0;
  /// Radial profile, direction of the shower from the generated particle
  TH1F* h_radialProfile_particle_x0;
  /// Long. profile, direction of the shower from the generated particle
  TH1F* h_longProfile_particle_x0;
  /// Radial profile, direction of the shower from hits in the first layer                                                                                                                    
  TH1F* h_radialProfile;
  /// Long. profile, direction of the shower from hits in the first layer                                                                                                                      
  TH1F* h_longProfile;
  /// Radial profile, direction of the shower from the generated particle                                                                                                                       
  TH1F* h_radialProfile_particle;
  /// Long. profile, direction of the shower from the generated particle                                                                                                                        
  TH1F* h_longProfile_particle;
  /// pt of the generated particle
  TH1F* h_ptGen;
  /// PDG code of the generated particle
  TH1F* h_pdgGen;

 private:
  void Initialize_histos();
  double m_energy;
  double m_sf;
  double SumE_hit;    // Total hit energy per event
  const double RcaloMin = 1920.;
  const double RcaloThickness = 650.;
  const double layerThickness = 80.;
  const double EtaMax = 10.0;
  const double X0 = 15.586; //average X0 for 4 mm Lar + 2 mm Pb

};

#endif /* SHOWERPROFILES_H */
