#ifndef PIONANALYSIS_H
#define PIONANALYSIS_H

#include "BaseAnalysis.h"
#include "Decoder.h"

#include "TObject.h"
#include "TH2F.h"

namespace podio {
  class EventStore;
  class ROOTReader;
}

class PionAnalysis: public BaseAnalysis {

 public:
  PionAnalysis(const std::string& aCluserCollName, const std::string& aPosHitCollName, double aEnergy, double aEtaMax, double aPhiMax);
  PionAnalysis(const std::string& aCluserCollName, const std::string& aPosHitCollName, double aEnergy, double aEtaMax, double aPhiMax,
    const std::string& aCellCollName,  const std::string& aBitfield, const std::string& aLayerField, int aLayerFirst, int aLayerLast, double aFirstLayerSF,
    double aP0p0, double aP0p1, double aP1p0, double aP1p1);
  ~PionAnalysis();

  void Initialize_histos();

  // energy
  TH1F* hPiEnergy;
  // energy, correction for material in front
  TH1F* hPiEnergyCorrected;

 private:
  virtual void processEvent(podio::EventStore& store, int aEventId, bool verbose) final;
  virtual void finishLoop(int aNumEvents, bool aVerbose) final;
  std::string m_clusterCollName;
  std::string m_particleCollName;
  std::string m_readout;
  double m_energy;
  double m_etaMax;
  double m_phiMax;
  bool m_ifCorrectForUpstream;
  std::string m_cellCollName;
  Decoder m_decoder;
  std::string m_layerFieldName;
  int m_firstLayerFirstId;
  int m_firstLayerLastId;
  double m_firstLayerSF;
  double m_P0p0;
  double m_P0p1;
  double m_P1p0;
  double m_P1p1;

};

#endif /* PIONANALYSIS_H */
