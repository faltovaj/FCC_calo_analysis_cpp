#include "PionAnalysis.h"

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHitCollection.h"

#include "TVector3.h"
#include "TLorentzVector.h"

// STL
#include <vector>
#include <iostream>
#include <bitset>
#include <cmath>

PionAnalysis::PionAnalysis(const std::string& aClusterCollName, const std::string& aParticleCollName,
			   double aEnergy, double aEtaMax, double aPhiMax):
  m_clusterCollName(aClusterCollName), m_particleCollName(aParticleCollName), m_energy(aEnergy), m_etaMax(aEtaMax), m_phiMax(aPhiMax),
  m_ifCorrectForUpstream(false), m_cellCollName(""), m_decoder("") , m_layerFieldName(""), m_firstLayerFirstId(0), m_firstLayerLastId(0), m_firstLayerSF(0),
  m_P0p0(0), m_P0p1(0), m_P1p0(0), m_P1p1(0){
  Initialize_histos();
}

PionAnalysis::PionAnalysis(const std::string& aClusterCollName, const std::string& aParticleCollName,
  double aEnergy, double aEtaMax, double aPhiMax,
  const std::string& aCellCollName, const std::string& aBitfield, const std::string& aLayerField, int aLayerFirst, int aLayerLast, double aFirstLayerSF,
  double aP0p0, double aP0p1, double aP1p0, double aP1p1):
  m_clusterCollName(aClusterCollName), m_particleCollName(aParticleCollName), m_energy(aEnergy),
  m_ifCorrectForUpstream(true), m_cellCollName(aCellCollName),  m_decoder(aBitfield), m_layerFieldName(aLayerField), m_firstLayerFirstId(aLayerFirst), m_firstLayerLastId(aLayerLast), m_firstLayerSF(aFirstLayerSF),
  m_P0p0(aP0p0), m_P0p1(aP0p1), m_P1p0(aP1p0), m_P1p1(aP1p1) {
  Initialize_histos();
}

PionAnalysis::~PionAnalysis(){}


void PionAnalysis::Initialize_histos() {

  hPiEnergy = new TH1F("piEnergy","Energy of pion", 99,0.6*m_energy,1.4*m_energy);
  hPiEnergyCorrected = new TH1F("piEnergyCorrected","Upstream corrected energy of pion", 99,0.2*m_energy,1.8*m_energy);

  m_histograms.push_back(hPiEnergy);
  m_histograms.push_back(hPiEnergyCorrected);

}

void PionAnalysis::processEvent(podio::EventStore& aStore, int aEventId, bool aVerbose) {
  // Get the collections
  const fcc::CaloClusterCollection* clusters(nullptr);
  const fcc::MCParticleCollection* particles(nullptr);

  bool testParticles = aStore.get(m_particleCollName, particles);       
  bool testClusters = aStore.get(m_clusterCollName, clusters);

  TVector3 momentum;

  // Get generated particles - assuming single particle events
  if (testParticles) {
    if (particles->size() > 1) {
      std::cout << "This is not a single particle event! Number of particles: " << particles->size() << std::endl;
    }
    //Loop through the collection
    for (const auto ipart = particles->begin(); ipart != particles->end(); ++ipart) {
      if (aVerbose) {
        std::cout << "Particle at " << ipart->core().vertex.x
                  << " , " <<  ipart->core().vertex.y
                  << " , " <<  ipart->core().vertex.z
                  << "  with momentum " << ipart->core().p4.px
                  << " , " <<  ipart->core().p4.py
                  << " , " <<  ipart->core().p4.pz
                  << "  and mass " <<  ipart->core().p4.mass << " GeV" << std::endl;
      }
      momentum = TVector3(ipart->core().p4.px, ipart->core().p4.py, ipart->core().p4.pz);
    }
  } else {
    std::cout << "No MC Particle Collection in the event." << std::endl;
    return;
  }

  const int nGamma = 2;
  // Get clusters reconstructed in an event
  if (testClusters) {
    if (aVerbose) {
      std::cout << "Number of clusters: " << clusters->size() << std::endl;
    }
    //Reconstructed photons (max. 2)
    TLorentzVector p_gamma[nGamma];
    for (int i = 0; i<nGamma; i++) {
      p_gamma[i].SetPxPyPzE(0,0,0,0);
    }
    int index = 0;
    for (const auto iclu = clusters->begin(); iclu != clusters->end(); ++iclu) {
      if (aVerbose) {
        std::cout << "Cluster reconstructed at " << iclu->core().position.x
                  << " , " <<  iclu->core().position.y
                  << " , " <<  iclu->core().position.z
                  << "  with energy " <<  iclu->core().energy << " GeV" << std::endl;
      }
      double E = iclu->core().energy;
      TVector3 pos (iclu->core().position.x, iclu->core().position.y, iclu->core().position.z);
      double phi = pos.Phi();
      double eta = pos.Eta();
      double pt = E / cosh(eta);
      if (index<nGamma) {
	p_gamma[index].SetE(E);
	p_gamma[index].SetPx(pt*cos(phi));
	p_gamma[index].SetPy(pt*sin(phi));
	p_gamma[index].SetPz(pt*sinh(eta));
      }
      else {
	break;
      }
      index += 1;
    }

    //Pion energy: one or two clusters 
    if (index>0) {
      TLorentzVector p_pi(0,0,0,0);
      for (int i = 0; i<nGamma; i++) {
	p_pi += p_gamma[i];
      }
   
      // fill pion energy into the histogram
      hPiEnergy->Fill(p_pi.E());
      //position of the pion
      double etaPiPos = p_pi.Eta();
      double phiPiPos = p_pi.Phi();
      // Get energy in the first layer, within the reconstructed cluster
      double EfirstLayer = 0.;
      if( m_ifCorrectForUpstream ) {
	// get cells to calculate energy deposited in first layer
	const fcc::CaloHitCollection* cells(nullptr);
	bool testCells = aStore.get(m_cellCollName, cells);
	if (testCells) {
	  if (aVerbose) {
	    std::cout << "Number of cells: " << cells->size() << std::endl;
	  }
	  uint verb=0;
	  //Setup for the cell matching to the pion
	  int m_noEta = 20;
	  int m_noPhi = 20;
	  double m_dEta = 0.01;
	  double m_dPhi = 0.01227;
	  for (const auto icell = cells->begin(); icell != cells->end(); ++icell) {
	    int layerId = m_decoder.value(m_layerFieldName,icell->core().cellId);
	    int etaId = m_decoder.value("eta",icell->core().cellId);
	    int phiId = m_decoder.value("phi",icell->core().cellId);
	    double etaPos = m_decoder.segmentationPosition(etaId, m_dEta, -fabs(m_etaMax));
	    double phiPos = m_decoder.segmentationPosition(phiId, m_dPhi, -fabs(m_phiMax));
	    if( layerId >= m_firstLayerFirstId && layerId <= m_firstLayerLastId ) {
	      // Check on eta & phi position: if within window
	      // Check is on position, but cluster position and cell position are centres of cells - no need to check ID
	      if ((fabs(etaPos - etaPiPos) < m_noEta * m_dEta) &&
		  ((fabs(phiPos - phiPiPos) < m_noPhi * m_dPhi) || (fabs(phiPos - phiPiPos) > 2 * M_PI - m_noPhi * m_dPhi)  )) {
		EfirstLayer += icell->core().energy;
	      }
	    }
	  }
	  // if cells were already calibrated to EM scale, scale them back
	  EfirstLayer *= m_firstLayerSF;
	} else {
	  std::cout << "No Cell Collection in the event." << std::endl;
	  return;
	}
	// correct for energy upstream (lost in tracker, cryostat...)
	// calculate parameters based on reconstructed energy
	double EupstreamP0 = m_P0p0 + m_P0p1 * p_pi.E();
	double EupstreamP1 = m_P1p0 + m_P1p1 / sqrt( p_pi.E() );
	double Eupstream = EupstreamP0 + EupstreamP1 * EfirstLayer;
	hPiEnergyCorrected->Fill(p_pi.E() + Eupstream);
      }
    }
  }
  else {
    std::cout << "No Cluster Collection in the event." << std::endl;
    return;
  }
}

void PionAnalysis::finishLoop(int aNumEvents, bool aVerbose) {
  int nentries = hPiEnergy->GetEntries();
  hPiEnergy->Scale(1./nentries);
  hPiEnergyCorrected->Scale(1./nentries);
}
