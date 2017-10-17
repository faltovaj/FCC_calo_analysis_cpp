#include "ShowerProfiles.h"

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/PositionedCaloHitCollection.h"

#include "TVector3.h"
#include "TMath.h"

// STL
#include <vector>
#include <iostream>

ShowerProfiles::ShowerProfiles(double aEnergy, double aSf):
  m_energy(aEnergy), m_sf(aSf) {
  Initialize_histos();
}

ShowerProfiles::~ShowerProfiles() {}


void ShowerProfiles::Initialize_histos()
{

  h_hitEnergy = new TH1F("h_hitEnergy","", 200, 0, m_energy);
  m_histograms.push_back(h_hitEnergy);

  h_cellEnergy = new TH1F("h_cellenergy","", 100, m_energy-0.2*m_energy, m_energy+0.2*m_energy);
  m_histograms.push_back(h_cellEnergy);

  h_radialProfile_x0 = new TH1F("h_radialProfile_x0","", 100, 0, 10);
  m_histograms.push_back(h_radialProfile_x0);

  h_longProfile_x0 = new TH1F("h_longProfile_x0","", 45, 0, 45);
  m_histograms.push_back(h_longProfile_x0);

  h_radialProfile_particle_x0 = new TH1F("h_radialProfile_particle_x0","", 100, 0, 10);
  m_histograms.push_back(h_radialProfile_particle_x0);

  h_longProfile_particle_x0 = new TH1F("h_longProfile_particle_x0","", 45, 0, 45);
  m_histograms.push_back(h_longProfile_particle_x0);

  h_radialProfile = new TH1F("h_radialProfile","", 100, 0, 100);
  m_histograms.push_back(h_radialProfile);

  h_longProfile = new TH1F("h_longProfile","", 45, 0, RcaloThickness);
  m_histograms.push_back(h_longProfile);

  h_radialProfile_particle = new TH1F("h_radialProfile_particle","", 100, 0, 100);
  m_histograms.push_back(h_radialProfile_particle);

  h_longProfile_particle = new TH1F("h_longProfile_particle","", 45, 0, RcaloThickness);
  m_histograms.push_back(h_longProfile_particle);

  h_ptGen = new TH1F("h_ptGen","", 100, m_energy-0.2*m_energy, m_energy+0.2*m_energy);
  m_histograms.push_back(h_ptGen);

  h_pdgGen = new TH1F("h_pdgGen","", 200,-100,100);
  m_histograms.push_back(h_pdgGen);
}


void ShowerProfiles::processEvent(podio::EventStore& aStore, int aEventId, bool aVerbose) {
  //Get the collections
  const fcc::MCParticleCollection*  colMCParticles(nullptr);
  const fcc::PositionedCaloHitCollection*     colECalPositionedHits(nullptr);

  bool colMCParticlesOK = aStore.get("GenParticles", colMCParticles);
  bool colECalPositionedHitsOK     = aStore.get("ECalPositions" , colECalPositionedHits);

  //Total hit energy per event
  SumE_hit = 0.;
  //EM shower axis - assuming single shower per event!!!!
  //Direction of gen. particle
  TVector3 directionParticle(0.,0.,0.);
  //Direction of hits in the first layer
  TVector3 directionHits(0.,0.,0.);

  //MCParticle and Vertices collection
  if (colMCParticlesOK) {
    if (aVerbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #MCTruthParticles:     " << colMCParticles->size()    << std::endl;
    }
    //Loop through the collection
    for (auto& iparticle=colMCParticles->begin(); iparticle!=colMCParticles->end(); ++iparticle) {
      TVector3 particle(iparticle->core().p4.px,iparticle->core().p4.py,iparticle->core().p4.pz);
      //unit vector
      directionParticle = particle.Unit();
      //Fill histograms
      h_ptGen->Fill( sqrt( pow(iparticle->core().p4.px,2)+pow(iparticle->core().p4.py,2) ) );
      h_pdgGen->Fill( iparticle->core().pdgId );
    }
  }
  else {
    if (aVerbose) {
      std::cout << "No MCTruth info available" << std::endl;
    }
  }


  //PositionedHits collection
  if (colECalPositionedHitsOK) {
    if (aVerbose) {
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #ECalPositionedHits:     " << colECalPositionedHits->size()    << std::endl;;
    }
    TH1F *h_hitR = new TH1F("h_hitR", "h_hitR", 1000, RcaloMin, RcaloMin+1.5*layerThickness);
    TH1F *h_hitPhi = new TH1F("h_hitPhi", "h_hitPhi", 100000, -TMath::Pi(), TMath::Pi());
    TH1F *h_hitEta = new TH1F("h_hitEta", "h_hitEta", 1000000, -EtaMax, EtaMax);

    //Loop through the collection
    //First: find the mean position in the first layer, calculate total energy
    for (auto& iecl=colECalPositionedHits->begin(); iecl!=colECalPositionedHits->end(); ++iecl)
      {
   double hitEnergy = iecl->core().energy;
   SumE_hit += hitEnergy;

   TVector3 hit_position(iecl->position().x,iecl->position().y,iecl->position().z);
   double R_firstLayer = RcaloMin+layerThickness;
   if (hit_position.Pt()<RcaloMin) {
     std::cout <<"Hit before the calorimeter????? Please check the value of RcaloMin!"<< std::endl;
   }
   if (hit_position.Pt()<R_firstLayer) {
     h_hitR->Fill(hit_position.Perp(), hitEnergy);
     h_hitPhi->Fill(hit_position.Phi(), hitEnergy);
     h_hitEta->Fill(hit_position.Eta(), hitEnergy);
   }

      }

    //Fill histograms
    h_hitEnergy->Fill(SumE_hit);
    h_cellEnergy->Fill(SumE_hit*m_sf);

    TVector3 meanFirstLayer_vector;
    /*
    //Mean value shifted from truth sometimes, maximum bin works better
    double meanR = h_hitR->GetMean();
    double meanPhi = h_hitPhi->GetMean();
    double meanEta = h_hitEta->GetMean();
    */
    //std::cout << "Low edge " << h_hitEta->GetBinLowEdge(1) << " maximum " << h_hitEta->GetMaximumBin() << " mean " << h_hitEta->GetMean() << std::endl;
    double meanR = h_hitR->GetBinLowEdge(1)+h_hitR->GetMaximumBin()*h_hitR->GetBinWidth(1)+h_hitR->GetBinWidth(1)*0.5;
    double meanPhi = h_hitPhi->GetBinLowEdge(1)+h_hitPhi->GetMaximumBin()*h_hitPhi->GetBinWidth(1)+h_hitPhi->GetBinWidth(1)*0.5;
    double meanEta = h_hitEta->GetBinLowEdge(1)+h_hitEta->GetMaximumBin()*h_hitEta->GetBinWidth(1)+h_hitEta->GetBinWidth(1)*0.5;

    /*
    if (i==1) {
      TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
      c1->Divide(2,2);
      c1->cd(1);
      h_hitR->Draw();
      c1->cd(2);
      h_hitPhi->Draw();
      c1->cd(3);
      h_hitEta->Draw();
      c1->SaveAs("hitFirstLayer.gif");
      delete c1;
    }
    */

    meanFirstLayer_vector.SetPtEtaPhi(meanR, meanEta, meanPhi);
    //direction from hits in the first calorimeter layer
    directionHits = meanFirstLayer_vector.Unit();

    delete h_hitR;
    delete h_hitPhi;
    delete h_hitEta;

    //Second loop: radial & longitudinal profiles
    for (auto& iecl=colECalPositionedHits->begin(); iecl!=colECalPositionedHits->end(); ++iecl)
     {
       double hitEnergy = iecl->core().energy;
       TVector3 hitPosition(iecl->position().x,iecl->position().y,iecl->position().z);

       //gen. particle direction
       TVector3 showerStart_particle = directionParticle*RcaloMin;
       TVector3 hitVector_particle = hitPosition-showerStart_particle;
       double hitLong_particle = directionParticle*hitVector_particle;
       double hitRadial_particle = hitVector_particle.Perp(directionParticle);

       //hits in first layer direction
       //Start of the hits in the first active layer
       TVector3 showerStart = directionHits*meanFirstLayer_vector.Mag();
       //TVector3 showerStart = directionHits*RcaloMin;
       TVector3 hitVector = hitPosition-showerStart;
       double hitLong = directionHits*hitVector;
       double hitRadial = hitVector.Perp(directionHits);

       /*
       if ((hitEnergy>1.2)&&((hitRadial-hitRadial_particle)/X0>1)) {
    std::cout << "Hit info: energy " << hitEnergy << std::endl;
    std::cout << "Position " << hitPosition.Pt() << " " << hitPosition.Phi() << " " << hitPosition.Eta() << std::endl;

    std::cout << "Direction " << directionHits.Pt() << " " << directionHits.Phi() << " " << directionHits.Eta() << std::endl;
    std::cout << "Direct p  " << directionParticle.Pt() << " " << directionParticle.Phi() << " " << directionParticle.Eta() << std::endl;
    std::cout << "Dir diff  " << (directionHits-directionParticle).Pt() << " "<< (directionHits-directionParticle).Phi() << " " << (directionHits-directionParticle).Eta() << std::endl;
       }
       */
       //Fill longitudinal and radial profile histograms
       h_radialProfile_particle_x0->Fill(hitRadial_particle/X0, hitEnergy*m_sf);
       h_radialProfile_particle->Fill(hitRadial_particle, hitEnergy*m_sf);
       h_longProfile_particle_x0->Fill(hitLong_particle/X0, hitEnergy*m_sf);
       h_longProfile_particle->Fill(hitLong_particle, hitEnergy*m_sf);
       h_radialProfile_x0->Fill(hitRadial/X0, hitEnergy*m_sf);
       h_radialProfile->Fill(hitRadial, hitEnergy*m_sf);
       h_longProfile_x0->Fill(hitLong/X0, hitEnergy*m_sf);
       h_longProfile->Fill(hitLong, hitEnergy*m_sf);

     }
   if (aVerbose) std::cout << "Total hit energy (GeV): " << SumE_hit << " total cell energy (GeV): " << SumE_hit*m_sf << " hit collection size: " << colECalPositionedHits->size() << std::endl;

  }
  else {
    if (aVerbose) {
      std::cout << "No CaloPositionedHits Collection!!!!!" << std::endl;
    }
  }
}

void ShowerProfiles::finishLoop(int aNumEvents, bool aVerbose) {
  h_radialProfile_x0->Scale(1./(double)aNumEvents);
  h_longProfile_x0->Scale(1./(double)aNumEvents);
  h_radialProfile_particle_x0->Scale(1./(double)aNumEvents);
  h_longProfile_particle_x0->Scale(1./(double)aNumEvents);

  h_radialProfile->Scale(1./(double)aNumEvents);
  h_longProfile->Scale(1./(double)aNumEvents);
  h_radialProfile_particle->Scale(1./(double)aNumEvents);
  h_longProfile_particle->Scale(1./(double)aNumEvents);

  //  std::cout << "Total energy: " << h_cellEnergy->GetMean() << std::endl;
  //  std::cout << "End of loop" << std::endl;
}
