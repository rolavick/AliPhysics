/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// c++ headers
#include <iostream>
#include <fstream>
#include <map>

// root headers
#include <TMath.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TObjString.h>
#include <TList.h>
#include <TChain.h>

// aliroot headers
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliMCEvent.h>
#include <AliMCParticle.h>
#include <AliAODInputHandler.h>
#include <AliMuonTrackCuts.h>

// my headers
#include "AliAnalysisTaskNanoFwdJpsiCentMult.h"
// ----------------------------------------------------------------------------------------------------------------------------------
class AliAnalysisTaskNanoFwdJpsiCentMult;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskNanoFwdJpsiCentMult) // classimp: necessary for root

// ----------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoFwdJpsiCentMult::AliAnalysisTaskNanoFwdJpsiCentMult() : AliAnalysisTaskSE(),
  fMuonTrackCuts(0x0), fAOD(0), fOutputList(0),fCounterH(0), fNumberMuonsH(0), fNtracksNoJpsiH(0), fNtracksCentralBarrelH(0),
  fAllTracksTree(0), fRunNum(0), fL0inputs(0),
  fZNCEnergy(-999), fZNAEnergy(-999),
  fV0ADecision(-10), fV0CDecision(-10),  fV0AFiredCells(-10), fV0CFiredCells(-10), fADADecision(-10), fADCDecision(-10), fIsZNAFired(-10), fIsZNCFired(-10),
  fJpsiPt(0), fJpsiY(0), fJpsiM(0),
  fCMUP6Decision(-10), fCMUP10Decision(-10), fCMUP11Decision(-10),
  fTrgTree(0), fTrgRunNum(0), 
  fCMUP6(-10), fCMUP10(-10), fCMUP11(-10),
  fCentralBarrelTracks(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
// ----------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoFwdJpsiCentMult::AliAnalysisTaskNanoFwdJpsiCentMult(const char* name) : AliAnalysisTaskSE(name),
  fMuonTrackCuts(0x0), fAOD(0), fOutputList(0),fCounterH(0), fNumberMuonsH(0), fNtracksNoJpsiH(0), fNtracksCentralBarrelH(0),
  fAllTracksTree(0), fRunNum(0), fL0inputs(0),
  fZNCEnergy(-999), fZNAEnergy(-999),
  fV0ADecision(-10), fV0CDecision(-10),  fV0AFiredCells(-10), fV0CFiredCells(-10), fADADecision(-10), fADCDecision(-10), fIsZNAFired(-10), fIsZNCFired(-10),
  fJpsiPt(0), fJpsiY(0), fJpsiM(0),
  fCMUP6Decision(-10), fCMUP10Decision(-10), fCMUP11Decision(-10),
  fTrgTree(0), fTrgRunNum(0), 
  fCMUP6(-10), fCMUP10(-10), fCMUP11(-10),
  fCentralBarrelTracks(0)
{
  // constructor
  DefineInput(0, TChain::Class());   
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());   
  DefineOutput(3, TTree::Class());
}
// ----------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoFwdJpsiCentMult::~AliAnalysisTaskNanoFwdJpsiCentMult()
{
  // destructor
  // liberate all allocated memory
  if(fOutputList) {delete fOutputList;}     	
  if(fMuonTrackCuts) {delete fMuonTrackCuts;}
  if(fAllTracksTree) {delete fAllTracksTree;}
  if(fTrgTree) {delete fTrgTree;}
  if(fCounterH) {delete fCounterH;}
  if(fNumberMuonsH) {delete fNumberMuonsH;}
  if(fNtracksNoJpsiH) {delete fNtracksNoJpsiH;}
  if(fNtracksCentralBarrelH) {delete fNtracksCentralBarrelH;}


}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoFwdJpsiCentMult::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)

  ////////////////////////////////////////
  //Muon track cuts
  ////////////////////////////////////////
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchLpt);	
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");


  ////////////////////////////////////////
  //All tracks tree
  ////////////////////////////////////////
  //tracks
  fCentralBarrelTracks = new TClonesArray("AliAODTrack", 1000);

  fAllTracksTree = new TTree("fRecTree", "fRecTree");
  fAllTracksTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fAllTracksTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  fAllTracksTree ->Branch("fZNCEnergy", &fZNCEnergy, "fZNCEnergy/F");
  fAllTracksTree ->Branch("fZNAEnergy", &fZNAEnergy, "fZNAEnergy/F");
  fAllTracksTree ->Branch("fZNATDC", &fZNATDC[0], "fZNATDC[4]/F");
  fAllTracksTree ->Branch("fZNCTDC", &fZNCTDC[0], "fZNCTDC[4]/F");
  fAllTracksTree ->Branch("fV0ADecision", &fV0ADecision, "fV0ADecision/I");
  fAllTracksTree ->Branch("fV0CDecision", &fV0CDecision, "fV0CDecision/I");
  fAllTracksTree ->Branch("fV0AFiredCells", &fV0AFiredCells, "fV0AFiredCells/I");
  fAllTracksTree ->Branch("fV0CFiredCells", &fV0CFiredCells, "fV0CFiredCells/I");
  fAllTracksTree ->Branch("fADADecision", &fADADecision, "fADADecision/I");
  fAllTracksTree ->Branch("fADCDecision", &fADCDecision, "fADCDecision/I");
  fAllTracksTree ->Branch("fIsZNAFired", &fIsZNAFired, "fIsZNAFired/I");
  fAllTracksTree ->Branch("fIsZNCFired", &fIsZNCFired, "fIsZNCFired/I");
  fAllTracksTree ->Branch("fJpsiPt", &fJpsiPt, "fJpsiPt/F");
  fAllTracksTree ->Branch("fJpsiY", &fJpsiY, "fJpsiY/F");
  fAllTracksTree ->Branch("fJpsiM", &fJpsiM, "fJpsiM/F");
  fAllTracksTree ->Branch("fCMUP6Decision", &fCMUP6Decision, "fCMUP6Decision/I");
  fAllTracksTree ->Branch("fCMUP10Decision", &fCMUP10Decision, "fCMUP10Decision/I");
  fAllTracksTree ->Branch("fCMUP11Decision", &fCMUP11Decision, "fCMUP11Decision/I");
  fAllTracksTree ->Branch("fCentralBarrelTracks", &fCentralBarrelTracks);

  // post data
  PostData(1, fAllTracksTree);


  ////////////////////////////////////////
  //Trigger information tree
  ////////////////////////////////////////
  fTrgTree = new TTree("fTrgTree", "fTrgTree");
  fTrgTree ->Branch("fTrgRunNum", &fTrgRunNum, "fTrgRunNum/I");
  fTrgTree ->Branch("fCMUP6", &fCMUP6, "fCMUP6/I");
  fTrgTree ->Branch("fCMUP10", &fCMUP10, "fCMUP10/I");
  fTrgTree ->Branch("fCMUP11", &fCMUP11, "fCMUP11/I");
  // post data
  PostData(3, fTrgTree);

  ////////////////////////////////////////
  //output histograms
  ////////////////////////////////////////
  fOutputList = new TList();          // this is a list which will contain all  histograms
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  //  counter for events passing each cut    
  fCounterH = new TH1F("fCounterH", "fCounterH", 25, 0., 25.);
  fOutputList->Add(fCounterH);
  // number of positive and negative muons passing the muon selection
  fNumberMuonsH = new TH2F("fNumberMuonsH", "fNumberMuonsH", 12, 0., 12.,12, 0., 12.); 
  fOutputList->Add(fNumberMuonsH);        // don't forget to add it to the list!
  // number of all other tracks except muons from Jpsi in MUON arm
  fNtracksNoJpsiH = new TH1F("fNtracksNoJpsiH", "fNtracksNoJpsiH", 100, 0., 100.);
  fOutputList->Add(fNtracksNoJpsiH);
  // number of all tracks in central barrel if Jpsi in MUON arm
  fNtracksCentralBarrelH = new TH1F("fNtracksCentralBarrelH", "fNtracksCentralBarrelH", 100, 0., 100.);
  fOutputList->Add(fNtracksCentralBarrelH);
  // post data
  PostData(2, fOutputList);
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoFwdJpsiCentMult::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoFwdJpsiCentMult::PostAllData()
{
  // Post data
  PostData(1, fAllTracksTree);
  PostData(2, fOutputList);
  PostData(3, fTrgTree);
}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoFwdJpsiCentMult::IfJpsiStoreAllTracks(Int_t *idxPosMuons, Int_t *idxNegMuons){

  // Get muon masss fromn PDG.
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  TParticlePDG *partMuon = pdgdat->GetParticle(13);
  Double_t MuonMass = partMuon->Mass();

  // Create all four vectors.
  // --  positive muon
  TLorentzVector PosMuon1;
  AliAODTrack *PosTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(idxPosMuons[0]));
  PosMuon1.SetPtEtaPhiM(PosTrack->Pt(), PosTrack->Eta(), PosTrack->Phi(), MuonMass);
  // --  negative muon
  TLorentzVector NegMuon1;
  AliAODTrack *NegTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(idxNegMuons[0]));
  NegMuon1.SetPtEtaPhiM(NegTrack->Pt(), NegTrack->Eta(), NegTrack->Phi(), MuonMass);

  // Create dimuon.
  TLorentzVector MuMu = NegMuon1+PosMuon1;

  // Is it Jpsi?
  if(MuMu.M() > 3.0 && 3.2 < MuMu.M()) return;

  // Is it in MUON arm rapidity?
  if(MuMu.Rapidity() > -2.5 && -4.0 < MuMu.Rapidity()) return;

  // Set tree variables.
  fJpsiPt = MuMu.Pt();
  fJpsiY = MuMu.Rapidity();
  fJpsiM = MuMu.M();

  // Loop over all tracks and save it.
  Int_t nCentralBarrelSelectedTracks(0);
  Int_t nTracks(fAOD->GetNumberOfTracks());
  if(nTracks<1) return;

  fCentralBarrelTracks->Clear("C");
  for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    // skip MUON arm muons
    if( (idxPosMuons[0] == iTrack) || (idxNegMuons[0] == iTrack) ) continue;
    // get track
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) {
      Printf("Track does not exist!");
      return;
    }
    // store only central barrel tracks
    if(TMath::Abs(track->Eta()) > 0.9) continue;
    new((*fCentralBarrelTracks)[nCentralBarrelSelectedTracks]) AliAODTrack(*track);
    nCentralBarrelSelectedTracks++;
  }

  fNtracksNoJpsiH->Fill(nTracks-2);
  fNtracksCentralBarrelH->Fill(nCentralBarrelSelectedTracks);

}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoFwdJpsiCentMult::UserExec(Option_t *)
{
  Int_t iSelectionCounter = 0; // no selection applied yet 
  fCounterH->Fill(iSelectionCounter); // entering UserExec 1/1 (data/MC)
  iSelectionCounter++;

  ////////////////////////////////////////////
  // Geting the AOD event
  ////////////////////////////////////////////
   // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
  if(!fAOD) {
    PostAllData();
    return;
  }                                  
  fCounterH->Fill(iSelectionCounter); // AOD event found 2/2
  iSelectionCounter++;

  ////////////////////////////////////////////
  //  Trigger information
  ////////////////////////////////////////////
  // in 2018 q,r : CMUP6-B-NOPF-MUFAST = *0VBA 0MUL ,  
  // in 2018 q,r and 2015 o:  CMUP11-B-NOPF-MUFAST = *0VBA *0UBA *0UBC 0MUL,
  // in 2015 o : CMUP10-B-NOPF-MUFAST = = *0VBA *0UBA *0UBC 0MSL , 
  TString trigger = fAOD->GetFiredTriggerClasses();
  
  Bool_t isTriggered = kFALSE;

  if (trigger.Contains("CMUP11-B-NOPF-MUFAST")) {
    isTriggered = kTRUE;
    fCMUP11Decision = 1;
    fCMUP11 = 1;
  } else {
    fCMUP11Decision = 0;
    fCMUP11 = 0;
  }
  if (trigger.Contains("CMUP10-B-NOPF-MUFAST")) {
    isTriggered = kTRUE;
    fCMUP10Decision = 1;
    fCMUP10 = 1;
  } else {
    fCMUP10Decision = 0;
    fCMUP10 = 0;
  }
  if (trigger.Contains("CMUP6-B-NOPF-MUFAST")) {
    isTriggered = kTRUE;
    fCMUP6Decision = 1;
    fCMUP6 = 1;
  } else {
    fCMUP6Decision = 0;
    fCMUP6 = 0;
  }

  if (!isTriggered) {
    PostAllData();
    return;
  }

  fTrgRunNum = fAOD->GetRunNumber();
  // Fill the trigger tree
  fTrgTree->Fill();

  fCounterH->Fill(iSelectionCounter); // right trigger found 4/7
  iSelectionCounter++;

  // get the run number and trigger inputs
  fRunNum = fAOD->GetRunNumber();
  fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  
  ////////////////////////////////////////////
  //  find muons
  ////////////////////////////////////////////
  //are there tracks at all?
  Int_t nTracks(fAOD->GetNumberOfTracks()); 
  if(nTracks<1) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // At least one track 5/8
  iSelectionCounter++;

  // loop over tracks and select good muons
  Int_t nGoodPosMuons = 0;
  Int_t nGoodNegMuons = 0;  
  Int_t *idxPosMuons = new Int_t[nTracks];
  Int_t *idxNegMuons = new Int_t[nTracks];  
  for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    // get track
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack)); 
    if(!track) {
      Printf("Track does not exist!");
      return;
    }

    // is it a good muon track?
    if(!track->IsMuonTrack()) continue;
    if(!fMuonTrackCuts->IsSelected(track)) continue;
    if( (track->GetRAtAbsorberEnd() < 17.5) || (track->GetRAtAbsorberEnd() > 89.5) ) continue;

    // increase counter and store indices
    if(track->Charge() > 0) {
      idxPosMuons[nGoodPosMuons] = iTrack;
      nGoodPosMuons++;
    }
    else if(track->Charge() < 0) {
      idxNegMuons[nGoodNegMuons] = iTrack;
      nGoodNegMuons++;
    } 
  }
  // store number of muons
  fNumberMuonsH->Fill(nGoodPosMuons,nGoodNegMuons);

  ////////////////////////////////////////////
  // two muon analysis
  ////////////////////////////////////////////
  if (!(nGoodPosMuons == 1 && nGoodNegMuons == 1)) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // exactly one positive and one negative muons 6/9
  iSelectionCounter++;
  IfJpsiStoreAllTracks(idxPosMuons,idxNegMuons);

  ////////////////////////////////////////////
  // info to determine exclusivity
  ////////////////////////////////////////////

  // ---ZDC 
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); // ZDC info is present 7/10
  iSelectionCounter++;

  fZNAEnergy = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy = dataZDC->GetZNCTowerEnergy()[0];
  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);

  // at least one ZDC hit in the timing window
  fIsZNAFired = 0;
  fIsZNCFired = 0;
  for (Int_t i=0;i<4;i++){
    if ( (fZNATDC[i]>-2.) && (fZNATDC[i]<2.) ) fIsZNAFired = 1;
    if ( (fZNCTDC[i]>-2.) && (fZNCTDC[i]<2.) ) fIsZNCFired = 1;  
  }

  // ---V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); //  V0 info 8/11
  iSelectionCounter++;

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  Int_t nV0CFiredCells = 0;
  Int_t nV0AFiredCells = 0;

  for(Int_t i = 0; i < 64; i++) {
    if(dataVZERO->GetBBFlag(i) == kTRUE) {
      if(i < 32) {
        nV0CFiredCells += 1;
      } else {
        nV0AFiredCells += 1;
      }
    }
  }

  fV0CFiredCells = nV0CFiredCells;
  fV0AFiredCells = nV0AFiredCells;

  // ---AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(!dataAD){
    PostAllData();
    return;
  } 
  fCounterH->Fill(iSelectionCounter); //  AD info 9/12
  iSelectionCounter++;

  fADADecision = dataAD->GetADADecision();
  fADCDecision = dataAD->GetADCDecision();


  // Fill the reconstruction tree
  fAllTracksTree->Fill();

  // post the data
  PostAllData();

  // clean up
  delete [] idxPosMuons;
  delete [] idxNegMuons;

}
// ----------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNanoFwdJpsiCentMult::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
// ----------------------------------------------------------------------------------------------------------------------------------


