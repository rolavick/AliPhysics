/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskNanoFwdJpsiCentMult_H
#define AliAnalysisTaskNanoFwdJpsiCentMult_H

#include <AliAnalysisTaskSE.h>

class AliMuonTrackCuts; 	// Include class for standard muon tack cuts
class TClonesArray;


class AliAnalysisTaskNanoFwdJpsiCentMult : public AliAnalysisTaskSE
{
public:
                            AliAnalysisTaskNanoFwdJpsiCentMult();
                            AliAnalysisTaskNanoFwdJpsiCentMult(const char *name);
    virtual                 ~AliAnalysisTaskNanoFwdJpsiCentMult();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    virtual void   			NotifyRun();								  // Implement the Notify run to search for the new parameters at each new runs
  void          IfJpsiStoreAllTracks(Int_t *pos, Int_t *neg); // Analyse two muons and if they have Jpsi mass, store all tracks in the event
	void 					PostAllData();	

    AliMuonTrackCuts* 		fMuonTrackCuts; 					// Use the class as a data member

private:

    AliAODEvent*            fAOD;       		//! input event

    TList*                  fOutputList; 		//! output list
    TH1F*                   fCounterH; 			//! counter for events passing each cut	
    TH2F*                   fNumberMuonsH; 	//! count good muons per event
    TH1F*                   fNtracksNoJpsiH; //! count additional tracks to Jpsi muons
    TH1F*                   fNtracksCentralBarrelH; //! count additional tracks to Jpsi muons

	TTree *fAllTracksTree;
	Int_t fRunNum;
	UInt_t fL0inputs;
	Float_t fZNCEnergy;
	Float_t fZNAEnergy;
	Float_t fZNATDC[4];
	Float_t fZNCTDC[4];
	Int_t fV0ADecision;
	Int_t fV0CDecision;
	Int_t fV0AFiredCells;
	Int_t fV0CFiredCells;
	Int_t fADADecision;
	Int_t fADCDecision;
	Int_t fIsZNAFired;
	Int_t fIsZNCFired;
	Float_t fJpsiPt;
	Float_t fJpsiY;
	Float_t fJpsiM;
	Int_t fCMUP6Decision;
	Int_t fCMUP10Decision;
	Int_t fCMUP11Decision;

	TTree *fTrgTree;
	Int_t fTrgRunNum;
	Int_t fCMUP6;
	Int_t fCMUP10;
	Int_t fCMUP11;

  TClonesArray *fCentralBarrelTracks;


    AliAnalysisTaskNanoFwdJpsiCentMult(const AliAnalysisTaskNanoFwdJpsiCentMult&); // not implemented
    AliAnalysisTaskNanoFwdJpsiCentMult& operator=(const AliAnalysisTaskNanoFwdJpsiCentMult&); // not implemented

    ClassDef(AliAnalysisTaskNanoFwdJpsiCentMult, 1);
};

#endif
