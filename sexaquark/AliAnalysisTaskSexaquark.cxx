#include "AliAnalysisTaskSexaquark.h"

ClassImp(AliAnalysisTaskSexaquark);

/*
 * Empty I/O constructor. Non-persistent members are initialized to their default values from here.
 */
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark()
    : AliAnalysisTaskSE(),
      /*  */
      fIsMC(0),
      fSourceOfV0s(),
      fReactionID(0),
      fDoQA(0),
      fStopAfter(),
      fReadSignalLogs(0),
      /*  */
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPrimaryVertex(0),
      fPIDResponse(0),
      fEventCuts(),
      fMagneticField(0.),
      /*  */
      fRunNumber(0),
      fDirNumber(0),
      fEventNumber(0),
      fCentrality(0.),
      /*  */
      fPDG(),
      fList_QA_Hists(0),
      fList_Tracks_Hists(0),
      fList_V0s_Hists(0),
      fList_Sexaquarks_Hists(0),
      fList_PosKaonPairs_Hists(0),
      /*  */
      fAliEnPath(),
      fIsFirstEvent(0),
      /*  */
      fProductsPDG(),
      fFinalStateProductsPDG(),
      fStruckNucleonPDG(0),
      getNegPdgCode_fromV0PdgCode(),
      getPosPdgCode_fromV0PdgCode(),
      /*  */
      fTree_Events(0),
      fTree_Injected(0),
      fTree_MC(0),
      fTree_Tracks(0),
      fTree_V0s(0),
      fTree_KaonPairs(0),
      fTree_Sexaquarks(0),
      /*  */
      tEvents_PV_TrueXv(0),
      tEvents_PV_TrueYv(0),
      tEvents_PV_TrueZv(0),
      tEvents_PV_RecXv(0),
      tEvents_PV_RecYv(0),
      tEvents_PV_RecZv(0),
      tEvents_NTracks(0),
      tEvents_IsMB(0),
      tEvents_IsHighMultV0(0),
      tEvents_IsHighMultSPD(0),
      tEvents_IsCentral(0),
      tEvents_IsSemiCentral(0),
      /*  */
      tInjected_RunNumber(0),
      tInjected_DirNumber(0),
      tInjected_EventNumber(0),
      tInjected_ReactionID(0),
      tInjected_Px(0),
      tInjected_Py(0),
      tInjected_Pz(0),
      tInjected_Mass(0),
      tInjected_Nucleon_PdgCode(0),
      tInjected_Nucleon_Px(0),
      tInjected_Nucleon_Py(0),
      tInjected_Nucleon_Pz(0),
      tInjected_ReactionChannel(0),
      /*  */
      tMC_Idx(0),
      tMC_PdgCode(0),
      tMC_Idx_Mother(0),
      tMC_NDaughters(0),
      tMC_Idx_FirstDau(0),
      tMC_Idx_LastDau(0),
      tMC_ReactionID(0),
      tMC_Px(0),
      tMC_Py(0),
      tMC_Pz(0),
      tMC_Xv_i(0),
      tMC_Yv_i(0),
      tMC_Zv_i(0),
      tMC_Status(0),
      tMC_IsSecondary(0),
      tMC_IsSignal(0),
      /*  */
      tTrack_Idx(0),
      tTrack_Idx_True(0),
      tTrack_Px(0),
      tTrack_Py(0),
      tTrack_Pz(0),
      tTrack_Charge(0),
      tTrack_NSigmaPion(0),
      tTrack_NSigmaKaon(0),
      tTrack_NSigmaProton(0),
      tTrack_IsTrue(0),
      tTrack_IsSecondary(0),
      tTrack_IsSignal(0),
      tTrack_ReactionID(0),
      /*  */
      tV0_Idx(0),
      tV0_Idx_Pos(0),
      tV0_Idx_Neg(0),
      tV0_Idx_True(0),
      tV0_PID(0),
      tV0_Px(0),
      tV0_Py(0),
      tV0_Pz(0),
      tV0_E(0),
      tV0_Xv_f(0),
      tV0_Yv_f(0),
      tV0_Zv_f(0),
      tV0_Neg_Px(0),
      tV0_Neg_Py(0),
      tV0_Neg_Pz(0),
      tV0_Pos_Px(0),
      tV0_Pos_Py(0),
      tV0_Pos_Pz(0),
      tV0_IsTrue(0),
      tV0_IsSecondary(0),
      tV0_IsSignal(0),
      tV0_ReactionID(0),
      /*  */
      tSexaquark_Idx_AL(0),
      tSexaquark_Idx_AL_Neg(0),
      tSexaquark_Idx_AL_Pos(0),
      tSexaquark_Px(0),
      tSexaquark_Py(0),
      tSexaquark_Pz(0),
      tSexaquark_E(0),
      tSexaquark_SV_Xv(0),
      tSexaquark_SV_Yv(0),
      tSexaquark_SV_Zv(0),
      tSexaquark_DecayLength(0),
      tSexaquark_CPAwrtPV(0),
      tSexaquark_DCAwrtPV(0),
      tSexaquark_Chi2ndf(0),
      tSexaquark_IsHybrid(0),
      tSexaquark_IsSignal(0),
      tSexaquark_ReactionID(0),
      /*  */
      tSexaquarkA_Idx_K0S(0),
      tSexaquarkA_Idx_K0S_Neg(0),
      tSexaquarkA_Idx_K0S_Pos(0),
      tSexaquarkA_DCAbtwV0s(0),
      tSexaquarkA_DCAV0aSV(0),
      tSexaquarkA_DCAV0bSV(0),
      tSexaquarkA_DCAV0anegSV(0),
      tSexaquarkA_DCAV0aposSV(0),
      tSexaquarkA_DCAV0bnegSV(0),
      tSexaquarkA_DCAV0bposSV(0),
      /*  */
      tSexaquarkD_Idx_PosKaon(0),
      tSexaquarkD_DCAV0Ba(0),
      tSexaquarkD_DCAV0negBa(0),
      tSexaquarkD_DCAV0posBa(0),
      tSexaquarkD_DCAV0SV(0),
      tSexaquarkD_DCABaSV(0),
      tSexaquarkD_DCAV0negSV(0),
      tSexaquarkD_DCAV0posSV(0),
      /*  */
      tSexaquarkE_Idx_PP(0),
      tSexaquarkE_Idx_PP_Neg(0),
      tSexaquarkE_Idx_PP_Pos(0),
      tSexaquarkE_Idx_PosKaon(0),
      tSexaquarkE_DCAV0SV(0),
      tSexaquarkE_DCAV0negSV(0),
      tSexaquarkE_DCAV0posSV(0),
      tSexaquarkE_DCApkSV(0),
      tSexaquarkE_DCApmSV(0),
      tSexaquarkE_DCAppSV(0),
      tSexaquarkE_DCApkV0(0),
      tSexaquarkE_DCApmV0(0),
      tSexaquarkE_DCAppV0(0),
      tSexaquarkE_DCApmPK(0),
      tSexaquarkE_DCAppPK(0),
      /*  */
      tKaonPair_Idx(0),
      tKaonPair_Idx_KaonA(0),
      tKaonPair_Idx_KaonB(0),
      tKaonPair_Px(0),
      tKaonPair_Py(0),
      tKaonPair_Pz(0),
      tKaonPair_E(0),
      tKaonPair_SV_Xv(0),
      tKaonPair_SV_Yv(0),
      tKaonPair_SV_Zv(0),
      tKaonPair_DecayLength(0),
      tKaonPair_DCAwrtPV(0),
      tKaonPair_DCAbtwKaons(0),
      tKaonPair_DCAkaSV(0),
      tKaonPair_DCAkbSV(0),
      tKaonPair_OpeningAngle(0),
      tKaonPair_Chi2ndf(0),
      tKaonPair_IsHybrid(0),
      tKaonPair_IsSignal(0),
      tKaonPair_ReactionID(0),
      /*  */
      fHist_QA(),
      f2DHist_QA(),
      fHist_Tracks_Bookkeep(),
      fHist_Tracks(),
      f2DHist_Tracks_TPCsignal(),
      fHist_V0s_Bookkeep(),
      fHist_V0s(),
      f2DHist_V0s_ArmenterosPodolanski(),
      fHist_AntiSexaquarks_Bookkeep(),
      fHist_AntiSexaquarks(),
      fHist_PosKaonPairs_Bookkeep(),
      fHist_PosKaonPairs(),
      /*  */
      getPdgCode_fromMcIdx(),
      isMcIdxSignal(),
      isMcIdxSecondary(),
      getReactionIdx_fromMcIdx(),
      getMcIdx_fromReactionIdx(),
      getPt_fromReactionIdx(),
      getRadius_fromReactionIdx(),
      doesMcIdxHaveMother(),
      getMotherMcIdx_fromMcIdx(),
      getNegDauMcIdx_fromMcIdx(),
      getPosDauMcIdx_fromMcIdx(),
      mcIndicesOfTrueV0s(),
      /*  */
      getMcIdx_fromEsdIdx(),
      esdIdxPassesTrackSelection_AsProton(),
      esdIdxPassesTrackSelection_AsAntiProton(),
      esdIdxPassesTrackSelection_AsPosKaon(),
      esdIdxPassesTrackSelection_AsNegKaon(),
      esdIdxPassesTrackSelection_AsPiPlus(),
      esdIdxPassesTrackSelection_AsPiMinus(),
      getEsdIdx_fromMcIdx(),
      isEsdIdxSignal(),
      getReactionIdx_fromEsdIdx(),
      getEsdIdx_fromReactionIdx(),
      doesEsdIdxHaveMother(),
      getMotherMcIdx_fromEsdIdx(),
      getNegDauEsdIdx_fromMcIdx(),
      getPosDauEsdIdx_fromMcIdx(),
      esdIndicesOfAntiProtonTracks(),
      esdIndicesOfPosKaonTracks(),
      esdIndicesOfPiMinusTracks(),
      esdIndicesOfPiPlusTracks(),
      /*  */
      getEsdIdxOfNegDau_fromAntiLambdaIdx(),
      getEsdIdxOfPosDau_fromAntiLambdaIdx(),
      getEsdIdxOfNegDau_fromKaonZeroShortIdx(),
      getEsdIdxOfPosDau_fromKaonZeroShortIdx(),
      getEsdIdxOfNegDau_fromPionPairIdx(),
      getEsdIdxOfPosDau_fromPionPairIdx(),
      esdAntiLambdas(),
      esdKaonsZeroShort(),
      kfAntiLambdas(),
      kfKaonsZeroShort(),
      kfPionPairs(),
      /*  */
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0),
      kTurnedOn_Track_RejectKinks(0),
      kMin_Track_DCA_wrtPV(0.),
      kMin_Track_DCAxy_wrtPV(0.),
      kMin_Track_DCAz_wrtPV(0.),
      kMin_Track_Pt(),
      /*  */
      kMin_V0_Mass(),
      kMax_V0_Mass(),
      kMax_V0_Eta(),
      kMax_V0_ArmPtOverAlpha(),
      kMin_V0_Pt(),
      kMin_V0_Radius(),
      kMin_V0_DecayLength(),
      kMax_V0_DecayLength(),
      kMin_V0_CPAwrtPV(),
      kMax_V0_CPAwrtPV(),
      kMin_V0_DCAwrtPV(),
      kMax_V0_DCAbtwDau(),
      kMax_V0_DCAnegV0(),
      kMax_V0_DCAposV0(),
      kMax_V0_ImprvDCAbtwDau(),
      kMax_V0_Chi2(),
      kMax_V0_Chi2ndf(),
      /*  */
      kMin_Sexa_Mass(0.),
      kMax_Sexa_Mass(0.),
      kMin_Sexa_Pt(0.),
      kMax_Sexa_Eta(0.),
      kMax_Sexa_Rapidity(0.),
      kMin_Sexa_DecayLength(0.),
      kMin_Sexa_Radius(0.),
      kMin_Sexa_CPAwrtPV(0.),
      kMax_Sexa_CPAwrtPV(0.),
      kMin_Sexa_DCAwrtPV(0.),
      kMax_Sexa_DCAwrtPV(0.),
      kMax_Sexa_Chi2ndf(0.),
      kMax_Sexa_DCAbtwV0s(0.),
      kMax_Sexa_DCAv0aSV(0.),
      kMax_Sexa_DCAv0bSV(0.),
      kMax_Sexa_DCAv0anegSV(0.),
      kMax_Sexa_DCAv0aposSV(0.),
      kMax_Sexa_DCAv0bnegSV(0.),
      kMax_Sexa_DCAv0bposSV(0.),
      kMax_Sexa_DCAbaSV(0.),
      kMax_Sexa_DCAv0ba(0.),
      kMax_Sexa_DCAv0negba(0.),
      kMax_Sexa_DCAv0posba(0.),
      kMax_Sexa_DCAv0SV(0.),
      kMax_Sexa_DCAv0negSV(0.),
      kMax_Sexa_DCAv0posSV(0.),
      kMax_Sexa_DCApkSV(0.),
      kMax_Sexa_DCApmSV(0.),
      kMax_Sexa_DCAppSV(0.),
      kMax_Sexa_DCAv0pk(0.),
      kMax_Sexa_DCAv0pm(0.),
      kMax_Sexa_DCAv0pp(0.),
      kMax_Sexa_DCApkpm(0.),
      kMax_Sexa_DCApkpp(0.),
      /*  */
      kMin_PosKaonPair_Radius(0.),
      kMin_PosKaonPair_Pt(0.),
      kMin_PosKaonPair_DCAwrtPV(0.),
      kMax_PosKaonPair_DCAbtwKaons(0.),
      kMax_PosKaonPair_DCApkaSV(0.),
      kMax_PosKaonPair_DCApkbSV(0.),
      kMax_PosKaonPair_Chi2ndf(0.) {}

/*
 * Constructor, called locally.
 */
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark(const char* name)
    : AliAnalysisTaskSE(name),
      /*  */
      fIsMC(0),
      fSourceOfV0s(),
      fReactionID(0),
      fDoQA(0),
      fStopAfter(),
      fReadSignalLogs(0),
      /*  */
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPrimaryVertex(0),
      fPIDResponse(0),
      fEventCuts(),
      fMagneticField(0.),
      fRunNumber(0),
      fDirNumber(0),
      fEventNumber(0),
      fCentrality(0.),
      /*  */
      fPDG(),
      fList_QA_Hists(0),
      fList_Tracks_Hists(0),
      fList_V0s_Hists(0),
      fList_Sexaquarks_Hists(0),
      fList_PosKaonPairs_Hists(0),
      /*  */
      fAliEnPath(),
      fIsFirstEvent(0),
      /*  */
      fProductsPDG(),
      fFinalStateProductsPDG(),
      fStruckNucleonPDG(0),
      getNegPdgCode_fromV0PdgCode(),
      getPosPdgCode_fromV0PdgCode(),
      /*  */
      fTree_Events(0),
      fTree_Injected(0),
      fTree_MC(0),
      fTree_Tracks(0),
      fTree_V0s(0),
      fTree_KaonPairs(0),
      fTree_Sexaquarks(0),
      /*  */
      tEvents_PV_TrueXv(0),
      tEvents_PV_TrueYv(0),
      tEvents_PV_TrueZv(0),
      tEvents_PV_RecXv(0),
      tEvents_PV_RecYv(0),
      tEvents_PV_RecZv(0),
      tEvents_NTracks(0),
      tEvents_IsMB(0),
      tEvents_IsHighMultV0(0),
      tEvents_IsHighMultSPD(0),
      tEvents_IsCentral(0),
      tEvents_IsSemiCentral(0),
      /*  */
      tInjected_RunNumber(0),
      tInjected_DirNumber(0),
      tInjected_EventNumber(0),
      tInjected_ReactionID(0),
      tInjected_Px(0),
      tInjected_Py(0),
      tInjected_Pz(0),
      tInjected_Mass(0),
      tInjected_Nucleon_PdgCode(0),
      tInjected_Nucleon_Px(0),
      tInjected_Nucleon_Py(0),
      tInjected_Nucleon_Pz(0),
      tInjected_ReactionChannel(0),
      /*  */
      tMC_Idx(0),
      tMC_PdgCode(0),
      tMC_Idx_Mother(0),
      tMC_NDaughters(0),
      tMC_Idx_FirstDau(0),
      tMC_Idx_LastDau(0),
      tMC_ReactionID(0),
      tMC_Px(0),
      tMC_Py(0),
      tMC_Pz(0),
      tMC_Xv_i(0),
      tMC_Yv_i(0),
      tMC_Zv_i(0),
      tMC_Status(0),
      tMC_IsSecondary(0),
      tMC_IsSignal(0),
      /*  */
      tTrack_Idx(0),
      tTrack_Idx_True(0),
      tTrack_Px(0),
      tTrack_Py(0),
      tTrack_Pz(0),
      tTrack_Charge(0),
      tTrack_NSigmaPion(0),
      tTrack_NSigmaKaon(0),
      tTrack_NSigmaProton(0),
      tTrack_IsTrue(0),
      tTrack_IsSecondary(0),
      tTrack_IsSignal(0),
      tTrack_ReactionID(0),
      /*  */
      tV0_Idx(0),
      tV0_Idx_Pos(0),
      tV0_Idx_Neg(0),
      tV0_Idx_True(0),
      tV0_PID(0),
      tV0_Px(0),
      tV0_Py(0),
      tV0_Pz(0),
      tV0_E(0),
      tV0_Xv_f(0),
      tV0_Yv_f(0),
      tV0_Zv_f(0),
      tV0_Neg_Px(0),
      tV0_Neg_Py(0),
      tV0_Neg_Pz(0),
      tV0_Pos_Px(0),
      tV0_Pos_Py(0),
      tV0_Pos_Pz(0),
      tV0_IsTrue(0),
      tV0_IsSecondary(0),
      tV0_IsSignal(0),
      tV0_ReactionID(0),
      /*  */
      tSexaquark_Idx_AL(0),
      tSexaquark_Idx_AL_Neg(0),
      tSexaquark_Idx_AL_Pos(0),
      tSexaquark_Px(0),
      tSexaquark_Py(0),
      tSexaquark_Pz(0),
      tSexaquark_E(0),
      tSexaquark_SV_Xv(0),
      tSexaquark_SV_Yv(0),
      tSexaquark_SV_Zv(0),
      tSexaquark_DecayLength(0),
      tSexaquark_CPAwrtPV(0),
      tSexaquark_DCAwrtPV(0),
      tSexaquark_Chi2ndf(0),
      tSexaquark_IsHybrid(0),
      tSexaquark_IsSignal(0),
      tSexaquark_ReactionID(0),
      /*  */
      tSexaquarkA_Idx_K0S(0),
      tSexaquarkA_Idx_K0S_Neg(0),
      tSexaquarkA_Idx_K0S_Pos(0),
      tSexaquarkA_DCAbtwV0s(0),
      tSexaquarkA_DCAV0aSV(0),
      tSexaquarkA_DCAV0bSV(0),
      tSexaquarkA_DCAV0anegSV(0),
      tSexaquarkA_DCAV0aposSV(0),
      tSexaquarkA_DCAV0bnegSV(0),
      tSexaquarkA_DCAV0bposSV(0),
      /*  */
      tSexaquarkD_Idx_PosKaon(0),
      tSexaquarkD_DCAV0Ba(0),
      tSexaquarkD_DCAV0negBa(0),
      tSexaquarkD_DCAV0posBa(0),
      tSexaquarkD_DCAV0SV(0),
      tSexaquarkD_DCABaSV(0),
      tSexaquarkD_DCAV0negSV(0),
      tSexaquarkD_DCAV0posSV(0),
      /*  */
      tSexaquarkE_Idx_PP(0),
      tSexaquarkE_Idx_PP_Neg(0),
      tSexaquarkE_Idx_PP_Pos(0),
      tSexaquarkE_Idx_PosKaon(0),
      tSexaquarkE_DCAV0SV(0),
      tSexaquarkE_DCAV0negSV(0),
      tSexaquarkE_DCAV0posSV(0),
      tSexaquarkE_DCApkSV(0),
      tSexaquarkE_DCApmSV(0),
      tSexaquarkE_DCAppSV(0),
      tSexaquarkE_DCApkV0(0),
      tSexaquarkE_DCApmV0(0),
      tSexaquarkE_DCAppV0(0),
      tSexaquarkE_DCApmPK(0),
      tSexaquarkE_DCAppPK(0),
      /*  */
      tKaonPair_Idx(0),
      tKaonPair_Idx_KaonA(0),
      tKaonPair_Idx_KaonB(0),
      tKaonPair_Px(0),
      tKaonPair_Py(0),
      tKaonPair_Pz(0),
      tKaonPair_E(0),
      tKaonPair_SV_Xv(0),
      tKaonPair_SV_Yv(0),
      tKaonPair_SV_Zv(0),
      tKaonPair_DecayLength(0),
      tKaonPair_DCAwrtPV(0),
      tKaonPair_DCAbtwKaons(0),
      tKaonPair_DCAkaSV(0),
      tKaonPair_DCAkbSV(0),
      tKaonPair_OpeningAngle(0),
      tKaonPair_Chi2ndf(0),
      tKaonPair_IsHybrid(0),
      tKaonPair_IsSignal(0),
      tKaonPair_ReactionID(0),
      /*  */
      fHist_QA(),
      f2DHist_QA(),
      fHist_Tracks_Bookkeep(),
      fHist_Tracks(),
      f2DHist_Tracks_TPCsignal(),
      fHist_V0s_Bookkeep(),
      fHist_V0s(),
      f2DHist_V0s_ArmenterosPodolanski(),
      fHist_AntiSexaquarks_Bookkeep(),
      fHist_AntiSexaquarks(),
      fHist_PosKaonPairs_Bookkeep(),
      fHist_PosKaonPairs(),
      /*  */
      getPdgCode_fromMcIdx(),
      isMcIdxSignal(),
      isMcIdxSecondary(),
      getReactionIdx_fromMcIdx(),
      getMcIdx_fromReactionIdx(),
      getPt_fromReactionIdx(),
      getRadius_fromReactionIdx(),
      doesMcIdxHaveMother(),
      getMotherMcIdx_fromMcIdx(),
      getNegDauMcIdx_fromMcIdx(),
      getPosDauMcIdx_fromMcIdx(),
      mcIndicesOfTrueV0s(),
      /*  */
      getMcIdx_fromEsdIdx(),
      esdIdxPassesTrackSelection_AsProton(),
      esdIdxPassesTrackSelection_AsAntiProton(),
      esdIdxPassesTrackSelection_AsPosKaon(),
      esdIdxPassesTrackSelection_AsNegKaon(),
      esdIdxPassesTrackSelection_AsPiPlus(),
      esdIdxPassesTrackSelection_AsPiMinus(),
      getEsdIdx_fromMcIdx(),
      isEsdIdxSignal(),
      getReactionIdx_fromEsdIdx(),
      getEsdIdx_fromReactionIdx(),
      doesEsdIdxHaveMother(),
      getMotherMcIdx_fromEsdIdx(),
      getNegDauEsdIdx_fromMcIdx(),
      getPosDauEsdIdx_fromMcIdx(),
      esdIndicesOfAntiProtonTracks(),
      esdIndicesOfPosKaonTracks(),
      esdIndicesOfPiMinusTracks(),
      esdIndicesOfPiPlusTracks(),
      /*  */
      getEsdIdxOfNegDau_fromAntiLambdaIdx(),
      getEsdIdxOfPosDau_fromAntiLambdaIdx(),
      getEsdIdxOfNegDau_fromKaonZeroShortIdx(),
      getEsdIdxOfPosDau_fromKaonZeroShortIdx(),
      getEsdIdxOfNegDau_fromPionPairIdx(),
      getEsdIdxOfPosDau_fromPionPairIdx(),
      esdAntiLambdas(),
      esdKaonsZeroShort(),
      kfAntiLambdas(),
      kfKaonsZeroShort(),
      kfPionPairs(),
      /*  */
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0),
      kTurnedOn_Track_RejectKinks(0),
      kMin_Track_DCA_wrtPV(0.),
      kMin_Track_DCAxy_wrtPV(0.),
      kMin_Track_DCAz_wrtPV(0.),
      kMin_Track_Pt(),
      /*  */
      kMin_V0_Mass(),
      kMax_V0_Mass(),
      kMax_V0_Eta(),
      kMax_V0_ArmPtOverAlpha(),
      kMin_V0_Pt(),
      kMin_V0_Radius(),
      kMin_V0_DecayLength(),
      kMax_V0_DecayLength(),
      kMin_V0_CPAwrtPV(),
      kMax_V0_CPAwrtPV(),
      kMin_V0_DCAwrtPV(),
      kMax_V0_DCAbtwDau(),
      kMax_V0_DCAnegV0(),
      kMax_V0_DCAposV0(),
      kMax_V0_ImprvDCAbtwDau(),
      kMax_V0_Chi2(),
      kMax_V0_Chi2ndf(),
      /*  */
      kMin_Sexa_Mass(0.),
      kMax_Sexa_Mass(0.),
      kMin_Sexa_Pt(0.),
      kMax_Sexa_Eta(0.),
      kMax_Sexa_Rapidity(0.),
      kMin_Sexa_DecayLength(0.),
      kMin_Sexa_Radius(0.),
      kMin_Sexa_CPAwrtPV(0.),
      kMax_Sexa_CPAwrtPV(0.),
      kMin_Sexa_DCAwrtPV(0.),
      kMax_Sexa_DCAwrtPV(0.),
      kMax_Sexa_Chi2ndf(0.),
      kMax_Sexa_DCAbtwV0s(0.),
      kMax_Sexa_DCAv0aSV(0.),
      kMax_Sexa_DCAv0bSV(0.),
      kMax_Sexa_DCAv0anegSV(0.),
      kMax_Sexa_DCAv0aposSV(0.),
      kMax_Sexa_DCAv0bnegSV(0.),
      kMax_Sexa_DCAv0bposSV(0.),
      kMax_Sexa_DCAbaSV(0.),
      kMax_Sexa_DCAv0ba(0.),
      kMax_Sexa_DCAv0negba(0.),
      kMax_Sexa_DCAv0posba(0.),
      kMax_Sexa_DCAv0SV(0.),
      kMax_Sexa_DCAv0negSV(0.),
      kMax_Sexa_DCAv0posSV(0.),
      kMax_Sexa_DCApkSV(0.),
      kMax_Sexa_DCApmSV(0.),
      kMax_Sexa_DCAppSV(0.),
      kMax_Sexa_DCAv0pk(0.),
      kMax_Sexa_DCAv0pm(0.),
      kMax_Sexa_DCAv0pp(0.),
      kMax_Sexa_DCApkpm(0.),
      kMax_Sexa_DCApkpp(0.),
      /*  */
      kMin_PosKaonPair_Radius(0.),
      kMin_PosKaonPair_Pt(0.),
      kMin_PosKaonPair_DCAwrtPV(0.),
      kMax_PosKaonPair_DCAbtwKaons(0.),
      kMax_PosKaonPair_DCApkaSV(0.),
      kMax_PosKaonPair_DCApkbSV(0.),
      kMax_PosKaonPair_Chi2ndf(0.) {
    DefineInput(0, TChain::Class());

    DefineOutput(1, TTree::Class());  // fTree_Events
    DefineOutput(2, TTree::Class());  // fTree_Injected
    DefineOutput(3, TTree::Class());  // fTree_MC
    DefineOutput(4, TTree::Class());  // fTree_Tracks
    DefineOutput(5, TTree::Class());  // fTree_V0s
    DefineOutput(6, TTree::Class());  // fTree_Sexaquarks
    DefineOutput(7, TTree::Class());  // fTree_KaonPairs

    DefineOutput(8, TList::Class());   // fList_QA_Hists
    DefineOutput(9, TList::Class());   // fList_Tracks_Hists
    DefineOutput(10, TList::Class());  // fList_V0s_Hists
    DefineOutput(11, TList::Class());  // fList_Sexaquarks_Hists
    DefineOutput(12, TList::Class());  // fList_PosKaonPairs_Hists
}

/*
 * Destructor.
 * - Note: if `TList::SetOwner(kTRUE)` was called, the TList destructor should delete all objects added to it
 */
AliAnalysisTaskSexaquark::~AliAnalysisTaskSexaquark() {
    if (fTree_Events) delete fTree_Events;
    if (fTree_Injected) delete fTree_Injected;
    if (fTree_MC) delete fTree_MC;
    if (fTree_Tracks) delete fTree_Tracks;
    if (fTree_V0s) delete fTree_V0s;
    if (fTree_Sexaquarks) delete fTree_Sexaquarks;
    if (fTree_KaonPairs) delete fTree_KaonPairs;

    if (fList_QA_Hists) delete fList_QA_Hists;
    if (fList_Tracks_Hists) delete fList_Tracks_Hists;
    if (fList_V0s_Hists) delete fList_V0s_Hists;
    if (fList_Sexaquarks_Hists) delete fList_Sexaquarks_Hists;
    if (fList_PosKaonPairs_Hists) delete fList_PosKaonPairs_Hists;
}

/*
 * Create output objects, called once at RUNTIME ~ execution on Grid
 */
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {

    AliInfo("!! It begins !!");

    /** Add mandatory routines **/

    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) AliFatal("ERROR: AliAnalysisManager couldn't be found.");

    AliESDInputHandler* inputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
    if (!inputHandler) AliFatal("ERROR: AliESDInputHandler couldn't be found.");

    fPIDResponse = inputHandler->GetPIDResponse();

    /** Prepare output  **/

    /* Trees */

    fTree_Events = new TTree("Events", "Events");
    PrepareEventsBranches();
    ClearEventsBranches();
    PostData(1, fTree_Events);

    fTree_Injected = new TTree("Injected", "Injected");
    PrepareInjectedBranches();
    ClearInjectedBranches();
    PostData(2, fTree_Injected);

    fTree_MC = new TTree("MCParticles", "MCParticles");
    PrepareMCParticlesBranches();
    ClearMCParticlesBranches();
    PostData(3, fTree_MC);

    fTree_Tracks = new TTree("Tracks", "Tracks");
    PrepareTracksBranches();
    ClearTracksBranches();
    PostData(4, fTree_Tracks);

    fTree_V0s = new TTree("V0s", "V0s");
    PrepareV0sBranches();
    ClearV0sBranches();
    PostData(5, fTree_V0s);

    fTree_Sexaquarks = new TTree("Sexaquarks", "Sexaquarks");
    PrepareSexaquarksBranches();
    ClearSexaquarksBranches();
    PostData(6, fTree_Sexaquarks);

    fTree_KaonPairs = new TTree("KaonPairs", "KaonPairs");
    PrepareKaonPairsBranches();
    ClearKaonPairsBranches();
    PostData(7, fTree_KaonPairs);

    /* Histograms */

    fList_QA_Hists = new TList();
    fList_QA_Hists->SetOwner(kTRUE);
    PrepareQAHistograms();
    PostData(8, fList_QA_Hists);

    fList_Tracks_Hists = new TList();
    fList_Tracks_Hists->SetOwner(kTRUE);
    PrepareTracksHistograms();
    PostData(9, fList_Tracks_Hists);

    fList_V0s_Hists = new TList();
    fList_V0s_Hists->SetOwner(kTRUE);
    PrepareV0Histograms();
    PostData(10, fList_V0s_Hists);

    fList_Sexaquarks_Hists = new TList();
    fList_Sexaquarks_Hists->SetOwner(kTRUE);
    PrepareSexaquarkHistograms();
    PostData(11, fList_Sexaquarks_Hists);

    fList_PosKaonPairs_Hists = new TList();
    fList_PosKaonPairs_Hists->SetOwner(kTRUE);
    PreparePosKaonPairHistograms();
    PostData(12, fList_PosKaonPairs_Hists);

    AliInfo("!! It ends !!");
}

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareEventsBranches() {

    fTree_Events->Branch("EventNumber", &fEventNumber);
    fTree_Events->Branch("DirNumber", &fDirNumber);
    fTree_Events->Branch("RunNumber", &fRunNumber);
    fTree_Events->Branch("Centrality", &fCentrality);

    fTree_Events->Branch("PV_TrueXv", &tEvents_PV_TrueXv);
    fTree_Events->Branch("PV_TrueYv", &tEvents_PV_TrueYv);
    fTree_Events->Branch("PV_TrueZv", &tEvents_PV_TrueZv);
    fTree_Events->Branch("PV_RecXv", &tEvents_PV_RecXv);
    fTree_Events->Branch("PV_RecYv", &tEvents_PV_RecYv);
    fTree_Events->Branch("PV_RecZv", &tEvents_PV_RecZv);
    fTree_Events->Branch("NTracks", &tEvents_NTracks);
    fTree_Events->Branch("IsMB", &tEvents_IsMB);
    fTree_Events->Branch("IsHighMultV0", &tEvents_IsHighMultV0);
    fTree_Events->Branch("tEvents_IsHighMultSPD", &tEvents_IsHighMultSPD);
    fTree_Events->Branch("IsCentral", &tEvents_IsCentral);
    fTree_Events->Branch("IsSemiCentral", &tEvents_IsSemiCentral);
}

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareMCParticlesBranches() {

    fTree_MC->Branch("EventNumber", &fEventNumber);
    fTree_MC->Branch("DirNumber", &fDirNumber);
    fTree_MC->Branch("RunNumber", &fRunNumber);
    fTree_MC->Branch("Idx", &tMC_Idx);
    fTree_MC->Branch("PdgCode", &tMC_PdgCode);
    fTree_MC->Branch("Idx_Mother", &tMC_Idx_Mother);
    fTree_MC->Branch("NDaughters", &tMC_NDaughters);
    fTree_MC->Branch("Idx_FirstDau", &tMC_Idx_FirstDau);
    fTree_MC->Branch("Idx_LastDau", &tMC_Idx_LastDau);
    fTree_MC->Branch("ReactionID", &tMC_ReactionID);
    fTree_MC->Branch("Px", &tMC_Px);
    fTree_MC->Branch("Py", &tMC_Py);
    fTree_MC->Branch("Pz", &tMC_Pz);
    fTree_MC->Branch("Xv_i", &tMC_Xv_i);
    fTree_MC->Branch("Yv_i", &tMC_Yv_i);
    fTree_MC->Branch("Zv_i", &tMC_Zv_i);
    fTree_MC->Branch("Status", &tMC_Status);
    fTree_MC->Branch("IsSecondary", &tMC_IsSecondary);
    fTree_MC->Branch("IsSignal", &tMC_IsSignal);
}

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareTracksBranches() {
    fTree_Tracks->Branch("EventNumber", &fEventNumber);
    fTree_Tracks->Branch("DirNumber", &fDirNumber);
    fTree_Tracks->Branch("RunNumber", &fRunNumber);
    fTree_Tracks->Branch("Idx", &tTrack_Idx);
    fTree_Tracks->Branch("Idx_True", &tTrack_Idx_True);
    fTree_Tracks->Branch("Px", &tTrack_Px);
    fTree_Tracks->Branch("Py", &tTrack_Py);
    fTree_Tracks->Branch("Pz", &tTrack_Pz);
    fTree_Tracks->Branch("Charge", &tTrack_Charge);
    fTree_Tracks->Branch("NSigmaPion", &tTrack_NSigmaPion);
    fTree_Tracks->Branch("NSigmaKaon", &tTrack_NSigmaKaon);
    fTree_Tracks->Branch("NSigmaProton", &tTrack_NSigmaProton);
    fTree_Tracks->Branch("IsTrue", &tTrack_IsTrue);
    fTree_Tracks->Branch("IsSecondary", &tTrack_IsSecondary);
    fTree_Tracks->Branch("IsSignal", &tTrack_IsSignal);
    fTree_Tracks->Branch("ReactionID", &tTrack_ReactionID);
}

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareV0sBranches() {

    fTree_V0s->Branch("EventNumber", &fEventNumber);
    fTree_V0s->Branch("DirNumber", &fDirNumber);
    fTree_V0s->Branch("RunNumber", &fRunNumber);
    fTree_V0s->Branch("Idx", &tV0_Idx);
    fTree_V0s->Branch("Idx_Pos", &tV0_Idx_Pos);
    fTree_V0s->Branch("Idx_Neg", &tV0_Idx_Neg);
    fTree_V0s->Branch("Idx_True", &tV0_Idx_True);
    fTree_V0s->Branch("PID", &tV0_PID);
    fTree_V0s->Branch("Px", &tV0_Px);
    fTree_V0s->Branch("Py", &tV0_Py);
    fTree_V0s->Branch("Pz", &tV0_Pz);
    fTree_V0s->Branch("E", &tV0_E);
    fTree_V0s->Branch("Xv_f", &tV0_Xv_f);
    fTree_V0s->Branch("Yv_f", &tV0_Yv_f);
    fTree_V0s->Branch("Zv_f", &tV0_Zv_f);
    fTree_V0s->Branch("Pos_Px", &tV0_Pos_Px);
    fTree_V0s->Branch("Pos_Py", &tV0_Pos_Py);
    fTree_V0s->Branch("Pos_Pz", &tV0_Pos_Pz);
    fTree_V0s->Branch("Neg_Px", &tV0_Neg_Px);
    fTree_V0s->Branch("Neg_Py", &tV0_Neg_Py);
    fTree_V0s->Branch("Neg_Pz", &tV0_Neg_Pz);
    fTree_V0s->Branch("IsTrue", &tV0_IsTrue);
    fTree_V0s->Branch("IsSecondary", &tV0_IsSecondary);
    fTree_V0s->Branch("IsSignal", &tV0_IsSignal);
    fTree_V0s->Branch("ReactionID", &tV0_ReactionID);
}

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareSexaquarksBranches() {

    fTree_Sexaquarks->Branch("EventNumber", &fEventNumber);
    fTree_Sexaquarks->Branch("DirNumber", &fDirNumber);
    fTree_Sexaquarks->Branch("RunNumber", &fRunNumber);
    fTree_Sexaquarks->Branch("Idx_AL", &tSexaquark_Idx_AL);
    fTree_Sexaquarks->Branch("Idx_AL_Neg", &tSexaquark_Idx_AL_Neg);
    fTree_Sexaquarks->Branch("Idx_AL_Pos", &tSexaquark_Idx_AL_Pos);
    fTree_Sexaquarks->Branch("Px", &tSexaquark_Px);
    fTree_Sexaquarks->Branch("Py", &tSexaquark_Py);
    fTree_Sexaquarks->Branch("Pz", &tSexaquark_Pz);
    fTree_Sexaquarks->Branch("E", &tSexaquark_E);
    fTree_Sexaquarks->Branch("SV_Xv", &tSexaquark_SV_Xv);
    fTree_Sexaquarks->Branch("SV_Yv", &tSexaquark_SV_Yv);
    fTree_Sexaquarks->Branch("SV_Zv", &tSexaquark_SV_Zv);
    fTree_Sexaquarks->Branch("DecayLength", &tSexaquark_DecayLength);
    fTree_Sexaquarks->Branch("CPAwrtPV", &tSexaquark_CPAwrtPV);
    fTree_Sexaquarks->Branch("DCAwrtPV", &tSexaquark_DCAwrtPV);
    fTree_Sexaquarks->Branch("Chi2ndf", &tSexaquark_Chi2ndf);
    fTree_Sexaquarks->Branch("IsHybrid", &tSexaquark_IsHybrid);
    fTree_Sexaquarks->Branch("IsSignal", &tSexaquark_IsSignal);
    fTree_Sexaquarks->Branch("ReactionID", &tSexaquark_ReactionID);

    if (fReactionID == 'A') {
        fTree_Sexaquarks->Branch("Idx_K0S", &tSexaquarkA_Idx_K0S);
        fTree_Sexaquarks->Branch("Idx_K0S_Neg", &tSexaquarkA_Idx_K0S_Neg);
        fTree_Sexaquarks->Branch("Idx_K0S_Pos", &tSexaquarkA_Idx_K0S_Pos);
        fTree_Sexaquarks->Branch("DCAbtwV0s", &tSexaquarkA_DCAbtwV0s);
        fTree_Sexaquarks->Branch("DCAV0aSV", &tSexaquarkA_DCAV0aSV);
        fTree_Sexaquarks->Branch("DCAV0bSV", &tSexaquarkA_DCAV0bSV);
        fTree_Sexaquarks->Branch("DCAV0anegSV", &tSexaquarkA_DCAV0anegSV);
        fTree_Sexaquarks->Branch("DCAV0aposSV", &tSexaquarkA_DCAV0aposSV);
        fTree_Sexaquarks->Branch("DCAV0bnegSV", &tSexaquarkA_DCAV0bnegSV);
        fTree_Sexaquarks->Branch("DCAV0bposSV", &tSexaquarkA_DCAV0bposSV);
    } else if (fReactionID == 'D') {
        fTree_Sexaquarks->Branch("Idx_PosKaon", &tSexaquarkD_Idx_PosKaon);
        fTree_Sexaquarks->Branch("DCAV0Ba", &tSexaquarkD_DCAV0Ba);
        fTree_Sexaquarks->Branch("DCAV0negBa", &tSexaquarkD_DCAV0negBa);
        fTree_Sexaquarks->Branch("DCAV0posBa", &tSexaquarkD_DCAV0posBa);
        fTree_Sexaquarks->Branch("DCAV0SV", &tSexaquarkD_DCAV0SV);
        fTree_Sexaquarks->Branch("DCAV0negSV", &tSexaquarkD_DCAV0negSV);
        fTree_Sexaquarks->Branch("DCAV0posSV", &tSexaquarkD_DCAV0posSV);
        fTree_Sexaquarks->Branch("DCABaSV", &tSexaquarkD_DCABaSV);
    } else if (fReactionID == 'E') {
        fTree_Sexaquarks->Branch("Idx_PP", &tSexaquarkE_Idx_PP);
        fTree_Sexaquarks->Branch("Idx_PP_Neg", &tSexaquarkE_Idx_PP_Neg);
        fTree_Sexaquarks->Branch("Idx_PP_Pos", &tSexaquarkE_Idx_PP_Pos);
        fTree_Sexaquarks->Branch("Idx_PosKaon", &tSexaquarkE_Idx_PosKaon);
        fTree_Sexaquarks->Branch("DCAV0SV", &tSexaquarkE_DCAV0SV);
        fTree_Sexaquarks->Branch("DCAV0negSV", &tSexaquarkE_DCAV0negSV);
        fTree_Sexaquarks->Branch("DCAV0posSV", &tSexaquarkE_DCAV0posSV);
        fTree_Sexaquarks->Branch("DCApkSV", &tSexaquarkE_DCApkSV);
        fTree_Sexaquarks->Branch("DCApmSV", &tSexaquarkE_DCApmSV);
        fTree_Sexaquarks->Branch("DCAppSV", &tSexaquarkE_DCAppSV);
        fTree_Sexaquarks->Branch("DCApkV0", &tSexaquarkE_DCApkV0);
        fTree_Sexaquarks->Branch("DCApmV0", &tSexaquarkE_DCApmV0);
        fTree_Sexaquarks->Branch("DCAppV0", &tSexaquarkE_DCAppV0);
        fTree_Sexaquarks->Branch("DCApmPK", &tSexaquarkE_DCApmPK);
        fTree_Sexaquarks->Branch("DCAppPK", &tSexaquarkE_DCAppPK);
    }
}

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareKaonPairsBranches() {

    fTree_KaonPairs->Branch("EventNumber", &fEventNumber);
    fTree_KaonPairs->Branch("DirNumber", &fDirNumber);
    fTree_KaonPairs->Branch("RunNumber", &fRunNumber);
    fTree_KaonPairs->Branch("Idx", &tKaonPair_Idx);
    fTree_KaonPairs->Branch("Idx_KaonA", &tKaonPair_Idx_KaonA);
    fTree_KaonPairs->Branch("Idx_KaonB", &tKaonPair_Idx_KaonB);
    fTree_KaonPairs->Branch("Px", &tKaonPair_Px);
    fTree_KaonPairs->Branch("Py", &tKaonPair_Py);
    fTree_KaonPairs->Branch("Pz", &tKaonPair_Pz);
    fTree_KaonPairs->Branch("E", &tKaonPair_E);
    fTree_KaonPairs->Branch("SV_Xv", &tKaonPair_SV_Xv);
    fTree_KaonPairs->Branch("SV_Yv", &tKaonPair_SV_Yv);
    fTree_KaonPairs->Branch("SV_Zv", &tKaonPair_SV_Zv);
    fTree_KaonPairs->Branch("DecayLength", &tKaonPair_DecayLength);
    fTree_KaonPairs->Branch("DCAwrtPV", &tKaonPair_DCAwrtPV);
    fTree_KaonPairs->Branch("DCAbtwKaons", &tKaonPair_DCAbtwKaons);
    fTree_KaonPairs->Branch("DCAkaSV", &tKaonPair_DCAkaSV);
    fTree_KaonPairs->Branch("DCAkbSV", &tKaonPair_DCAkbSV);
    fTree_KaonPairs->Branch("OpeningAngle", &tKaonPair_OpeningAngle);
    fTree_KaonPairs->Branch("Chi2ndf", &tKaonPair_Chi2ndf);
    fTree_KaonPairs->Branch("IsHybrid", &tKaonPair_IsHybrid);
    fTree_KaonPairs->Branch("IsSignal", &tKaonPair_IsSignal);
    fTree_KaonPairs->Branch("ReactionID", &tKaonPair_ReactionID);
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearEventsBranches() {
    // not clearing fRunNumber nor fDirNumber
    fEventNumber = 0;
    fCentrality = 0;

    tEvents_PV_TrueXv = -999.;
    tEvents_PV_TrueYv = -999.;
    tEvents_PV_TrueZv = -999.;
    tEvents_PV_RecXv = 0.;
    tEvents_PV_RecYv = 0.;
    tEvents_PV_RecZv = 0.;
    tEvents_NTracks = 0;
    tEvents_IsMB = kFALSE;
    tEvents_IsHighMultV0 = kFALSE;
    tEvents_IsHighMultSPD = kFALSE;
    tEvents_IsCentral = kFALSE;
    tEvents_IsSemiCentral = kFALSE;
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearMCParticlesBranches() {
    tMC_Idx = 0;
    tMC_PdgCode = 0;
    tMC_Idx_Mother = 0;
    tMC_NDaughters = 0;
    tMC_Idx_FirstDau = 0;
    tMC_Idx_LastDau = 0;
    tMC_ReactionID = 0;
    tMC_Px = 0.;
    tMC_Py = 0.;
    tMC_Pz = 0.;
    tMC_Xv_i = 0.;
    tMC_Yv_i = 0.;
    tMC_Zv_i = 0.;
    tMC_Status = 0;
    tMC_IsSecondary = kFALSE;
    tMC_IsSignal = kFALSE;
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearTracksBranches() {
    tTrack_Idx = 0;
    tTrack_Idx_True = 0;
    tTrack_Px = 0.;
    tTrack_Py = 0.;
    tTrack_Pz = 0.;
    tTrack_Charge = 0;
    tTrack_NSigmaPion = 0.;
    tTrack_NSigmaKaon = 0.;
    tTrack_NSigmaProton = 0.;
    tTrack_IsTrue = kFALSE;
    tTrack_IsSecondary = kFALSE;
    tTrack_IsSignal = kFALSE;
    tTrack_ReactionID = 0;
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearV0sBranches() {
    tV0_Idx = 0;
    tV0_Idx_Pos = 0;
    tV0_Idx_Neg = 0;
    tV0_Idx_True = 0;
    tV0_PID = 0;
    tV0_Px = 0.;
    tV0_Py = 0.;
    tV0_Pz = 0.;
    tV0_E = 0.;
    tV0_Xv_f = 0.;
    tV0_Yv_f = 0.;
    tV0_Zv_f = 0.;
    tV0_Pos_Px = 0.;
    tV0_Pos_Py = 0.;
    tV0_Pos_Pz = 0.;
    tV0_Neg_Px = 0.;
    tV0_Neg_Py = 0.;
    tV0_Neg_Pz = 0.;
    tV0_IsTrue = kFALSE;
    tV0_IsSecondary = kFALSE;
    tV0_IsSignal = kFALSE;
    tV0_ReactionID = 0;
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearSexaquarksBranches() {
    tSexaquark_Idx_AL = 0;
    tSexaquark_Idx_AL_Neg = 0;
    tSexaquark_Idx_AL_Pos = 0;
    tSexaquark_Px = 0.;
    tSexaquark_Py = 0.;
    tSexaquark_Pz = 0.;
    tSexaquark_E = 0.;
    tSexaquark_SV_Xv = 0.;
    tSexaquark_SV_Yv = 0.;
    tSexaquark_SV_Zv = 0.;
    tSexaquark_DecayLength = 0.;
    tSexaquark_CPAwrtPV = 0.;
    tSexaquark_DCAwrtPV = 0.;
    tSexaquark_Chi2ndf = 0.;
    tSexaquark_IsHybrid = kFALSE;
    tSexaquark_IsSignal = kFALSE;
    tSexaquark_ReactionID = 0;

    if (fReactionID == 'A') {
        tSexaquarkA_Idx_K0S = 0;
        tSexaquarkA_Idx_K0S_Neg = 0;
        tSexaquarkA_Idx_K0S_Pos = 0;
        tSexaquarkA_DCAbtwV0s = 0.;
        tSexaquarkA_DCAV0aSV = 0.;
        tSexaquarkA_DCAV0bSV = 0.;
        tSexaquarkA_DCAV0anegSV = 0.;
        tSexaquarkA_DCAV0aposSV = 0.;
        tSexaquarkA_DCAV0bnegSV = 0.;
        tSexaquarkA_DCAV0bposSV = 0.;
    } else if (fReactionID == 'D') {
        tSexaquarkD_Idx_PosKaon = 0;
        tSexaquarkD_DCAV0Ba = 0.;
        tSexaquarkD_DCAV0negBa = 0.;
        tSexaquarkD_DCAV0posBa = 0.;
        tSexaquarkD_DCAV0SV = 0.;
        tSexaquarkD_DCAV0negSV = 0.;
        tSexaquarkD_DCAV0posSV = 0.;
        tSexaquarkD_DCABaSV = 0.;
    } else if (fReactionID == 'E') {
        tSexaquarkE_Idx_PP = 0;
        tSexaquarkE_Idx_PP_Neg = 0;
        tSexaquarkE_Idx_PP_Pos = 0;
        tSexaquarkE_Idx_PosKaon = 0;
        tSexaquarkE_DCAV0SV = 0.;
        tSexaquarkE_DCAV0negSV = 0.;
        tSexaquarkE_DCAV0posSV = 0.;
        tSexaquarkE_DCApkSV = 0.;
        tSexaquarkE_DCApmSV = 0.;
        tSexaquarkE_DCAppSV = 0.;
        tSexaquarkE_DCApkV0 = 0.;
        tSexaquarkE_DCApmV0 = 0.;
        tSexaquarkE_DCAppV0 = 0.;
        tSexaquarkE_DCApmPK = 0.;
        tSexaquarkE_DCAppPK = 0.;
    }
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearKaonPairsBranches() {
    tKaonPair_Idx = 0;
    tKaonPair_Idx_KaonA = 0;
    tKaonPair_Idx_KaonB = 0;
    tKaonPair_Px = 0.;
    tKaonPair_Py = 0.;
    tKaonPair_Pz = 0.;
    tKaonPair_E = 0.;
    tKaonPair_SV_Xv = 0.;
    tKaonPair_SV_Yv = 0.;
    tKaonPair_SV_Zv = 0.;
    tKaonPair_DecayLength = 0.;
    tKaonPair_DCAwrtPV = 0.;
    tKaonPair_DCAbtwKaons = 0.;
    tKaonPair_DCAkaSV = 0.;
    tKaonPair_DCAkbSV = 0.;
    tKaonPair_OpeningAngle = 0.;
    tKaonPair_Chi2ndf = 0.;
    tKaonPair_IsHybrid = kFALSE;
    tKaonPair_IsSignal = kFALSE;
    tKaonPair_ReactionID = 0;
}

/*
 * Define and add QA histograms to the histograms output list.
 * - Uses: `fHist_QA`, `f2DHist_QA`
 */
void AliAnalysisTaskSexaquark::PrepareQAHistograms() {

    /* 1D Hists */

    // clang-format off
    std::vector<TString> QA_props = {// MC. Gen. Event
                                     "MCGen_VertexX", "MCGen_VertexY", "MCGen_VertexZ",
                                     // Event
                                     "Events_Bookkeep", "VertexX", "VertexY", "VertexZ",
                                     "NTracks", "NSelectedTracks",
                                     // Tracks
                                     "Tracks_DCAxy", "Tracks_DCAz", "Tracks_Pt",
                                     // SPD
                                     "SPD_NTracklets",
                                     // TPC
                                     "TPC_NClusters"};
    std::vector<Int_t> QA_nbins = {100, 100, 200,
                                   5, 100, 100, 200,
                                   2000, 2000,
                                   100, 100, 200,
                                   60,
                                   200};
    std::vector<Float_t> QA_min = {-1, -1, -50,
                                   0., -1, -1, -50,
                                   0., 0.,
                                   -5., -5., 0.,
                                   0.,
                                   0.};
    std::vector<Float_t> QA_max = {1, 1, 50,
                                   5., 1, 1, 50,
                                   20000., 20000.,
                                   5., 5., 100.,
                                   6000,
                                   5000000.};
    // clang-format on

    TString histKey;
    for (Int_t prop_idx = 0; prop_idx < (Int_t)QA_props.size(); prop_idx++) {
        histKey = QA_props[prop_idx];
        fHist_QA[histKey] = new TH1F("QA_" + histKey, "", QA_nbins[prop_idx], QA_min[prop_idx], QA_max[prop_idx]);
        fList_QA_Hists->Add(fHist_QA[histKey]);
    }

    /* 2D Hists */

    // clang-format off
    QA_props = {// Tracks
                "Tracks_PhiVsEta",
                // SPD
                "SPD_MultiplicityVsVertexZ", "SPD_NTrackletsClu0VsNTracklets", "SPD_NTrackletsClu1VsNTracklets", "SPD_PhiVsEta",
                // ITS
                "ITS_LayerNoVsPhi", "ITS_LayerNoVsEta",
                // TPC
                "TPC_NSigmaProtonVsInnerParamP", "TPC_NSigmaKaonVsInnerParamP", "TPC_NSigmaPionVsInnerParamP",
                "TPC_NSigmaProtonVsEta", "TPC_NSigmaKaonVsEta", "TPC_NSigmaPionVsEta"};
    std::vector<Int_t> QA_nbinsx = {100,
                                   10, 200, 200, 400,
                                   200, 40,
                                   200, 200, 200, 200, 200, 200};
    std::vector<Float_t> QA_xlow = {-1.,
                                   -10., 0., 0., -1.,
                                   0., -1.,
                                   0., 0., 0., -1., -1., -1.};
    std::vector<Float_t> QA_xup = {1.,
                                   10., 3000., 3000., 1.,
                                   2 * TMath::Pi(), 1.,
                                   10., 10., 10., 1., 1., 1.};
    std::vector<Int_t> QA_nbinsy = {200,
                                   30, 200, 200, 400,
                                   6, 6,
                                   200, 200, 200, 200, 200, 200};
    std::vector<Float_t> QA_ylow = {0.,
                                    0., 0., 0., 0.,
                                    0., 0.,
                                    -5., -5., -5., -5., -5., -5.};
    std::vector<Float_t> QA_yup = {2 * TMath::Pi(),
                                   6000., 6000., 6000., 2 * TMath::Pi(),
                                   6., 6.,
                                   5., 5., 5., 5., 5., 5.};
    // clang-format on

    for (Int_t prop_idx = 0; prop_idx < (Int_t)QA_props.size(); prop_idx++) {
        histKey = QA_props[prop_idx];
        if (fReactionID == 'A' && histKey.Contains("Kaon")) continue;
        f2DHist_QA[histKey] = new TH2F("QA_" + histKey, "", QA_nbinsx[prop_idx], QA_xlow[prop_idx], QA_xup[prop_idx], QA_nbinsy[prop_idx],
                                       QA_ylow[prop_idx], QA_yup[prop_idx]);
        fList_QA_Hists->Add(f2DHist_QA[histKey]);
    }
}

/*
 * Define and add tracks' histograms to the histograms output list.
 * - Uses: `fHist_Tracks_Bookkeep`, `fHist_Tracks`, `fOutputListOfHists`
 */
void AliAnalysisTaskSexaquark::PrepareTracksHistograms() {

    TString tracks_stages[2] = {"MCGen", "Found"};
    std::vector<TString> tracks_sets = {"All", "True", "Secondary", "Signal"};
    std::set<Int_t> tracks_species(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end());
    std::map<Int_t, TString> tracks_name = {{-2212, "AntiProton"}, {321, "PosKaon"}, {211, "PiPlus"}, {-211, "PiMinus"}};

    // clang-format off
    std::vector<TString> tracks_props = {"Px",           "Py",          "Pt",         "Pz",           "Eta",
                                         "DCA_wrtPV",    "DCAxy_wrtPV", "DCAz_wrtPV", "NTPCClusters", "NCrossedRows",
                                         "NFindableCls", "NSharedCls",  "Chi2/NCls",
                                         "NSigmaProton", "NSigmaKaon",  "NSigmaPion", "Status"};
    std::vector<Int_t> tracks_nbins = {100, 100, 100, 100, 100,
                                       100, 100, 100, 160, 160,
                                       100, 100, 140,
                                       100, 100, 100, 16};
    std::vector<Double_t> tracks_min = {-20., -20., 0.,  -20., -1.2,
                                        0.,   0.,   0.,  0.,   0.,
                                        0.,   0.,   0.,
                                        -10., -10,  -10, 0.};
    std::vector<Double_t> tracks_max = {20.,  20.,  10.,  20.,  1.2,
                                        200., 200., 200., 160., 160.,
                                        100., 100., 7.,
                                        10,   10,   10,   16};
    // clang-format on

    for (Int_t species : tracks_species) {
        fHist_Tracks_Bookkeep[species] = new TH1F(Form("%s_Bookkeep", tracks_name[species].Data()), "", 150, 0., 150.);
        fList_Tracks_Hists->Add(fHist_Tracks_Bookkeep[species]);
    }

    for (TString& stage : tracks_stages) {
        for (TString& set : tracks_sets) {
            for (Int_t species : tracks_species) {
                for (Int_t prop_idx = 0; prop_idx < (Int_t)tracks_props.size(); prop_idx++) {

                    if (!fIsMC && (stage == "MCGen" || set == "True" || set == "Secondary" || set == "Signal")) continue;

                    if (stage == "MCGen" && set == "True") continue;
                    if (stage == "MCGen" && !(tracks_props[prop_idx] == "Px" || tracks_props[prop_idx] == "Py" || tracks_props[prop_idx] == "Pt" ||
                                              tracks_props[prop_idx] == "Pz" || tracks_props[prop_idx] == "Eta")) {
                        continue;
                    }

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), tracks_name[species].Data(), tracks_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, species, tracks_props[prop_idx]);
                    fHist_Tracks[histKey] = new TH1F(histName, "", tracks_nbins[prop_idx], tracks_min[prop_idx], tracks_max[prop_idx]);
                    fList_Tracks_Hists->Add(fHist_Tracks[histKey]);
                }
            }
        }
    }

    tracks_sets.insert(tracks_sets.begin() + 1, "Selected");
    for (TString& set : tracks_sets) {
        TString histName = Form("Tracks_%s_TPCsignal", set.Data());
        f2DHist_Tracks_TPCsignal[set] = new TH2F(histName, "", 200, 0., 10., 200, 0., 400.);
        fList_Tracks_Hists->Add(f2DHist_Tracks_TPCsignal[set]);
    }
}

/*
 * Define and add V0s' histograms to the histograms output list.
 * - Uses: `fSourceOfV0s`, `fOutputListOfHists`
 */
void AliAnalysisTaskSexaquark::PrepareV0Histograms() {

    TString V0_stages[3] = {"MCGen", "Findable", "Found"};
    TString V0_sets[4] = {"All", "True", "Secondary", "Signal"};

    std::vector<Int_t> V0_species = {-3122, 310, 422};
    std::vector<TString> V0_name = {"AntiLambda", "KaonZeroShort", "PionPair"};

    // clang-format off
    std::vector<TString> V0_props = {"Mass",        "Pt",           "Px",          "Py",       "Pz",        "Phi",
                                     "OriginRadius", "DecayRadius", "Zv",       "Eta",       "Rapidity",
                                     "DecayLength", "OpeningAngle", "ArmQt",    "ArmAlpha",       "CPAwrtPV",    "DCAwrtPV",
                                     "DCAbtwDau",   "DCAnegV0",     "DCAposV0", "ImprvDCAbtwDau", "Chi2",        "Chi2ndf"};
    std::vector<Int_t> V0_nbins = {100, 100, 100, 100, 100, 100,
                                   100, 100, 100, 100, 100,
                                   100, 100, 100, 100, 100, 100,
                                   100, 100, 100, 100, 100, 100};
    std::vector<Double_t> V0_min = {0., 0., -20., -20., -20., 0.,
                                    0., 0.,   -50., -4.,  -4.,
                                    0., 0., 0.,   -1.,  0.,   0.,
                                    0., 0., 0.,   0.,   0.,   0.};
    std::vector<Double_t> V0_max = {10.,  10.,         20.,  20., 20., 2. * TMath::Pi(),
                                    250.,        250., 50., 4.,  4.,
                                    500., TMath::Pi(), 0.3,  1.,  1.,  200.,
                                    2.,   2.,          2.,   1.,  2.,  0.5};
    // clang-format on

    for (Int_t sp = 0; sp < (Int_t)V0_species.size(); sp++) {
        if (fReactionID != 'A' && V0_species[sp] == 310) continue;
        if (fReactionID != 'E' && V0_species[sp] == 422) continue;
        fHist_V0s_Bookkeep[V0_species[sp]] = new TH1F(Form("%s_Bookkeep", V0_name[sp].Data()), "", 200, 0., 200.);
        fList_V0s_Hists->Add(fHist_V0s_Bookkeep[V0_species[sp]]);
    }

    for (TString& stage : V0_stages) {
        for (TString& set : V0_sets) {
            for (Int_t sp = 0; sp < (Int_t)V0_species.size(); sp++) {
                for (Int_t prop_idx = 0; prop_idx < (Int_t)V0_props.size(); prop_idx++) {

                    if (!fIsMC && (stage == "MCGen" || stage == "Findable" || set == "True" || set == "Secondary" || set == "Signal")) continue;
                    if (fReactionID != 'A' && V0_species[sp] == 310) continue;
                    if (fReactionID != 'E' && V0_species[sp] == 422) continue;

                    if (fSourceOfV0s == "kalman" && (V0_props[prop_idx] == "ImprvDCAbtwDau" || V0_props[prop_idx] == "Chi2")) continue;
                    if ((fSourceOfV0s == "offline" || fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "custom") && V0_props[prop_idx] == "Chi2ndf") {
                        continue;
                    }

                    if (stage == "Findable" && set == "All") continue;
                    if (stage == "MCGen" && set == "True") continue;

                    if (V0_species[sp] == 422) {
                        if (stage == "MCGen") continue;
                        if (stage == "Findable" || stage == "Found") {
                            if (set == "True" || set == "Secondary") continue;
                        }
                    }

                    if ((stage == "Findable" || stage == "MCGen") &&
                        (V0_props[prop_idx] == "DCAbtwDau" || V0_props[prop_idx] == "DCAnegV0" || V0_props[prop_idx] == "DCAposV0" ||
                         V0_props[prop_idx] == "ImprvDCAbtwDau" || V0_props[prop_idx] == "Chi2" || V0_props[prop_idx] == "Chi2ndf")) {
                        continue;
                    }

                    if (stage == "Found" && V0_props[prop_idx] == "OriginRadius") continue;

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), V0_name[sp].Data(), V0_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, V0_species[sp], V0_props[prop_idx]);
                    Float_t minEdge = V0_props[prop_idx] == "Mass" ? kMin_V0_Mass[V0_species[sp]] : V0_min[prop_idx];
                    Float_t maxEdge = V0_props[prop_idx] == "Mass" ? kMax_V0_Mass[V0_species[sp]] : V0_max[prop_idx];
                    fHist_V0s[histKey] = new TH1F(histName, "", V0_nbins[prop_idx], minEdge, maxEdge);
                    fList_V0s_Hists->Add(fHist_V0s[histKey]);
                }
            }
        }
    }

    for (TString& set : V0_sets) {
        TString histName = Form("Found_%s_V0s_ArmenterosPodolanski", set.Data());
        f2DHist_V0s_ArmenterosPodolanski[set] = new TH2F(histName, "", 200, -1., 1., 200, 0., 0.3);
        fList_V0s_Hists->Add(f2DHist_V0s_ArmenterosPodolanski[set]);
    }
}

/*
 * Define and add anti-sexaquark histograms to the histograms output list.
 * - Uses: `fHist_AntiSexaquarks_Bookkeep`, `fList_Sexaquarks_Hists`, `fHist_AntiSexaquarks`, `fOutputListOfHists`
 */
void AliAnalysisTaskSexaquark::PrepareSexaquarkHistograms() {

    TString SQ_stages[3] = {"MCGen", "Findable", "Found"};
    TString SQ_sets[2] = {"All", "Signal"};

    std::vector<TString> SQ_props;
    std::vector<Int_t> SQ_nbins;
    std::vector<Double_t> SQ_min;
    std::vector<Double_t> SQ_max;

    // clang-format off
    if (fReactionID == 'A') {
        SQ_props = {"Mass",     "Pt",       "Px",          "Py",          "Pz",          "Phi",         "Radius",
                    "Zv",       "Eta",      "Rapidity",    "DecayLength", "CPAwrtPV",    "DCAwrtPV",    "DCAbtwV0s",
                    "DCAV0aSV", "DCAV0bSV", "DCAV0anegSV", "DCAV0aposSV", "DCAV0bnegSV", "DCAV0bposSV", "Chi2ndf"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100};
        SQ_min = {-10., -10., -10., -10., -50., -2 * TMath::Pi(), 0.,
                  -50., -2.,  -5.,  0.,   0.,   0.,               0.,
                  0.,   0.,   0.,   0.,   0.,   0.,               0.};
        SQ_max = {10., 10., 10., 10.,  50., 2 * TMath::Pi(), 200.,
                  50., 2.,  5.,  500., 1.,  10.,             2.,
                  2.,  2.,  10., 10.,  10., 10.,             1.};
    } else if (fReactionID == 'D') {
        SQ_props = {"Mass",       "Pt",         "Px",       "Py",          "Pz",         "Phi",        "Radius",
                    "Zv",         "Eta",        "Rapidity", "DecayLength", "CPAwrtPV",   "DCAwrtPV",   "DCAV0Ba",
                    "DCAV0negBa", "DCAV0posBa", "DCAV0SV",  "DCABaSV",     "DCAV0negSV", "DCAV0posSV", "Chi2ndf"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100};
        SQ_min = {-10., -10., -10., -10., -50., 0., 0.,
                  -50., -2.,  -5.,  0.,   0.,   0., 0.,
                  0.,   0.,   0.,   0.,   0.,   0., 0.};
        SQ_max = {10., 10., 10., 10.,  50., 2 * TMath::Pi(), 200.,
                  50., 2.,  5.,  500., 1.,  10.,             10.,
                  10., 10., 10., 10.,  10., 10.,             1.};
    } else if (fReactionID == 'E') {
        SQ_props = {"Mass",    "Pt",       "Px",          "Py",       "Pz",       "Phi",     "Radius",     "Zv",
                    "Eta",     "Rapidity", "DecayLength", "CPAwrtPV", "DCAwrtPV", "DCAV0SV", "DCAV0negSV", "DCAV0posSV",
                    "DCApkSV", "DCApmSV",  "DCAppSV",     "DCApkV0",  "DCApmV0",  "DCAppV0", "DCApmPK",    "DCAppPK",
                    "Chi2ndf"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100, 100,
                    100};
        SQ_min = {-10., -10., -10., -10., -50., 0., 0., -50.,
                  -2.,  -5.,  0.,   0.,   0.,   0., 0., 0.,
                  0.,   0.,   0.,   0.,   0.,   0., 0., 0.,
                  0.};
        SQ_max = {10., 10., 10.,  10., 50., 2 * TMath::Pi(), 200., 50.,
                  2.,  5.,  500., 1.,  10., 10.,             10.,  10.,
                  10., 10., 10.,  10., 10., 10.,             10.,  1.,
                  30.};
    } else if (fReactionID == 'H') {
        SQ_props = {"Mass", "Pt",       "Px",          "Py",       "Pz",       "Phi",   "Radius",   "Zv",
                    "Eta",  "Rapidity", "DecayLength", "CPAwrtPV", "DCAwrtPV"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100};
        SQ_min = {-10., -10., -10., -10., -50., 0., 0., -50.,
                  -2.,  -5.,  0.,   0.,   0.};
        SQ_max = {10., 10., 10.,  10., 50., 2 * TMath::Pi(), 200., 50.,
                  2.,  5.,  500., 1.,  10.};
    }
    // clang-format on

    fHist_AntiSexaquarks_Bookkeep = new TH1D("AntiSexaquarks_Bookkeep", "", 100, 0., 100.);
    fList_Sexaquarks_Hists->Add(fHist_AntiSexaquarks_Bookkeep);

    for (TString& stage : SQ_stages) {
        for (TString& set : SQ_sets) {
            for (Int_t prop_idx = 0; prop_idx < (Int_t)SQ_props.size(); prop_idx++) {
                if (!fIsMC && (stage.Contains("MCGen") || stage == "Findable" || set == "Signal")) continue;
                if (fReactionID == 'H' && (stage == "Findable" || stage == "Found")) continue;

                if (set == "Signal" && !stage.Contains("Found")) continue;

                if (fSourceOfV0s != "kalman" && SQ_props[prop_idx] == "Chi2ndf") continue;

                if ((stage.Contains("MCGen") || stage == "Findable") &&
                    ((SQ_props[prop_idx].Contains("DCA") && SQ_props[prop_idx] != "DCAwrtPV") || SQ_props[prop_idx].Contains("Chi2"))) {
                    continue;
                }

                TString histName = Form("%s_%s_AntiSexaquark_%s", stage.Data(), set.Data(), SQ_props[prop_idx].Data());
                std::tuple<TString, TString, TString> histKey = std::make_tuple(stage, set, SQ_props[prop_idx]);
                fHist_AntiSexaquarks[histKey] = new TH1D(histName, "", SQ_nbins[prop_idx], SQ_min[prop_idx], SQ_max[prop_idx]);

                fList_Sexaquarks_Hists->Add(fHist_AntiSexaquarks[histKey]);
            }
        }
    }
}

/*
 * Define and add histograms for the K+K+ pairs to the histograms output list.
 * - Uses: `fSourceOfV0s`, `fOutputListOfHists`
 */
void AliAnalysisTaskSexaquark::PreparePosKaonPairHistograms() {

    TString stages[2] = {"Findable", "Found"};
    TString sets[2] = {"All", "Signal"};

    // clang-format off
    std::vector<TString> props = {// true + rec. vars.
                                  "Mass", "OriginRadius", "DecayLength", "Xv", "Yv", "Zv", "Eta", "Rapidity", "Px", "Py", "Pz", "Pt", "OpeningAngle",
                                  // rec. vars.
                                  "DCAwrtPV", "DCAbtwKaons", "DCAkaSV", "DCAkbSV", "Chi2ndf"};
    std::vector<Int_t> nbins = {// true + rec. vars.
                                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                                 // rec. vars.
                                100, 100, 100, 100, 150};
    std::vector<Double_t> min = {// true + rec. vars
                                 0., 0., 0., -50., -300., -300., -4., -1., -10., -10., -10., 0., 0.,
                                 // rec. vars.
                                 0., 0., 0., 0., 0.};
    std::vector<Double_t> max = {// true + rec. vars
                                 10., 250., 500., 50., 300., 300., 4., 1., 10., 10., 10., 15., TMath::Pi(),
                                 // rec. vars.
                                 200., 1., 1., 1., 1.5};
    // clang-format on

    fHist_PosKaonPairs_Bookkeep = new TH1D("PosKaonPairs_Bookkeep", "", 100, 0., 100.);
    fList_PosKaonPairs_Hists->Add(fHist_PosKaonPairs_Bookkeep);

    for (TString& stage : stages) {
        for (TString& set : sets) {
            for (Int_t prop_idx = 0; prop_idx < (Int_t)props.size(); prop_idx++) {

                if (!fIsMC && (stage == "Findable" || set == "Signal")) continue;

                if (stage == "Findable" && set == "Signal") continue;

                if (stage == "Findable" && prop_idx > 12) continue;

                TString histName = Form("%s_%s_PosKaonPair_%s", stage.Data(), set.Data(), props[prop_idx].Data());
                std::tuple<TString, TString, TString> histKey = std::make_tuple(stage, set, props[prop_idx]);
                fHist_PosKaonPairs[histKey] = new TH1D(histName, "", nbins[prop_idx], min[prop_idx], max[prop_idx]);
                fList_PosKaonPairs_Hists->Add(fHist_PosKaonPairs[histKey]);
            }
        }
    }
}

/*
 * Main function, called per each event at RUNTIME ~ execution on Grid
 * - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fPrimaryVertex`, `fMagneticField`, `fSourceOfV0s`, `fOutputListOfHists`
 */
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {

    /* Handle external logs */

    if (fIsFirstEvent) {
        if (fAliEnPath == "") AliInfo("!! No luck finding fAliEnPath !!");
        AliInfoF("!! fAliEnPath: %s !!", fAliEnPath.Data());
        if (fReadSignalLogs) LoadLogsIntoTree();
        fIsFirstEvent = kFALSE;
    }

    /* Load MC Gen. Event and PV */

    if (fIsMC) {
        fMC = MCEvent();
        if (!fMC) AliFatal("ERROR: AliMCEvent couldn't be found.");
        fMC_PrimaryVertex = const_cast<AliVVertex*>(fMC->GetPrimaryVertex());
    }

    /* Load Reconstructed Event, PV and Magnetic Field */

    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) AliFatal("ERROR: AliESDEvent couldn't be found.");
    fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());

    fMagneticField = fESD->GetMagneticField();

    fRunNumber = fESD->GetRunNumber();
    fEventNumber = fESD->GetEventNumberInFile();

    // print event centrality class
    AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*>(fESD->FindListObject("MultSelection"));
    if (!MultSelection) AliFatal("ERROR: AliMultSelection couldn't be found.");
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    // AliInfoF("!! fCentrality = %f !!", fCentrality); // DEBUG

    // init KFparticle
    KFParticle::SetField(fMagneticField);

    /* Assign branches and fill tree */

    if (fMC_PrimaryVertex) {
        tEvents_PV_TrueXv = fMC_PrimaryVertex->GetX();
        tEvents_PV_TrueYv = fMC_PrimaryVertex->GetY();
        tEvents_PV_TrueZv = fMC_PrimaryVertex->GetZ();
    }
    tEvents_PV_RecXv = fPrimaryVertex->GetX();
    tEvents_PV_RecYv = fPrimaryVertex->GetY();
    tEvents_PV_RecZv = fPrimaryVertex->GetZ();
    tEvents_NTracks = fESD->GetNumberOfTracks();
    tEvents_IsMB = fInputHandler->IsEventSelected() & AliVEvent::kINT7;
    tEvents_IsHighMultV0 = fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0;
    tEvents_IsHighMultSPD = fInputHandler->IsEventSelected() & AliVEvent::AliVEvent::kHighMultSPD;
    tEvents_IsCentral = fInputHandler->IsEventSelected() & AliVEvent::kCentral;
    tEvents_IsSemiCentral = fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral;

    fTree_Events->Fill();

    /* MC */

    if (fIsMC) {
        ProcessMCGen();
        ProcessSignalInteractions();
    }

    if (fStopAfter == "MC") {
        EndOfEvent();
        return;
    }

    /* Event */

    if (!PassesEventSelection()) {
        EndOfEvent();
        return;
    }

    /* Tracks */

    ProcessTracks();

    if (fStopAfter == "Tracks") {
        EndOfEvent();
        return;
    }

    /* Findables */

    if (fIsMC) {
        if (fReactionID != 'H') {
            ProcessFindableV0s();
            if (fReactionID == 'E') ProcessFindablePionPairs();
            ProcessFindableSexaquarks();
        } else {
            ProcessFindablePosKaonPairs();
        }
    }

    if (fStopAfter == "Findables") {
        EndOfEvent();
        return;
    }

    /* V0s */

    if (fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "offline") {
        GetV0sFromESD(fSourceOfV0s == "on-the-fly");
    }

    if (fSourceOfV0s == "custom") {
        if (fReactionID == 'A') {  // "AntiSexaquark,N->AntiLambda,K0S"
            CustomV0Finder(-3122);
            CustomV0Finder(310);
        }
    }

    if (fSourceOfV0s == "kalman") {
        if (fReactionID == 'A') {  // "AntiSexaquark,N->AntiLambda,K0S"
            KalmanV0Finder(-3122, -2212, 211);
            KalmanV0Finder(310, -211, 211);
        } else if (fReactionID == 'D') {  // "AntiSexaquark,P->AntiLambda,K+"
            KalmanV0Finder(-3122, -2212, 211);
        } else if (fReactionID == 'E') {  // "AntiSexaquark,P->AntiLambda,K+,pi-,pi+"
            KalmanV0Finder(-3122, -2212, 211);
            KalmanV0Finder(422, -211, 211);
        } else if (fReactionID == 'H') {  // "AntiSexaquark,P->AntiProton,K+,K+,pi0"
            KalmanPosKaonPairFinder();
        }
    }

    if (fStopAfter == "V0s") {
        EndOfEvent();
        return;
    }

    /* Sexaquarks */

    // clang-format off
    if (fSourceOfV0s != "kalman") {
        if (fReactionID == 'A') GeoSexaquarkFinder_ChannelA();
        else if (fReactionID == 'D') GeoSexaquarkFinder_ChannelD();
        else if (fReactionID == 'E') GeoSexaquarkFinder_ChannelE();
    } else {
        if (fReactionID == 'A') KalmanSexaquarkFinder_ChannelA();
        else if (fReactionID == 'D') KalmanSexaquarkFinder_ChannelD();
        else if (fReactionID == 'E') KalmanSexaquarkFinder_ChannelE();
    }
    // clang-format on

    EndOfEvent();
}

/*
 * End of event.
 */
void AliAnalysisTaskSexaquark::EndOfEvent() {

    ClearEventsBranches();
    ClearContainers();

    PostData(1, fTree_Events);
    PostData(2, fTree_Injected);
    PostData(3, fTree_MC);
    PostData(4, fTree_Tracks);
    PostData(5, fTree_V0s);
    PostData(6, fTree_Sexaquarks);
    PostData(7, fTree_KaonPairs);

    PostData(8, fList_QA_Hists);
    PostData(9, fList_Tracks_Hists);
    PostData(10, fList_V0s_Hists);
    PostData(11, fList_Sexaquarks_Hists);
    PostData(12, fList_PosKaonPairs_Hists);
}

/*
 * User implementation of Notify(). Needed for reading the AliEn path.
 * This function is loaded during AliAnalysisManager::Notify().
 * It's called after UserCreateOutputObjects(), for each new file, and before the first UserExec().
 */
Bool_t AliAnalysisTaskSexaquark::UserNotify() {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) AliFatal("!! Analysis Manager not found !!");
    TTree* man_tree = man->GetTree();
    if (!man_tree) AliFatal("!! Analysis Manager Tree not found !!");
    TFile* man_file = man_tree->GetCurrentFile();
    if (!man_file) AliFatal("!! Analysis Manager File not found !!");
    fAliEnPath = man_file->GetName();
    AliInfoF("!! Loaded AliEn Path: %s !!", fAliEnPath.Data());

    TObjArray* tokens = fAliEnPath.Tokenize("/");
    fDirNumber = ((TObjString*)tokens->At(tokens->GetEntries() - 2))->GetString().Atoi();
    AliInfoF("!! Dir Number: %i !!", fDirNumber);

    fIsFirstEvent = kTRUE;

    return kTRUE;
}

/*
 * Initialize analysis task.
 */
void AliAnalysisTaskSexaquark::Initialize() {

    std::array<TString, 4> validSourcesOfV0s = {"offline", "on-the-fly", "custom", "kalman"};

    if (std::find(validSourcesOfV0s.begin(), validSourcesOfV0s.end(), fSourceOfV0s) == validSourcesOfV0s.end()) {
        AliWarning("Source of V0s must be \"offline\", \"on-the-fly\", \"custom\" or \"kalman\".");
        AliWarning("=> Setting it to \"kalman\".");
        fSourceOfV0s = "kalman";
    }

    if (fReactionID == 'H' && (fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "offline")) {
        AliWarning("Source of V0s must be \"custom\" or \"kalman\" for the reaction channel H.");
        AliWarning("=> Setting it to \"kalman\".");
        fSourceOfV0s = "kalman";
    }

    std::array<Char_t, 4> validReactionIDs = {'A', 'D', 'E', 'H'};
    if (std::find(validReactionIDs.begin(), validReactionIDs.end(), fReactionID) == validReactionIDs.end()) {
        AliWarning("Reaction ID must be set.");
        AliWarning("=> Setting it to 'A'.");
        fReactionID = 'A';
    }

    /* Define cuts */

    DefineTracksCuts("standard");
    // clang-format off
    if (fReactionID != 'H') DefineV0Cuts("standard");
    else DefinePosKaonPairCuts("standard");
    // clang-format on
    DefineSexaquarkCuts("standard");

    getNegPdgCode_fromV0PdgCode[-3122] = -2212;
    getPosPdgCode_fromV0PdgCode[-3122] = 211;
    getNegPdgCode_fromV0PdgCode[310] = -211;
    getPosPdgCode_fromV0PdgCode[310] = 211;

    /* Assign the rest of properties related to the respective Reaction Channel */

    TString ReactionChannel;
    if (fReactionID == 'A') {
        ReactionChannel = "AntiSexaquark,N->AntiLambda,K0S";
        fProductsPDG = {-3122, 310};
        fFinalStateProductsPDG = {-2212, 211, -211, 211};
        fStruckNucleonPDG = 2112;
    } else if (fReactionID == 'D') {
        ReactionChannel = "AntiSexaquark,P->AntiLambda,K+";
        fProductsPDG = {-3122, 321};
        fFinalStateProductsPDG = {-2212, 211, 321};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'E') {
        ReactionChannel = "AntiSexaquark,P->AntiLambda,K+,pi-,pi+";
        fProductsPDG = {-3122, 321, -211, 211};
        fFinalStateProductsPDG = {-2212, 211, 321, -211, 211};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'H') {
        ReactionChannel = "AntiSexaquark,P->AntiProton,K+,K+,pi0";
        fProductsPDG = {-2212, 321, 321, 111};
        fFinalStateProductsPDG = {321, 321};  // focused only on K+K+
        fStruckNucleonPDG = 2212;
    }

    AliInfo("Initializing...");
    AliInfo("Settings:");
    AliInfo("========");
    AliInfoF(">> IsMC                = %i", (Int_t)fIsMC);
    AliInfoF(">> Source of V0s       = \"%s\"", fSourceOfV0s.Data());
    AliInfoF(">> Reaction ID         = \'%c\'", fReactionID);
    AliInfoF(">> Reaction Channel    = \"%s\"", ReactionChannel.Data());
    AliInfoF(">> Stop After          = \"%s\"", fStopAfter.Data());
    AliInfoF(">> DoQA                = %i", (Int_t)fDoQA);
    AliInfoF(">> Read Signal Logs    = %i", (Int_t)fReadSignalLogs);
}

/*
 * Define track selection cuts.
 * - Input: `cuts_option`
 */
void AliAnalysisTaskSexaquark::DefineTracksCuts(TString cuts_option) {

    kMax_NSigma_Pion = 3.;
    kMax_NSigma_Kaon = 3.;
    kMax_NSigma_Proton = 3.;

    kMax_Track_Eta = 1.;
    kMin_Track_NTPCClusters = 50;
    kMax_Track_Chi2PerNTPCClusters = 2.;
    kTurnedOn_Track_StatusCuts = kTRUE;
    kTurnedOn_Track_RejectKinks = kFALSE;  // PENDING: need to study further
    kMin_Track_DCAxy_wrtPV = 2.;           // PENDING: could be improved for pions with DCAxy > 4 ??

    kMin_Track_Pt[2212] = 0.2;
    kMin_Track_Pt[321] = 0.0;  // PENDING: to check when analyzing the next reaction channels
    kMin_Track_Pt[211] = 0.2;
}

/*
 * Define V0 selection cuts.
 * - Uses: `fSourceOfV0s`
 * - Input: `cuts_option`
 */
void AliAnalysisTaskSexaquark::DefineV0Cuts(TString cuts_option) {

    /* Common */

    kMin_V0_Mass[-3122] = 1.08;
    kMax_V0_Mass[-3122] = 1.16;
    kMin_V0_Pt[-3122] = 1.0;

    kMin_V0_Mass[310] = 0.475;
    kMax_V0_Mass[310] = 0.525;
    kMin_V0_Pt[310] = 1.0;

    kMin_V0_Mass[422] = 2 * fPDG.GetParticle(211)->Mass();
    kMax_V0_Mass[422] = 2.;
    kMin_V0_Pt[422] = 0.1;

    if (fSourceOfV0s == "offline") {

        // none so far
    }

    if (fSourceOfV0s == "on-the-fly") {

        kMax_V0_ImprvDCAbtwDau[-3122] = 0.6;

        kMax_V0_ImprvDCAbtwDau[310] = 0.6;
    }

    if (fSourceOfV0s == "custom") {

        kMax_V0_CPAwrtPV[-3122] = 0.99;
        kMin_V0_DCAwrtPV[-3122] = 4.;
        kMax_V0_Chi2[-3122] = 0.3;

        kMax_V0_Eta[310] = 0.6;
        kMax_V0_CPAwrtPV[310] = 0.99;
        kMax_V0_DCAposV0[310] = 0.15;
        kMin_V0_DCAwrtPV[310] = 4.;
        kMax_V0_Chi2[310] = 0.3;
        kMax_V0_ImprvDCAbtwDau[310] = 0.3;
    }

    if (fSourceOfV0s == "kalman") {

        kMax_V0_Eta[-3122] = 0.9;
        kMin_V0_Radius[-3122] = 20.;
        kMin_V0_DecayLength[-3122] = 40.;
        kMin_V0_CPAwrtPV[-3122] = 0.1;
        kMax_V0_CPAwrtPV[-3122] = 0.99;
        kMin_V0_DCAwrtPV[-3122] = 4.;
        kMax_V0_DCAbtwDau[-3122] = 2.;
        kMax_V0_DCAnegV0[-3122] = 2.;
        kMax_V0_DCAposV0[-3122] = 2.;
        kMax_V0_ArmPtOverAlpha[-3122] = 0.2;
        // kMax_V0_Chi2ndf[-3122] = 0.4;

        kMax_V0_Eta[310] = 0.8;
        kMin_V0_Radius[310] = 15.;
        kMin_V0_DecayLength[310] = 30.;
        kMax_V0_DecayLength[310] = 175.;
        kMin_V0_CPAwrtPV[310] = 0.2;
        kMax_V0_CPAwrtPV[310] = 0.99;
        kMin_V0_DCAwrtPV[310] = 4.;
        kMax_V0_DCAbtwDau[310] = 0.25;
        kMax_V0_DCAnegV0[310] = 0.25;
        kMax_V0_DCAposV0[310] = 0.25;
        // kMax_V0_Chi2ndf[310] = 0.4;

        kMax_V0_Eta[422] = 0.9;
        kMin_V0_Radius[422] = 15.;
        kMin_V0_CPAwrtPV[422] = 0.8;
        kMax_V0_CPAwrtPV[422] = 0.98;
        kMin_V0_DCAwrtPV[422] = 4.;
        kMax_V0_DCAbtwDau[422] = 0.3;
        kMax_V0_DCAnegV0[422] = 0.2;
        kMax_V0_DCAposV0[422] = 0.2;
        // kMax_V0_Chi2ndf[422] = 0.1;
    }
}

/*
 * Define anti-sexaquark candidates selection cuts.
 * - Uses: `fSourceOfV0s`
 * - Input: `cuts_option`
 */
void AliAnalysisTaskSexaquark::DefineSexaquarkCuts(TString cuts_option) {

    /* Common */

    // none so far

    if (fSourceOfV0s == "offline") {
        // none so far
    }
    if (fSourceOfV0s == "on-the-fly") {
        // none so far
    }
    if (fSourceOfV0s == "custom") {
        kMin_Sexa_Radius = 40;
        kMin_Sexa_CPAwrtPV = 0.99;
        kMax_Sexa_CPAwrtPV = 1.;
        kMax_Sexa_DCAv0aSV = 0.1;
        kMax_Sexa_DCAv0bSV = 0.1;
    }
    if (fSourceOfV0s == "kalman") {
        if (fReactionID == 'A') {
            kMin_Sexa_Radius = 40;
            kMin_Sexa_CPAwrtPV = 0.99;
            kMax_Sexa_CPAwrtPV = 1.;
            kMax_Sexa_DCAbtwV0s = 2.;
            kMax_Sexa_DCAv0aSV = 2.;
            kMax_Sexa_DCAv0bSV = 2.;
            kMax_Sexa_DCAv0anegSV = 10.;
            kMax_Sexa_DCAv0aposSV = 10.;
            kMax_Sexa_DCAv0bnegSV = 10.;
            kMax_Sexa_DCAv0bposSV = 10.;
            // kMax_Sexa_Chi2ndf = 0.5;
        } else if (fReactionID == 'D') {
            kMin_Sexa_Radius = 30;
            kMax_Sexa_Rapidity = 0.8;
            kMin_Sexa_CPAwrtPV = 0.99;
            kMax_Sexa_CPAwrtPV = 1.;
            kMax_Sexa_DCAv0ba = 0.3;
            kMax_Sexa_DCAbaSV = 0.2;
            kMax_Sexa_DCAv0SV = 0.2;
            // kMax_Sexa_Chi2ndf = 1.;
        } else if (fReactionID == 'E') {
            kMin_Sexa_Radius = 30;
            kMin_Sexa_CPAwrtPV = 0.99;
            kMax_Sexa_CPAwrtPV = 1.;
            kMax_Sexa_DCAv0SV = 0.1;
            kMax_Sexa_DCApkSV = 0.1;
            kMax_Sexa_DCApmSV = 0.1;
            kMax_Sexa_DCAppSV = 0.1;
            kMax_Sexa_DCAv0pk = 0.1;
            kMax_Sexa_DCAv0pm = 0.1;
            kMax_Sexa_DCAv0pp = 0.1;
            kMax_Sexa_DCApkpm = 0.1;
            kMax_Sexa_DCApkpp = 0.1;
            // kMax_Sexa_Chi2ndf = 1.;
        }
    }
}

/*
 * Define K+K+ pair selection cuts.
 * - Uses: `fSourceOfV0s`
 * - Input: `cuts_option`
 */
void AliAnalysisTaskSexaquark::DefinePosKaonPairCuts(TString cuts_option) {

    if (fSourceOfV0s == "kalman") {
        kMin_PosKaonPair_Radius = 20;
        kMin_PosKaonPair_Pt = 0.2;
        kMax_PosKaonPair_DCAbtwKaons = 0.2;
        kMax_PosKaonPair_DCApkaSV = 0.2;
        kMax_PosKaonPair_DCApkbSV = 0.2;
        // kMax_PosKaonPair_Chi2ndf = 0.1;
    }
}

/*                  */
/**  MC Generated  **/
/*** ============ ***/

/*
 * Loop over MC particles in a single event. Store the indices of the signal particles.
 * - Uses: `fMC`, `fPDG`, `fMC_PrimaryVertex`
 */
void AliAnalysisTaskSexaquark::ProcessMCGen() {

    Int_t mcIdx;
    AliMCParticle* mcPart;
    Int_t pdg_mc;
    Int_t mother_idx;

    TString generator_name;

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    Int_t idx_first_dau;
    Int_t idx_last_dau;
    AliMCParticle* mcNegDau;
    AliMCParticle* mcPosDau;

    Int_t n_daughters = 0;
    Int_t mcIdxNegDaughter = 0;
    Int_t mcIdxPosDaughter = 0;

    TVector3 originVertex;
    TVector3 decayVertex;

    Float_t decay_length;
    Float_t opening_angle;
    Float_t dca_wrt_pv, cpa_wrt_pv;
    Float_t armenteros_alpha, armenteros_qt;
    Bool_t reconstructable_v0;

    if (fDoQA) {
        fHist_QA["MCGen_VertexX"]->Fill(fMC_PrimaryVertex->GetX());
        fHist_QA["MCGen_VertexY"]->Fill(fMC_PrimaryVertex->GetY());
        fHist_QA["MCGen_VertexZ"]->Fill(fMC_PrimaryVertex->GetZ());
    }

    /* Loop over MC gen. particles in this event */

    for (mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        /* Protection: exclude particles from the anti-neutron injector */

        // `GetCocktailGenerator` not only considers primaries, but subsequent daughters/secondaries all the way down, as well
        fMC->GetCocktailGenerator(mcIdx, generator_name);
        if (generator_name.Contains("Injector (Anti-Neutron)")) continue;

        /* Get MC particle info */

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);

        pdg_mc = mcPart->PdgCode();
        getPdgCode_fromMcIdx[mcIdx] = pdg_mc;

        mother_idx = mcPart->GetMother();
        n_daughters = mcPart->GetNDaughters();
        idx_first_dau = mcPart->GetDaughterFirst();
        idx_last_dau = mcPart->GetDaughterLast();

        originVertex.SetXYZ(mcPart->Xv(), mcPart->Yv(), mcPart->Zv());

        /* Check if particle is signal */
        /* - It should be produced by the "Anti-Sexaquark Interaction (Channel %c)" generator */
        /* - Then, the unique reaction ID is provided by MCStatusCode of the injected signal particles */

        getReactionIdx_fromMcIdx[mcIdx] = 0;
        isMcIdxSignal[mcIdx] = generator_name.Contains("Anti-Sexaquark Interaction");
        if (isMcIdxSignal[mcIdx]) {
            getReactionIdx_fromMcIdx[mcIdx] = mother_idx < 0 ? mcPart->MCStatusCode() : getReactionIdx_fromMcIdx[mother_idx];
            getMcIdx_fromReactionIdx[getReactionIdx_fromMcIdx[mcIdx]].push_back(mcIdx);
        }

        /* Check if particle is secondary */
        /* - It includes signal particles too, because despite they were injected as primaries, they are really not */
        isMcIdxSecondary[mcIdx] = mcPart->IsSecondaryFromMaterial() || mcPart->IsSecondaryFromWeakDecay() || isMcIdxSignal[mcIdx];

        /* Assign branches */

        tMC_Idx = mcIdx;
        tMC_PdgCode = pdg_mc;
        tMC_Idx_Mother = mother_idx;
        tMC_NDaughters = n_daughters;
        tMC_Idx_FirstDau = idx_first_dau;
        tMC_Idx_LastDau = idx_last_dau;
        tMC_Px = (Float_t)mcPart->Px();
        tMC_Py = (Float_t)mcPart->Py();
        tMC_Pz = (Float_t)mcPart->Pz();
        tMC_Xv_i = (Float_t)originVertex.X();
        tMC_Yv_i = (Float_t)originVertex.Y();
        tMC_Zv_i = (Float_t)originVertex.Z();
        tMC_Status = mcPart->MCStatusCode();
        tMC_IsSecondary = isMcIdxSecondary[mcIdx];
        tMC_IsSignal = isMcIdxSignal[mcIdx];
        tMC_ReactionID = getReactionIdx_fromMcIdx[mcIdx];

        /** Charged particles: anti-protons, pos. kaons, pions **/

        if (std::find(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end(), pdg_mc) != fFinalStateProductsPDG.end()) {

            /* Fill MC tree */

            fTree_MC->Fill();
            ClearMCParticlesBranches();

            /* Fill histograms */

            fHist_Tracks_Bookkeep[pdg_mc]->Fill(0);
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Px")]->Fill(mcPart->Px());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Py")]->Fill(mcPart->Py());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());

            if (isMcIdxSecondary[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(10);
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Px")]->Fill(mcPart->Px());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Py")]->Fill(mcPart->Py());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }

            if (isMcIdxSignal[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(20);
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Px")]->Fill(mcPart->Px());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Py")]->Fill(mcPart->Py());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }
        }  // end of anti-protons, kaons, pions condition

        /** Neutral particles: anti-lambdas and K0S **/

        if (fReactionID != 'H' && (pdg_mc == -3122 || (fReactionID == 'A' && pdg_mc == 310))) {

            /* Fill MC tree */

            fTree_MC->Fill();
            ClearMCParticlesBranches();

            /* Fill histograms (reconstructable or not) */

            fHist_V0s_Bookkeep[pdg_mc]->Fill(0);
            if (isMcIdxSecondary[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(10);
            if (isMcIdxSignal[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(20);

            /* Check if this V0 can be reconstructed */

            // conditions:
            // (1) has daughters
            // (2) and they correspond to:
            //     -- anti-lambdas -> anti-proton pi+
            //     -- K0S -> pi+ pi-

            n_daughters = mcPart->GetNDaughters();
            mcIdxNegDaughter = 0;
            mcIdxPosDaughter = 0;

            // loop over daughters and fill mcIdxNegDaughter and mcIdxPosDaughter in case of relevant decays
            if (n_daughters) {
                for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {

                    mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
                    pdg_dau = mcDaughter->PdgCode();

                    if (pdg_mc == -3122) {
                        if (pdg_dau == -2212) mcIdxNegDaughter = mcIdxDaughter;
                        if (pdg_dau == 211) mcIdxPosDaughter = mcIdxDaughter;
                    }
                    if (pdg_mc == 310) {
                        if (pdg_dau == -211) mcIdxNegDaughter = mcIdxDaughter;
                        if (pdg_dau == 211) mcIdxPosDaughter = mcIdxDaughter;
                    }
                }
            }  // end of n_daughters condition

            /* In any case, get the decay vertex of the V0, which is the origin vertex of its daughters */

            reconstructable_v0 = n_daughters && mcIdxNegDaughter && mcIdxPosDaughter;
            if (!reconstructable_v0) continue;

            /* Store true V0 and mother-daughter info */

            mcIndicesOfTrueV0s.push_back(mcIdx);

            doesMcIdxHaveMother[mcIdxNegDaughter] = kTRUE;
            getMotherMcIdx_fromMcIdx[mcIdxNegDaughter] = mcIdx;
            getNegDauMcIdx_fromMcIdx[mcIdx] = mcIdxNegDaughter;

            doesMcIdxHaveMother[mcIdxPosDaughter] = kTRUE;
            getMotherMcIdx_fromMcIdx[mcIdxPosDaughter] = mcIdx;
            getPosDauMcIdx_fromMcIdx[mcIdx] = mcIdxPosDaughter;

            /* Get particle's properties */

            mcNegDau = (AliMCParticle*)fMC->GetTrack(mcIdxNegDaughter);
            mcPosDau = (AliMCParticle*)fMC->GetTrack(mcIdxPosDaughter);

            armenteros_qt = ArmenterosQt(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz());
            armenteros_alpha = ArmenterosAlpha(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz(),
                                               mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz());

            decayVertex.SetXYZ(mcNegDau->Xv(), mcNegDau->Yv(), mcNegDau->Zv());
            cpa_wrt_pv = CosinePointingAngle(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                             fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                      fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

            originVertex.SetXYZ(mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
            decay_length = (originVertex - decayVertex).Mag();

            TVector3 neg(mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz());
            TVector3 pos(mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz());
            opening_angle = neg.Angle(pos);

            /* Fill histograms (reconstructables only) */

            fHist_V0s_Bookkeep[pdg_mc]->Fill(1);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Mass")]->Fill(mcPart->M());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Px")]->Fill(mcPart->Px());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Py")]->Fill(mcPart->Py());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Phi")]->Fill(mcPart->Phi());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "OriginRadius")]->Fill(originVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DecayRadius")]->Fill(decayVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Zv")]->Fill(decayVertex.Z());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Rapidity")]->Fill(mcPart->Y());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (!isMcIdxSecondary[mcIdx]) continue;

            fHist_V0s_Bookkeep[pdg_mc]->Fill(11);
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Mass")]->Fill(mcPart->M());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Px")]->Fill(mcPart->Px());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Py")]->Fill(mcPart->Py());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Phi")]->Fill(mcPart->Phi());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "OriginRadius")]->Fill(originVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DecayRadius")]->Fill(decayVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Zv")]->Fill(decayVertex.Z());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Rapidity")]->Fill(mcPart->Y());
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (!isMcIdxSignal[mcIdx]) continue;

            fHist_V0s_Bookkeep[pdg_mc]->Fill(21);
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Mass")]->Fill(mcPart->M());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Px")]->Fill(mcPart->Px());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Py")]->Fill(mcPart->Py());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Phi")]->Fill(mcPart->Phi());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "OriginRadius")]->Fill(originVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DecayRadius")]->Fill(decayVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Zv")]->Fill(decayVertex.Z());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Rapidity")]->Fill(mcPart->Y());
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
        }  // end of anti-lambda && K0s condition
    }      // end of loop over MC particles
}

/*
 * Loop over generated antisexaquark-nucleon reactions.
 * - Uses: `fMC`, `fPDG`, `fStruckNucleonPDG`
 * - Containers: `getMcIdx_fromReactionIdx`
 */
void AliAnalysisTaskSexaquark::ProcessSignalInteractions() {

    AliMCParticle* mcProduct;

    /** Declare TLorentzVectors **/

    TLorentzVector lvProduct;
    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /** Declare variables **/

    TVector3 secondary_vertex(0., 0., 0.);
    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over reactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /** Protection in case of general-purpose simulations **/

        if (!getMcIdx_fromReactionIdx[reactionIdx].size()) continue;

        /** Reset vectors **/

        secondary_vertex.SetXYZ(0., 0., 0.);
        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        /** Loop over reaction products **/

        for (Int_t& mcIdx : getMcIdx_fromReactionIdx[reactionIdx]) {

            mcProduct = (AliMCParticle*)fMC->GetTrack(mcIdx);

            /* Get first generation products, the direct daughters of the reaction */

            if (mcProduct->MCStatusCode() != reactionIdx) continue;

            lvProduct.SetPxPyPzE(mcProduct->Px(), mcProduct->Py(), mcProduct->Pz(), mcProduct->E());

            lvAntiSexaquark = lvAntiSexaquark + lvProduct;

            /* Get sec. vertex (if it hasn't been defined yet)*/

            if (secondary_vertex.Mag() < 1E-6) secondary_vertex.SetXYZ(mcProduct->Xv(), mcProduct->Yv(), mcProduct->Zv());
        }

        lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest
        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /** Get properties **/

        radius = secondary_vertex.Perp();
        getRadius_fromReactionIdx[reactionIdx] = radius;

        decay_length = TMath::Sqrt(TMath::Power(secondary_vertex.X() - fMC_PrimaryVertex->GetX(), 2) +
                                   TMath::Power(secondary_vertex.Y() - fMC_PrimaryVertex->GetY(), 2) +
                                   TMath::Power(secondary_vertex.Z() - fMC_PrimaryVertex->GetZ(), 2));
        cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), secondary_vertex.X(), secondary_vertex.Y(),
                                  secondary_vertex.Z(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /** Fill histograms **/

        /* With information A.I. -- After Interaction  */

        if (fReactionID == 'H') fHist_PosKaonPairs_Bookkeep->Fill(0);

        fHist_AntiSexaquarks_Bookkeep->Fill(0);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Mass")]->Fill(lvAntiSexaquark.M());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Px")]->Fill(lvAntiSexaquark.Px());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Py")]->Fill(lvAntiSexaquark.Py());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Zv")]->Fill(secondary_vertex.Z());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Radius")]->Fill(radius);  // interaction radius
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "DecayLength")]->Fill(decay_length);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
    }  // end of loop over reactions
}

/*                   */
/**  Reconstructed  **/
/*** ============= ***/

/*
 * Apply event selection
 */
Bool_t AliAnalysisTaskSexaquark::PassesEventSelection() {

    if (fDoQA) fHist_QA["Events_Bookkeep"]->Fill(0);

    // ! REFERENCE: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGRunList18r ! //
    if (!fIsMC && (fRunNumber == 296749 || fRunNumber == 296750 || fRunNumber == 296849 || fRunNumber == 296890 || fRunNumber == 297029 ||
                   fRunNumber == 297194 || fRunNumber == 297219 || fRunNumber == 297481)) {
        fEventCuts.UseTimeRangeCut();
        fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny);  // PENDING: remove kAny!
        if (!fEventCuts.AcceptEvent(fESD)) return kFALSE;
        AliInfo("!! MC Event Passed z-vertex range on PV !!");  // DEBUG
    }

    /* Trigger Selection */

    Bool_t IsMB = kFALSE;
    Bool_t IsHighMultV0 = kFALSE;
    Bool_t IsHighMultSPD = kFALSE;
    Bool_t IsCentral = kFALSE;
    Bool_t IsSemiCentral = kFALSE;

    if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) IsMB = kTRUE;
    if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0) IsHighMultV0 = kTRUE;
    if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) IsHighMultSPD = kTRUE;
    if (fInputHandler->IsEventSelected() & AliVEvent::kCentral) IsCentral = kTRUE;
    if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) IsSemiCentral = kTRUE;

    if (!IsMB && !IsHighMultV0 && !IsHighMultSPD && !IsCentral && !IsSemiCentral) {
        AliInfoF("!! Event Rejected -- %u %u %u %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, AliVEvent::kCentral,
                 AliVEvent::kSemiCentral);  // DEBUG
        return kFALSE;
    }
    AliInfoF("!! Event Selected -- %u %u %u %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, AliVEvent::kCentral,
             AliVEvent::kSemiCentral);  // DEBUG

    if (fDoQA) fHist_QA["Events_Bookkeep"]->Fill(1);

    /* rec. PV z-vertex range */

    if (TMath::Abs(fPrimaryVertex->GetZ()) > 10.) return kFALSE;
    AliInfo("!! Event Passed z-vertex range on PV !!");  // DEBUG

    if (fDoQA) {
        fHist_QA["Events_Bookkeep"]->Fill(2);

        /* Vertex */

        fHist_QA["VertexX"]->Fill(fPrimaryVertex->GetX());
        fHist_QA["VertexY"]->Fill(fPrimaryVertex->GetY());
        fHist_QA["VertexZ"]->Fill(fPrimaryVertex->GetZ());

        /* N Tracks */

        Int_t ntracks = fESD->GetNumberOfTracks();
        fHist_QA["NTracks"]->Fill(ntracks);

        /* SPD */

        const AliESDVertex* spd_vertex = fESD->GetPrimaryVertexSPD();
        const AliMultiplicity* mult = fESD->GetMultiplicity();

        fHist_QA["SPD_NTracklets"]->Fill(mult->GetNumberOfTracklets());

        f2DHist_QA["SPD_MultiplicityVsVertexZ"]->Fill(spd_vertex->GetZ(), mult->GetNumberOfTracklets());
        f2DHist_QA["SPD_NTrackletsClu0VsNTracklets"]->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(0));
        f2DHist_QA["SPD_NTrackletsClu1VsNTracklets"]->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(1));

        for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets(); iTracklet++) {
            Float_t phiTr = mult->GetPhi(iTracklet);
            Float_t etaTr = mult->GetEta(iTracklet);
            f2DHist_QA["SPD_PhiVsEta"]->Fill(etaTr, phiTr);
        }

        /* TPC */

        fHist_QA["TPC_NClusters"]->Fill(fESD->GetNumberOfTPCClusters());
    }

    return kTRUE;
}

/*
 * Loop over the reconstructed tracks in a single event.
 */
void AliAnalysisTaskSexaquark::ProcessTracks() {

    Int_t esdIdxTrack;
    AliESDtrack* track;

    Int_t mcIdx;
    Int_t mcPdgCode;
    Int_t motherMcIdx;

    TString generator_name;

    Float_t dca_wrt_pv, dca_xy_wrt_pv, dca_z_wrt_pv;
    Float_t n_tpc_clusters;
    Float_t n_crossed_rows;
    Float_t chi2_over_nclusters;
    Float_t n_findable_clusters, n_shared_clusters;
    Float_t chi2_over_ncls;

    Float_t n_sigma_pion;
    Float_t n_sigma_kaon;
    Float_t n_sigma_proton;

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Bool_t already_selected;  // as another species?
    Bool_t already_in_tree;
    UInt_t nSelectedTracks = 0;

    /** Loop over tracks in a single event **/

    for (esdIdxTrack = 0; esdIdxTrack < fESD->GetNumberOfTracks(); esdIdxTrack++) {

        /* Get track */

        track = fESD->GetTrack(esdIdxTrack);
        mcIdx = TMath::Abs(track->GetLabel());

        /* Sanity check */

        if (track->Pt() < 1E-3 || track->Pt() > 1E3) continue;

        /* Protection: exclude particles from the anti-neutron injector */

        // `GetCocktailGenerator` not only considers primaries, but subsequent daughters/secondaries all the way down, as well
        fMC->GetCocktailGenerator(mcIdx, generator_name);
        if (generator_name.Contains("Injector (Anti-Neutron)")) continue;

        /* Get MC info */

        getMcIdx_fromEsdIdx[esdIdxTrack] = mcIdx;
        mcPdgCode = getPdgCode_fromMcIdx[mcIdx];

        /* Fill histograms before track selection */

        f2DHist_Tracks_TPCsignal["All"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

        /* QA : Tracks */

        if (fDoQA) {
            const Int_t NITSLayers = 6;
            for (Int_t iLayer = 0; iLayer < NITSLayers; iLayer++) {
                if (track->HasPointOnITSLayer(iLayer)) {
                    f2DHist_QA["ITS_LayerNoVsPhi"]->Fill(track->Phi(), iLayer);
                    f2DHist_QA["ITS_LayerNoVsEta"]->Fill(track->Eta(), iLayer);
                }
            }
            if (TMath::Abs(track->Eta()) < kMax_Track_Eta) {
                track->GetImpactParameters(dca_xy_wrt_pv, dca_z_wrt_pv);  // pre-calculated DCA w.r.t. PV
                fHist_QA["Tracks_DCAxy"]->Fill(dca_xy_wrt_pv);
                fHist_QA["Tracks_DCAz"]->Fill(dca_z_wrt_pv);
                fHist_QA["Tracks_Pt"]->Fill(track->Pt());
                f2DHist_QA["Tracks_PhiVsEta"]->Fill(track->Eta(), track->Phi());
            }
        }  // end of QA

        /* Fill containers */

        isEsdIdxSignal[esdIdxTrack] = isMcIdxSignal[mcIdx];
        if (isEsdIdxSignal[esdIdxTrack]) {
            getReactionIdx_fromEsdIdx[esdIdxTrack] = getReactionIdx_fromMcIdx[mcIdx];
            getEsdIdx_fromReactionIdx[getReactionIdx_fromEsdIdx[esdIdxTrack]].push_back(esdIdxTrack);
        }

        doesEsdIdxHaveMother[esdIdxTrack] = doesMcIdxHaveMother[mcIdx];
        if (doesEsdIdxHaveMother[esdIdxTrack]) {
            motherMcIdx = getMotherMcIdx_fromMcIdx[mcIdx];
            getMotherMcIdx_fromEsdIdx[esdIdxTrack] = motherMcIdx;
            if (getNegDauMcIdx_fromMcIdx[motherMcIdx] == mcIdx) getNegDauEsdIdx_fromMcIdx[motherMcIdx] = esdIdxTrack;
            if (getPosDauMcIdx_fromMcIdx[motherMcIdx] == mcIdx) getPosDauEsdIdx_fromMcIdx[motherMcIdx] = esdIdxTrack;
        }

        /* Get track properties */

        n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

        track->GetImpactParameters(dca_xy_wrt_pv, dca_z_wrt_pv);  // pre-calculated DCA w.r.t. PV
        dca_xy_wrt_pv = TMath::Abs(dca_xy_wrt_pv);
        dca_z_wrt_pv = TMath::Abs(dca_z_wrt_pv);
        dca_wrt_pv = TMath::Sqrt(dca_xy_wrt_pv * dca_xy_wrt_pv + dca_z_wrt_pv * dca_z_wrt_pv);
        n_tpc_clusters = track->GetTPCNcls();  // NOTE: capital N
        n_crossed_rows = track->GetTPCCrossedRows();
        n_findable_clusters = (Int_t)track->GetTPCNclsF();
        n_shared_clusters = (Int_t)track->GetTPCnclsS();
        chi2_over_ncls = n_tpc_clusters > 1E-4 ? track->GetTPCchi2() / (Double_t)n_tpc_clusters : 999.;

        /* Loop over PID hypotheses */

        // get only the tracks that have a relevant PID for the chosen channel
        // (note: fFinalStateProductsPDG is copied into a set to remove duplicates)
        std::set<Int_t> possiblePID(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end());
        already_selected = kFALSE;
        already_in_tree = kFALSE;

        for (Int_t esdPdgCode : possiblePID) {

            /* Apply track selection */

            is_true = mcPdgCode == esdPdgCode;
            is_secondary = isMcIdxSecondary[mcIdx] && is_true;
            is_signal = isEsdIdxSignal[esdIdxTrack] && is_secondary;

            if (!PassesTrackSelectionAs(esdPdgCode, is_true, is_secondary, is_signal, track)) continue;

            if (esdPdgCode == 2212) esdIdxPassesTrackSelection_AsProton[esdIdxTrack] = kTRUE;
            if (esdPdgCode == -2212) esdIdxPassesTrackSelection_AsAntiProton[esdIdxTrack] = kTRUE;
            if (esdPdgCode == 321) esdIdxPassesTrackSelection_AsPosKaon[esdIdxTrack] = kTRUE;
            if (esdPdgCode == -321) esdIdxPassesTrackSelection_AsNegKaon[esdIdxTrack] = kTRUE;
            if (esdPdgCode == 211) esdIdxPassesTrackSelection_AsPiPlus[esdIdxTrack] = kTRUE;
            if (esdPdgCode == -211) esdIdxPassesTrackSelection_AsPiMinus[esdIdxTrack] = kTRUE;

            if (!already_selected) {
                getEsdIdx_fromMcIdx[mcIdx] = esdIdxTrack;  // purposefully done after selection
                f2DHist_Tracks_TPCsignal["Selected"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
                if (fDoQA) nSelectedTracks++;
                already_selected = kTRUE;
            }

            if (fDoQA) {
                if (esdPdgCode == -2212) {
                    f2DHist_QA["TPC_NSigmaProtonVsInnerParamP"]->Fill(track->GetInnerParam()->P(), n_sigma_proton);
                    f2DHist_QA["TPC_NSigmaProtonVsEta"]->Fill(track->GetInnerParam()->Eta(), n_sigma_proton);
                }
                if (esdPdgCode == 321) {
                    f2DHist_QA["TPC_NSigmaKaonVsInnerParamP"]->Fill(track->GetInnerParam()->P(), n_sigma_kaon);
                    f2DHist_QA["TPC_NSigmaKaonVsEta"]->Fill(track->GetInnerParam()->Eta(), n_sigma_kaon);
                }
                if (TMath::Abs(esdPdgCode) == 211) {
                    f2DHist_QA["TPC_NSigmaPionVsInnerParamP"]->Fill(track->GetInnerParam()->P(), n_sigma_pion);
                    f2DHist_QA["TPC_NSigmaPionVsEta"]->Fill(track->GetInnerParam()->Eta(), n_sigma_pion);
                }
            }

            if (esdPdgCode == -2212) esdIndicesOfAntiProtonTracks.push_back(esdIdxTrack);
            if (esdPdgCode == 321) esdIndicesOfPosKaonTracks.push_back(esdIdxTrack);
            if (esdPdgCode == -211) esdIndicesOfPiMinusTracks.push_back(esdIdxTrack);
            if (esdPdgCode == 211) esdIndicesOfPiPlusTracks.push_back(esdIdxTrack);

            /* Assign branches and fill tree */

            if (!already_in_tree) {

                tTrack_Idx = esdIdxTrack;
                tTrack_Idx_True = mcIdx;
                tTrack_Px = (Float_t)track->Px();
                tTrack_Py = (Float_t)track->Py();
                tTrack_Pz = (Float_t)track->Pz();
                tTrack_Charge = track->Charge();
                tTrack_NSigmaPion = n_sigma_pion;
                tTrack_NSigmaKaon = n_sigma_kaon;
                tTrack_NSigmaProton = n_sigma_proton;
                tTrack_IsTrue = is_true;
                tTrack_IsSecondary = is_secondary;
                tTrack_IsSignal = is_signal;
                tTrack_ReactionID = getReactionIdx_fromEsdIdx[esdIdxTrack];

                fTree_Tracks->Fill();
                ClearTracksBranches();

                already_in_tree = kTRUE;
            }

            /* Fill histograms */

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(50);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "All", esdPdgCode);

            if (!is_true) continue;

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(80);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "True", esdPdgCode);
            f2DHist_Tracks_TPCsignal["True"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

            if (!is_secondary) continue;

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(110);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "Secondary", esdPdgCode);
            f2DHist_Tracks_TPCsignal["Secondary"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

            if (!is_signal) continue;

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(140);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "Signal", esdPdgCode);
            f2DHist_Tracks_TPCsignal["Signal"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
        }  // end of loop over PID hypotheses
    }      // end of loop over tracks

    if (fDoQA) fHist_QA["NSelectedTracks"]->Fill(nSelectedTracks);
}

/*
 * Determine if current AliESDtrack passes track selection,
 * and fill the bookkeeping histograms to measure the effect of the cuts.
 * - Input: `track`
 * - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesTrackSelectionAs(Int_t HypothesisPID, Bool_t is_true, Bool_t is_secondary, Bool_t is_signal,
                                                        AliESDtrack* track) {

    AliPID::EParticleType ParticleType;
    Double_t kMax_NSigma;
    Bool_t chargeShouldBePositive = HypothesisPID > 0;  // not valid for electrons/positrons, for example
    if (TMath::Abs(HypothesisPID) == 2212) {
        ParticleType = AliPID::kProton;
        kMax_NSigma = kMax_NSigma_Proton;
    }
    if (TMath::Abs(HypothesisPID) == 321) {
        ParticleType = AliPID::kKaon;
        kMax_NSigma = kMax_NSigma_Kaon;
    }
    if (TMath::Abs(HypothesisPID) == 211) {
        ParticleType = AliPID::kPion;
        kMax_NSigma = kMax_NSigma_Pion;
    }

    Double_t charge = track->Charge();
    if (chargeShouldBePositive && charge < 0.) return kFALSE;
    if (!chargeShouldBePositive && charge > 0.) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(30);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(60);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(90);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(120);

    Double_t n_sigma = fPIDResponse->NumberOfSigmasTPC(track, ParticleType);
    if (TMath::Abs(n_sigma) > kMax_NSigma) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(31);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(61);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(91);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(121);

    Double_t Eta = TMath::Abs(track->Eta());
    if (kMax_Track_Eta && Eta > kMax_Track_Eta) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(32);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(62);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(92);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(122);

    Double_t NTPCClusters = track->GetTPCNcls();
    if (kMin_Track_NTPCClusters && NTPCClusters < kMin_Track_NTPCClusters) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(33);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(63);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(93);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(123);

    Double_t Chi2PerNTPCClusters = NTPCClusters > 1E-4 ? track->GetTPCchi2() / (Double_t)track->GetTPCNcls() : 999;
    if (kMax_Track_Chi2PerNTPCClusters && Chi2PerNTPCClusters > kMax_Track_Chi2PerNTPCClusters) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(34);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(64);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(94);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(124);

    // >> TPC and ITS status
    Bool_t tpc_status =
        (track->GetStatus() & AliESDtrack::kTPCin) && (track->GetStatus() & AliESDtrack::kTPCout) && (track->GetStatus() & AliESDtrack::kTPCrefit);
    if (kTurnedOn_Track_StatusCuts && !tpc_status) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(35);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(65);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(95);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(125);

    Bool_t its_status =
        !(track->GetStatus() & AliESDtrack::kITSin) && !(track->GetStatus() & AliESDtrack::kITSout) && !(track->GetStatus() & AliESDtrack::kITSrefit);
    if (kTurnedOn_Track_StatusCuts && !its_status) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(36);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(66);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(96);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(126);

    // >> reject kinks // !! PENDING !!
    if (kTurnedOn_Track_RejectKinks && track->GetKinkIndex(0) != 0) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(37);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(67);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(97);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(127);

    Float_t DCAxy_wrtPV, DCAz_wrtPV;
    track->GetImpactParameters(DCAxy_wrtPV, DCAz_wrtPV);
    Float_t DCA_wrtPV = TMath::Sqrt((Double_t)DCAxy_wrtPV * (Double_t)DCAxy_wrtPV + (Double_t)DCAz_wrtPV * (Double_t)DCAz_wrtPV);
    if (kMin_Track_DCA_wrtPV && DCA_wrtPV < kMin_Track_DCA_wrtPV) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(38);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(68);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(98);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(128);

    if (kMin_Track_DCAxy_wrtPV && TMath::Abs(DCAxy_wrtPV) < kMin_Track_DCAxy_wrtPV) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(39);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(69);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(99);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(129);

    if (kMin_Track_DCAz_wrtPV && TMath::Abs(DCAz_wrtPV) < kMin_Track_DCAz_wrtPV) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(40);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(70);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(100);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(130);

    if (kMin_Track_Pt[HypothesisPID] && track->Pt() < kMin_Track_Pt[HypothesisPID]) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(41);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(71);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(101);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(131);

    return kTRUE;
}

/*
 * Plot the status of the track.
 * - Input: `track`, `stage`, `esdPdgCode`
 */
void AliAnalysisTaskSexaquark::PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode) {

    ULong64_t StatusCollection[16] = {AliESDtrack::kITSin,   AliESDtrack::kITSout,     AliESDtrack::kITSrefit,    AliESDtrack::kITSpid,
                                      AliESDtrack::kTPCin,   AliESDtrack::kTPCout,     AliESDtrack::kTPCrefit,    AliESDtrack::kTPCpid,
                                      AliESDtrack::kITSupg,  AliESDtrack::kSkipFriend, AliESDtrack::kGlobalMerge, AliESDtrack::kMultInV0,
                                      AliESDtrack::kMultSec, AliESDtrack::kEmbedded,   AliESDtrack::kITSpureSA,   AliESDtrack::kESDpid};

    for (Int_t i = 0; i < 16; i++) {
        if ((track->GetStatus() & StatusCollection[i])) fHist_Tracks[std::make_tuple("Found", set, esdPdgCode, "Status")]->Fill(i);
    }
}

/*                                                                   */
/**  V0s, Anti-Sexaquarks and K+K+ Pairs -- That Require True Info  **/
/*** ============================================================= ***/

/*
 * Find all true V0s that both of their daughters passed track selection.
 * Variables are true MC vertexing & kinematics.
 * - Uses: `fMC`, `fMC_PrimaryVertex`, `fHist_V0s`, `fHist_V0s_Bookkeep`
 */
void AliAnalysisTaskSexaquark::ProcessFindableV0s() {

    AliMCParticle* mcV0;
    AliMCParticle* mcNegDaughter;
    AliMCParticle* mcPosDaughter;

    Bool_t is_signal;
    Bool_t is_secondary;

    TVector3 neg, pos;

    Int_t pdg_code;
    Float_t origin_radius, decay_radius;
    Float_t decay_length;
    Float_t opening_angle;
    Float_t arm_qt, arm_alpha;
    Float_t cpa_wrt_pv, dca_wrt_pv;

    /* Loop over true V0s */

    for (Int_t& mcIdxOfTrueV0 : mcIndicesOfTrueV0s) {

        mcV0 = (AliMCParticle*)fMC->GetTrack(mcIdxOfTrueV0);

        /* Check if both daughters were reconstructed */

        if (!getNegDauEsdIdx_fromMcIdx[mcIdxOfTrueV0] || !getPosDauEsdIdx_fromMcIdx[mcIdxOfTrueV0]) {
            continue;
        }

        is_signal = mcV0->MCStatusCode() >= 600 && mcV0->MCStatusCode() < 620;
        is_secondary = mcV0->IsSecondaryFromMaterial() || mcV0->IsSecondaryFromWeakDecay() || is_signal;

        /* Get daughters' properties */

        mcNegDaughter = (AliMCParticle*)fMC->GetTrack(getNegDauMcIdx_fromMcIdx[mcIdxOfTrueV0]);
        mcPosDaughter = (AliMCParticle*)fMC->GetTrack(getPosDauMcIdx_fromMcIdx[mcIdxOfTrueV0]);

        neg.SetXYZ(mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz());
        pos.SetXYZ(mcPosDaughter->Px(), mcPosDaughter->Py(), mcPosDaughter->Pz());

        /* Get V0 properties */

        pdg_code = mcV0->PdgCode();
        origin_radius = TMath::Sqrt(mcV0->Xv() * mcV0->Xv() + mcV0->Yv() * mcV0->Yv());
        decay_radius = TMath::Sqrt(mcNegDaughter->Xv() * mcNegDaughter->Xv() + mcNegDaughter->Yv() * mcNegDaughter->Yv());
        decay_length = TMath::Sqrt(TMath::Power(mcNegDaughter->Xv() - mcV0->Xv(), 2) + TMath::Power(mcNegDaughter->Yv() - mcV0->Yv(), 2) +
                                   TMath::Power(mcNegDaughter->Zv() - mcV0->Zv(), 2));
        opening_angle = neg.Angle(pos);
        arm_qt = ArmenterosQt(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz());
        arm_alpha = ArmenterosAlpha(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz(),
                                    mcPosDaughter->Px(), mcPosDaughter->Py(), mcPosDaughter->Pz());
        cpa_wrt_pv = CosinePointingAngle(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                         fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill bookkeeping histograms */

        if (pdg_code == -3122 || pdg_code == 310) {

            fHist_V0s_Bookkeep[pdg_code]->Fill(30);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Mass")]->Fill(mcV0->M());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Px")]->Fill(mcV0->Px());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Py")]->Fill(mcV0->Py());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Phi")]->Fill(mcV0->Phi());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "OriginRadius")]->Fill(origin_radius);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Rapidity")]->Fill(mcV0->Y());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[pdg_code]->Fill(40);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Mass")]->Fill(mcV0->M());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Px")]->Fill(mcV0->Px());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Py")]->Fill(mcV0->Py());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Phi")]->Fill(mcV0->Phi());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "OriginRadius")]->Fill(origin_radius);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Rapidity")]->Fill(mcV0->Y());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdg_code]->Fill(50);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Mass")]->Fill(mcV0->M());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Px")]->Fill(mcV0->Px());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Py")]->Fill(mcV0->Py());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Phi")]->Fill(mcV0->Phi());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "OriginRadius")]->Fill(origin_radius);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Rapidity")]->Fill(mcV0->Y());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);
        }
    }  // end of loop over true V0s
}

/*
 * Find the reconstructed signal pi+ (that don't come from the anti-lambda) and pi-
 * (only valid for Reaction Channel E: anti-lambda, K+, pi+, pi-)
 */
void AliAnalysisTaskSexaquark::ProcessFindablePionPairs() {

    Int_t single_pion_pdg_code;
    Int_t single_pion_mother_pdg_code;

    Int_t mcIdxPiPlus = -1;
    Int_t mcIdxPiMinus = -1;
    AliMCParticle* mcPiPlus;
    AliMCParticle* mcPiMinus;

    /* Declare TLorentzVectors and TVector3 */

    TLorentzVector lvPiPlus, lvPiMinus;
    TLorentzVector lvPionPair;
    TVector3 secondary_vertex(0., 0., 0.);

    /* Declare properties */

    Float_t origin_radius;
    Float_t decay_length;
    Float_t opening_angle;
    Float_t arm_qt;
    Float_t arm_alpha;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over interactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection against general-purpose simulations */
        // There should be tracks that belong to the interaction //

        if (!getEsdIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Collect true info from rec. signal tracks */

        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {

            single_pion_pdg_code = getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[esdIdx]];
            single_pion_mother_pdg_code = 0;
            if (doesEsdIdxHaveMother[esdIdx]) single_pion_mother_pdg_code = getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdx]];

            if (TMath::Abs(single_pion_pdg_code) != 211) continue;
            if (single_pion_mother_pdg_code == -3122) continue;

            // if it hasn't been set yet, set it
            if (mcIdxPiPlus == -1) mcIdxPiPlus = getMcIdx_fromEsdIdx[esdIdx];
            if (mcIdxPiMinus == -1) mcIdxPiMinus = getMcIdx_fromEsdIdx[esdIdx];
        }

        /* Findable condition */
        // At least 1 pi+ and 1 pi- should have been reconstructed and passed track selection

        if (mcIdxPiPlus < 0 || mcIdxPiMinus < 0) continue;

        /* Get MC particles */

        mcPiPlus = (AliMCParticle*)fMC->GetTrack(mcIdxPiPlus);
        mcPiMinus = (AliMCParticle*)fMC->GetTrack(mcIdxPiMinus);

        /* Get true momentum */

        lvPiPlus.SetXYZM(mcPiPlus->Px(), mcPiPlus->Py(), mcPiPlus->Pz(), fPDG.GetParticle(211)->Mass());
        lvPiMinus.SetXYZM(mcPiMinus->Px(), mcPiMinus->Py(), mcPiMinus->Pz(), fPDG.GetParticle(211)->Mass());

        lvPionPair = lvPiPlus + lvPiMinus;

        /* Get sec. vertex */

        secondary_vertex.SetXYZ(mcPiPlus->Xv(), mcPiPlus->Yv(), mcPiPlus->Zv());

        /* Get properties */

        origin_radius = TMath::Sqrt(TMath::Power(secondary_vertex.X(), 2) + TMath::Power(secondary_vertex.Y(), 2));
        decay_length = TMath::Sqrt(TMath::Power(secondary_vertex.X() - fMC_PrimaryVertex->GetX(), 2) +
                                   TMath::Power(secondary_vertex.Y() - fMC_PrimaryVertex->GetY(), 2) +
                                   TMath::Power(secondary_vertex.Z() - fMC_PrimaryVertex->GetZ(), 2));
        opening_angle = lvPionPair.Vect().Angle(lvPiMinus.Vect());
        arm_qt = ArmenterosQt(lvPionPair.Px(), lvPionPair.Py(), lvPionPair.Pz(), lvPiMinus.Px(), lvPiMinus.Py(), lvPiMinus.Pz());
        arm_alpha = ArmenterosAlpha(lvPionPair.Px(), lvPionPair.Py(), lvPionPair.Pz(), lvPiMinus.Px(), lvPiMinus.Py(), lvPiMinus.Pz(), lvPiPlus.Px(),
                                    lvPiPlus.Py(), lvPiPlus.Pz());
        cpa_wrt_pv = CosinePointingAngle(lvPiPlus, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvPiPlus.Px(), lvPiPlus.Py(), lvPiPlus.Pz(), secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill histograms  */

        fHist_V0s_Bookkeep[422]->Fill(50);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Mass")]->Fill(lvPionPair.M());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Pt")]->Fill(lvPionPair.Pt());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Px")]->Fill(lvPionPair.Px());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Py")]->Fill(lvPionPair.Py());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Pz")]->Fill(lvPionPair.Pz());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Phi")]->Fill(lvPionPair.Phi());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "OriginRadius")]->Fill(origin_radius);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "DecayRadius")]->Fill(origin_radius);  // NOTE: is the same on purpose
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Zv")]->Fill(secondary_vertex.Z());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Eta")]->Fill(lvPionPair.Eta());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Rapidity")]->Fill(lvPionPair.Rapidity());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "DecayLength")]->Fill(decay_length);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "OpeningAngle")]->Fill(opening_angle);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "ArmQt")]->Fill(arm_qt);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "ArmAlpha")]->Fill(arm_alpha);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "DCAwrtPV")]->Fill(dca_wrt_pv);
    }
}

/*
 * Find the anti-sexaquarks that have all of their *final-state particles* reconstructed. Valid for:
 * - Reaction Channel A: 1 anti-proton, 2 pi+, 1 pi-
 * - Reaction Channel D: 1 anti-proton, 1 pi+, 1 K+
 * - Reaction Channel E: 1 anti-proton, 2 pi+, 1 pi-, 1 K+
 */
void AliAnalysisTaskSexaquark::ProcessFindableSexaquarks() {

    AliMCParticle* mcFSPart;  // FS: final-state
    Int_t mcIdx;
    AliMCParticle* mcFGPart;  // FG: first gen.
    Int_t motherMcIdx;

    /* Declare TLorentzVectors */

    TLorentzVector lvFinalStateParticle;
    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Declare variables */

    TVector3 secondary_vertex(0., 0., 0.);
    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over interactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection in case of general-purpose simulations */

        if (!getEsdIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Findable condition */
        // All final state products should have been reconstructed and passed track selection

        std::unordered_multiset<Int_t> pdg_reconstructed_particles;
        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {
            pdg_reconstructed_particles.insert(getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[esdIdx]]);
        }

        std::unordered_multiset<Int_t> expected_pdg_fs_particles(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end());

        if (expected_pdg_fs_particles != pdg_reconstructed_particles) continue;

        /* Reset vectors */

        secondary_vertex.SetXYZ(0., 0., 0.);
        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        /* Loop, to get the reconstructed final-state particles info */

        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {

            mcIdx = getMcIdx_fromEsdIdx[esdIdx];
            mcFSPart = (AliMCParticle*)fMC->GetTrack(mcIdx);

            lvFinalStateParticle.SetPxPyPzE(mcFSPart->Px(), mcFSPart->Py(), mcFSPart->Pz(), mcFSPart->E());

            lvAntiSexaquark = lvAntiSexaquark + lvFinalStateParticle;

            /* Get sec. vertex from the origin of the first-generation products of the reaction */

            if (secondary_vertex.Mag() < 1E-6) {
                if (mcFSPart->MCStatusCode() == reactionIdx) {
                    secondary_vertex.SetXYZ(mcFSPart->Xv(), mcFSPart->Yv(), mcFSPart->Zv());
                } else if (doesMcIdxHaveMother[mcIdx]) {
                    motherMcIdx = getMotherMcIdx_fromMcIdx[mcIdx];
                    mcFGPart = (AliMCParticle*)fMC->GetTrack(motherMcIdx);
                    if (mcFGPart->MCStatusCode() == reactionIdx) secondary_vertex.SetXYZ(mcFGPart->Xv(), mcFGPart->Yv(), mcFGPart->Zv());
                }
            }
        }

        lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest
        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /* Get properties */

        radius = TMath::Sqrt(TMath::Power(secondary_vertex.X(), 2) + TMath::Power(secondary_vertex.Y(), 2));
        decay_length = TMath::Sqrt(TMath::Power(secondary_vertex.X() - fMC_PrimaryVertex->GetX(), 2) +
                                   TMath::Power(secondary_vertex.Y() - fMC_PrimaryVertex->GetY(), 2) +
                                   TMath::Power(secondary_vertex.Z() - fMC_PrimaryVertex->GetZ(), 2));
        cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), secondary_vertex.X(), secondary_vertex.Y(),
                                  secondary_vertex.Z(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill histograms  */

        fHist_AntiSexaquarks_Bookkeep->Fill(10);
        if (fReactionID == 'H') fHist_PosKaonPairs_Bookkeep->Fill(1);  // PENDING: should be ~ 0, because it's a channel so difficult to rec.
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Mass")]->Fill(lvAntiSexaquark.M());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Px")]->Fill(lvAntiSexaquark.Px());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Py")]->Fill(lvAntiSexaquark.Py());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Radius")]->Fill(radius);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Zv")]->Fill(secondary_vertex.Z());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "DecayLength")]->Fill(decay_length);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
    }
}

/*
 * Find the anti-sexaquarks that have two K+s reconstructed.
 * Valid for Reaction Channel H: 1 anti-proton, 2 K+, 1 pi0
 */
void AliAnalysisTaskSexaquark::ProcessFindablePosKaonPairs() {

    AliMCParticle* mcPosKaonA;
    AliMCParticle* mcPosKaonB;

    /* Declare TLorentzVectors */

    TLorentzVector lvPosKaonA, lvPosKaonB;
    TLorentzVector lvPosKaonPair;

    /* Declare vertices */

    TVector3 primary_vertex(fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
    TVector3 secondary_vertex(0., 0., 0.);

    /* Auxiliary variables */

    Int_t pdgCode;
    Int_t mcIdxPosKaonA = -1;
    Int_t mcIdxPosKaonB = -1;
    Int_t motherMcIdx;

    /* Loop over interactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection in case of general-purpose simulations */

        if (!getEsdIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Collect true info from rec. signal tracks */

        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {

            pdgCode = getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[esdIdx]];

            if (pdgCode != 321) continue;

            if (mcIdxPosKaonA == -1)
                mcIdxPosKaonA = getMcIdx_fromEsdIdx[esdIdx];
            else if (mcIdxPosKaonB == -1)
                mcIdxPosKaonB = getMcIdx_fromEsdIdx[esdIdx];
        }

        /* Findable condition */
        // Two positive kaons should have been reconstructed and passed track selection

        if (mcIdxPosKaonA < 0 || mcIdxPosKaonB < 0) continue;

        /* Get MC particles */

        mcPosKaonA = (AliMCParticle*)fMC->GetTrack(mcIdxPosKaonA);
        mcPosKaonB = (AliMCParticle*)fMC->GetTrack(mcIdxPosKaonB);

        /* Get true momentum */

        lvPosKaonA.SetXYZM(mcPosKaonA->Px(), mcPosKaonA->Py(), mcPosKaonA->Pz(), fPDG.GetParticle(321)->Mass());
        lvPosKaonB.SetXYZM(mcPosKaonB->Px(), mcPosKaonB->Py(), mcPosKaonB->Pz(), fPDG.GetParticle(321)->Mass());

        lvPosKaonPair = lvPosKaonA + lvPosKaonB;

        /* Get sec. vertex */

        secondary_vertex.SetXYZ(mcPosKaonA->Xv(), mcPosKaonA->Yv(), mcPosKaonA->Zv());
    }

    /* Fill histograms  */

    fHist_PosKaonPairs_Bookkeep->Fill(10);
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Mass")]->Fill(lvPosKaonPair.M());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "OriginRadius")]->Fill(secondary_vertex.Perp());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "DecayLength")]->Fill((secondary_vertex - primary_vertex).Mag());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Xv")]->Fill(secondary_vertex.X());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Yv")]->Fill(secondary_vertex.Y());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Zv")]->Fill(secondary_vertex.Z());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Eta")]->Fill(lvPosKaonPair.Eta());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Rapidity")]->Fill(lvPosKaonPair.Rapidity());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Px")]->Fill(lvPosKaonPair.Px());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Py")]->Fill(lvPosKaonPair.Py());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Pz")]->Fill(lvPosKaonPair.Pz());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Pt")]->Fill(lvPosKaonPair.Pt());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "OpeningAngle")]->Fill(lvPosKaonA.Vect().Angle(lvPosKaonB.Vect()));
}

/*                             */
/**  V0s -- ALICE V0 Finders  **/
/*** ======================= ***/

/*
 * Loop over all the V0s that were found and stored in the ESD file by the Official V0 Finders.
 * - Uses: `fESD`, `fMC`, `fPDG`, `fPrimaryVertex`
 * - Containers: `esdAntiLambdas`, `getEsdIdxOfNegDau_fromAntiLambdaIdx`, `getEsdIdxOfPosDau_fromAntiLambdaIdx`, `isAntiLambdaIdxATrueV0`,
 * `isEsdIdxDaughterOfTrueV0`, `getMcIdxOfTrueV0_fromEsdIdx`, `isMcIdxSignal`, `esdKaonsZeroShort`, `getEsdIdxOfNegDau_fromKaonZeroShortIdx`,
 * `getEsdIdxOfPosDau_fromKaonZeroShortIdx`, `isKaonZeroShortIdxATrueV0`
 */
void AliAnalysisTaskSexaquark::GetV0sFromESD(Bool_t onTheFly) {

    AliESDv0* V0;

    Int_t esdIdxPos;
    Int_t esdIdxNeg;

    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    Float_t decay_radius;
    Float_t decay_length;  // decay vertex - primary vertex, could be useful
    Float_t opening_angle;
    Float_t cpa_wrt_pv, dca_wrt_pv;
    Float_t dca_btw_dau, dca_neg_v0, dca_pos_v0;

    AliESDv0 outputV0;

    Double_t xthiss, xpp;
    Float_t impar_neg_v0[2], impar_pos_v0[2];

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Auxiliary containers */

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Bool_t passes_cuts_as_antilambda;
    Bool_t passes_cuts_as_kaonzero;

    std::vector<AliESDv0>* esdFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    /* Loop over V0s that are stored in the ESD */

    for (Int_t esdIdxV0 = 0; esdIdxV0 < fESD->GetNumberOfV0s(); esdIdxV0++) {

        V0 = fESD->GetV0(esdIdxV0);

        /* Choose between offline and on-the-fly V0s */

        if (!onTheFly && V0->GetOnFlyStatus()) {
            continue;
        }

        if (onTheFly && !V0->GetOnFlyStatus()) {
            continue;
        }

        esdIdxNeg = V0->GetNindex();
        esdIdxPos = V0->GetPindex();

        /* Check if both daughters passed track selection */

        if (!esdIdxPassesTrackSelection_AsAntiProton[esdIdxNeg] && !esdIdxPassesTrackSelection_AsNegKaon[esdIdxNeg] &&
            !esdIdxPassesTrackSelection_AsPiMinus[esdIdxNeg]) {
            continue;
        }

        if (!esdIdxPassesTrackSelection_AsProton[esdIdxPos] && !esdIdxPassesTrackSelection_AsPosKaon[esdIdxPos] &&
            !esdIdxPassesTrackSelection_AsPiPlus[esdIdxPos]) {
            continue;
        }

        /* Get tracks */

        neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
        pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

        /* Apply cuts */

        is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                  getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                  getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == -3122;
        is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
        is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;

        passes_cuts_as_antilambda = PassesV0CutsAs(-3122, is_true, is_secondary, is_signal, V0, neg_track, pos_track);

        is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                  getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                  getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == 310;
        is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
        is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;

        passes_cuts_as_kaonzero = PassesV0CutsAs(310, is_true, is_secondary, is_signal, V0, neg_track, pos_track);

        if (!passes_cuts_as_antilambda && !passes_cuts_as_kaonzero) {
            continue;
        }

        /* Get properties */

        V0->GetXYZ(V0_X, V0_Y, V0_Z);
        V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
        V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

        neg_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_neg_v0);
        pos_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_pos_v0);

        decay_radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);
        decay_length = TMath::Sqrt((V0_X - fPrimaryVertex->GetX()) * (V0_X - fPrimaryVertex->GetX()) +
                                   (V0_Y - fPrimaryVertex->GetY()) * (V0_Y - fPrimaryVertex->GetY()) +
                                   (V0_Z - fPrimaryVertex->GetZ()) * (V0_Z - fPrimaryVertex->GetZ()));

        cpa_wrt_pv = V0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_wrt_pv =
            LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_btw_dau = TMath::Abs(neg_track->GetDCA(pos_track, fMagneticField, xthiss, xpp));              // DCAxyz
        dca_neg_v0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);  // DCAxyz
        dca_pos_v0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);  // DCAxyz

        /* Fill histograms */

        std::vector<Int_t> possiblePID;
        if (passes_cuts_as_antilambda) possiblePID.push_back(-3122);
        if (passes_cuts_as_kaonzero) possiblePID.push_back(310);

        for (Int_t& v0PdgCode : possiblePID) {

            /* Update TLorentzVectors to get reconstructed mass */

            lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[v0PdgCode])->Mass());
            lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[v0PdgCode])->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Determine containers to use */

            if (v0PdgCode == -3122) {
                esdFoundV0s = &esdAntiLambdas;
                getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
                getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
            } else if (v0PdgCode == 310) {
                esdFoundV0s = &esdKaonsZeroShort;
                getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
                getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
            }

            /* Store V0 */

            outputV0 = *V0;

            esdFoundV0s->push_back(outputV0);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxPos;

            /* Fill histograms */

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(80);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(V0->AlphaV0(), V0->PtArmV0());

            if (!is_true) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(110);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(V0->AlphaV0(), V0->PtArmV0());

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(140);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(V0->AlphaV0(), V0->PtArmV0());

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(170);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(V0->AlphaV0(), V0->PtArmV0());
        }
    }  // end of loop over V0s
}

/*
 * Custom V0 Finder.
 * (Adapted from `AliRoot/STEER/ESD/AliV0vertexer.cxx`)
 * - Uses: `fESD`, `fMC`, `fPDG`, `fPrimaryVertex`, `fMagneticField`
 * - Input: `pdgV0`
 */
void AliAnalysisTaskSexaquark::CustomV0Finder(Int_t pdgV0) {

    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Auxiliary variables */

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Double_t neg_d_wrt_pv, pos_d_wrt_pv;
    Double_t decay_radius;
    Double_t decay_length;  // decay vertex - primary vertex, could be useful
    Double_t opening_angle;
    Double_t cpa_wrt_pv, dca_wrt_pv;
    Double_t dca_btw_dau, dca_neg_v0, dca_pos_v0;

    Double_t xthiss, xpp;
    Float_t impar_neg_v0[2], impar_pos_v0[2];

    /* Custom V0 Finder's Technical Cuts */

    Double_t COV0F_DCAmax = 0.6;  // (d=1.) max DCA between the daughter tracks
    Double_t COV0F_Rmin = 50.;    // (d=0.9) min radius of the V0 radius
    Double_t COV0F_Rmax = 200.;   // (d=100) max radius of the V0 radius
    Bool_t COV0F_UseImprovedFinding = kTRUE;

    /*  Choose between anti-lambda and K0S */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;

    std::vector<AliESDv0>* esdFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        esdFoundV0s = &esdAntiLambdas;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
    } else if (pdgV0 == 310) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        esdFoundV0s = &esdKaonsZeroShort;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
    }

    /* Loop over all possible pairs of selected tracks */

    for (Int_t& esdIdxNeg : esdIndicesNegTracks) {
        for (Int_t& esdIdxPos : esdIndicesPosTracks) {

            /* Sanity check */

            if (esdIdxNeg == esdIdxPos) continue;

            /* Get tracks */

            neg_track = fESD->GetTrack(esdIdxNeg);
            pos_track = fESD->GetTrack(esdIdxPos);

            /* --------------------------------------------- */
            /* --- Beginning of Core of Custom V0 Finder --- */

            neg_d_wrt_pv = neg_track->GetD(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fMagneticField);
            if (TMath::Abs(neg_d_wrt_pv) > COV0F_Rmax) continue;
            pos_d_wrt_pv = pos_track->GetD(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fMagneticField);
            if (TMath::Abs(pos_d_wrt_pv) > COV0F_Rmax) continue;

            AliExternalTrackParam nt(*neg_track), pt(*pos_track), *ntp = &nt, *ptp = &pt;

            Double_t xn, xp;
            if (COV0F_UseImprovedFinding && Preoptimize(ntp, ptp, &xn, &xp, fMagneticField)) {
                nt.PropagateTo(xn, fMagneticField);
                pt.PropagateTo(xp, fMagneticField);
            }

            Double_t dca_after_preopt = ntp->GetDCA(ptp, fMagneticField, xn, xp);
            if (dca_after_preopt > COV0F_DCAmax) continue;
            if ((xn + xp) > 2 * COV0F_Rmax) continue;
            if ((xn + xp) < 2 * COV0F_Rmin) continue;

            nt.PropagateTo(xn, fMagneticField);
            pt.PropagateTo(xp, fMagneticField);

            /* Note: declaration of AliESDv0 within the loop, because it has important functions on its constructor */
            AliESDv0 vertex(nt, esdIdxNeg, pt, esdIdxPos);
            if (COV0F_UseImprovedFinding) vertex.Refit();

            vertex.SetDcaV0Daughters(dca_after_preopt);

            /* --- End of Core of Custom V0 Finder --- */
            /* --------------------------------------- */

            /* Apply cuts */

            is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                      getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                      getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == pdgV0;
            is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
            is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;

            if (!PassesV0CutsAs(pdgV0, is_true, is_secondary, is_signal, &vertex, neg_track, pos_track)) continue;

            /* Get V0 properties */

            vertex.GetXYZ(V0_X, V0_Y, V0_Z);
            vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
            vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

            lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[pdgV0])->Mass());
            lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[pdgV0])->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            decay_radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);

            cpa_wrt_pv = vertex.GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());

            dca_wrt_pv = LinePointDCA(vertex.Px(), vertex.Py(), vertex.Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                      fPrimaryVertex->GetZ());

            dca_btw_dau = TMath::Abs(neg_track->GetDCA(pos_track, fMagneticField, xthiss, xpp));  // DCAxyz

            neg_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_neg_v0);
            dca_neg_v0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);  // DCAxyz

            pos_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_pos_v0);
            dca_pos_v0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);  // DCAxyz

            decay_length = TMath::Sqrt((V0_X - fPrimaryVertex->GetX()) * (V0_X - fPrimaryVertex->GetX()) +
                                       (V0_Y - fPrimaryVertex->GetY()) * (V0_Y - fPrimaryVertex->GetY()) +
                                       (V0_Z - fPrimaryVertex->GetZ()) * (V0_Z - fPrimaryVertex->GetZ()));
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Store V0 */

            esdFoundV0s->push_back(vertex);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxPos;

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(80);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());

            if (!is_true) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(110);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(140);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(170);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*
 * Apply cuts to a V0 candidate.
 * - Input: `pdgV0`, `v0`, `neg_track`, `pos_track`
 * - Uses: `fPrimaryVertex`
 * - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesV0CutsAs(Int_t pdgV0, Bool_t IsTrue, Bool_t IsSecondary, Bool_t IsSignal, AliESDv0* v0, AliESDtrack* neg_track,
                                                AliESDtrack* pos_track) {

    fHist_V0s_Bookkeep[pdgV0]->Fill(60);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(90);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(120);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(150);

    Double_t V0_X, V0_Y, V0_Z;
    v0->GetXYZ(V0_X, V0_Y, V0_Z);

    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;
    v0->GetNPxPyPz(N_Px, N_Py, N_Pz);
    v0->GetPPxPyPz(P_Px, P_Py, P_Pz);
    TLorentzVector lvTrackNeg, lvTrackPos;
    lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[pdgV0])->Mass());
    lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[pdgV0])->Mass());
    TLorentzVector lvV0 = lvTrackNeg + lvTrackPos;

    Double_t Radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);
    if (kMin_V0_Radius[pdgV0] && Radius < kMin_V0_Radius[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(61);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(91);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(121);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(151);

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(62);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(92);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(122);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(152);

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(63);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(93);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(123);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(153);

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(64);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(94);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(124);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(154);

    Double_t DecayLength = TMath::Sqrt(TMath::Power(v0->Xv() - fPrimaryVertex->GetX(), 2) + TMath::Power(v0->Yv() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(v0->Zv() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(65);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(95);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(125);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(155);

    Double_t CPAwrtPV = v0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(66);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(96);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(126);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(156);

    Double_t DCAwrtPV = LinePointDCA(v0->Px(), v0->Py(), v0->Pz(), v0->Xv(), v0->Yv(), v0->Zv(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                     fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(67);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(97);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(127);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(157);

    Double_t DCAbtwDau = TMath::Abs(v0->GetDcaV0Daughters());
    if (kMax_V0_DCAbtwDau[pdgV0] && DCAbtwDau > kMax_V0_DCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(68);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(98);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(128);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(158);

    Float_t impar_neg_v0[2];
    neg_track->GetDZ(v0->Xv(), v0->Yv(), v0->Zv(), fMagneticField, impar_neg_v0);
    Double_t DCAnegV0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(69);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(99);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(129);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(159);

    Float_t impar_pos_v0[2];
    pos_track->GetDZ(v0->Xv(), v0->Yv(), v0->Zv(), fMagneticField, impar_pos_v0);
    Double_t DCAposV0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(70);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(100);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(130);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(160);

    Double_t ArmPt = ArmenterosQt(v0->Px(), v0->Py(), v0->Pz(), N_Px, N_Py, N_Pz);
    Double_t ArmAlpha = ArmenterosAlpha(v0->Px(), v0->Py(), v0->Pz(), N_Px, N_Py, N_Pz, P_Px, P_Py, P_Pz);
    Double_t ArmPtOverAlpha = ArmPt / TMath::Abs(ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(71);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(101);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(131);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(161);

    /* Exclusive for Offline, On-The-Fly, and Custom V0 Finders */

    Double_t Chi2 = v0->GetChi2V0();
    if (kMax_V0_Chi2[pdgV0] && Chi2 > kMax_V0_Chi2[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(72);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(102);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(132);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(162);

    Double_t ImprvDCAbtwDau = v0->GetDcaV0Daughters();
    if (kMax_V0_ImprvDCAbtwDau[pdgV0] && ImprvDCAbtwDau > kMax_V0_ImprvDCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(73);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(103);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(133);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(163);

    return kTRUE;
}

/*
 * Apply cuts to a V0 candidate.
 * - Input: `pdgV0`, `kfV0`, `kfDaughterNeg`, `kfDaughterPos`, `lvV0`, `lvTrackNeg`, `lvTrackPos`
 * - Uses: `fPrimaryVertex`
 * - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesV0CutsAs(Int_t pdgV0, Bool_t IsTrue, Bool_t IsSecondary, Bool_t IsSignal, KFParticleMother kfV0,
                                                KFParticle kfDaughterNeg, KFParticle kfDaughterPos, TLorentzVector lvV0, TLorentzVector lvTrackNeg,
                                                TLorentzVector lvTrackPos) {

    fHist_V0s_Bookkeep[pdgV0]->Fill(60);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(90);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(120);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(150);

    Double_t Radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
    if (kMin_V0_Radius[pdgV0] && Radius < kMin_V0_Radius[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(61);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(91);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(121);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(151);

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(62);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(92);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(122);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(152);

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(63);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(93);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(123);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(153);

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(64);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(94);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(124);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(154);

    Double_t DecayLength = TMath::Sqrt(TMath::Power(kfV0.GetX() - fPrimaryVertex->GetX(), 2) + TMath::Power(kfV0.GetY() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(kfV0.GetZ() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(65);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(95);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(125);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(155);

    Double_t CPAwrtPV =
        CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(66);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(96);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(126);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(156);

    Double_t DCAwrtPV = LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(), kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(),
                                     fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(67);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(97);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(127);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(157);

    Double_t DCAbtwDau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
    if (kMax_V0_DCAbtwDau[pdgV0] && DCAbtwDau > kMax_V0_DCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(68);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(98);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(128);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(158);

    Double_t DCAnegV0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(69);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(99);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(129);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(159);

    Double_t DCAposV0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(70);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(100);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(130);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(160);

    Double_t ArmPt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
    Double_t ArmAlpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());
    Double_t ArmPtOverAlpha = TMath::Abs(ArmPt / ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(71);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(101);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(131);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(161);

    /* Exclusive for KF method */

    Double_t Chi2ndf = (Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF();
    if (kMax_V0_Chi2ndf[pdgV0] && Chi2ndf > kMax_V0_Chi2ndf[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(72);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(102);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(132);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(162);

    return kTRUE;
}

/*                                     */
/**  Sexaquark -- Geometrical Finder  **/
/*** =============================== ***/

/*
 * Using two AliESDv0 objects, reconstruct the anti-sexaquark candidate.
 * - Containers: `esdAntiLambdas`, `esdKaonsZeroShort`, `getEsdIdxOfNegDau_fromAntiLambdaIdx`, `getEsdIdxOfNegDau_fromKaonZeroShortIdx`,
 * `getEsdIdxOfPosDau_fromKaonZeroShortIdx`, `getEsdIdxOfPosDau_fromAntiLambdaIdx`
 */
void AliAnalysisTaskSexaquark::GeoSexaquarkFinder_ChannelA() {

    /* Declare AliRoot objects */

    AliESDv0 AntiLambda;
    AliESDv0 KaonZeroShort;

    AliESDtrack* AntiLambda_NegDau;
    AliESDtrack* AntiLambda_PosDau;
    AliESDtrack* KaonZeroShort_NegDau;
    AliESDtrack* KaonZeroShort_PosDau;

    /* Declare TVector3s  */

    TVector3 AntiLambda_momentum, KaonZeroShort_momentum;
    TVector3 AntiLambda_vertex, KaonZeroShort_vertex;

    TVector3 SecondaryVertex;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Declare auxiliar variables */

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxNeg_KaonZeroShort, esdIdxPos_KaonZeroShort;

    TVector3 pca1, pca2;
    Float_t radius;
    Float_t decay_length;
    Float_t dca_wrt_pv;
    Float_t cpa_wrt_pv;
    Float_t dca_btw_v0s;
    Float_t dca_v0a_sv;
    Float_t dca_v0b_sv;
    Float_t impar_v0a_neg_sv[2], impar_v0a_pos_sv[2];
    Float_t impar_v0b_neg_sv[2], impar_v0b_pos_sv[2];
    Float_t dca_v0a_neg_sv;
    Float_t dca_v0a_pos_sv;
    Float_t dca_v0b_neg_sv;
    Float_t dca_v0b_pos_sv;

    Bool_t is_signal;

    /* Loop over all possible pairs of V0s */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)esdAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)esdKaonsZeroShort.size(); idxKaonZeroShort++) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxNeg_KaonZeroShort = getEsdIdxOfNegDau_fromKaonZeroShortIdx[idxKaonZeroShort];
            esdIdxPos_KaonZeroShort = getEsdIdxOfPosDau_fromKaonZeroShortIdx[idxKaonZeroShort];

            if (esdIdxNeg_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxNeg_AntiLambda == esdIdxPos_KaonZeroShort ||
                esdIdxPos_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxPos_AntiLambda == esdIdxPos_KaonZeroShort) {
                continue;
            }

            /* Get AliRoot objects */

            AntiLambda = esdAntiLambdas[idxAntiLambda];
            KaonZeroShort = esdKaonsZeroShort[idxKaonZeroShort];

            AntiLambda_NegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            AntiLambda_PosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));
            KaonZeroShort_NegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_KaonZeroShort));
            KaonZeroShort_PosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_KaonZeroShort));

            /* Transport V0s as if they were straight lines */

            AntiLambda_momentum.SetXYZ(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz());
            KaonZeroShort_momentum.SetXYZ(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz());

            AntiLambda_vertex.SetXYZ(AntiLambda.Xv(), AntiLambda.Yv(), AntiLambda.Zv());
            KaonZeroShort_vertex.SetXYZ(KaonZeroShort.Xv(), KaonZeroShort.Yv(), KaonZeroShort.Zv());

            dca_btw_v0s = Calculate_TwoLinesDCA_v2(AntiLambda_vertex, AntiLambda_momentum, KaonZeroShort_vertex, KaonZeroShort_momentum, pca1, pca2);

            SecondaryVertex = 0.5 * (pca1 + pca2);  // middle point between the two points of closest approach

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda.SetPxPyPzE(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz(), AntiLambda.E());
            lvKaonZeroShort.SetPxPyPzE(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz(), KaonZeroShort.E());

            lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assuming struck nucleon at rest

            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvStruckNucleon;

            /* Apply cuts */

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            if (!PassesSexaquarkCuts_ChannelA(is_signal, SecondaryVertex, lvAntiSexaquark, AntiLambda, KaonZeroShort, AntiLambda_NegDau,
                                              AntiLambda_PosDau, KaonZeroShort_NegDau, KaonZeroShort_PosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(SecondaryVertex.X() * SecondaryVertex.X() + SecondaryVertex.Y() * SecondaryVertex.Y());
            decay_length = TMath::Sqrt(TMath::Power(SecondaryVertex.X() - fPrimaryVertex->GetX(), 2) +
                                       TMath::Power(SecondaryVertex.Y() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(SecondaryVertex.Z() - fPrimaryVertex->GetZ(), 2));
            cpa_wrt_pv =
                CosinePointingAngle(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                    SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                      SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_v0a_sv = LinePointDCA(lvAntiLambda.Px(), lvAntiLambda.Py(), lvAntiLambda.Pz(), AntiLambda_vertex.X(), AntiLambda_vertex.Y(),
                                      AntiLambda_vertex.Z(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
            dca_v0b_sv =
                LinePointDCA(lvKaonZeroShort.Px(), lvKaonZeroShort.Py(), lvKaonZeroShort.Pz(), KaonZeroShort_vertex.X(), KaonZeroShort_vertex.Y(),
                             KaonZeroShort_vertex.Z(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());

            AntiLambda_NegDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_neg_sv);
            dca_v0a_neg_sv = TMath::Sqrt(impar_v0a_neg_sv[0] * impar_v0a_neg_sv[0] + impar_v0a_neg_sv[1] * impar_v0a_neg_sv[1]);

            AntiLambda_PosDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_pos_sv);
            dca_v0a_pos_sv = TMath::Sqrt(impar_v0a_pos_sv[0] * impar_v0a_pos_sv[0] + impar_v0a_pos_sv[1] * impar_v0a_pos_sv[1]);

            KaonZeroShort_NegDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_neg_sv);
            dca_v0b_neg_sv = TMath::Sqrt(impar_v0b_neg_sv[0] * impar_v0b_neg_sv[0] + impar_v0b_neg_sv[1] * impar_v0b_neg_sv[1]);

            KaonZeroShort_PosDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_pos_sv);
            dca_v0b_pos_sv = TMath::Sqrt(impar_v0b_pos_sv[0] * impar_v0b_pos_sv[0] + impar_v0b_pos_sv[1] * impar_v0b_pos_sv[1]);

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(50);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(SecondaryVertex.Z());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0bposSV")]->Fill(dca_v0b_pos_sv);

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(90);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(SecondaryVertex.Z());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0bposSV")]->Fill(dca_v0b_pos_sv);
        }  // end of loop over K0S
    }      // end of loop over anti-lambdas

    return;
}

/*
 * Apply anti-sexaquark candidate cuts. Exclusive to reaction channel A: `AntiSexaquark,N->AntiLambda,K0S`.
 * - Use: `fPrimaryVertex`, `fMagneticField`
 * - Input: `SecondaryVertex`, `lvAntiSexaquark`, `esdAntiLambda`, `esdKaonZeroShort`, `esdAntiLambdaNeg`, `esdAntiLambdaPos`, `esdKaonZeroShortNeg`,
 * `esdKaonZeroShortPos`
 * - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA(Bool_t isSignal, TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark,
                                                              AliESDv0 esdAntiLambda, AliESDv0 esdKaonZeroShort, AliESDtrack* esdAntiLambdaNeg,
                                                              AliESDtrack* esdAntiLambdaPos, AliESDtrack* esdKaonZeroShortNeg,
                                                              AliESDtrack* esdKaonZeroShortPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(SecondaryVertex.X() * SecondaryVertex.X() + SecondaryVertex.Y() * SecondaryVertex.Y());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((SecondaryVertex.X() - fPrimaryVertex->GetX()) * (SecondaryVertex.X() - fPrimaryVertex->GetX()) +
                                        (SecondaryVertex.Y() - fPrimaryVertex->GetY()) * (SecondaryVertex.Y() - fPrimaryVertex->GetY()) +
                                        (SecondaryVertex.Z() - fPrimaryVertex->GetZ()) * (SecondaryVertex.Z() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fPrimaryVertex->GetX(),
                                              fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    Double_t dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                       SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    TVector3 AntiLambda_vertex(esdAntiLambda.Xv(), esdAntiLambda.Yv(), esdAntiLambda.Zv());
    TVector3 AntiLambda_momentum(esdAntiLambda.Px(), esdAntiLambda.Py(), esdAntiLambda.Pz());
    TVector3 KaonZeroShort_vertex(esdKaonZeroShort.Xv(), esdKaonZeroShort.Yv(), esdKaonZeroShort.Zv());
    TVector3 KaonZeroShort_momentum(esdKaonZeroShort.Px(), esdKaonZeroShort.Py(), esdKaonZeroShort.Pz());
    TVector3 pca1, pca2;
    Double_t dca_btw_v0s = Calculate_TwoLinesDCA_v2(AntiLambda_vertex, AntiLambda_momentum, KaonZeroShort_vertex, KaonZeroShort_momentum, pca1, pca2);
    if (kMax_Sexa_DCAbtwV0s && dca_btw_v0s > kMax_Sexa_DCAbtwV0s) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0a_sv = LinePointDCA(esdAntiLambda.Px(), esdAntiLambda.Py(), esdAntiLambda.Pz(), esdAntiLambda.Xv(), esdAntiLambda.Yv(),
                                       esdAntiLambda.Zv(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
    if (kMax_Sexa_DCAv0aSV && dca_v0a_sv > kMax_Sexa_DCAv0aSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_v0b_sv = LinePointDCA(esdKaonZeroShort.Px(), esdKaonZeroShort.Py(), esdKaonZeroShort.Pz(), esdKaonZeroShort.Xv(),
                                       esdKaonZeroShort.Yv(), esdKaonZeroShort.Zv(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
    if (kMax_Sexa_DCAv0bSV && dca_v0b_sv > kMax_Sexa_DCAv0bSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Float_t impar_v0a_neg_sv[2];
    esdAntiLambdaNeg->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_neg_sv);
    Double_t dca_v0a_neg_sv = TMath::Sqrt(impar_v0a_neg_sv[0] * impar_v0a_neg_sv[0] + impar_v0a_neg_sv[1] * impar_v0a_neg_sv[1]);
    if (kMax_Sexa_DCAv0anegSV && dca_v0a_neg_sv > kMax_Sexa_DCAv0anegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Float_t impar_v0a_pos_sv[2];
    esdAntiLambdaPos->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_pos_sv);
    Double_t dca_v0a_pos_sv = TMath::Sqrt(impar_v0a_pos_sv[0] * impar_v0a_pos_sv[0] + impar_v0a_pos_sv[1] * impar_v0a_pos_sv[1]);
    if (kMax_Sexa_DCAv0aposSV && dca_v0a_pos_sv > kMax_Sexa_DCAv0aposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Float_t impar_v0b_neg_sv[2];
    esdKaonZeroShortNeg->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_neg_sv);
    Double_t dca_v0b_neg_sv = TMath::Sqrt(impar_v0b_neg_sv[0] * impar_v0b_neg_sv[0] + impar_v0b_neg_sv[1] * impar_v0b_neg_sv[1]);
    if (kMax_Sexa_DCAv0bnegSV && dca_v0b_neg_sv > kMax_Sexa_DCAv0bnegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Float_t impar_v0b_pos_sv[2];
    esdKaonZeroShortPos->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_pos_sv);
    Double_t dca_v0b_pos_sv = TMath::Sqrt(impar_v0b_pos_sv[0] * impar_v0b_pos_sv[0] + impar_v0b_pos_sv[1] * impar_v0b_pos_sv[1]);
    if (kMax_Sexa_DCAv0bposSV && dca_v0b_pos_sv > kMax_Sexa_DCAv0bposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    return kTRUE;
}

/*
 */
void AliAnalysisTaskSexaquark::GeoSexaquarkFinder_ChannelD() {
    // PENDING
}

/*
 */
void AliAnalysisTaskSexaquark::GeoSexaquarkFinder_ChannelE() {
    // PENDING
}

/*                          */
/**  V0s -- Kalman Filter  **/
/*** ==================== ***/

/*
 * Find all V0s via Kalman Filter.
 * - Uses: `fESD`, `fMC`
 * - Input: `pdgV0`, `pdgTrackNeg`, `pdgTrackPos`
 * - Output: `kfV0s`, `idxDaughters`
 */
void AliAnalysisTaskSexaquark::KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos) {

    AliESDtrack* esdTrackNeg;
    AliESDtrack* esdTrackPos;

    Double_t impar_neg[2], impar_pos[2];

    Double_t decay_radius;
    Double_t decay_length;  // decay radius - primary vertex, could be useful
    Double_t opening_angle;
    Double_t arm_qt, arm_alpha;
    Double_t cpa_wrt_pv, dca_wrt_pv;
    Double_t dca_btw_dau, dca_neg_v0, dca_pos_v0;
    Double_t chi2_ndf;

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Int_t mc_idx_neg, mc_idx_pos;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Choose between anti-lambda, K0S or just pion pairs */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;
    std::vector<KFParticleMother>* kfFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        kfFoundV0s = &kfAntiLambdas;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
    } else {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        if (pdgV0 == 310) {
            kfFoundV0s = &kfKaonsZeroShort;
            getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
            getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
        } else if (pdgV0 == 422) {
            kfFoundV0s = &kfPionPairs;
            getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromPionPairIdx;
            getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromPionPairIdx;
        }
    }

    /* Loop over all possible pairs of tracks */

    for (Int_t& esdIdxNeg : esdIndicesNegTracks) {
        for (Int_t& esdIdxPos : esdIndicesPosTracks) {

            /* Sanity check */

            if (esdIdxNeg == esdIdxPos) continue;

            /* Get tracks */

            esdTrackNeg = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
            esdTrackPos = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

            /* Kalman Filter */

            KFParticle kfDaughterNeg = CreateKFParticle(*esdTrackNeg, fPDG.GetParticle(pdgTrackNeg)->Mass(), (Int_t)esdTrackNeg->Charge());
            KFParticle kfDaughterPos = CreateKFParticle(*esdTrackPos, fPDG.GetParticle(pdgTrackPos)->Mass(), (Int_t)esdTrackPos->Charge());

            KFParticleMother kfV0;
            kfV0.AddDaughter(kfDaughterNeg);
            kfV0.AddDaughter(kfDaughterPos);

            /* Transport V0 and daughters */

            kfV0.TransportToDecayVertex();

            KFParticle kfTransportedNeg = TransportKFParticle(kfDaughterNeg, kfDaughterPos, pdgTrackNeg, (Int_t)esdTrackNeg->Charge());
            KFParticle kfTransportedPos = TransportKFParticle(kfDaughterPos, kfDaughterNeg, pdgTrackPos, (Int_t)esdTrackPos->Charge());

            /* Reconstruct V0 */

            lvTrackNeg.SetXYZM(kfTransportedNeg.Px(), kfTransportedNeg.Py(), kfTransportedNeg.Pz(), fPDG.GetParticle(pdgTrackNeg)->Mass());
            lvTrackPos.SetXYZM(kfTransportedPos.Px(), kfTransportedPos.Py(), kfTransportedPos.Pz(), fPDG.GetParticle(pdgTrackPos)->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            /* Apply cuts */

            is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                      getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                      getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == pdgV0;
            is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
            is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;
            if (pdgV0 == 422) {
                mc_idx_neg = getMcIdx_fromEsdIdx[esdIdxNeg];
                mc_idx_pos = getMcIdx_fromEsdIdx[esdIdxPos];
                is_signal = isMcIdxSignal[mc_idx_neg] && isMcIdxSignal[mc_idx_pos] && getPdgCode_fromMcIdx[mc_idx_neg] == -211 &&
                            getPdgCode_fromMcIdx[mc_idx_pos] == 211 && !doesMcIdxHaveMother[mc_idx_neg] && !doesMcIdxHaveMother[mc_idx_pos] &&
                            getReactionIdx_fromMcIdx[mc_idx_neg] == getReactionIdx_fromMcIdx[mc_idx_pos];
            }

            if (!PassesV0CutsAs(pdgV0, is_true, is_secondary, is_signal, kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos)) continue;

            /* Calculate V0 properties */

            decay_radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
            decay_length = TMath::Sqrt((kfV0.GetX() - fPrimaryVertex->GetX()) * (kfV0.GetX() - fPrimaryVertex->GetX()) +
                                       (kfV0.GetY() - fPrimaryVertex->GetY()) * (kfV0.GetY() - fPrimaryVertex->GetY()) +
                                       (kfV0.GetZ() - fPrimaryVertex->GetZ()) * (kfV0.GetZ() - fPrimaryVertex->GetZ()));
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            arm_qt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
            arm_alpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());

            cpa_wrt_pv = CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                             fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfV0.GetDistanceFromVertex(kfPrimaryVertex));

            dca_btw_dau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
            dca_neg_v0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
            dca_pos_v0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));

            chi2_ndf = (Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF();

            /* Store V0 */

            kfFoundV0s->push_back(kfV0);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxPos;

            /* Assign branches and fill tree */

            tV0_Idx = (Int_t)kfFoundV0s->size() - 1;
            tV0_Idx_Pos = esdIdxPos;
            tV0_Idx_Neg = esdIdxNeg;
            tV0_Idx_True = getMotherMcIdx_fromEsdIdx[esdIdxNeg];
            tV0_PID = pdgV0;
            tV0_Px = (Float_t)lvV0.Px();
            tV0_Py = (Float_t)lvV0.Py();
            tV0_Pz = (Float_t)lvV0.Pz();
            tV0_E = (Float_t)lvV0.E();
            tV0_Xv_f = (Float_t)kfV0.GetX();
            tV0_Yv_f = (Float_t)kfV0.GetY();
            tV0_Zv_f = (Float_t)kfV0.GetZ();
            tV0_Neg_Px = (Float_t)lvTrackNeg.Px();
            tV0_Neg_Py = (Float_t)lvTrackNeg.Py();
            tV0_Neg_Pz = (Float_t)lvTrackNeg.Pz();
            tV0_Pos_Px = (Float_t)lvTrackPos.Px();
            tV0_Pos_Py = (Float_t)lvTrackPos.Py();
            tV0_Pos_Pz = (Float_t)lvTrackPos.Pz();
            tV0_IsTrue = is_true;
            tV0_IsSecondary = is_secondary;
            tV0_IsSignal = is_signal;
            tV0_ReactionID = getReactionIdx_fromEsdIdx[esdIdxNeg];

            fTree_V0s->Fill();
            ClearV0sBranches();

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(80);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(arm_alpha, arm_qt);

            if (pdgV0 != 422) {

                if (!is_true) continue;

                fHist_V0s_Bookkeep[pdgV0]->Fill(110);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pt")]->Fill(lvV0.Pt());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Px")]->Fill(lvV0.Px());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Py")]->Fill(lvV0.Py());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pz")]->Fill(lvV0.Pz());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Phi")]->Fill(lvV0.Phi());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayRadius")]->Fill(decay_radius);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Zv")]->Fill(kfV0.GetZ());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Eta")]->Fill(lvV0.Eta());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayLength")]->Fill(decay_length);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "OpeningAngle")]->Fill(opening_angle);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmQt")]->Fill(arm_qt);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
                f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(arm_alpha, arm_qt);

                if (!is_secondary) continue;

                fHist_V0s_Bookkeep[pdgV0]->Fill(140);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Mass")]->Fill(lvV0.M());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pt")]->Fill(lvV0.Pt());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Px")]->Fill(lvV0.Px());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Py")]->Fill(lvV0.Py());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pz")]->Fill(lvV0.Pz());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Phi")]->Fill(lvV0.Phi());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayRadius")]->Fill(decay_radius);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Zv")]->Fill(kfV0.GetZ());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Eta")]->Fill(lvV0.Eta());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayLength")]->Fill(decay_length);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "OpeningAngle")]->Fill(opening_angle);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmQt")]->Fill(arm_qt);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
                f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(arm_alpha, arm_qt);
            }

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(170);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(arm_alpha, arm_qt);
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/*
 * Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel A: `AntiSexaquark,N->AntiLambda,K0S`
 * - Uses: `fStruckNucleonPDG`, `fPDG`, `fPrimaryVertex`, `fESD`
 * - Containers: `kf[AntiLambdas|KaonsZeroShort]`, `getEsdIdxOf[Neg|Pos]Dau_from[AntiLambda|KaonZeroShort]Idx`
 */
void AliAnalysisTaskSexaquark::KalmanSexaquarkFinder_ChannelA() {

    AliESDtrack *esdAntiLambdaNegDau, *esdAntiLambdaPosDau;
    AliESDtrack *esdKaonZeroShortNegDau, *esdKaonZeroShortPosDau;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambdaNegDau, lvAntiLambdaPosDau;
    TLorentzVector lvKaonZeroShortNegDau, lvKaonZeroShortPosDau;

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv, dca_wrt_pv;
    Float_t dca_btw_v0s;
    Float_t dca_v0a_sv, dca_v0b_sv;
    Float_t dca_v0a_neg_sv, dca_v0a_pos_sv;
    Float_t dca_v0b_neg_sv, dca_v0b_pos_sv;

    Float_t chi2_ndf;

    Bool_t is_signal;
    Bool_t is_hybrid;

    Int_t idxAntiLambda;
    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t idxKaonZeroShort;
    Int_t esdIdxNeg_KaonZeroShort, esdIdxPos_KaonZeroShort;

    /* Loop over all possible pairs of anti-lambdas and K0S */

    for (idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)kfKaonsZeroShort.size(); idxKaonZeroShort++) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxNeg_KaonZeroShort = getEsdIdxOfNegDau_fromKaonZeroShortIdx[idxKaonZeroShort];
            esdIdxPos_KaonZeroShort = getEsdIdxOfPosDau_fromKaonZeroShortIdx[idxKaonZeroShort];

            if (esdIdxNeg_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxNeg_AntiLambda == esdIdxPos_KaonZeroShort ||
                esdIdxPos_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxPos_AntiLambda == esdIdxPos_KaonZeroShort) {
                continue;
            }

            /* Kalman Filter */

            KFParticleMother kfAntiLambda = kfAntiLambdas[idxAntiLambda];
            KFParticleMother kfKaonZeroShort = kfKaonsZeroShort[idxKaonZeroShort];

            KFParticleMother kfAntiSexaquark;
            kfAntiSexaquark.AddDaughter(kfAntiLambda);
            kfAntiSexaquark.AddDaughter(kfKaonZeroShort);
            kfAntiSexaquark.TransportToDecayVertex();

            /* Transport V0s to the secondary vertex */

            kfAntiLambda.SetProductionVertex(kfAntiSexaquark);
            kfKaonZeroShort.SetProductionVertex(kfAntiSexaquark);

            kfAntiLambda.TransportToProductionVertex();
            kfKaonZeroShort.TransportToProductionVertex();

            /* Get and transport the final-state AliESDtracks */

            esdAntiLambdaNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            esdAntiLambdaPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));
            esdKaonZeroShortNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_KaonZeroShort));
            esdKaonZeroShortPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_KaonZeroShort));

            KFParticle kfAntiLambdaNegDau =
                CreateKFParticle(*esdAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass(), (Int_t)esdAntiLambdaNegDau->Charge());
            KFParticle kfAntiLambdaPosDau =
                CreateKFParticle(*esdAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass(), (Int_t)esdAntiLambdaPosDau->Charge());
            KFParticle kfKaonZeroShortNegDau = CreateKFParticle(*esdKaonZeroShortNegDau, fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass(),
                                                                (Int_t)esdKaonZeroShortNegDau->Charge());
            KFParticle kfKaonZeroShortPosDau = CreateKFParticle(*esdKaonZeroShortPosDau, fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass(),
                                                                (Int_t)esdKaonZeroShortPosDau->Charge());

            kfAntiLambdaNegDau = TransportKFParticle(kfAntiLambdaNegDau, kfAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->PdgCode(),
                                                     (Int_t)esdAntiLambdaNegDau->Charge());
            kfAntiLambdaPosDau = TransportKFParticle(kfAntiLambdaPosDau, kfAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->PdgCode(),
                                                     (Int_t)esdAntiLambdaPosDau->Charge());
            kfKaonZeroShortNegDau =
                TransportKFParticle(kfKaonZeroShortNegDau, kfKaonZeroShortPosDau, fPDG.GetParticle(fFinalStateProductsPDG[2])->PdgCode(),
                                    (Int_t)esdKaonZeroShortNegDau->Charge());
            kfKaonZeroShortPosDau =
                TransportKFParticle(kfKaonZeroShortPosDau, kfKaonZeroShortNegDau, fPDG.GetParticle(fFinalStateProductsPDG[3])->PdgCode(),
                                    (Int_t)esdKaonZeroShortPosDau->Charge());

            lvAntiLambdaNegDau.SetXYZM(kfAntiLambdaNegDau.Px(), kfAntiLambdaNegDau.Py(), kfAntiLambdaNegDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
            lvAntiLambdaPosDau.SetXYZM(kfAntiLambdaPosDau.Px(), kfAntiLambdaPosDau.Py(), kfAntiLambdaPosDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
            lvKaonZeroShortNegDau.SetXYZM(kfKaonZeroShortNegDau.Px(), kfKaonZeroShortNegDau.Py(), kfKaonZeroShortNegDau.Pz(),
                                          fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());
            lvKaonZeroShortPosDau.SetXYZM(kfKaonZeroShortPosDau.Px(), kfKaonZeroShortPosDau.Py(), kfKaonZeroShortPosDau.Pz(),
                                          fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass());

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda = lvAntiLambdaNegDau + lvAntiLambdaPosDau;
            lvKaonZeroShort = lvKaonZeroShortNegDau + lvKaonZeroShortPosDau;

            lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assuming struck nucleon at rest

            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvStruckNucleon;

            /* Apply cuts */

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda] && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] &&
                        isEsdIdxSignal[esdIdxPos_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];
            is_hybrid = (isEsdIdxSignal[esdIdxNeg_AntiLambda] || isEsdIdxSignal[esdIdxPos_AntiLambda] || isEsdIdxSignal[esdIdxNeg_KaonZeroShort] ||
                         isEsdIdxSignal[esdIdxPos_KaonZeroShort]) &&
                        !is_signal;

            if (!PassesSexaquarkCuts_ChannelA(is_signal, kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfKaonZeroShort, kfAntiLambdaNegDau,
                                              kfAntiLambdaPosDau, kfKaonZeroShortNegDau, kfKaonZeroShortPosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
            decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                       (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                       (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                             kfAntiSexaquark.GetZ(),  //
                                             fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
            dca_btw_v0s = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfKaonZeroShort));
            dca_v0a_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_sv = TMath::Abs(kfKaonZeroShort.GetDistanceFromVertex(kfAntiSexaquark));

            dca_v0a_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0a_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_neg_sv = TMath::Abs(kfKaonZeroShortNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_pos_sv = TMath::Abs(kfKaonZeroShortPosDau.GetDistanceFromVertex(kfAntiSexaquark));

            chi2_ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();

            /* Assign branches and fill tree */

            tSexaquark_Idx_AL = idxAntiLambda;
            tSexaquark_Idx_AL_Neg = esdIdxNeg_AntiLambda;
            tSexaquark_Idx_AL_Pos = esdIdxPos_AntiLambda;
            tSexaquark_Px = (Float_t)lvAntiSexaquark.Px();
            tSexaquark_Py = (Float_t)lvAntiSexaquark.Py();
            tSexaquark_Pz = (Float_t)lvAntiSexaquark.Pz();
            tSexaquark_E = (Float_t)lvAntiSexaquark.E();
            tSexaquark_SV_Xv = (Float_t)kfAntiSexaquark.GetX();
            tSexaquark_SV_Yv = (Float_t)kfAntiSexaquark.GetY();
            tSexaquark_SV_Zv = (Float_t)kfAntiSexaquark.GetZ();
            tSexaquark_DecayLength = decay_length;
            tSexaquark_CPAwrtPV = cpa_wrt_pv;
            tSexaquark_DCAwrtPV = dca_wrt_pv;
            tSexaquark_Chi2ndf = chi2_ndf;
            tSexaquark_IsHybrid = is_hybrid;
            tSexaquark_IsSignal = is_signal;
            tSexaquark_ReactionID = getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda];

            tSexaquarkA_Idx_K0S = idxKaonZeroShort;
            tSexaquarkA_Idx_K0S_Neg = esdIdxNeg_KaonZeroShort;
            tSexaquarkA_Idx_K0S_Pos = esdIdxPos_KaonZeroShort;
            tSexaquarkA_DCAbtwV0s = dca_btw_v0s;
            tSexaquarkA_DCAV0aSV = dca_v0a_sv;
            tSexaquarkA_DCAV0bSV = dca_v0b_sv;
            tSexaquarkA_DCAV0anegSV = dca_v0a_neg_sv;
            tSexaquarkA_DCAV0aposSV = dca_v0a_pos_sv;
            tSexaquarkA_DCAV0bnegSV = dca_v0b_neg_sv;
            tSexaquarkA_DCAV0bposSV = dca_v0b_pos_sv;

            fTree_Sexaquarks->Fill();
            ClearSexaquarksBranches();

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(50);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0bposSV")]->Fill(dca_v0b_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(90);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0bposSV")]->Fill(dca_v0b_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);

            // -- DEBUG -- //
            /*
            AliInfo("!! Anti-Sexaquark Found !!");
            AliInfoF("RN: %i, EV: %i, DIR: %i", fRunNumber, fEventNumber, fDirNumber);
            AliInfoF("AL: %i, K0S: %i, ESDs: %i %i %i %i, MC: %i %i %i %i, MC Mothers: %i %i %i %i, Reaction: %i %i %i %i",  //
                     idxAntiLambda, idxKaonZeroShort,                                                                        //
                     esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda, esdIdxNeg_KaonZeroShort,
                     esdIdxPos_KaonZeroShort,  //
                     getMcIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getMcIdx_fromEsdIdx[esdIdxPos_AntiLambda],
                     getMcIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort],
                     getMcIdx_fromEsdIdx[esdIdxPos_KaonZeroShort],  //
                     getMotherMcIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getMotherMcIdx_fromEsdIdx[esdIdxPos_AntiLambda],
                     getMotherMcIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort], getMotherMcIdx_fromEsdIdx[esdIdxPos_KaonZeroShort],  //
                     getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda],
                     getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort], getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort]);
            */
            // -- END OF DEBUG -- //
        }  // end of loop over K0S
    }      // end of loop over anti-lambdas
}

/*
 * Apply anti-sexaquark candidate cuts. Exclusive to reaction channel A: `AntiSexaquark,N->AntiLambda,K0S`.
 * - Use: `fPrimaryVertex`
 * - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort|AntiLambdaNeg|AntiLambdaPos|KaonZeroShortNeg|KaonZeroShortPos]`, `lvAntiSexaquark`
 * - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg,
                                                              KFParticle kfAntiLambdaPos, KFParticle kfKaonZeroShortNeg,
                                                              KFParticle kfKaonZeroShortPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    Double_t dca_btw_v0s = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfKaonZeroShort));
    if (kMax_Sexa_DCAbtwV0s && dca_btw_v0s > kMax_Sexa_DCAbtwV0s) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0a_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0aSV && dca_v0a_sv > kMax_Sexa_DCAv0aSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_v0b_sv = TMath::Abs(kfKaonZeroShort.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bSV && dca_v0b_sv > kMax_Sexa_DCAv0bSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Double_t dca_v0a_neg_sv = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0anegSV && dca_v0a_neg_sv > kMax_Sexa_DCAv0anegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Double_t dca_v0a_pos_sv = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0aposSV && dca_v0a_pos_sv > kMax_Sexa_DCAv0aposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Double_t dca_v0b_neg_sv = TMath::Abs(kfKaonZeroShortNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bnegSV && dca_v0b_neg_sv > kMax_Sexa_DCAv0bnegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Double_t dca_v0b_pos_sv = TMath::Abs(kfKaonZeroShortPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bposSV && dca_v0b_pos_sv > kMax_Sexa_DCAv0bposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(36);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(76);

    return kTRUE;
}

/*
 * Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel D: `AntiSexaquark,P->AntiLambda,K+`
 * - Uses: `fStruckNucleonPDG`, `fPDG`, `fPrimaryVertex`, `fESD`
 * - Containers: `kfAntiLambdas`, `esdIndicesOfPosKaonTracks`, `getEsdIdxOf[Neg|Pos]Dau_fromAntiLambdaIdx`
 */
void AliAnalysisTaskSexaquark::KalmanSexaquarkFinder_ChannelD() {

    AliESDtrack *esdAntiLambdaNegDau, *esdAntiLambdaPosDau;
    AliESDtrack* esdPosKaon;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvAntiLambdaNegDau, lvAntiLambdaPosDau;
    TLorentzVector lvPosKaon;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv, dca_wrt_pv;
    Float_t dca_v0_ba;
    Float_t dca_v0_neg_ba, dca_v0_pos_ba;
    Float_t dca_v0_sv;
    Float_t dca_v0_neg_sv, dca_v0_pos_sv;
    Float_t dca_ba_sv;

    Float_t chi2_ndf;

    Bool_t is_signal;
    Bool_t is_hybrid;

    Int_t idxAntiLambda;
    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxPosKaon;

    /* Loop over all possible pairs of anti-lambdas and positive kaons */

    for (idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxPosKaon = 0; idxPosKaon < (Int_t)esdIndicesOfPosKaonTracks.size(); idxPosKaon++) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxPosKaon = esdIndicesOfPosKaonTracks[idxPosKaon];
            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];

            if (esdIdxNeg_AntiLambda == esdIdxPosKaon || esdIdxPosKaon == esdIdxPos_AntiLambda) {
                continue;
            }

            /* Kalman Filter */

            KFParticleMother kfAntiLambda = kfAntiLambdas[idxAntiLambda];

            esdPosKaon = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaon));
            KFParticle kfPosKaon = CreateKFParticle(*esdPosKaon, fPDG.GetParticle(fProductsPDG[1])->Mass(), (Int_t)esdPosKaon->Charge());

            KFParticleMother kfAntiSexaquark;
            kfAntiSexaquark.AddDaughter(kfAntiLambda);
            kfAntiSexaquark.AddDaughter(kfPosKaon);
            kfAntiSexaquark.TransportToDecayVertex();

            /* Transport V0s to the secondary vertex */

            kfAntiLambda.SetProductionVertex(kfAntiSexaquark);
            kfPosKaon.SetProductionVertex(kfAntiSexaquark);

            kfAntiLambda.TransportToProductionVertex();
            kfPosKaon.TransportToProductionVertex();

            /* Get and transport AliESDtracks of anti-lambda daughters */

            esdAntiLambdaNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            esdAntiLambdaPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));

            KFParticle kfAntiLambdaNegDau =
                CreateKFParticle(*esdAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass(), (Int_t)esdAntiLambdaNegDau->Charge());
            KFParticle kfAntiLambdaPosDau =
                CreateKFParticle(*esdAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass(), (Int_t)esdAntiLambdaPosDau->Charge());

            kfAntiLambdaNegDau = TransportKFParticle(kfAntiLambdaNegDau, kfAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->PdgCode(),
                                                     (Int_t)esdAntiLambdaNegDau->Charge());
            kfAntiLambdaPosDau = TransportKFParticle(kfAntiLambdaPosDau, kfAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->PdgCode(),
                                                     (Int_t)esdAntiLambdaPosDau->Charge());
            kfPosKaon =
                TransportKFParticle(kfPosKaon, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[2])->PdgCode(), (Int_t)esdPosKaon->Charge());

            lvAntiLambdaNegDau.SetXYZM(kfAntiLambdaNegDau.Px(), kfAntiLambdaNegDau.Py(), kfAntiLambdaNegDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
            lvAntiLambdaPosDau.SetXYZM(kfAntiLambdaPosDau.Px(), kfAntiLambdaPosDau.Py(), kfAntiLambdaPosDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
            lvPosKaon.SetXYZM(kfPosKaon.Px(), kfPosKaon.Py(), kfPosKaon.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda = lvAntiLambdaNegDau + lvAntiLambdaPosDau;

            lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assuming struck nucleon at rest

            lvAntiSexaquark = lvAntiLambda + lvPosKaon - lvStruckNucleon;

            /* Apply cuts */

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda] && isEsdIdxSignal[esdIdxPosKaon] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPosKaon];
            is_hybrid = (isEsdIdxSignal[esdIdxNeg_AntiLambda] || isEsdIdxSignal[esdIdxPos_AntiLambda] || isEsdIdxSignal[esdIdxPosKaon]) && !is_signal;

            if (!PassesSexaquarkCuts_ChannelD(is_signal, kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfPosKaon, kfAntiLambdaNegDau,
                                              kfAntiLambdaPosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
            decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                       (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                       (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                             kfAntiSexaquark.GetZ(),  //
                                             fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));

            dca_v0_ba = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPosKaon));
            dca_v0_neg_ba = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfPosKaon));
            dca_v0_pos_ba = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfPosKaon));

            dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
            dca_ba_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));

            chi2_ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();  // exclusive to KF

            /* Assign branches and fill trees */

            tSexaquark_Idx_AL = idxAntiLambda;
            tSexaquark_Idx_AL_Neg = esdIdxNeg_AntiLambda;
            tSexaquark_Idx_AL_Pos = esdIdxPos_AntiLambda;
            tSexaquark_Px = (Float_t)lvAntiSexaquark.Px();
            tSexaquark_Py = (Float_t)lvAntiSexaquark.Py();
            tSexaquark_Pz = (Float_t)lvAntiSexaquark.Pz();
            tSexaquark_E = (Float_t)lvAntiSexaquark.E();
            tSexaquark_SV_Xv = (Float_t)kfAntiSexaquark.GetX();
            tSexaquark_SV_Yv = (Float_t)kfAntiSexaquark.GetY();
            tSexaquark_SV_Zv = (Float_t)kfAntiSexaquark.GetZ();
            tSexaquark_DecayLength = decay_length;
            tSexaquark_CPAwrtPV = cpa_wrt_pv;
            tSexaquark_DCAwrtPV = dca_wrt_pv;
            tSexaquark_Chi2ndf = chi2_ndf;
            tSexaquark_IsHybrid = is_hybrid;
            tSexaquark_IsSignal = is_signal;
            tSexaquark_ReactionID = getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda];

            tSexaquarkD_Idx_PosKaon = esdIdxPosKaon;
            tSexaquarkD_DCAV0Ba = dca_v0_ba;
            tSexaquarkD_DCAV0negBa = dca_v0_neg_ba;
            tSexaquarkD_DCAV0posBa = dca_v0_pos_ba;
            tSexaquarkD_DCAV0SV = dca_v0_sv;
            tSexaquarkD_DCABaSV = dca_ba_sv;
            tSexaquarkD_DCAV0negSV = dca_v0_neg_sv;
            tSexaquarkD_DCAV0posSV = dca_v0_pos_sv;

            fTree_Sexaquarks->Fill();
            ClearSexaquarksBranches();

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(50);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0Ba")]->Fill(dca_v0_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0negBa")]->Fill(dca_v0_neg_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0posBa")]->Fill(dca_v0_pos_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0SV")]->Fill(dca_v0_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCABaSV")]->Fill(dca_ba_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0negSV")]->Fill(dca_v0_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0posSV")]->Fill(dca_v0_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(90);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0Ba")]->Fill(dca_v0_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0negBa")]->Fill(dca_v0_neg_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0posBa")]->Fill(dca_v0_pos_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0SV")]->Fill(dca_v0_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCABaSV")]->Fill(dca_ba_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0negSV")]->Fill(dca_v0_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0posSV")]->Fill(dca_v0_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);
        }  // end of loop over pos. kaons
    }      // end of loop over anti-lambdas
}

/*
 * Apply anti-sexaquark candidate cuts.
 * - Use: `fPrimaryVertex`
 * - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort|AntiLambdaNeg|AntiLambdaPos|KaonZeroShortNeg|KaonZeroShortPos]`, `lvAntiSexaquark`
 * - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelD(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfPosKaon, KFParticle kfAntiLambdaNeg,
                                                              KFParticle kfAntiLambdaPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    Double_t dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0SV && dca_v0_sv > kMax_Sexa_DCAv0SV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0_ba = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfPosKaon));
    if (kMax_Sexa_DCAv0ba && dca_v0_ba > kMax_Sexa_DCAv0ba) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_ba_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAbaSV && dca_ba_sv > kMax_Sexa_DCAbaSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Double_t dca_v0_neg_sv = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0negSV && dca_v0_neg_sv > kMax_Sexa_DCAv0negSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Double_t dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0posSV && dca_v0_pos_sv > kMax_Sexa_DCAv0posSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Double_t dca_v0_neg_ba = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfPosKaon));
    if (kMax_Sexa_DCAv0negba && dca_v0_neg_ba > kMax_Sexa_DCAv0negba) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Double_t dca_v0_pos_ba = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfPosKaon));
    if (kMax_Sexa_DCAv0posba && dca_v0_pos_ba > kMax_Sexa_DCAv0posba) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(36);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(76);

    return kTRUE;
}

/*
 * Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel E: `AntiSexaquark,P->AntiLambda,K+,pi-,pi+`
 */
void AliAnalysisTaskSexaquark::KalmanSexaquarkFinder_ChannelE() {

    AliESDtrack *esdAntiLambdaNegDau, *esdAntiLambdaPosDau;
    AliESDtrack *esdPiMinus, *esdPiPlus;
    AliESDtrack* esdPosKaon;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvAntiLambdaNegDau, lvAntiLambdaPosDau;

    TLorentzVector lvPiMinus;
    TLorentzVector lvPiPlus;
    TLorentzVector lvPosKaon;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv, dca_wrt_pv;
    Float_t dca_v0_sv, dca_v0_neg_sv, dca_v0_pos_sv;
    Float_t dca_pk_sv, dca_pm_sv, dca_pp_sv;
    Float_t dca_pk_v0, dca_pm_v0, dca_pp_v0;
    Float_t dca_pm_pk, dca_pp_pk;

    Float_t chi2_ndf;

    Bool_t is_signal;
    Bool_t is_hybrid;

    Int_t idxAntiLambda;
    Int_t idxPionPair;

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxPos_PionPair, esdIdxNeg_PionPair;
    Int_t esdIdxPosKaon;

    /* Loop over all possible groups of anti-lambdas, pion pairs and pos. kaons */

    for (idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (idxPionPair = 0; idxPionPair < (Int_t)kfPionPairs.size(); idxPionPair++) {
            for (Int_t idxPosKaon = 0; idxPosKaon < (Int_t)esdIndicesOfPosKaonTracks.size(); idxPosKaon++) {

                /* Check that no daughter is repeated */

                esdIdxPosKaon = esdIndicesOfPosKaonTracks[idxPosKaon];
                esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
                esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
                esdIdxNeg_PionPair = getEsdIdxOfNegDau_fromPionPairIdx[idxPionPair];
                esdIdxPos_PionPair = getEsdIdxOfPosDau_fromPionPairIdx[idxPionPair];

                std::set<Int_t> uniqueIndices = {esdIdxPosKaon, esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda, esdIdxNeg_PionPair, esdIdxPos_PionPair};
                if (uniqueIndices.size() != 5) continue;

                /* Kalman Filter */

                KFParticleMother kfAntiLambda = kfAntiLambdas[idxAntiLambda];

                esdPosKaon = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaon));
                KFParticle kfPosKaon = CreateKFParticle(*esdPosKaon, fPDG.GetParticle(fProductsPDG[1])->Mass(), (Int_t)esdPosKaon->Charge());

                esdPiMinus = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_PionPair));
                KFParticle kfPiMinus = CreateKFParticle(*esdPiMinus, fPDG.GetParticle(fProductsPDG[2])->Mass(), (Int_t)esdPiMinus->Charge());

                esdPiPlus = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_PionPair));
                KFParticle kfPiPlus = CreateKFParticle(*esdPiPlus, fPDG.GetParticle(fProductsPDG[3])->Mass(), (Int_t)esdPiPlus->Charge());

                KFParticleMother kfAntiSexaquark;
                kfAntiSexaquark.AddDaughter(kfAntiLambda);
                kfAntiSexaquark.AddDaughter(kfPosKaon);
                kfAntiSexaquark.AddDaughter(kfPiMinus);
                kfAntiSexaquark.AddDaughter(kfPiPlus);
                kfAntiSexaquark.TransportToDecayVertex();

                /* Transport V0s to the secondary vertex */

                kfAntiLambda.SetProductionVertex(kfAntiSexaquark);
                kfPosKaon.SetProductionVertex(kfAntiSexaquark);
                kfPiMinus.SetProductionVertex(kfAntiSexaquark);
                kfPiPlus.SetProductionVertex(kfAntiSexaquark);

                kfAntiLambda.TransportToProductionVertex();
                kfPosKaon.TransportToProductionVertex();
                kfPiMinus.TransportToProductionVertex();
                kfPiPlus.TransportToProductionVertex();

                /* Transport KFParticles of final-state particles */

                esdAntiLambdaNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
                esdAntiLambdaPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));

                KFParticle kfAntiLambdaNegDau =
                    CreateKFParticle(*esdAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass(), (Int_t)esdAntiLambdaNegDau->Charge());
                KFParticle kfAntiLambdaPosDau =
                    CreateKFParticle(*esdAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass(), (Int_t)esdAntiLambdaPosDau->Charge());

                kfAntiLambdaNegDau =
                    TransportKFParticle(kfAntiLambdaNegDau, kfAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->PdgCode(),
                                        (Int_t)esdAntiLambdaNegDau->Charge());
                kfAntiLambdaPosDau =
                    TransportKFParticle(kfAntiLambdaPosDau, kfAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->PdgCode(),
                                        (Int_t)esdAntiLambdaPosDau->Charge());
                kfPosKaon =
                    TransportKFParticle(kfPosKaon, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[2])->PdgCode(), (Int_t)esdPosKaon->Charge());
                kfPiMinus =
                    TransportKFParticle(kfPiMinus, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[3])->PdgCode(), (Int_t)esdPiMinus->Charge());
                kfPiPlus =
                    TransportKFParticle(kfPiPlus, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[4])->PdgCode(), (Int_t)esdPiPlus->Charge());

                lvAntiLambdaNegDau.SetXYZM(kfAntiLambdaNegDau.Px(), kfAntiLambdaNegDau.Py(), kfAntiLambdaNegDau.Pz(),
                                           fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
                lvAntiLambdaPosDau.SetXYZM(kfAntiLambdaPosDau.Px(), kfAntiLambdaPosDau.Py(), kfAntiLambdaPosDau.Pz(),
                                           fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
                lvPosKaon.SetXYZM(kfPosKaon.Px(), kfPosKaon.Py(), kfPosKaon.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());
                lvPiMinus.SetXYZM(kfPiMinus.Px(), kfPiMinus.Py(), kfPiMinus.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass());
                lvPiPlus.SetXYZM(kfPiPlus.Px(), kfPiPlus.Py(), kfPiPlus.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[4])->Mass());

                /* Reconstruct anti-sexaquark candidate */

                lvAntiLambda = lvAntiLambdaNegDau + lvAntiLambdaPosDau;

                lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assuming struck nucleon at rest

                lvAntiSexaquark = lvAntiLambda + lvPosKaon + lvPiMinus + lvPiPlus - lvStruckNucleon;

                /* Apply cuts */

                std::set<UInt_t> uniqueReactionID = {getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda],
                                                     getReactionIdx_fromEsdIdx[esdIdxPosKaon], getReactionIdx_fromEsdIdx[esdIdxNeg_PionPair],
                                                     getReactionIdx_fromEsdIdx[esdIdxPos_PionPair]};
                is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda] && isEsdIdxSignal[esdIdxPosKaon] &&
                            isEsdIdxSignal[esdIdxNeg_PionPair] && isEsdIdxSignal[esdIdxPos_PionPair] && (Int_t)uniqueReactionID.size() == 1;
                is_hybrid = (isEsdIdxSignal[esdIdxNeg_AntiLambda] || isEsdIdxSignal[esdIdxPos_AntiLambda] || isEsdIdxSignal[esdIdxPosKaon] ||
                             isEsdIdxSignal[esdIdxNeg_PionPair] || isEsdIdxSignal[esdIdxPos_PionPair]) &&
                            !is_signal;

                if (!PassesSexaquarkCuts_ChannelE(is_signal, kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfAntiLambdaNegDau, kfAntiLambdaPosDau,
                                                  kfPosKaon, kfPiMinus, kfPiPlus)) {
                    continue;
                }

                /* Get anti-sexaquark's properties */

                radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
                decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                           (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                           (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
                cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                                 kfAntiSexaquark.GetZ(),  //
                                                 fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
                dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));

                dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
                dca_v0_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
                dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
                dca_pk_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
                dca_pm_sv = TMath::Abs(kfPiMinus.GetDistanceFromVertex(kfAntiSexaquark));
                dca_pp_sv = TMath::Abs(kfPiPlus.GetDistanceFromVertex(kfAntiSexaquark));

                dca_pk_v0 = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPosKaon));
                dca_pm_v0 = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiMinus));
                dca_pp_v0 = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiPlus));

                dca_pm_pk = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiMinus));
                dca_pp_pk = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiPlus));

                chi2_ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();  // exclusive to KF method

                /* Assign branches and fill trees */

                tSexaquark_Idx_AL = idxAntiLambda;
                tSexaquark_Idx_AL_Neg = esdIdxNeg_AntiLambda;
                tSexaquark_Idx_AL_Pos = esdIdxPos_AntiLambda;
                tSexaquark_Px = (Float_t)lvAntiSexaquark.Px();
                tSexaquark_Py = (Float_t)lvAntiSexaquark.Py();
                tSexaquark_Pz = (Float_t)lvAntiSexaquark.Pz();
                tSexaquark_E = (Float_t)lvAntiSexaquark.E();
                tSexaquark_SV_Xv = (Float_t)kfAntiSexaquark.GetX();
                tSexaquark_SV_Yv = (Float_t)kfAntiSexaquark.GetY();
                tSexaquark_SV_Zv = (Float_t)kfAntiSexaquark.GetZ();
                tSexaquark_DecayLength = decay_length;
                tSexaquark_CPAwrtPV = cpa_wrt_pv;
                tSexaquark_DCAwrtPV = dca_wrt_pv;
                tSexaquark_Chi2ndf = chi2_ndf;
                tSexaquark_IsHybrid = is_hybrid;
                tSexaquark_IsSignal = is_signal;
                tSexaquark_ReactionID = getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda];

                tSexaquarkE_Idx_PP = idxPionPair;
                tSexaquarkE_Idx_PP_Neg = esdIdxNeg_PionPair;
                tSexaquarkE_Idx_PP_Pos = esdIdxPos_PionPair;
                tSexaquarkE_Idx_PosKaon = esdIdxPosKaon;
                tSexaquarkE_DCAV0SV = dca_v0_sv;
                tSexaquarkE_DCAV0negSV = dca_v0_neg_sv;
                tSexaquarkE_DCAV0posSV = dca_v0_pos_sv;
                tSexaquarkE_DCApkSV = dca_pk_sv;
                tSexaquarkE_DCApmSV = dca_pm_sv;
                tSexaquarkE_DCAppSV = dca_pp_sv;
                tSexaquarkE_DCApkV0 = dca_pk_v0;
                tSexaquarkE_DCApmV0 = dca_pm_v0;
                tSexaquarkE_DCAppV0 = dca_pp_v0;
                tSexaquarkE_DCApmPK = dca_pm_pk;
                tSexaquarkE_DCAppPK = dca_pp_pk;

                fTree_Sexaquarks->Fill();
                ClearSexaquarksBranches();

                /* Fill histograms */

                fHist_AntiSexaquarks_Bookkeep->Fill(50);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(kfAntiSexaquark.GetZ());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0SV")]->Fill(dca_v0_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0negSV")]->Fill(dca_v0_neg_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAV0posSV")]->Fill(dca_v0_pos_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApkSV")]->Fill(dca_pk_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApmSV")]->Fill(dca_pm_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAppSV")]->Fill(dca_pp_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApkV0")]->Fill(dca_pk_v0);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApmV0")]->Fill(dca_pm_v0);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAppV0")]->Fill(dca_pp_v0);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApmPK")]->Fill(dca_pm_pk);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAppPK")]->Fill(dca_pp_pk);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

                if (!is_signal) continue;

                fHist_AntiSexaquarks_Bookkeep->Fill(90);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0SV")]->Fill(dca_v0_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0negSV")]->Fill(dca_v0_neg_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAV0posSV")]->Fill(dca_v0_pos_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApkSV")]->Fill(dca_pk_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApmSV")]->Fill(dca_pm_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAppSV")]->Fill(dca_pp_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApkV0")]->Fill(dca_pk_v0);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApmV0")]->Fill(dca_pm_v0);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAppV0")]->Fill(dca_pp_v0);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApmPK")]->Fill(dca_pm_pk);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAppPK")]->Fill(dca_pp_pk);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);

                // -- DEBUG -- //
                // clang-format off
                /*
                AliInfo("!! Anti-Sexaquark Found !!");
                AliInfoF("RN: %i, EV: %i, DIR: %i", fRunNumber, fEventNumber, fDirNumber);
                AliInfoF("AL: %i, PK: %i, PP: %i, ESDs: %i %i %i %i %i, MC: %i %i %i %i %i, MC Mothers: %i %i %i %i %i, Reaction: %i %i %i %i %i",
                         idxAntiLambda, esdIdxPosKaon, idxPionPair,
                         esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda, esdIdxPosKaon, esdIdxNeg_PionPair, esdIdxPos_PionPair,
                         getMcIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getMcIdx_fromEsdIdx[esdIdxPos_AntiLambda], getMcIdx_fromEsdIdx[esdIdxPosKaon],
                         getMcIdx_fromEsdIdx[esdIdxNeg_PionPair], getMcIdx_fromEsdIdx[esdIdxPos_PionPair], //
                         getMotherMcIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getMotherMcIdx_fromEsdIdx[esdIdxPos_AntiLambda], getMotherMcIdx_fromEsdIdx[esdIdxPosKaon],
                         getMotherMcIdx_fromEsdIdx[esdIdxNeg_PionPair], getMotherMcIdx_fromEsdIdx[esdIdxPos_PionPair],  //
                         getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda], getReactionIdx_fromEsdIdx[esdIdxPosKaon],
                         getReactionIdx_fromEsdIdx[esdIdxNeg_PionPair], getReactionIdx_fromEsdIdx[esdIdxPos_PionPair]);
                */
                // clang-format on
                // -- END OF DEBUG -- //
            }  // end of loop over K+
        }      // end of loop over pi+pi- pairs
    }          // end of loop over anti-lambdas
}

/*
 * Apply anti-sexaquark candidate cuts via channel E reconstruction.
 * - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelE(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfAntiLambdaNegDau, KFParticle kfAntiLambdaPosDau,
                                                              KFParticle kfPosKaon, KFParticle kfPiMinus, KFParticle kfPiPlus) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    Double_t dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0SV && dca_v0_sv > kMax_Sexa_DCAv0SV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0posSV && dca_v0_pos_sv > kMax_Sexa_DCAv0posSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_pk_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCApkSV && dca_pk_sv > kMax_Sexa_DCApkSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Double_t dca_pm_sv = TMath::Abs(kfPiMinus.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCApmSV && dca_pm_sv > kMax_Sexa_DCApmSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Double_t dca_pp_sv = TMath::Abs(kfPiPlus.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAppSV && dca_pp_sv > kMax_Sexa_DCAppSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Double_t dca_v0_pk = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPosKaon));
    if (kMax_Sexa_DCAv0pk && dca_v0_pk > kMax_Sexa_DCAv0pk) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Double_t dca_v0_pm = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiMinus));
    if (kMax_Sexa_DCAv0pm && dca_v0_pm > kMax_Sexa_DCAv0pm) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    Double_t dca_v0_pp = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiPlus));
    if (kMax_Sexa_DCAv0pp && dca_v0_pp > kMax_Sexa_DCAv0pp) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(36);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(76);

    Double_t dca_pk_pm = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiMinus));
    if (kMax_Sexa_DCApkpm && dca_pk_pm > kMax_Sexa_DCApkpm) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(37);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(77);

    Double_t dca_pk_pp = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiPlus));
    if (kMax_Sexa_DCApkpp && dca_pk_pp > kMax_Sexa_DCApkpp) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(38);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(78);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(39);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(79);

    return kTRUE;
}

/*                                */
/**  Channel H -- Kalman Filter  **/
/*** ========================== ***/

/*
 * Find all K+K+ pairs via Kalman Filter.
 * - Uses: `fESD`, `esdIndicesOfPosKaonTracks`, `fHist_PosKaonPairs`
 */
void AliAnalysisTaskSexaquark::KalmanPosKaonPairFinder() {

    const Int_t kKaonPDGCode = 321;
    const Float_t kKaonMass = fPDG.GetParticle(kKaonPDGCode)->Mass();

    AliESDtrack* esdPosKaonA;
    AliESDtrack* esdPosKaonB;

    /* Declare TLorentzVectors */

    TLorentzVector lvPosKaonA, lvPosKaonB;
    TLorentzVector lvPosKaonPair;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t origin_radius;
    Float_t decay_length;  // origin radius - primary vertex, i.e., the decay length of the mother particle
    Float_t dca_wrt_pv;
    Float_t dca_btw_kaons;
    Float_t dca_ka_sv, dca_kb_sv;
    Float_t opening_angle;

    Float_t chi2_ndf;

    Bool_t is_signal;
    Bool_t is_hybrid;

    Int_t idxKaonPair;
    Int_t esdIdxPosKaonA, esdIdxPosKaonB;

    /* Loop over all possible pairs of tracks */

    idxKaonPair = 0;
    for (Int_t i = 0; i < (Int_t)esdIndicesOfPosKaonTracks.size() - 1; i++) {
        for (Int_t j = i + 1; j < (Int_t)esdIndicesOfPosKaonTracks.size(); j++) {

            esdIdxPosKaonA = esdIndicesOfPosKaonTracks[i];
            esdIdxPosKaonB = esdIndicesOfPosKaonTracks[j];

            /* Sanity check */

            if (esdIdxPosKaonA == esdIdxPosKaonB) continue;

            /* Get tracks */

            esdPosKaonA = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaonA));
            esdPosKaonB = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaonB));

            /* Kalman Filter */

            KFParticle kfPosKaonA = CreateKFParticle(*esdPosKaonA, kKaonMass, (Int_t)esdPosKaonA->Charge());
            KFParticle kfPosKaonB = CreateKFParticle(*esdPosKaonB, kKaonMass, (Int_t)esdPosKaonB->Charge());

            KFParticleMother kfPosKaonPair;
            kfPosKaonPair.AddDaughter(kfPosKaonA);
            kfPosKaonPair.AddDaughter(kfPosKaonB);

            /* Transport V0 and daughters */

            kfPosKaonPair.TransportToDecayVertex();

            KFParticle kfTransportedPosKaonA = TransportKFParticle(kfPosKaonA, kfPosKaonB, kKaonPDGCode, (Int_t)esdPosKaonA->Charge());
            KFParticle kfTransportedPosKaonB = TransportKFParticle(kfPosKaonB, kfPosKaonA, kKaonPDGCode, (Int_t)esdPosKaonB->Charge());

            /* Reconstruct V0 */

            lvPosKaonA.SetXYZM(kfTransportedPosKaonA.Px(), kfTransportedPosKaonA.Py(), kfTransportedPosKaonA.Pz(), kKaonMass);
            lvPosKaonB.SetXYZM(kfTransportedPosKaonB.Px(), kfTransportedPosKaonB.Py(), kfTransportedPosKaonB.Pz(), kKaonMass);
            lvPosKaonPair = lvPosKaonA + lvPosKaonB;

            /* Apply cuts */

            is_signal = isEsdIdxSignal[esdIdxPosKaonA] && isEsdIdxSignal[esdIdxPosKaonB] &&
                        getReactionIdx_fromEsdIdx[esdIdxPosKaonA] == getReactionIdx_fromEsdIdx[esdIdxPosKaonB];
            is_hybrid = (isEsdIdxSignal[esdIdxPosKaonA] || isEsdIdxSignal[esdIdxPosKaonB]) && !is_signal;

            if (!PassesPosKaonPairCuts(is_signal, kfPosKaonPair, kfPosKaonA, kfPosKaonB, lvPosKaonPair, lvPosKaonA, lvPosKaonB)) {
                continue;
            }

            /* Calculate variables */

            origin_radius = TMath::Sqrt(kfPosKaonPair.GetX() * kfPosKaonPair.GetX() + kfPosKaonPair.GetY() * kfPosKaonPair.GetY());

            dca_wrt_pv = TMath::Abs(kfPosKaonPair.GetDistanceFromVertex(kfPrimaryVertex));
            dca_btw_kaons = TMath::Abs(kfPosKaonA.GetDistanceFromParticle(kfPosKaonB));
            dca_ka_sv = TMath::Abs(kfPosKaonA.GetDistanceFromVertex(kfPosKaonPair));
            dca_kb_sv = TMath::Abs(kfPosKaonB.GetDistanceFromVertex(kfPosKaonPair));
            decay_length = TMath::Sqrt((kfPosKaonPair.GetX() - fPrimaryVertex->GetX()) * (kfPosKaonPair.GetX() - fPrimaryVertex->GetX()) +
                                       (kfPosKaonPair.GetY() - fPrimaryVertex->GetY()) * (kfPosKaonPair.GetY() - fPrimaryVertex->GetY()) +
                                       (kfPosKaonPair.GetZ() - fPrimaryVertex->GetZ()) * (kfPosKaonPair.GetZ() - fPrimaryVertex->GetZ()));
            opening_angle = lvPosKaonA.Angle(lvPosKaonB.Vect());

            chi2_ndf = (Double_t)kfPosKaonPair.GetChi2() / (Double_t)kfPosKaonPair.GetNDF();

            /* Assign branches and fill tree */

            tKaonPair_Idx = idxKaonPair;
            tKaonPair_Idx_KaonA = esdIdxPosKaonA;
            tKaonPair_Idx_KaonB = esdIdxPosKaonB;
            tKaonPair_Px = (Float_t)lvPosKaonPair.Px();
            tKaonPair_Py = (Float_t)lvPosKaonPair.Py();
            tKaonPair_Pz = (Float_t)lvPosKaonPair.Pz();
            tKaonPair_E = (Float_t)lvPosKaonPair.E();
            tKaonPair_SV_Xv = (Float_t)kfPosKaonPair.GetX();
            tKaonPair_SV_Yv = (Float_t)kfPosKaonPair.GetY();
            tKaonPair_SV_Zv = (Float_t)kfPosKaonPair.GetZ();
            tKaonPair_DecayLength = decay_length;
            tKaonPair_DCAwrtPV = dca_wrt_pv;
            tKaonPair_DCAbtwKaons = dca_btw_kaons;
            tKaonPair_DCAkaSV = dca_ka_sv;
            tKaonPair_DCAkbSV = dca_kb_sv;
            tKaonPair_OpeningAngle = opening_angle;
            tKaonPair_Chi2ndf = chi2_ndf;
            tSexaquark_IsHybrid = is_hybrid;
            tSexaquark_IsSignal = is_signal;
            tKaonPair_ReactionID = getReactionIdx_fromEsdIdx[esdIdxPosKaonA];

            fTree_KaonPairs->Fill();
            ClearKaonPairsBranches();

            /* Fill histograms */

            fHist_PosKaonPairs_Bookkeep->Fill(50);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Mass")]->Fill(lvPosKaonPair.M());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "OriginRadius")]->Fill(origin_radius);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Xv")]->Fill(kfPosKaonPair.GetX());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Yv")]->Fill(kfPosKaonPair.GetY());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Zv")]->Fill(kfPosKaonPair.GetZ());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Eta")]->Fill(lvPosKaonPair.Eta());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvPosKaonPair.Rapidity());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Px")]->Fill(lvPosKaonPair.Px());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Py")]->Fill(lvPosKaonPair.Py());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Pz")]->Fill(lvPosKaonPair.Pz());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Pt")]->Fill(lvPosKaonPair.Pt());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "OpeningAngle")]->Fill(opening_angle);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCAbtwKaons")]->Fill(dca_btw_kaons);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCAkaSV")]->Fill(dca_ka_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCAkbSV")]->Fill(dca_kb_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

            if (!is_signal) continue;

            fHist_PosKaonPairs_Bookkeep->Fill(90);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvPosKaonPair.M());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "OriginRadius")]->Fill(origin_radius);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Xv")]->Fill(kfPosKaonPair.GetX());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Yv")]->Fill(kfPosKaonPair.GetY());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfPosKaonPair.GetZ());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvPosKaonPair.Eta());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvPosKaonPair.Rapidity());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Px")]->Fill(lvPosKaonPair.Px());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Py")]->Fill(lvPosKaonPair.Py());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvPosKaonPair.Pz());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvPosKaonPair.Pt());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "OpeningAngle")]->Fill(opening_angle);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCAbtwKaons")]->Fill(dca_btw_kaons);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCAkaSV")]->Fill(dca_ka_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCAkbSV")]->Fill(dca_kb_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);

            idxKaonPair++;
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*
 * Apply cuts to a K+K+ pair
 * - Input: `kfPosKaonPair`, `kfPosKaonA`, `kfPosKaonB`, `lvPosKaonPair`, `lvPosKaonA`, `lvPosKaonB`
 * - Uses: `fPrimaryVertex`
 * - Return: `kTRUE` if the K+K+ pair passes the cuts, `kFALSE` otherwise
 */
Bool_t AliAnalysisTaskSexaquark::PassesPosKaonPairCuts(Bool_t isSignal, KFParticleMother kfPosKaonPair, KFParticle kfPosKaonA, KFParticle kfPosKaonB,
                                                       TLorentzVector lvPosKaonPair, TLorentzVector lvPosKaonA, TLorentzVector lvPosKaonB) {

    fHist_PosKaonPairs_Bookkeep->Fill(20);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(60);

    Double_t Radius = TMath::Sqrt(kfPosKaonPair.GetX() * kfPosKaonPair.GetX() + kfPosKaonPair.GetY() * kfPosKaonPair.GetY());
    if (kMin_PosKaonPair_Radius && Radius < kMin_PosKaonPair_Radius) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(21);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(61);

    Double_t Pt = lvPosKaonPair.Pt();
    if (kMin_PosKaonPair_Pt && Pt < kMin_PosKaonPair_Pt) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(22);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(62);

    Double_t DCAwrtPV = LinePointDCA(lvPosKaonPair.Px(), lvPosKaonPair.Py(), lvPosKaonPair.Pz(), kfPosKaonPair.GetX(), kfPosKaonPair.GetY(),
                                     kfPosKaonPair.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_PosKaonPair_DCAwrtPV && DCAwrtPV < kMin_PosKaonPair_DCAwrtPV) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(23);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(63);

    Double_t DCAbtwKaons = TMath::Abs(kfPosKaonA.GetDistanceFromParticle(kfPosKaonB));
    if (kMax_PosKaonPair_DCAbtwKaons && DCAbtwKaons > kMax_PosKaonPair_DCAbtwKaons) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(24);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(64);

    Double_t DCApkaSV = TMath::Abs(kfPosKaonA.GetDistanceFromVertex(kfPosKaonPair));
    if (kMax_PosKaonPair_DCApkaSV && DCApkaSV > kMax_PosKaonPair_DCApkaSV) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(25);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(65);

    Double_t DCApkbSV = TMath::Abs(kfPosKaonB.GetDistanceFromVertex(kfPosKaonPair));
    if (kMax_PosKaonPair_DCApkbSV && DCApkbSV > kMax_PosKaonPair_DCApkbSV) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(26);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(66);

    Double_t Chi2ndf = (Double_t)kfPosKaonPair.GetChi2() / (Double_t)kfPosKaonPair.GetNDF();
    if (kMax_PosKaonPair_Chi2ndf && Chi2ndf > kMax_PosKaonPair_Chi2ndf) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(27);
    if (isSignal) fHist_PosKaonPairs_Bookkeep->Fill(67);

    return kTRUE;
}

/*                            */
/**  Mathematical Functions  **/
/*** ====================== ***/

/*
 * Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 * - Input: `lvParticle`, `X`, `Y`, `Z`, `refPointX`, `refPointY`, `refPointZ`
 * - Return: `cos(theta)`
 * (Based on `AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()`)
 */
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX,
                                                       Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    return TMath::Cos(lvParticle.Angle(posRelativeToRef));
}

/*
 * Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's TLorentzVector.
 */
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX,
                                                       Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    TVector3 momParticle(Px, Py, Pz);
    return TMath::Cos(momParticle.Angle(posRelativeToRef));
}

/*
 * This gives the Armenteros-Podolanski alpha.
 * (Based on `AliRoot/STEER/ESD/AliESDv0::AlphaV0()`)
 */
Double_t AliAnalysisTaskSexaquark::ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz,
                                                   Double_t Pos_Px, Double_t Pos_Py, Double_t Pos_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);
    TVector3 momPos(Pos_Px, Pos_Py, Pos_Pz);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    // (protection)
    if (lQlPos + lQlNeg == 0.) {
        return 2;
    }  // closure
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/*
 * This gives the Armenteros-Podolanski qT
 * (Based on `AliRoot/STEER/ESD/AliESDv0::PtArmV0()`)
 */
Double_t AliAnalysisTaskSexaquark::ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);

    return momNeg.Perp(momTot);
}

/*
 * Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 * This function stores the point of closest approach, and returns its distance to the PV.
 * (PENDING: are we sure of that?)
 * - Input: `V0_Px`, `V0_Py`, `V0_Pz`, `V0_X`, `V0_Y`, `V0_Z`, `PV_X`, `PV_Y`, `PV_Z`
 * - Return: the distance of closest approach to the PV
 */
Double_t AliAnalysisTaskSexaquark::LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z,
                                                Double_t refPointX, Double_t refPointY, Double_t refPointZ) {

    TVector3 V0Momentum(V0_Px, V0_Py, V0_Pz);
    TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
    TVector3 RefVertex(refPointX, refPointY, refPointZ);

    TVector3 CrossProduct = (RefVertex - V0Vertex).Cross(V0Momentum);

    return CrossProduct.Mag() / V0Momentum.Mag();

    /*
    // convert input vectors to arrays
    Double_t ref[3] = {refPointX, refPointY, refPointZ};
    Double_t pos0[3] = {V0_X, V0_Y, V0_Z};
    Double_t dir0[3] = {V0_Px, V0_Py, V0_Pz};

    // lambda function
    auto func = [this, &ref, &pos0, &dir0](const Double_t* t) { return SquaredDistancePointToLine(t, ref, pos0, dir0); };
    ROOT::Math::Functor f(func, 1);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetVariable(0, "t", 0., 0.01);
    min->Minimize();

    // print results
    const Double_t* xs = min->X();

    return TMath::Sqrt(SquaredDistancePointToLine(xs, ref, pos0, dir0));
    */
}

/*
 * Calculate transverse impact parameter w.r.t. Primary Vertex.
 */
Float_t AliAnalysisTaskSexaquark::MCRec_GetImpactParameter(AliESDtrack* track) {
    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

    return (Float_t)track->GetD(PV[0], PV[1], fESD->GetMagneticField());
}

/*
 * Find a common vertex for two neutral particles, given their position and momentum, assuming they propagate in straight lines.
 * If a common vertex cannot be found, then return the DCA between the two lines.
 * (Based on Marty Cohen's answer in https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d)
 */
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v1(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0,
                                                            TVector3& PCA1) {

    if (v3_dir0.Mag() == 0. || v3_dir1.Mag() == 0.) {
        AliWarning("One of the directions is null!");
        return -1.;
    }

    /* Require perpendicularity, that means both lines must intersect eventually */

    if (v3_dir0.Cross(v3_dir1).Mag() <= 0.) {
        AliWarning("TwoLinesDCA :: Not possible to intersect!");
        return -1.;
    }

    TVector3 diff = v3_pos0 - v3_pos1;
    Double_t determinant = v3_dir0.Dot(v3_dir1) * v3_dir0.Dot(v3_dir1) - v3_dir0.Mag2() * v3_dir1.Mag2();
    if (determinant == 0.) {
        AliWarning("TwoLinesDCA :: Lines are parallel!");
        return -1.;
    }
    Double_t t0 = (v3_dir1.Mag2() * v3_dir0.Dot(diff) - v3_dir1.Dot(diff) * v3_dir1.Dot(v3_dir0)) / determinant;
    Double_t t1 = (-1 * v3_dir0.Mag2() * v3_dir1.Dot(diff) + v3_dir0.Dot(diff) * v3_dir1.Dot(v3_dir0)) / determinant;

    TString aux = "";
    if (t0 >= 0. || t1 >= 0.) {
        aux = "*";
    }

    Double_t dca = (diff + v3_dir0 * t0 - v3_dir1 * t1).Mag();
    PCA0 = v3_pos0 + t0 * v3_dir0;
    PCA1 = v3_pos1 + t1 * v3_dir1;
    TVector3 middle_point((PCA0.X() + PCA1.X()) * 0.5, (PCA0.Y() + PCA1.Y()) * 0.5, (PCA0.Z() + PCA1.Z()) * 0.5);

    AliInfoF("dca pca1 pca2 radius = %f (%f, %f, %f) (%f, %f, %f) %f %s", dca, PCA0.X(), PCA0.Y(), PCA0.Z(), PCA1.X(), PCA1.Y(), PCA1.Z(),
             middle_point.Perp(), aux.Data());

    return dca;
}

/*
 * Given the parametric form of a line, find the distance between it and a point.
 * - Variable: `t`, the parameter of the line.
 * - Parameters: `point`, the point to which we want to calculate the distance.
 * - Parameters: `linePoint`, `lineDir`, the position and direction of the line.
 * - Return: the square of the distance between the point and the line.
 */
Double_t AliAnalysisTaskSexaquark::SquaredDistancePointToLine(const Double_t* t, Double_t point[], Double_t linePoint[], Double_t lineDir[]) {
    return TMath::Power(point[0] - linePoint[0] - t[0] * lineDir[0], 2) +  //
           TMath::Power(point[1] - linePoint[1] - t[0] * lineDir[1], 2) +  //
           TMath::Power(point[2] - linePoint[2] - t[0] * lineDir[2], 2);
}

/*
 * Given the parametric form of two lines, find their distance.
 * - Variable: `t`, array of size 2, where `t[0]` and `t[1] are the parameters for the first and second line, respectively.
 * - Parameters: `pos0`, `dir0`, `pos1`, `dir1`, the position and direction of the two lines.
 * - Return: the square of the distance between the two lines.
 */
Double_t AliAnalysisTaskSexaquark::SquaredDistanceBetweenLines(const Double_t* t, Double_t pos0[], Double_t dir0[], Double_t pos1[],
                                                               Double_t dir1[]) {
    return TMath::Power(dir1[0] * t[1] + pos1[0] - (dir0[0] * t[0] + pos0[0]), 2) +
           TMath::Power(dir1[1] * t[1] + pos1[1] - (dir0[1] * t[0] + pos0[1]), 2) +
           TMath::Power(dir1[2] * t[1] + pos1[2] - (dir0[2] * t[0] + pos0[2]), 2);
}

/*
 * Calculate the distance of closest approach between two lines, given their position and momentum, assuming they propagate in straight lines.
 * - Input: `pos0`, `dir0`, `pos1`, `dir1`, the position and direction of the two lines.
 * - Output: `CPA0`, `CPA1`, the closest point of approach for the two lines.
 * - Return: the distance of closest approach between the two lines.
 */
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v2(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0,
                                                            TVector3& PCA1) {

    // convert input vectors to arrays
    Double_t pos0[3] = {v3_pos0.X(), v3_pos0.Y(), v3_pos0.Z()};
    Double_t dir0[3] = {v3_dir0.X(), v3_dir0.Y(), v3_dir0.Z()};

    Double_t pos1[3] = {v3_pos1.X(), v3_pos1.Y(), v3_pos1.Z()};
    Double_t dir1[3] = {v3_dir1.X(), v3_dir1.Y(), v3_dir1.Z()};

    // lambda function
    auto func = [this, &pos0, &dir0, &pos1, &dir1](const Double_t* t) { return SquaredDistanceBetweenLines(t, pos0, dir0, pos1, dir1); };
    ROOT::Math::Functor f(func, 2);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetLimitedVariable(0, "t0", 0.0, 0.1, -999., 0.);
    min->SetLimitedVariable(1, "t1", 0.0, 0.1, -999., 0.);
    min->Minimize();

    // print results
    const Double_t* xs = min->X();
    Double_t t0 = xs[0];
    Double_t t1 = xs[1];

    Double_t dca = TMath::Sqrt(SquaredDistanceBetweenLines(xs, pos0, dir0, pos1, dir1));
    PCA0.SetXYZ(dir0[0] * t0 + pos0[0], dir0[1] * t0 + pos0[1], dir0[2] * t0 + pos0[2]);
    PCA1.SetXYZ(dir1[0] * t1 + pos1[0], dir1[1] * t1 + pos1[1], dir1[2] * t1 + pos1[2]);

    return dca;
}

/*                                */
/**  Custom V0 Finder Utilities  **/
/*** ========================== ***/

/*
 * This function pre-optimizes a two-track pair in the XY plane
 * and provides two X values for the tracks if successful
 * [copied from AliRoot/STEER/ESD/AliV0vertexer.cxx]
 */
Bool_t AliAnalysisTaskSexaquark::Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn,
                                             Double_t* lPreprocessxp, const Double_t b) {
    Double_t lMinimumX = -3.0;
    Double_t lMaximumX = 300.0;

    Double_t nhelix[6], phelix[6];
    nt->GetHelixParameters(nhelix, b);
    pt->GetHelixParameters(phelix, b);
    Double_t lNegCenterR[2], lPosCenterR[2];

    // Negative track parameters in XY
    GetHelixCenter(nt, lNegCenterR, b);
    Double_t xNegCenter = lNegCenterR[0];
    Double_t yNegCenter = lNegCenterR[1];
    Double_t NegRadius = TMath::Abs(1. / nhelix[4]);

    // Positive track parameters in XY
    GetHelixCenter(pt, lPosCenterR, b);
    Double_t xPosCenter = lPosCenterR[0];
    Double_t yPosCenter = lPosCenterR[1];
    Double_t PosRadius = TMath::Abs(1. / phelix[4]);

    // Define convenient coordinate system
    // Logical zero: position of negative center
    Double_t ux = xPosCenter - xNegCenter;
    Double_t uy = yPosCenter - yNegCenter;

    // Check center-to-center distance
    Double_t lDist = TMath::Sqrt(TMath::Power(xNegCenter - xPosCenter, 2) + TMath::Power(yNegCenter - yPosCenter, 2));
    // Normalize ux, uz to unit vector
    ux /= lDist;
    uy /= lDist;

    // Calculate perpendicular vector (normalized)
    Double_t vx = -uy;
    Double_t vy = +ux;

    Double_t lPreprocessDCAxy = 1e+3;  // define outside scope
    *lPreprocessxp = pt->GetX();       // start at current location
    *lPreprocessxn = nt->GetX();       // start at current location

    //============================================================
    // Pre-optimization in the XY plane: cases considered here
    //============================================================
    //
    //  Case 1: Circles do not touch, centers far away
    //          (D > R1 + R2)
    //
    //  Case 2: Circles touch, centers at reasonable distance wrt D
    //          (D < R1 + R2) && (D > |R1-R2|)
    //
    //  Case 3: Circles do not touch, one inside the other
    //          (D < |R1-R2|)
    //
    //  Cases 1 and 2 are treated. Case 3 is not treated (unlikely
    //  to be a problem with unlike-sign charged tracks): brute
    //  force minimization takes place in any case
    //
    //============================================================

    //______________________
    // CASE 1
    if ((lDist > NegRadius + PosRadius)) {
        //================================================================
        // Case 1: distance bigger than sum of radii ("gamma-like")
        //        re-position tracks along the center-to-center axis
        // Re-position negative track
        Double_t xNegOptPosition = xNegCenter + NegRadius * ux;
        Double_t yNegOptPosition = yNegCenter + NegRadius * uy;
        Double_t csNeg = TMath::Cos(nt->GetAlpha());
        Double_t snNeg = TMath::Sin(nt->GetAlpha());
        Double_t xThisNeg = xNegOptPosition * csNeg + yNegOptPosition * snNeg;

        // Re-position positive track
        Double_t xPosOptPosition = xPosCenter - PosRadius * ux;
        Double_t yPosOptPosition = yPosCenter - PosRadius * uy;
        Double_t csPos = TMath::Cos(pt->GetAlpha());
        Double_t snPos = TMath::Sin(pt->GetAlpha());
        Double_t xThisPos = xPosOptPosition * csPos + yPosOptPosition * snPos;

        if (xThisNeg < lMaximumX && xThisPos < lMaximumX && xThisNeg > lMinimumX && xThisPos > lMinimumX) {
            Double_t lCase1NegR[3] = {0.}, lCase1PosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg, b, lCase1NegR) && pt->GetXYZAt(xThisPos, b, lCase1PosR)) {
                lPreprocessDCAxy = TMath::Sqrt(TMath::Power(lCase1NegR[0] - lCase1PosR[0], 2) + TMath::Power(lCase1NegR[1] - lCase1PosR[1], 2) +
                                               TMath::Power(lCase1NegR[2] - lCase1PosR[2], 2));
                // Pass coordinates
                if (lPreprocessDCAxy < 999) {
                    *lPreprocessxp = xThisPos;
                    *lPreprocessxn = xThisNeg;
                }
            }
        }
        //================================================================
    }

    //______________________
    // CASE 2
    if ((lDist > TMath::Abs(NegRadius - PosRadius)) && (lDist < NegRadius + PosRadius)) {
        //================================================================
        // Case 2: distance smaller than sum of radii (cowboy/sailor configs)

        // Calculate coordinate for radical line
        Double_t lRadical = (lDist * lDist - PosRadius * PosRadius + NegRadius * NegRadius) / (2 * lDist);

        // Calculate absolute displacement from center-to-center axis
        Double_t lDisplace = (0.5 / lDist) * TMath::Sqrt((-lDist + PosRadius - NegRadius) * (-lDist - PosRadius + NegRadius) *
                                                         (-lDist + PosRadius + NegRadius) * (lDist + PosRadius + NegRadius));

        // 3D distances in the two cases studied (prefer smallest)
        Double_t lCase2aDCA = 1e+3;
        Double_t lCase2bDCA = 1e+3;

        // 2 cases: positive and negative displacement
        Double_t xNegOptPosition[2], yNegOptPosition[2], xPosOptPosition[2], yPosOptPosition[2];
        Double_t csNeg, snNeg, csPos, snPos;
        Double_t xThisNeg[2], xThisPos[2];

        csNeg = TMath::Cos(nt->GetAlpha());
        snNeg = TMath::Sin(nt->GetAlpha());
        csPos = TMath::Cos(pt->GetAlpha());
        snPos = TMath::Sin(pt->GetAlpha());

        // Case 2a: Positive displacement along v vector
        // Re-position negative track
        xNegOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
        yNegOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
        xThisNeg[0] = xNegOptPosition[0] * csNeg + yNegOptPosition[0] * snNeg;
        // Re-position positive track
        xPosOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
        yPosOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
        xThisPos[0] = xPosOptPosition[0] * csPos + yPosOptPosition[0] * snPos;

        // Case 2b: Negative displacement along v vector
        // Re-position negative track
        xNegOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
        yNegOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
        xThisNeg[1] = xNegOptPosition[1] * csNeg + yNegOptPosition[1] * snNeg;
        // Re-position positive track
        xPosOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
        yPosOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
        xThisPos[1] = xPosOptPosition[1] * csPos + yPosOptPosition[1] * snPos;

        // Test the two cases, please

        // Case 2a
        if (xThisNeg[0] < lMaximumX && xThisPos[0] < lMaximumX && xThisNeg[0] > lMinimumX && xThisPos[0] > lMinimumX) {
            Double_t lCase2aNegR[3] = {0.}, lCase2aPosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg[0], b, lCase2aNegR) && pt->GetXYZAt(xThisPos[0], b, lCase2aPosR)) {
                lCase2aDCA = TMath::Sqrt(TMath::Power(lCase2aNegR[0] - lCase2aPosR[0], 2) + TMath::Power(lCase2aNegR[1] - lCase2aPosR[1], 2) +
                                         TMath::Power(lCase2aNegR[2] - lCase2aPosR[2], 2));
            }
        }

        // Case 2b
        if (xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX) {
            Double_t lCase2bNegR[3] = {0.}, lCase2bPosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg[1], b, lCase2bNegR) && pt->GetXYZAt(xThisPos[1], b, lCase2bPosR)) {
                lCase2bDCA = TMath::Sqrt(TMath::Power(lCase2bNegR[0] - lCase2bPosR[0], 2) + TMath::Power(lCase2bNegR[1] - lCase2bPosR[1], 2) +
                                         TMath::Power(lCase2bNegR[2] - lCase2bPosR[2], 2));
            }
        }
        // Minor detail: all things being equal, prefer closest X
        Double_t lCase2aSumX = xThisPos[0] + xThisNeg[0];
        Double_t lCase2bSumX = xThisPos[1] + xThisNeg[1];

        Double_t lDCAxySmallestR = lCase2aDCA;
        Double_t lxpSmallestR = xThisPos[0];
        Double_t lxnSmallestR = xThisNeg[0];

        Double_t lDCAxyLargestR = lCase2bDCA;
        Double_t lxpLargestR = xThisPos[1];
        Double_t lxnLargestR = xThisNeg[1];

        if (lCase2bSumX + 1e-6 < lCase2aSumX) {
            lDCAxySmallestR = lCase2bDCA;
            lxpSmallestR = xThisPos[1];
            lxnSmallestR = xThisNeg[1];
            lDCAxyLargestR = lCase2aDCA;
            lxpLargestR = xThisPos[0];
            lxnLargestR = xThisNeg[0];
        }

        // Pass conclusion to lPreprocess variables, please
        lPreprocessDCAxy = lDCAxySmallestR;
        *lPreprocessxp = lxpSmallestR;
        *lPreprocessxn = lxnSmallestR;
        if (lDCAxyLargestR + 1e-6 < lDCAxySmallestR) {  // beware epsilon: numerical calculations are unstable here
            lPreprocessDCAxy = lDCAxyLargestR;
            *lPreprocessxp = lxpLargestR;
            *lPreprocessxn = lxnLargestR;
        }
        // Protection against something too crazy, please
        if (lPreprocessDCAxy > 999) {
            *lPreprocessxp = pt->GetX();  // start at current location
            *lPreprocessxn = nt->GetX();  // start at current location
        }
    }
    // End of preprocessing stage!
    // at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
    Bool_t lWorked = kFALSE;
    if (lPreprocessDCAxy < 999) {
        lWorked = kTRUE;
    }
    return lWorked;
}

/*
 * Copied from AliV0ReaderV1::GetHelixCenter
 * Get Center of the helix track parametrization
 * (based on AliRoot/STEER/ESD/AliV0vertexer.cxx)
 */
void AliAnalysisTaskSexaquark::GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b) {

    Int_t charge = track->Charge();

    Double_t helix[6];
    track->GetHelixParameters(helix, b);

    Double_t xpos = helix[5];
    Double_t ypos = helix[0];
    Double_t radius = TMath::Abs(1. / helix[4]);
    Double_t phi = helix[2];
    if (phi < 0) {
        phi = phi + 2 * TMath::Pi();
    }

    phi -= TMath::Pi() / 2.;
    Double_t xpoint = radius * TMath::Cos(phi);
    Double_t ypoint = radius * TMath::Sin(phi);

    if (b < 0 && charge > 0) {
        xpoint = -xpoint;
        ypoint = -ypoint;
    }

    if (b > 0 && charge < 0) {
        xpoint = -xpoint;
        ypoint = -ypoint;
    }

    center[0] = xpos + xpoint;
    center[1] = ypos + ypoint;
    return;
}

/*
 * Calculate position of a point on a track and some derivatives
 * (Copied from `AliRoot/STEER/STEERBase/AliExternalTrackParam.cxx`)
 */
void AliAnalysisTaskSexaquark::Evaluate(const Double_t* h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]) {
    Double_t phase = h[4] * t + h[2];
    Double_t sn = TMath::Sin(phase), cs = TMath::Cos(phase);

    r[0] = h[5];
    r[1] = h[0];
    if (TMath::Abs(h[4]) > kAlmost0) {
        r[0] += (sn - h[6]) / h[4];
        r[1] -= (cs - h[7]) / h[4];
    } else {
        r[0] += t * cs;
        r[1] -= -t * sn;
    }
    r[2] = h[1] + h[3] * t;

    g[0] = cs;
    g[1] = sn;
    g[2] = h[3];

    gg[0] = -h[4] * sn;
    gg[1] = h[4] * cs;
    gg[2] = 0.;
}

/*                             */
/**  Kalman Filter Utilities  **/
/*** ======================= ***/

/*
 * Correct initialization of a KFParticle.
 * (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFParticle AliAnalysisTaskSexaquark::CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge) {

    Double_t fP[6];
    track.GetXYZ(fP);
    track.PxPyPz(fP + 3);

    Int_t fQ = track.Charge() * TMath::Abs(charge);
    fP[3] *= TMath::Abs(charge);
    fP[4] *= TMath::Abs(charge);
    fP[5] *= TMath::Abs(charge);

    Double_t pt = 1. / TMath::Abs(track.GetParameter()[4]) * TMath::Abs(charge);
    Double_t cs = TMath::Cos(track.GetAlpha());
    Double_t sn = TMath::Sin(track.GetAlpha());
    Double_t r = TMath::Sqrt((1. - track.GetParameter()[2]) * (1. + track.GetParameter()[2]));

    Double_t m00 = -sn;
    Double_t m10 = cs;
    Double_t m23 = -pt * (sn + track.GetParameter()[2] * cs / r);
    Double_t m43 = -pt * pt * (r * cs - track.GetParameter()[2] * sn);
    Double_t m24 = pt * (cs - track.GetParameter()[2] * sn / r);
    Double_t m44 = -pt * pt * (r * sn + track.GetParameter()[2] * cs);
    Double_t m35 = pt;
    Double_t m45 = -pt * pt * track.GetParameter()[3];

    m43 *= track.GetSign();
    m44 *= track.GetSign();
    m45 *= track.GetSign();

    const Double_t* cTr = track.GetCovariance();
    Double_t fC[21];
    fC[0] = cTr[0] * m00 * m00;
    fC[1] = cTr[0] * m00 * m10;
    fC[2] = cTr[0] * m10 * m10;
    fC[3] = cTr[1] * m00;
    fC[4] = cTr[1] * m10;
    fC[5] = cTr[2];
    fC[6] = m00 * (cTr[3] * m23 + cTr[10] * m43);
    fC[7] = m10 * (cTr[3] * m23 + cTr[10] * m43);
    fC[8] = cTr[4] * m23 + cTr[11] * m43;
    fC[9] = m23 * (cTr[5] * m23 + cTr[12] * m43) + m43 * (cTr[12] * m23 + cTr[14] * m43);
    fC[10] = m00 * (cTr[3] * m24 + cTr[10] * m44);
    fC[11] = m10 * (cTr[3] * m24 + cTr[10] * m44);
    fC[12] = cTr[4] * m24 + cTr[11] * m44;
    fC[13] = m23 * (cTr[5] * m24 + cTr[12] * m44) + m43 * (cTr[12] * m24 + cTr[14] * m44);
    fC[14] = m24 * (cTr[5] * m24 + cTr[12] * m44) + m44 * (cTr[12] * m24 + cTr[14] * m44);
    fC[15] = m00 * (cTr[6] * m35 + cTr[10] * m45);
    fC[16] = m10 * (cTr[6] * m35 + cTr[10] * m45);
    fC[17] = cTr[7] * m35 + cTr[11] * m45;
    fC[18] = m23 * (cTr[8] * m35 + cTr[12] * m45) + m43 * (cTr[13] * m35 + cTr[14] * m45);
    fC[19] = m24 * (cTr[8] * m35 + cTr[12] * m45) + m44 * (cTr[13] * m35 + cTr[14] * m45);
    fC[20] = m35 * (cTr[9] * m35 + cTr[13] * m45) + m45 * (cTr[13] * m35 + cTr[14] * m45);

    KFParticle part;
    part.Create(fP, fC, fQ, mass);

    return part;
}

/*
 * Correct initialization of a KFVertex.
 * (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFVertex AliAnalysisTaskSexaquark::CreateKFVertex(const AliVVertex& vertex) {

    Double_t param[6];
    vertex.GetXYZ(param);

    Double_t cov[6];
    vertex.GetCovarianceMatrix(cov);

    KFPVertex kfpVtx;
    Float_t paramF[3] = {(Float_t)param[0], (Float_t)param[1], (Float_t)param[2]};
    kfpVtx.SetXYZ(paramF);
    Float_t covF[6] = {(Float_t)cov[0], (Float_t)cov[1], (Float_t)cov[2], (Float_t)cov[3], (Float_t)cov[4], (Float_t)cov[5]};
    kfpVtx.SetCovarianceMatrix(covF);

    KFVertex KFVtx(kfpVtx);

    return KFVtx;
}

/*
 * Transport a KFParticle to the point of closest approach w.r.t. another KFParticle.
 * - Uses: `fPDG`
 * - Input: `kfThis`, `kfOther`, `pdgThis`, `chargeThis`
 * - Return: `kfTransported`
 */
KFParticle AliAnalysisTaskSexaquark::TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis) {

    float dS[2];
    float dsdr[4][6];
    // kfThis.GetDStoParticle(kfOther, dS, dsdr);
    GetDStoParticleBz(fMagneticField, kfThis, kfOther, dS, dsdr);

    float mP[8], mC[36];
    kfThis.Transport(dS[0], dsdr[0], mP, mC);

    float mM = fPDG.GetParticle(pdgThis)->Mass();
    float mQ = chargeThis;  // only valid for charged particles with Q = +/- 1

    KFParticle kfTransported;
    kfTransported.Create(mP, mC, mQ, mM);

    return kfTransported;
}

/*
TEST.
*/
// clang-format off
void AliAnalysisTaskSexaquark::GetDStoParticleBz( float Bz, const KFParticleBase &p, const KFParticleBase &q, float dS[2], float dsdr[4][6], const float* param1, const float* param2 )  const
{
  /** Calculates dS = l/p parameters for two particles, where \n
   ** 1) l - signed distance to the DCA point with the other particle;\n
   ** 2) p - momentum of the particle; \n
   ** under the assumption of the constant homogeneous field Bz. dS[0] is the transport parameter for the current particle,
   ** dS[1] - for the particle "p".
   ** Also calculates partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:\n
   ** 1) dsdr[0][6] = d(dS[0])/d(param1);\n
   ** 2) dsdr[1][6] = d(dS[0])/d(param2);\n
   ** 3) dsdr[2][6] = d(dS[1])/d(param1);\n
   ** 4) dsdr[3][6] = d(dS[1])/d(param2);\n
   ** where param1 are parameters of the current particle (if the pointer is not provided it is initialised with fP) and
   ** param2 are parameters of the second particle "q" (if the pointer is not provided it is initialised with q.fP). Parameters
   ** param1 and param2 should be either provided both or both set to null pointers.
   ** \param[in] Bz - magnetic field Bz
   ** \param[in] p - second particle
   ** \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
   ** \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
   ** \param[in] param1 - optional parameter, is used in case if the parameters of the current particles are rotated
   ** to other coordinate system (see GetDStoParticleBy() function), otherwise p.fP are used
   ** \param[in] param2 - optional parameter, is used in case if the parameters of the second particles are rotated
   ** to other coordinate system (see GetDStoParticleBy() function), otherwise q.fP are used
   **/

  float p_fP[8];           ///< Particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}.
  float q_fP[8];           ///< Particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}.
    for (Int_t i = 0; i < 8; i++) {
        p_fP[i] = p.GetParameter(i);
        q_fP[i] = q.GetParameter(i);
    }

  if(!param1)
  {
    param1 = p_fP;
    param2 = q_fP;
  }

  //* Get dS to another particle for Bz field
  const float kOvSqr6 = 1.f/sqrt(float(6.f));
  const float kCLight = 0.000299792458f;

  //in XY plane
  //first root
  const float& bq1 = Bz*p.GetQ()*kCLight;
  const float& bq2 = Bz*q.GetQ()*kCLight;

  const bool& isStraight1 = fabs(bq1) < 1.e-8f;
  const bool& isStraight2 = fabs(bq2) < 1.e-8f;

  if( isStraight1 && isStraight2 )
  {
    p.GetDStoParticleLine(q, dS, dsdr);
    return;
  }

  const float& px1 = param1[3];
  const float& py1 = param1[4];
  const float& pz1 = param1[5];

  const float& px2 = param2[3];
  const float& py2 = param2[4];
  const float& pz2 = param2[5];

  const float& pt12 = px1*px1 + py1*py1;
  const float& pt22 = px2*px2 + py2*py2;

  const float& x01 = param1[0];
  const float& y01 = param1[1];
  const float& z01 = param1[2];

  const float& x02 = param2[0];
  const float& y02 = param2[1];
  const float& z02 = param2[2];

  float dS1[2] = {0.f}, dS2[2]={0.f};

  const float& dx0 = (x01 - x02);
  const float& dy0 = (y01 - y02);
  const float& dr02 = dx0*dx0 + dy0*dy0;
  const float& drp1  = dx0*px1 + dy0*py1;
  const float& dxyp1 = dx0*py1 - dy0*px1;
  const float& drp2  = dx0*px2 + dy0*py2;
  const float& dxyp2 = dx0*py2 - dy0*px2;
  const float& p1p2 = px1*px2 + py1*py2;
  const float& dp1p2 = px1*py2 - px2*py1;

  const float& k11 = (bq2*drp1 - dp1p2);
  const float& k21 = (bq1*(bq2*dxyp1 - p1p2) + bq2*pt12);
  const float& k12 = ((bq1*drp2 - dp1p2));
  const float& k22 = (bq2*(bq1*dxyp2 + p1p2) - bq1*pt22);

  const float& kp = (dxyp1*bq2 - dxyp2*bq1 - p1p2);
  const float& kd = dr02/2.f*bq1*bq2 + kp;
  const float& c1 = -(bq1*kd + pt12*bq2);
  const float& c2 = bq2*kd + pt22*bq1;

  float d1 = pt12*pt22 - kd*kd;
  if(d1<0.f)
    d1 = float(0.f);
  d1 = sqrt( d1 );
  float d2 = pt12*pt22 - kd*kd;
  if(d2<0.f)
    d2 = float(0.f);
  d2 = sqrt( d2 );

  // find two points of closest approach in XY plane

  float dS1dR1[2][6];
  float dS2dR2[2][6];

  float dS1dR2[2][6];
  float dS2dR1[2][6];

  float dk11dr1[6] = {bq2*px1, bq2*py1, 0, bq2*dx0 - py2, bq2*dy0 + px2, 0};
  float dk11dr2[6] = {-bq2*px1, -bq2*py1, 0, py1, -px1, 0};
  float dk12dr1[6] = {bq1*px2, bq1*py2, 0, -py2, px2, 0};
  float dk12dr2[6] = {-bq1*px2, -bq1*py2, 0, bq1*dx0 + py1, bq1*dy0 - px1, 0};
  float dk21dr1[6] = {bq1*bq2*py1, -bq1*bq2*px1, 0, 2*bq2*px1 + bq1*(-(bq2*dy0) - px2), 2*bq2*py1 + bq1*(bq2*dx0 - py2), 0};
  float dk21dr2[6] = {-(bq1*bq2*py1), bq1*bq2*px1, 0, -(bq1*px1), -(bq1*py1), 0};
  float dk22dr1[6] = {bq1*bq2*py2, -(bq1*bq2*px2), 0, bq2*px2, bq2*py2, 0};
  float dk22dr2[6] = {-(bq1*bq2*py2), bq1*bq2*px2, 0, bq2*(-(bq1*dy0) + px1) - 2*bq1*px2, bq2*(bq1*dx0 + py1) - 2*bq1*py2, 0};

  float dkddr1[6] = {bq1*bq2*dx0 + bq2*py1 - bq1*py2, bq1*bq2*dy0 - bq2*px1 + bq1*px2, 0, -bq2*dy0 - px2, bq2*dx0 - py2, 0};
  float dkddr2[6] = {-bq1*bq2*dx0 - bq2*py1 + bq1*py2, -bq1*bq2*dy0 + bq2*px1 - bq1*px2, 0, bq1*dy0 - px1, -bq1*dx0 - py1, 0};

  float dc1dr1[6] = {-(bq1*(bq1*bq2*dx0 + bq2*py1 - bq1*py2)), -(bq1*(bq1*bq2*dy0 - bq2*px1 + bq1*px2)), 0, -2*bq2*px1 - bq1*(-(bq2*dy0) - px2), -2*bq2*py1 - bq1*(bq2*dx0 - py2), 0};
  float dc1dr2[6] = {-(bq1*(-(bq1*bq2*dx0) - bq2*py1 + bq1*py2)), -(bq1*(-(bq1*bq2*dy0) + bq2*px1 - bq1*px2)), 0, -(bq1*(bq1*dy0 - px1)), -(bq1*(-(bq1*dx0) - py1)), 0};

  float dc2dr1[6] = {bq2*(bq1*bq2*dx0 + bq2*py1 - bq1*py2), bq2*(bq1*bq2*dy0 - bq2*px1 + bq1*px2), 0, bq2*(-(bq2*dy0) - px2), bq2*(bq2*dx0 - py2), 0};
  float dc2dr2[6] = {bq2*(-(bq1*bq2*dx0) - bq2*py1 + bq1*py2), bq2*(-(bq1*bq2*dy0) + bq2*px1 - bq1*px2), 0, bq2*(bq1*dy0 - px1) + 2*bq1*px2, bq2*(-(bq1*dx0) - py1) + 2*bq1*py2, 0};

  float dd1dr1[6] = {0,0,0,0,0,0};
  float dd1dr2[6] = {0,0,0,0,0,0};
  if(d1>0)
  {
    for(int i=0; i<6; i++)
    {
      dd1dr1[i] = -kd/d1*dkddr1[i];
      dd1dr2[i] = -kd/d1*dkddr2[i];
    }
    dd1dr1[3] += px1/d1*pt22; dd1dr1[4] += py1/d1*pt22;
    dd1dr2[3] += px2/d1*pt12; dd1dr2[4] += py2/d1*pt12;
  }

  if(!isStraight1)
  {
    dS1[0] = atan2( bq1*(k11*c1 + k21*d1), (bq1*k11*d1*bq1 - k21*c1) )/bq1;
    dS1[1] = atan2( bq1*(k11*c1 - k21*d1), (-bq1*k11*d1*bq1 - k21*c1) )/bq1;

    float a = bq1*(k11*c1 + k21*d1);
    float b = bq1*k11*d1*bq1 - k21*c1;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq1*( dk11dr1[iP]*c1 + k11*dc1dr1[iP] + dk21dr1[iP]*d1 + k21*dd1dr1[iP] );
        const float dadr2 = bq1*( dk11dr2[iP]*c1 + k11*dc1dr2[iP] + dk21dr2[iP]*d1 + k21*dd1dr2[iP] );
        const float dbdr1 = bq1*bq1*( dk11dr1[iP]*d1 + k11*dd1dr1[iP] ) - ( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = bq1*bq1*( dk11dr2[iP]*d1 + k11*dd1dr2[iP] ) - ( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[0][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS1dR2[0][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS1dR1[0][iP] = 0;
        dS1dR2[0][iP] = 0;
      }
    }

    a = bq1*(k11*c1 - k21*d1);
    b = -bq1*k11*d1*bq1 - k21*c1;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq1*( dk11dr1[iP]*c1 + k11*dc1dr1[iP] - (dk21dr1[iP]*d1 + k21*dd1dr1[iP]) );
        const float dadr2 = bq1*( dk11dr2[iP]*c1 + k11*dc1dr2[iP] - (dk21dr2[iP]*d1 + k21*dd1dr2[iP]) );
        const float dbdr1 = -bq1*bq1*( dk11dr1[iP]*d1 + k11*dd1dr1[iP] ) - ( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = -bq1*bq1*( dk11dr2[iP]*d1 + k11*dd1dr2[iP] ) - ( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[1][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS1dR2[1][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS1dR1[1][iP] = 0;
        dS1dR2[1][iP] = 0;
      }
    }
  }
  if(!isStraight2)
  {
    dS2[0] = atan2( (bq2*k12*c2 + k22*d2*bq2), (bq2*k12*d2*bq2 - k22*c2) )/bq2;
    dS2[1] = atan2( (bq2*k12*c2 - k22*d2*bq2), (-bq2*k12*d2*bq2 - k22*c2) )/bq2;

    float a = bq2*(k12*c2 + k22*d2);
    float b = bq2*k12*d2*bq2 - k22*c2;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq2*( dk12dr1[iP]*c2 + k12*dc2dr1[iP] + dk22dr1[iP]*d1 + k22*dd1dr1[iP] );
        const float dadr2 = bq2*( dk12dr2[iP]*c2 + k12*dc2dr2[iP] + dk22dr2[iP]*d1 + k22*dd1dr2[iP] );
        const float dbdr1 = bq2*bq2*( dk12dr1[iP]*d1 + k12*dd1dr1[iP] ) - (dk22dr1[iP]*c2 + k22*dc2dr1[iP]);
        const float dbdr2 = bq2*bq2*( dk12dr2[iP]*d1 + k12*dd1dr2[iP] ) - (dk22dr2[iP]*c2 + k22*dc2dr2[iP]);

        dS2dR1[0][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS2dR2[0][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS2dR1[0][iP] = 0;
        dS2dR2[0][iP] = 0;
      }
    }

    a = bq2*(k12*c2 - k22*d2);
    b = -bq2*k12*d2*bq2 - k22*c2;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq2*( dk12dr1[iP]*c2 + k12*dc2dr1[iP] - (dk22dr1[iP]*d1 + k22*dd1dr1[iP]) );
        const float dadr2 = bq2*( dk12dr2[iP]*c2 + k12*dc2dr2[iP] - (dk22dr2[iP]*d1 + k22*dd1dr2[iP]) );
        const float dbdr1 = -bq2*bq2*( dk12dr1[iP]*d1 + k12*dd1dr1[iP] ) - (dk22dr1[iP]*c2 + k22*dc2dr1[iP]);
        const float dbdr2 = -bq2*bq2*( dk12dr2[iP]*d1 + k12*dd1dr2[iP] ) - (dk22dr2[iP]*c2 + k22*dc2dr2[iP]);

        dS2dR1[1][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS2dR2[1][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS2dR1[1][iP] = 0;
        dS2dR2[1][iP] = 0;
      }
    }
  }
  if(isStraight1 && (pt12>0.f) )
  {
    dS1[0] = (k11*c1 + k21*d1)/(- k21*c1);
    dS1[1] = (k11*c1 - k21*d1)/(- k21*c1);

    float a = k11*c1 + k21*d1;
    float b = -k21*c1;

    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk11dr1[iP]*c1 + k11*dc1dr1[iP] + dk21dr1[iP]*d1 + k21*dd1dr1[iP] );
        const float dadr2 = ( dk11dr2[iP]*c1 + k11*dc1dr2[iP] + dk21dr2[iP]*d1 + k21*dd1dr2[iP] );
        const float dbdr1 = -( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = -( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[0][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS1dR2[0][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS1dR1[0][iP] = 0;
        dS1dR2[0][iP] = 0;
      }
    }

    a = k11*c1 - k21*d1;
    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk11dr1[iP]*c1 + k11*dc1dr1[iP] - dk21dr1[iP]*d1 - k21*dd1dr1[iP] );
        const float dadr2 = ( dk11dr2[iP]*c1 + k11*dc1dr2[iP] - dk21dr2[iP]*d1 - k21*dd1dr2[iP] );
        const float dbdr1 = -( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = -( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[1][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS1dR2[1][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS1dR1[1][iP] = 0;
        dS1dR2[1][iP] = 0;
      }
    }
  }
  if(isStraight2 && (pt22>0.f) )
  {
    dS2[0] = (k12*c2 + k22*d2)/(- k22*c2);
    dS2[1] = (k12*c2 - k22*d2)/(- k22*c2);

    float a = k12*c2 + k22*d1;
    float b = -k22*c2;

    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk12dr1[iP]*c2 + k12*dc2dr1[iP] + dk22dr1[iP]*d1 + k22*dd1dr1[iP] );
        const float dadr2 = ( dk12dr2[iP]*c2 + k12*dc2dr2[iP] + dk22dr2[iP]*d1 + k22*dd1dr2[iP] );
        const float dbdr1 = -( dk22dr1[iP]*c2 + k22*dc2dr1[iP] );
        const float dbdr2 = -( dk22dr2[iP]*c2 + k22*dc2dr2[iP] );

        dS2dR1[0][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS2dR2[0][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS2dR1[0][iP] = 0;
        dS2dR2[0][iP] = 0;
      }
    }

    a = k12*c2 - k22*d1;
    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk12dr1[iP]*c2 + k12*dc2dr1[iP] - dk22dr1[iP]*d1 - k22*dd1dr1[iP] );
        const float dadr2 = ( dk12dr2[iP]*c2 + k12*dc2dr2[iP] - dk22dr2[iP]*d1 - k22*dd1dr2[iP] );
        const float dbdr1 = -( dk22dr1[iP]*c2 + k22*dc2dr1[iP] );
        const float dbdr2 = -( dk22dr2[iP]*c2 + k22*dc2dr2[iP] );

        dS2dR1[1][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS2dR2[1][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS2dR1[1][iP] = 0;
        dS2dR2[1][iP] = 0;
      }
    }
  }

  //select a point which is close to the primary vertex (with the smallest r)

  float dr2[2];
  for(int iP = 0; iP<2; iP++)
  {
    const float& bs1 = bq1*dS1[iP];
    const float& bs2 = bq2*dS2[iP];
    float sss = sin(bs1), ccc = cos(bs1);

    const bool& bs1Big = fabs(bs1) > 1.e-8f;
    const bool& bs2Big = fabs(bs2) > 1.e-8f;

    float sB(0.f), cB(0.f);
    if(bs1Big)
    {
      sB = sss/bq1;
      cB = (1.f-ccc)/bq1;
    }
    else
    {
      sB = ((1.f-bs1*kOvSqr6)*(1.f+bs1*kOvSqr6)*dS1[iP]);
      cB = .5f*sB*bs1;
    }

    const float& x1 = param1[0] + sB*px1 + cB*py1;
    const float& y1 = param1[1] - cB*px1 + sB*py1;
    const float& z1 = param1[2] + dS1[iP]*param1[5];

    sss = sin(bs2), ccc = cos(bs2);

    if(bs2Big)
    {
      sB = sss/bq2;
      cB = (1.f-ccc)/bq2;
    }
    else
    {
      sB = ((1.f-bs2*kOvSqr6)*(1.f+bs2*kOvSqr6)*dS2[iP]);
      cB = .5f*sB*bs2;
    }

    const float& x2 = param2[0] + sB*px2 + cB*py2;
    const float& y2 = param2[1] - cB*px2 + sB*py2;
    const float& z2 = param2[2] + dS2[iP]*param2[5];

    float dx = (x1-x2);
    float dy = (y1-y2);
    float dz = (z1-z2);

    dr2[iP] = dx*dx + dy*dy + dz*dz;
  }

  const bool isFirstRoot = dr2[0] < dr2[1];
  if(isFirstRoot)
  {
    dS[0]  = dS1[0];
    dS[1] = dS2[0];

    for(int iP=0; iP<6; iP++)
    {
      dsdr[0][iP] = dS1dR1[0][iP];
      dsdr[1][iP] = dS1dR2[0][iP];
      dsdr[2][iP] = dS2dR1[0][iP];
      dsdr[3][iP] = dS2dR2[0][iP];
    }
  }
  else
  {
    dS[0]  = dS1[1];
    dS[1] = dS2[1];

    for(int iP=0; iP<6; iP++)
    {
      dsdr[0][iP] = dS1dR1[1][iP];
      dsdr[1][iP] = dS1dR2[1][iP];
      dsdr[2][iP] = dS2dR1[1][iP];
      dsdr[3][iP] = dS2dR2[1][iP];
    }
  }

  //find correct parts of helices
//   int n1(0);
//   int n2(0);
//   float dzMin = fabs( (z01-z02) + dS[0]*pz1 - dS[1]*pz2 );
//   const float pi2(6.283185307f);

  //TODO optimise for loops for neutral particles
//   const float& i1Float = -bq1/pi2*(z01/pz1+dS[0]);
//   for(int di1=-1; di1<=1; di1++)
//   {
//     int i1(0);
//     if(!isStraight1)
//       i1 = int(i1Float) + di1;
//
//     const float& i2Float = ( ((z01-z02) + (dS[0]+pi2*i1/bq1)*pz1)/pz2 - dS[1]) * bq2/pi2;
//     for(int di2 = -1; di2<=1; di2++)
//     {
//       int i2(0);
//       if(!isStraight2)
//         i2 = int(i2Float) + di2;
//
//       const float& z1 = z01 + (dS[0]+pi2*i1/bq1)*pz1;
//       const float& z2 = z02 + (dS[1]+pi2*i2/bq2)*pz2;
//       const float& dz = fabs( z1-z2 );
//
//       if(dz < dzMin)
//       {
//         n1 = i1;
//         n2 = i2;
//         dzMin = dz;
//       }
//     }
//   }
//
//   if(!isStraight1)
//     dS[0] += float(n1)*pi2/bq1;
//   if(!isStraight2)
//     dS[1] += float(n2)*pi2/bq2;

  //Line correction
  {
    const float& bs1 = bq1*dS[0];
    const float& bs2 = bq2*dS[1];
    float sss = sin(bs1), ccc = cos(bs1);

    const bool& bs1Big = fabs(bs1) > 1.e-8f;
    const bool& bs2Big = fabs(bs2) > 1.e-8f;

    float sB(0.f), cB(0.f);
    if(bs1Big)
    {
      sB = sss/bq1;
      cB = (1.f-ccc)/bq1;
    }
    else
    {
      sB = ((1.f-bs1*kOvSqr6)*(1.f+bs1*kOvSqr6)*dS[0]);
      cB = .5f*sB*bs1;
    }

    const float& x1 = x01 + sB*px1 + cB*py1;
    const float& y1 = y01 - cB*px1 + sB*py1;
    const float& z1 = z01 + dS[0]*pz1;
    const float& ppx1 =  ccc*px1 + sss*py1;
    const float& ppy1 = -sss*px1 + ccc*py1;
    const float& ppz1 = pz1;

    float sss1 = sin(bs2), ccc1 = cos(bs2);

    float sB1(0.f), cB1(0.f);
    if(bs2Big)
    {
      sB1 = sss1/bq2;
      cB1 = (1.f-ccc1)/bq2;
    }
    else
    {
      sB1 = ((1.f-bs2*kOvSqr6)*(1.f+bs2*kOvSqr6)*dS[1]);
      cB1 = .5f*sB1*bs2;
    }

    const float& x2 = x02 + sB1*px2 + cB1*py2;
    const float& y2 = y02 - cB1*px2 + sB1*py2;
    const float& z2 = z02 + dS[1]*pz2;
    const float& ppx2 =  ccc1*px2 + sss1*py2;
    const float& ppy2 = -sss1*px2 + ccc1*py2;
    const float& ppz2 = pz2;

    const float& p12  = ppx1*ppx1 + ppy1*ppy1 + ppz1*ppz1;
    const float& p22  = ppx2*ppx2 + ppy2*ppy2 + ppz2*ppz2;
    const float& lp1p2 = ppx1*ppx2 + ppy1*ppy2 + ppz1*ppz2;

    const float& dx = (x2 - x1);
    const float& dy = (y2 - y1);
    const float& dz = (z2 - z1);

    const float& ldrp1 = ppx1*dx + ppy1*dy + ppz1*dz;
    const float& ldrp2 = ppx2*dx + ppy2*dy + ppz2*dz;

    float detp =  lp1p2*lp1p2 - p12*p22;
    if( fabs(detp)<1.e-4 ) detp = 1; //TODO correct!!!

    //dsdr calculation
    const float a1 = ldrp2*lp1p2 - ldrp1*p22;
    const float a2 = ldrp2*p12 - ldrp1*lp1p2;
    const float lp1p2_ds0 = bq1*( ppx2*ppy1 - ppy2*ppx1);
    const float lp1p2_ds1 = bq2*( ppx1*ppy2 - ppy1*ppx2);
    const float ldrp1_ds0 = -p12 + bq1*(ppy1*dx - ppx1*dy);
    const float ldrp1_ds1 =  lp1p2;
    const float ldrp2_ds0 = -lp1p2;
    const float ldrp2_ds1 =  p22 + bq2*(ppy2*dx - ppx2*dy);
    const float detp_ds0 = 2*lp1p2*lp1p2_ds0;
    const float detp_ds1 = 2*lp1p2*lp1p2_ds1;
    const float a1_ds0 = ldrp2_ds0*lp1p2 + ldrp2*lp1p2_ds0 - ldrp1_ds0*p22;
    const float a1_ds1 = ldrp2_ds1*lp1p2 + ldrp2*lp1p2_ds1 - ldrp1_ds1*p22;
    const float a2_ds0 = ldrp2_ds0*p12 - ldrp1_ds0*lp1p2 - ldrp1*lp1p2_ds0;
    const float a2_ds1 = ldrp2_ds1*p12 - ldrp1_ds1*lp1p2 - ldrp1*lp1p2_ds1;

    // AliInfoF("a1_ds0 a1 detp_ds0 detp a1_ds1 detp_ds1 a2_ds0 a2 a2_ds1 = %f %f %f %f %f %f %f %f %f", a1_ds0, a1, detp_ds0, detp, a1_ds1, detp_ds1, a2_ds0, a2, a2_ds1);

    const float dsl1ds0 = a1_ds0/detp - a1*detp_ds0/(detp*detp);
    const float dsl1ds1 = a1_ds1/detp - a1*detp_ds1/(detp*detp);
    const float dsl2ds0 = a2_ds0/detp - a2*detp_ds0/(detp*detp);
    const float dsl2ds1 = a2_ds1/detp - a2*detp_ds1/(detp*detp);

    float dsldr[4][6];
    for(int iP=0; iP<6; iP++)
    {
      dsldr[0][iP] = dsl1ds0*dsdr[0][iP] + dsl1ds1*dsdr[2][iP];
      dsldr[1][iP] = dsl1ds0*dsdr[1][iP] + dsl1ds1*dsdr[3][iP];
      dsldr[2][iP] = dsl2ds0*dsdr[0][iP] + dsl2ds1*dsdr[2][iP];
      dsldr[3][iP] = dsl2ds0*dsdr[1][iP] + dsl2ds1*dsdr[3][iP];
    }

    for(int iDS=0; iDS<4; iDS++)
      for(int iP=0; iP<6; iP++)
        dsdr[iDS][iP] += dsldr[iDS][iP];

    const float lp1p2_dr0[6] = {0, 0, 0, ccc*ppx2 - ppy2*sss, ccc*ppy2 + ppx2*sss, pz2};
    const float lp1p2_dr1[6] = {0, 0, 0, ccc1*ppx1 - ppy1*sss1, ccc1*ppy1 + ppx1*sss1, pz1};
    const float ldrp1_dr0[6] = {-ppx1, -ppy1, -pz1,  cB*ppy1 - ppx1*sB + ccc*dx - sss*dy, -cB*ppx1-ppy1*sB + sss*dx + ccc*dy, -dS[0]*pz1 + dz};
    const float ldrp1_dr1[6] = { ppx1,  ppy1,  pz1, -cB1*ppy1 + ppx1*sB1, cB1*ppx1 + ppy1*sB1, dS[1]*pz1};
    const float ldrp2_dr0[6] = {-ppx2, -ppy2, -pz2, cB*ppy2 - ppx2*sB, -cB*ppx2-ppy2*sB, -dS[0]*pz2};
    const float ldrp2_dr1[6] = {ppx2, ppy2, pz2, -cB1*ppy2 + ppx2*sB1 + ccc1*dx- sss1*dy, cB1*ppx2 + ppy2*sB1 + sss1*dx + ccc1*dy, dz + dS[1]*pz2};
    const float p12_dr0[6] = {0, 0, 0, 2*px1, 2*py1, 2*pz1};
    const float p22_dr1[6] = {0, 0, 0, 2*px2, 2*py2, 2*pz2};
    float a1_dr0[6], a1_dr1[6], a2_dr0[6], a2_dr1[6], detp_dr0[6], detp_dr1[6];
    for(int iP=0; iP<6; iP++)
    {
      a1_dr0[iP] = ldrp2_dr0[iP]*lp1p2 + ldrp2*lp1p2_dr0[iP] - ldrp1_dr0[iP]*p22;
      a1_dr1[iP] = ldrp2_dr1[iP]*lp1p2 + ldrp2*lp1p2_dr1[iP] - ldrp1_dr1[iP]*p22 - ldrp1*p22_dr1[iP];
      a2_dr0[iP] = ldrp2_dr0[iP]*p12 + ldrp2*p12_dr0[iP] - ldrp1_dr0[iP]*lp1p2 - ldrp1*lp1p2_dr0[iP];
      a2_dr1[iP] = ldrp2_dr1[iP]*p12 - ldrp1_dr1[iP]*lp1p2 - ldrp1*lp1p2_dr1[iP];
      detp_dr0[iP] = 2*lp1p2*lp1p2_dr0[iP] - p12_dr0[iP]*p22;
      detp_dr1[iP] = 2*lp1p2*lp1p2_dr1[iP] - p12*p22_dr1[iP];

      dsdr[0][iP] += a1_dr0[iP]/detp - a1*detp_dr0[iP]/(detp*detp);
      dsdr[1][iP] += a1_dr1[iP]/detp - a1*detp_dr1[iP]/(detp*detp);
      dsdr[2][iP] += a2_dr0[iP]/detp - a2*detp_dr0[iP]/(detp*detp);
      dsdr[3][iP] += a2_dr1[iP]/detp - a2*detp_dr1[iP]/(detp*detp);
    }

    dS[0] += (ldrp2*lp1p2 - ldrp1*p22) /detp;
    dS[1] += (ldrp2*p12 - ldrp1*lp1p2)/detp;
  }
}
// clang-format on

/*
 * Clear all containers.
 */
void AliAnalysisTaskSexaquark::ClearContainers() {

    /*  */
    getPdgCode_fromMcIdx.clear();

    isMcIdxSignal.clear();
    isMcIdxSecondary.clear();

    getReactionIdx_fromMcIdx.clear();
    getMcIdx_fromReactionIdx.clear();
    getPt_fromReactionIdx.clear();
    getRadius_fromReactionIdx.clear();

    doesMcIdxHaveMother.clear();
    getMotherMcIdx_fromMcIdx.clear();
    getNegDauMcIdx_fromMcIdx.clear();
    getPosDauMcIdx_fromMcIdx.clear();

    mcIndicesOfTrueV0s.clear();

    /*  */
    getMcIdx_fromEsdIdx.clear();
    esdIdxPassesTrackSelection_AsProton.clear();
    esdIdxPassesTrackSelection_AsAntiProton.clear();
    esdIdxPassesTrackSelection_AsPosKaon.clear();
    esdIdxPassesTrackSelection_AsNegKaon.clear();
    esdIdxPassesTrackSelection_AsPiPlus.clear();
    esdIdxPassesTrackSelection_AsPiMinus.clear();
    getEsdIdx_fromMcIdx.clear();

    isEsdIdxSignal.clear();
    getReactionIdx_fromEsdIdx.clear();
    getEsdIdx_fromReactionIdx.clear();

    doesEsdIdxHaveMother.clear();
    getMotherMcIdx_fromEsdIdx.clear();
    getNegDauEsdIdx_fromMcIdx.clear();
    getPosDauEsdIdx_fromMcIdx.clear();

    esdIndicesOfAntiProtonTracks.clear();
    esdIndicesOfPosKaonTracks.clear();
    esdIndicesOfPiMinusTracks.clear();
    esdIndicesOfPiPlusTracks.clear();

    /*  */
    getEsdIdxOfNegDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfPosDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfNegDau_fromKaonZeroShortIdx.clear();
    getEsdIdxOfPosDau_fromKaonZeroShortIdx.clear();
    getEsdIdxOfNegDau_fromPionPairIdx.clear();
    getEsdIdxOfPosDau_fromPionPairIdx.clear();

    esdAntiLambdas.clear();
    esdKaonsZeroShort.clear();

    kfAntiLambdas.clear();
    kfKaonsZeroShort.clear();
    kfPionPairs.clear();
}

/*                    */
/**  External Files  **/
/*** ============== ***/

/*
 *
 */
void AliAnalysisTaskSexaquark::PrepareInjectedBranches() {
    fTree_Injected->Branch("RunNumber", &tInjected_RunNumber);
    fTree_Injected->Branch("DirNumber", &tInjected_DirNumber);
    fTree_Injected->Branch("EventNumber", &tInjected_EventNumber);
    fTree_Injected->Branch("ReactionID", &tInjected_ReactionID);
    fTree_Injected->Branch("Px", &tInjected_Px);
    fTree_Injected->Branch("Py", &tInjected_Py);
    fTree_Injected->Branch("Pz", &tInjected_Pz);
    fTree_Injected->Branch("Mass", &tInjected_Mass);
    fTree_Injected->Branch("Nucleon_PdgCode", &tInjected_Nucleon_PdgCode);
    fTree_Injected->Branch("Nucleon_Px", &tInjected_Nucleon_Px);
    fTree_Injected->Branch("Nucleon_Py", &tInjected_Nucleon_Py);
    fTree_Injected->Branch("Nucleon_Pz", &tInjected_Nucleon_Pz);
    fTree_Injected->Branch("ReactionChannel", &tInjected_ReactionChannel);
}

/*
 * Open the respective `sim.log` that corresponds to the `RunNumber`+`DirNumber` that's being analyzed.
 * From it, read the injected anti-sexaquark and struck nucleon kinematics and store them into a tree.
 * NOTE: it must be executed ONLY when analyzing SIGNAL SIMs, protection PENDING!
 */
Bool_t AliAnalysisTaskSexaquark::LoadLogsIntoTree() {

    TGrid* alien = nullptr;
    if (!gGrid) {
        alien = TGrid::Connect("alien://");
        if (!alien) return kFALSE;
    }

    TString AliEn_Dir = fAliEnPath(0, fAliEnPath.Last('/'));

    TString orig_path = Form("%s/sim.log", AliEn_Dir.Data());
    AliInfoF("!! Copying file %s ... !!", orig_path.Data());

    // ! assuming path ends with format .../LHC23l1a3/A1.73/297595/001/sim.log ! //
    Int_t AliEn_DirNumber = fDirNumber;
    TObjArray* tokens = fAliEnPath.Tokenize("/");
    Int_t AliEn_RunNumber = ((TObjString*)tokens->At(tokens->GetEntries() - 3))->GetString().Atoi();
    TString AliEn_SimSubSet = ((TObjString*)tokens->At(tokens->GetEntries() - 4))->GetString();

    TString Log_NewBasename = Form("sim_%s_%i_%03i.log", AliEn_SimSubSet.Data(), AliEn_RunNumber, AliEn_DirNumber);

    if (AliEn_Dir.BeginsWith("alien://")) {
        gSystem->Exec(Form("alien.py cp %s file://./%s", orig_path.Data(), Log_NewBasename.Data()));
    } else {
        gSystem->Exec(Form("cp %s ./%s", orig_path.Data(), Log_NewBasename.Data()));
    }

    TString new_path = Form("%s/%s", gSystem->pwd(), Log_NewBasename.Data());
    AliInfoF("!! Reading file %s ... !!", new_path.Data());

    std::ifstream SimLog(new_path);
    if (!SimLog.is_open()) {
        AliInfo("!! Unable to open file !!");
        return kFALSE;
    }

    Int_t CurrentEventNumber = -1;

    // auxiliary variables
    std::string cstr_line;
    TString tstr_line, csv;
    TObjArray* csv_arr = nullptr;

    Char_t aux_char;

    /* Read lines */

    while (std::getline(SimLog, cstr_line)) {

        tstr_line = cstr_line;

        // a new event has appeared
        if (tstr_line.Contains("I-AliGenCocktail::Generate: Generator 1: AliGenHijing")) CurrentEventNumber++;

        if (!tstr_line.Contains("I-AliGenSexaquarkReaction::GenerateN: 6")) continue;

        csv = static_cast<TString>(tstr_line(38, tstr_line.Length() - 1));
        csv_arr = csv.Tokenize(",");

        /* Assign branches and fill tree */

        tInjected_RunNumber = AliEn_RunNumber;
        tInjected_DirNumber = AliEn_DirNumber;  // PENDING: could be another var to simplify things
        tInjected_EventNumber = CurrentEventNumber;
        tInjected_ReactionID = dynamic_cast<TObjString*>(csv_arr->At(0))->String().Atoi();
        tInjected_Px = dynamic_cast<TObjString*>(csv_arr->At(1))->String().Atof();
        tInjected_Py = dynamic_cast<TObjString*>(csv_arr->At(2))->String().Atof();
        tInjected_Pz = dynamic_cast<TObjString*>(csv_arr->At(3))->String().Atof();
        tInjected_Mass = TString(AliEn_SimSubSet(1, 4)).Atof();
        tInjected_Nucleon_PdgCode = AliEn_SimSubSet[0] == 'A' ? 2112 : 2212;  // considering only ADEH
        tInjected_Nucleon_Px = dynamic_cast<TObjString*>(csv_arr->At(4))->String().Atof();
        tInjected_Nucleon_Py = dynamic_cast<TObjString*>(csv_arr->At(5))->String().Atof();
        tInjected_Nucleon_Pz = dynamic_cast<TObjString*>(csv_arr->At(6))->String().Atof();
        tInjected_ReactionChannel = (UInt_t)AliEn_SimSubSet[0];

        fTree_Injected->Fill();
        ClearInjectedBranches();
    }  // end of loop over lines

    AliInfo("!! Closing file ... !!");
    SimLog.close();

    return kTRUE;
}

/*
 *
 */
void AliAnalysisTaskSexaquark::ClearInjectedBranches() {
    tInjected_RunNumber = 0;
    tInjected_DirNumber = 0;
    tInjected_EventNumber = 0;
    tInjected_ReactionID = 0;
    tInjected_Px = 0.;
    tInjected_Py = 0.;
    tInjected_Pz = 0.;
    tInjected_Mass = 0.;
    tInjected_Nucleon_PdgCode = 0;
    tInjected_Nucleon_Px = 0.;
    tInjected_Nucleon_Py = 0.;
    tInjected_Nucleon_Pz = 0.;
    tInjected_ReactionChannel = 0;
}
