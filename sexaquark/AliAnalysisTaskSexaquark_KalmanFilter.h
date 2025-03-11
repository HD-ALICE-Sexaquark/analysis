#ifndef TASKSEXAQUARK_KF_H
#define TASKSEXAQUARK_KF_H

#include "RtypesCore.h"

#include "AliExternalTrackParam.h"
#include "AliVVertex.h"

#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#include "KFParticle.h"
#include "KFVertex.h"

/*
 * Correct initialization of a KFParticle.
 * (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFParticle CreateKFParticle(const AliExternalTrackParam* track_param, Double_t mass, Int_t charge) {

    Double_t fP[6];
    track_param->GetXYZ(fP);
    track_param->PxPyPz(fP + 3);

    Int_t fQ = track_param->Charge() * TMath::Abs(charge);
    fP[3] *= TMath::Abs(charge);
    fP[4] *= TMath::Abs(charge);
    fP[5] *= TMath::Abs(charge);

    Double_t pt = 1. / TMath::Abs(track_param->GetParameter()[4]) * TMath::Abs(charge);
    Double_t cs = TMath::Cos(track_param->GetAlpha());
    Double_t sn = TMath::Sin(track_param->GetAlpha());
    Double_t r = TMath::Sqrt((1. - track_param->GetParameter()[2]) * (1. + track_param->GetParameter()[2]));

    Double_t m00 = -sn;
    Double_t m10 = cs;
    Double_t m23 = -pt * (sn + track_param->GetParameter()[2] * cs / r);
    Double_t m43 = -pt * pt * (r * cs - track_param->GetParameter()[2] * sn);
    Double_t m24 = pt * (cs - track_param->GetParameter()[2] * sn / r);
    Double_t m44 = -pt * pt * (r * sn + track_param->GetParameter()[2] * cs);
    Double_t m35 = pt;
    Double_t m45 = -pt * pt * track_param->GetParameter()[3];

    m43 *= track_param->GetSign();
    m44 *= track_param->GetSign();
    m45 *= track_param->GetSign();

    const Double_t* cTr = track_param->GetCovariance();
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
KFVertex CreateKFVertex(const AliVVertex* vertex) {

    Double_t param[6];
    vertex->GetXYZ(param);

    Double_t cov[6];
    vertex->GetCovarianceMatrix(cov);

    KFPVertex kfpVtx;
    Float_t paramF[3] = {(Float_t)param[0], (Float_t)param[1], (Float_t)param[2]};
    kfpVtx.SetXYZ(paramF);
    Float_t covF[6] = {(Float_t)cov[0], (Float_t)cov[1], (Float_t)cov[2], (Float_t)cov[3], (Float_t)cov[4], (Float_t)cov[5]};
    kfpVtx.SetCovarianceMatrix(covF);

    KFVertex KFVtx(kfpVtx);

    return KFVtx;
}

#endif  // TASKSEXAQUARK_KF_H
