#ifndef MATH_HXX
#define MATH_HXX

//________________________________________________________________________
Double_t TwoLinesDCA(TVector3 pos1, TVector3 dir1,  //
                     TVector3 pos2, TVector3 dir2,  //
                     TVector3& P1, TVector3& P2) {
  //
  // find a common vertex for two neutral particles, given their
  // this function stores the origin vertices for each particle, and returns the DCA between them
  //

  // require perpendicularity <-> lines must intersect
  if (dir1.Cross(dir2).Mag() > 0.) {
    // based on https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d
    // (Marty Cohen's answer)

    TVector3 diff = pos1 - pos2;
    Double_t determinant = dir1.Dot(dir2) * dir1.Dot(dir2) - dir1.Mag2() * dir2.Mag2();
    Double_t sol1 = (dir2.Mag2() * dir1.Dot(diff) - dir2.Dot(diff) * dir2.Dot(dir1)) / determinant;
    Double_t sol2 = (-1 * dir1.Mag2() * dir2.Dot(diff) + dir1.Dot(diff) * dir2.Dot(dir1)) / determinant;

    // Sven's solution
    if (sol1 >= 0. || sol2 >= 0.) {
      printf("TwoLinesDCA :: Not possible to intersect in the right direction\n");
      return -1.;
    }

    // endpoints of the closest line
    P1 = pos1 + sol1 * dir1;
    P2 = pos2 + sol2 * dir2;

    // the distance of the closest line
    return (diff + dir1 * sol1 - dir2 * sol2).Mag();
  }
  printf("TwoLinesDCA :: Not possible to intersect!\n");
  return -1.;
}

//________________________________________________________________________
Double_t LinePointDCA(TVector3 part_momentum, TVector3 part_vertex, TVector3 ref_vertex) {
  //
  // Find the distance of closest approach to a reference point, after backtracking a particle.
  // [copied from AliAnalysisTaskSexaquark.cxx]
  //
  TVector3 cross_product = (ref_vertex - part_vertex).Cross(part_momentum);

  return cross_product.Mag() / part_momentum.Mag();
}

//________________________________________________________________________
Double_t CPA(TVector3 part_momentum, TVector3 part_vertex, TVector3 ref_vertex) {
  //
  // calculates the cosine of the pointing angle
  // of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
  // [based on AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()]
  //

  // vector between the reference point and the particle's vertex
  TVector3 delta = part_vertex - ref_vertex;

  // (protection)
  if (delta.Mag() == 0. || part_momentum.Mag() == 0.) {
    return -2.;
  }

  return part_momentum.Dot(delta) / (part_momentum.Mag() * delta.Mag());
}

//________________________________________________________________________
void EvaluateHelix(Double_t aliHelix_params[6], Double_t t, Double_t r[3]) {
  //
  // Calculate position of a point on a track and some derivatives at a given phase
  //
  Double_t phase = aliHelix_params[4] * t + aliHelix_params[2];
  Double_t sn = TMath::Sin(phase);
  Double_t cs = TMath::Cos(phase);

  r[0] = aliHelix_params[5] + sn / aliHelix_params[4];
  r[1] = aliHelix_params[0] - cs / aliHelix_params[4];
  r[2] = aliHelix_params[1] + aliHelix_params[3] * t;
}

//________________________________________________________________________
void HelixPointDCA(Double_t aliHelix_params[6],  //
                   TVector3 space_vec,           //
                   Float_t path_initA,           // 0. ?
                   Float_t path_initB,           // 30. ?
                   Float_t& pathA, Float_t& dcaAB) {
  //
  // << under construction >>
  //
  Float_t pA[2] = {path_initA, path_initB};  // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
  Float_t distarray[2];

  TVector3 testA;
  for (Int_t r = 0; r < 2; r++) {
    Double_t helix_point[3];
    EvaluateHelix(aliHelix_params, pA[r], helix_point);
    testA.SetXYZ(helix_point[0], helix_point[1], helix_point[2]);  // 3D-vector of helixA at path pA[r]
    distarray[r] = (testA - space_vec).Mag();                      // dca between helixA and helixB
  }

  Int_t loopcounter = 0;
  Float_t scale = 1.0;
  Float_t flip = 1.0;  // checks if the minimization direction changed
  Float_t scale_length = 30.0;

  while (TMath::Abs(scale_length) > 0.1 && loopcounter < 100)  // stops when the length is too small
  {
    if (distarray[0] > distarray[1]) {
      if (loopcounter != 0) {
        if (flip == 1.0) {
          // if minimization direction changes -> go back, but only the way * 0.4
          scale = 0.4;
        } else {
          // go on in this direction but only by the way * 0.7
          scale = 0.7;
        }
      }
      scale_length = (pA[1] - pA[0]) * scale;  // the next length interval
      pA[0] = pA[1] + scale_length;            // the new path

      Double_t helix_point[3];
      EvaluateHelix(aliHelix_params, pA[0], helix_point);
      testA.SetXYZ(helix_point[0], helix_point[1], helix_point[2]);  // 3D-vector of helixA at path pA[0]
      distarray[0] = (testA - space_vec).Mag();                      // new dca
      flip = 1.0;

    } else {
      if (loopcounter != 0) {
        if (flip == -1.0) {
          scale = 0.4;
        } else {
          scale = 0.7;
        }
      }
      scale_length = (pA[0] - pA[1]) * scale;
      pA[1] = pA[0] + scale_length;

      Double_t helix_point[3];
      EvaluateHelix(aliHelix_params, pA[1], helix_point);
      testA.SetXYZ(helix_point[0], helix_point[1], helix_point[2]);  // 3D-vector of helixA at path pA[1]
      distarray[1] = (testA - space_vec).Mag();
      flip = -1.0;
    }
    loopcounter++;
  }

  if (loopcounter >= 100) {
    std::cout << "WARNING: FindDCAHelixPoint exceeded maximum of 100 loops" << std::endl;
  }

  if (distarray[0] < distarray[1]) {
    pathA = pA[0];
    dcaAB = distarray[0];
  } else {
    pathA = pA[1];
    dcaAB = distarray[1];
  }
}

#endif
