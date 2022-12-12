CUTS
====

(date: 7/December/2022)

This document corresponds to a summary of all the cuts that were applied on the different stages of this analysis.

## I. AliAnalysisTaskSexaquark

1) Cut on the PID of **MC Particles**. Only keep: positive pion, negative pion, positive kaons, neutral kaon short, anti-lambda, anti-proton, Xi plus. These are the **relevant MC particles**.

2) Cut on the PID of **Reconstructed Tracks**. Only keep the ones that are pions, kaons, and protons within 3 sigma of the expected TPC signal. These are the **relevant rec. particles**.

3) Also, keep the **Reconstructed Tracks** that come from a true positive pion, negative pion, positive kaon or anti-proton. Regardless if they were misidentified as something else. These are the **relevant MC particles -> non-relevant rec. particles**.

4) Also, keep the **MC Particles** that were reconstructed as a pion, kaon or proton. Regardless of their true PID. For example, a true positron that was reconstructed as a proton. These are the **non-relevant MC particles -> relevant rec. particles**.

* **True V0s**:

  While looping over **MC Particles**,

  1) Cut on PID, particle must be a true neutral kaon short or anti-lambda.

  2) Cut on number of daughters. It must be exactly 2.

  3) Cut on PID of daughters to make sure it's a **relevant decay channel**. That means, neutral kaon short into positive and negative pion, and anti-lambda into positive pion and anti-proton.

  4) Both daughters must have been reconstructed and belong to the **Reconstructed Tracks** selection (see above).

  Finally, store all possible positive-negative tracks combinations that survived this cuts.

* **Official V0s**:

  While looping over the V0s stored in the ESD file:

  1) Both daughters must have been reconstructed and belong to the **Reconstructed Tracks** selection.

  2) Cut on fiducial radius of the V0 vertex, it must be between 5 and 180 cm.

  Finally, store V0s that obbey these cuts.

* **Custom V0s**:

  First, loop over all **Reconstructed Tracks**, and store them according to its charge. Then, loop over every positive-negative possible. While doing so,

  1) Replicate the **Offline V0 Finder**, only excluding the CPA cut.

  2) Both daughters must have been reconstructed and belong to the **Reconstructed Tracks** selection.

  Finally, store, the found V0s that obey those cuts.

## II. SexaquarkFinder

  Loop over events. Then, loop over all the possible V0 pairs in one event. While doing so,

  1) Check that the positive daughter of the V0 "A" is different than the positive daughter of the V0 "B".

  2) One of the V0s must be a neutral kaon short and the other one must be an anti-lambda.

  Finally, store V0 pair as an anti-sexaquark candidate.

## III. Visualization Macros


