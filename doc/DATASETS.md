# Datasets

## Run Numbers

Choosing the run numbers with `Central Barrel Tracking_hadron PID`. In detail:

- `LHC15o_pass2_rn.txt` (total: 144)
- `LHC18q_pass3_rn.txt` (total: 135)
- `LHC18r_pass3_rn.txt` (total: 90)
- `LHC18qr_pass3_rn.txt` (total: 225)
- `LHC18qr_extra_pass3_rn.txt` (total: 232)
  - Added 7 extra run numbers! Only valid in data and gen. purpose MC
  - Reference: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGRunList18r

## Simulations

- **Dedicated for anti-sexaquark and anti-neutron studies**

  - `LHC23l1`

    - Anchored to: `LHC18qr_pass3`, `LHC15o_pass2`
    - Transport: Geant3
    - Jira Ticket: [8862](https://alice.its.cern.ch/jira/browse/ALIROOT-8862)

- **Dedicated for $\Lambda(1520)$ studies**

  1. `LHC20g14[a,b,c]`

     - Anchored to: `LHC18qr_pass3`
     - Jira Ticket: [8525](https://alice.its.cern.ch/jira/browse/ALIROOT-8525)

  2. `LHC21l7[a,b,c]`

     - Anchored to: `LHC15o_pass2`
     - Jira Ticket: [8773](https://alice.its.cern.ch/jira/browse/ALIROOT-8773)

- **General purpose**

  1. `LHC20e3a`

     - Anchored to: `LHC18qr_pass3`
     - Generator: HIJING
     - Transport: Geant3
     - Jira Ticket: [8462](https://alice.its.cern.ch/jira/browse/ALIROOT-8462)

  2. `LHC20j6a`

     - Anchored to: `LHC15o_pass2`
     - Generator: HIJING
     - Transport: Geant3
     - Jira Ticket: [8568](https://alice.its.cern.ch/jira/browse/ALIROOT-8568)

## Data

- `LHC18qr`

  - Collision system: Pb-Pb, 5.02 TeV
  - Pass: `pass3`

- `LHC15o`

  - Collision system: Pb-Pb, 5.02 TeV
  - Pass: `pass2`
