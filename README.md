# Overview

Flat tree for electron trigger efficiency study in HWW channels.

The definition of **efficiency** is under the condition of some specific electron ID.
The study is based on tag-and-probe(TnP) method, the input source can be `SingleElectron` or `DoubleEG` `MINIAOD` data. In the consideration of testing various IDs with a single ntuple, the effect electron IDs are separated from the other cuts. If a tag/probe passed all the cut under some PID, then `tagPID`/`probePID` will be true. `passProbeXXLegs` is kept when all probe condition are passed under some PID.

For the PID, some cut other than the official electron IDs are applied, they're summarized in `plugins/MyPatElectronMVAIdSelector.cc`.


## Setup

```bash
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv

# Retrieve CutBased, MVA xml files, EffectiveArea, and etc..
# Ref: https://twiki.cern.ch/twiki/bin/view/CMSPublic/ElectronTagAndProbe
git cms-init
git cms-merge-topic fcouderc:tnp_egm_80X_Moriond17_v1.0
git clone -b egm_id_80X_v1 https://github.com/ikrav/RecoEgamma-ElectronIdentification.git ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data
git clone -b egm_id_80X_v1 https://github.com/ikrav/RecoEgamma-PhotonIdentification.git ../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data
rm -f ../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V0.weights.xml
rm -f ../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V0.weights.xml

# Get ntupler
git clone ssh://git@github.com/pageyu/ElectronEfficiencyHWW EleEfficiencyHWW

scram b -j 8
```

## Run

```bash

cd ${CMSSW_BASE}/src/EleEfficiencyHWW/runEffAnaMiniAOD
cmsRun runEffAnaMiniAOD_cfg.py

```
