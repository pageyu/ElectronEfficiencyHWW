# Flat tree for trigger efficiency study in HWW channels.
**Why a new recipe is needed?**
Quite a simple reason, it's not straitforward to apply some cuts such as `dz`, `d0` in the official recipe, while in Arun's recipe, the MVA value and valueMaps are inaccessible in MiniAOD. The only way is running under full CMSSW frame. It's nothing more than porting Arun's recipe from fwlite to the full framework.

```bash
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv

git cms-init
git cms-merge-topic fcouderc:tnp_egm_80X_Moriond17_v1.0

git clone -b egm_id_80X_v1 https://github.com/ikrav/RecoEgamma-ElectronIdentification.git ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data
git clone -b egm_id_80X_v1 https://github.com/ikrav/RecoEgamma-PhotonIdentification.git ../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data

rm -f ../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V0.weights.xml
rm -f ../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V0.weights.xml

git clone ssh://git@github.com/pageyu/ElectronEfficiencyHWW EleEfficiencyHWW

scram b -j 4
```
