## Contributors 
```Alibordi, Artur```
-The code is reading the same data format produced by the set-up 
-https://github.com/mumuhamm/UserCode/tree/devel_cmssw13_algosKB
-This set-up takes latest improvements by implementing OMTF algorithms for L1 Trigger 


### Installation instructions:

```
git clone -b devel_ANYDATAFORMAT https://github.com/mumuhamm/RootAnalysis.git 
cd RootAnalysis
mkdir build; cd build
cmake ../
make install -j 4
```
### Remember to 
```git checkout relevant_tag```


### Run instructions:

Before running update path to the data files in the config/omtf_emulator.ini file

```
cd RootAnalysis/build
./bin/omtfAnalysis config/omtf_emulator.ini
```
The location of the sample produced by the *UserCode* is mentioned in th *.ini* file 
Please change the tag 

```processName = PRIVATE_omtfTree``` : while running on the CMSSW data format (AOD,MINIAOD,RAW-RECO)
```processName = PRIVATE_nanoAOD```   : while runnning on the nano-AOD

The resulting plots are stored in the fig_png directory.

### Statements 
- Comming soon , requires a bit of explanation on the : displaced muons, rate analysis , nano - no subsystem information 
