**L1T Kalman Filter Upgrade**

Code related to Kalman Filter algorithm upgrades of the CMS Level 1 Global Muon Trigger.

C++ algorithm found at cmssw: L1Trigger/Phase2L1GMT in this link: https://github.com/cms-sw/cmssw/tree/master/L1Trigger/Phase2L1GMT

Step 1 is to setup an up-to-date cmsenv, and run a cmsDriver command that runs the emulator and will produce a python file (--fileout argument), and a processed data file. The example below processes 10 events with -n for a DY sample through --filein argument and runs L1 and L1P2GT processes for the emulator through -s argument. Per these insturctions for the emulator recipe: https://cmsl1tofflinesw.docs.cern.ch/Phase2/phase2recipe/ look for a recent version of CMSSW with scram then setup an env with it. Since this was developing mostly the Kalman Filter code, I use cms-addpkg L1Trigger/Phase2L1GMT:

```text
scram list CMSSW_16_1_X_*
cmsrel CMSSW_16_1_X_
cd CMSSW_16_1_X_*/src
cmsenv
git cms-init
git cms-addpkg L1Trigger/Phase2L1GMT
scram b -j 8
git checkout -b mybranch
```
```text
cmsDriver.py -s L1,L1P2GT --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2.addHcalTriggerPrimitives --filein /store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/5bed6cad-6cdb-4a5e-87af-08cc8c8e0ff8.root --fileout file:output_Phase2_L1T.root --python_filename run.py --inputCommands="keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT" --mc -n 10 --nThreads 1
```

That will create your run.py file with many processes and algorithms which make processing the events VERY slow. For Kalman Filter Algorithm development, changes to the C++ code happen regularly and waiting long durations to simply process 10 events is not feasible. To solve this problem, build a "reduced" run.py called `run_reduced.py` which is a fraction of the processes before. Only the necessary algorithms, processes and modules/collections are kept to study the KMTF quickly. One must run the previous cmsDriver.py command first since the new `run_reduced.py` will take the first output as input.

The idea is that the driver command is run once with all processes/collections, etc, then  C++ algorithm modifications can be made Edits to the C++ algorithm can be made, then events reprocessed quickly with:

```text
cmsRun run_reduced.py
```

If that works fine, the next step is to process a full dataset which uses `submitChain.py` together with `run.py`. The output data files which will ultimately be sent to your eos area configured inside `run.py`. The code is currently setup to process a Drell-Yan full dataset and points to my eos area in `run.py`. change these and also change the settings in `submitChain.py` to keep more collection modules or remove some depending on the goals: 

```text
voms-proxy-init --voms cms --valid 192:00
python3 submitChain.py
```

If a new plot needs to be made, add the function in `interactive_plotters.py`.
Each function takes `data` and `station`, where `data` comes from the main
`event_loop` function in `main.py`.

Example that plots Δz vs pT for station 1:

```text
$ python3 -i interactive_plotter.py
  successful event loop. events: 100
  plot_delta_z_vs_pT(data, 1)
```
To add more events and choose to not convert k, z out of digital units then plot propagation from MB1 to vertex using the non-converted data:

```text
$ python3 -i interactive_plotter.py
  successful event loop. events: 100
  data=event_loop(10000, False, False)
  plot_st1_to_vtx_vs_k1(data, False, False)
```

Or choose to do the same except with converted z and k:

```text
$ python3 -i interactive_plotter.py
  successful event loop. events: 100
  data=event_loop(10000)
  plot_st1_to_vtx_vs_k1(data)
```
