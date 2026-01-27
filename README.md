**L1T Kalman Filter Upgrade**

Code related to algorithm upgrades of the CMS Level 1 Muon Trigger.

Step 1 is to process data which uses `submitChain.py` together with `run.py`. `run.py` specifies settings of the driver command which produces the output data files which will ultimately be sent to your eos area configured inside `run.py`. The code is currently setup to process a Drell-Yan full dataset and points to my eos area in `run.py`. change these and also change the settings in `submitChain.py` to keep more collection modules or remove some depending on the goals: 

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
