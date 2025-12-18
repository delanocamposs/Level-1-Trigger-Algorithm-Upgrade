**L1T Kalman Filter Upgrade**

Code related to algorithm upgrades of the CMS Level 1 Muon Trigger.

If a new plot needs to be made, add the function in `interactive_plotters.py`.
Each function takes `data` and `station`, where `data` comes from the main
`event_loop` function in `dummy_z.py`.

Example that plots Δz vs pT for station 1:

```text
$ python3 -i interactive_plotters.py
  successful event loop. events: 100
  plot_delta_z_vs_pT(data, 1)
```
To add more events and choose to not convert k, z out of digital units then plot propagation from MB1 to vertex using the non-converted data:

```text
$ python3 -i interactive_plotters.py
  successful event loop. events: 100
  data=event_loop(10000, False, False)
  plot_st1_to_vtx_vs_k1(data, False, False)
```

Or choose to do the same except with converted z and k:

```text
$ python3 -i interactive_plotters.py
  successful event loop. events: 100
  data=event_loop(10000)
  plot_st1_to_vtx_vs_k1(data)
```
