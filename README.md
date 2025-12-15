**L1T Kalman Filter Upgrade**

Code related to algorithm upgrades of the CMS Level 1 Muon Trigger.

If a new plot needs to be made, add the function in `interactive_plotters.py`.
Each function takes `data` and `station`, where `data` comes from the main
`event_loop` function in `dummy_z.py`.

Example that plots Δz vs pT for station 1:

```text
$ python3 -i interactive_plotters.py
  successful event loop. events: 1000
  plot_delta_z_vs_pT(data, 1)
```
