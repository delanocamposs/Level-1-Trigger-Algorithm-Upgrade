import ROOT 
import numpy as np 
from dummy_z import event_loop

store_plots={"canvas":{}, "histos":{}}

def plot_rho(data, station):
    ROOT.gROOT.SetBatch(False)
    c=ROOT.TCanvas(f"c_{station}", "", 800, 600)
    h=ROOT.TH1F(f"h_{station}", f"Station {station};#rho (cm);Events",100, 0, 700)
    h.SetDirectory(0)
    for vr, eta_gen, z_stub, vz in zip(data["gen_vr"][station], data["gen_eta"][station], data["stub_z"][station], data["gen_vz"][station]):
        theta_gen=2*np.arctan(np.exp(-1*eta_gen))
        rho=vr+(z_stub-vz)*np.tan(theta_gen)
        h.Fill(rho)
    max_bin=h.GetMaximumBin()
    peak_rho=h.GetXaxis().GetBinCenter(max_bin)
    h.GetXaxis().SetRangeUser(peak_rho-(peak_rho)/2, peak_rho+(peak_rho)/2)
    h.SetStats(0)
    h.Draw("")
    store_plots['canvas'][f'rho_{station}']=c
    store_plots['histos'][f'rho_{station}']=h
    c.Update()
    c.SaveAs(f"rho_{station}.png")
    return c, h

def plot_delta_z_vs_pT(data, station):
    ROOT.gROOT.SetBatch(False)
    c=ROOT.TCanvas(f"c_{station}", "", 800, 600)
    h=ROOT.TH2F(f"h_{station}", f"Station {station};pT (GeV);#Delta z (cm)",50, 0, 100, 50, -100, 100)
    h.SetDirectory(0)
    for pt, dz in zip(data["gen_pt"][station], data["delta_z"][station]):
        h.Fill(pt, dz)
    h.SetStats(0)
    h.Draw("COLZ")
    store_plots['canvas'][f'delta_z_vs_pT_{station}']=c
    store_plots['histos'][f'delta_z_vs_pT_{station}']=h
    c.Update()
    c.SaveAs(f"delta_z_vs_pT_{station}.png")
    return c, h

def plot_sigma_z_vs_pT(data, station):
    c=ROOT.TCanvas(f"c_{station}", "", 800, 600)
    fit_gaus=ROOT.TF1(f"fit_gauss_{station}","gaus", -100, 100)
    if f'delta_z_vs_pT_{station}' not in store_plots['histos'].keys():
        raise RuntimeError("ERROR: run plot_delta_z_vz_pT function first. the sigma slice relies on delta_z data")
    h=store_plots['histos'][f'delta_z_vs_pT_{station}']
    h.SetDirectory(0)
    h.FitSlicesY(fit_gaus, 1, -1, 0, "QNR")
    h.SetTitle(f"Station {station};pT (GeV);#sigma [cm]")
    h.SetStats(0)
    h.SetMarkerStyle(8)
    h.SetMarkerSize(1.2)
    h.SetLineWidth(2)
    h.SetStats(0)
    h.Draw("PE")
    store_plots['canvas'][f'sigma_z_vs_pT_{station}']=c
    store_plots['histos'][f'sigma_z_vs_pT_{station}']=h
    c.Update()
    c.SaveAs(f"sigma_z_vs_pT_{station}.png")
    return 

if __name__=="__main__":
    data=event_loop()
    

