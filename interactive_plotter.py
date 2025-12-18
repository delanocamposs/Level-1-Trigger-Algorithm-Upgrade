import ROOT 
import numpy as np 
from dummy_z import event_loop
from trigger_helpers import *

store_plots={"canvas":{}, "histos":{}, "fits":{}, "profiles":{}}

def plot_rho(data,station):
    ROOT.gROOT.SetBatch(False)
    c=ROOT.TCanvas(f"c_rho_{station}","",800,600)
    rhos=[]
    for vr,eta_gen,z_stub,vz in zip(data["gen_vr"][station],data["gen_eta"][station],data["stub_z"][station],data["gen_vz"][station]):
        theta_gen=2*np.arctan(np.exp(-1*eta_gen))
        rho=vr+(z_stub-vz)*np.tan(theta_gen)
        if rho>300 and rho<800:
            rhos.append(rho) #the tails of the distributiosn were messing up the plotting
    rhos=np.array(rhos)
    mu=rhos.mean()
    print(mu)
    rms=rhos.std()
    lo=mu-5*rms
    hi=mu+5*rms
    h=ROOT.TH1F(f"h_rho_{station}",f"Station {station};#rho (cm);Events",100, lo, hi)
    h.SetDirectory(0)
    for rho in rhos:
        h.Fill(rho)
    h.SetStats(0)
    h.Draw("HIST")
    fit=ROOT.TF1(f"fit_rho_{station}","gaus",lo,hi)
    fit.SetLineWidth(2)
    res=h.Fit(fit,"RQS")
    fit.Draw("SAME")
    mean=fit.GetParameter(1)
    mean_err=fit.GetParError(1)
    sigma=fit.GetParameter(2)
    sigma_err=fit.GetParError(2)
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextSize(0.035)
    txt.DrawLatex(0.15,0.80,f"#mu={mean:.3g} #pm {mean_err:.2g} cm")
    txt.DrawLatex(0.15,0.75,f"#sigma={sigma:.3g} #pm {sigma_err:.2g} cm")
    store_plots["canvas"][f"rho_{station}"]=c
    store_plots["histos"][f"rho_{station}"]=h
    store_plots["fits"][f"rho_{station}"]=fit
    c.Update()
    c.SaveAs(f"rho_{station}.png")
    return c,h,fit,(mean,mean_err,sigma,sigma_err)

def plot_delta_z_vs_pT(data,station):
    ROOT.gROOT.SetBatch(False)
    c=ROOT.TCanvas(f"c_delta_{station}","",800,600)
    h=ROOT.TH2F(f"h_{station}",f"Station {station};pT (GeV);#Delta z (cm)",50,0,100,50,-100,100)
    h.SetDirectory(0)
    for pt,dz in zip(data["gen_pt"][station],data["delta_z"][station]):
        h.Fill(pt,dz)
    h.SetStats(0)
    h.Draw("COLZ")
    store_plots["canvas"][f"delta_z_vs_pT_{station}"]=c
    store_plots["histos"][f"delta_z_vs_pT_{station}"]=h
    c.Update()
    c.SaveAs(f"delta_z_vs_pT_{station}.png")
    return c,h

def plot_sigma_z_vs_pT(data,station):
    key=f"delta_z_vs_pT_{station}"
    if key not in store_plots["histos"]:
        raise RuntimeError("run plot_delta_z_vs_pT first")
    h=store_plots["histos"][key]
    h.SetDirectory(0)
    c=ROOT.TCanvas(f"c_sigma_{station}","",800,600)
    fit_gaus=ROOT.TF1(f"fit_gauss_{station}","gaus",-100,100)
    h.FitSlicesY(fit_gaus,1,-1,0,"QNR")
    h_name=f"{h.GetName()}_2"
    h_sigma=ROOT.gDirectory.Get(h_name)
    if not h_sigma:
        raise RuntimeError(f"missing FitSlicesY output {h_name}")
    h_sigma=h_sigma.Clone(f"sigma_z_vs_pT_{station}")
    h_sigma.SetDirectory(0)
    h_sigma.SetTitle(f"Station {station};pT (GeV);#sigma [cm]")
    h_sigma.SetStats(0)
    h_sigma.SetMarkerStyle(8)
    h_sigma.SetMarkerSize(1.2)
    h_sigma.SetLineWidth(2)
    h_sigma.Draw("PE")
    store_plots["canvas"][f"sigma_z_vs_pT_{station}"]=c
    store_plots["histos"][f"sigma_z_vs_pT_{station}"]=h_sigma
    c.Update()
    c.SaveAs(f"sigma_z_vs_pT_{station}.png")
    return

def plot_Yslices_of_delta_z(data, station):
    ROOT.gROOT.SetBatch(True)
    if f'delta_z_vs_pT_{station}' not in store_plots['histos'].keys():
        raise RuntimeError("ERROR: run plot_delta_z_vz_pT function first. the sigma slice relies on delta_z data")
    h=store_plots['histos'][f'delta_z_vs_pT_{station}']
    for i in range(1,h.GetNbinsX()+1):
        c_sl,h_sl,g_sl=plot_slice_and_gaussian(h,i,station,-100,100)
    print("all y slice plots saved.")
    return 

def plot_deltak_and_profile_vs_eta(data,st1,st2):
    h=ROOT.TH2F(f"h_dk_eta_{st1}_{st2}",f"k_{st2}-k_{st1} vs stub #eta_{st1};stub #eta;#Delta k",60,-1.3,1.3,80,-1.0,1.0)
    h.SetDirectory(0)
    m1={}
    for muid,k,eta in zip(data["mu_id"][st1],data["stub_k"][st1],data["stub_eta"][st1]):
        if muid not in m1:
            m1[muid]=(k,eta)
    for muid,k2 in zip(data["mu_id"][st2],data["stub_k"][st2]):
        if muid not in m1:
            continue
        k1,eta1=m1[muid]
        h.Fill(eta1,k2-k1)
    c1=ROOT.TCanvas(f"c_dk_eta_{st1}_{st2}","",900,700)
    h.SetStats(0)
    h.Draw("COLZ")
    c1.Update()
    c1.SaveAs(f"dk_vs_eta1_MB{st2}_MB{st1}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"<k_{st2}-k_{st1}> vs stub #eta_{st1};stub #eta;<#Delta k>")
    c2=ROOT.TCanvas(f"c_dkprof_eta_{st1}_{st2}","",900,700)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"dk_profile_vs_eta1_MB{st2}_MB{st1}.png")
    store_plots["canvas"][f"dk_vs_eta1_{st1}_{st2}"]=c1
    store_plots["canvas"][f"dkprof_vs_eta1_{st1}_{st2}"]=c2
    store_plots["histos"][f"dk_vs_eta1_{st1}_{st2}"]=h
    store_plots["profiles"][f"dkprof_vs_eta1_{st1}_{st2}"]=p
    return c1,c2,h,p

def plot_dz_vs_1_over_k1(data,st1,st2):
    h=ROOT.TH2F(f"h_dz_invk_{st1}_{st2}",f"z{st1}-z{st2} vs 1/k{st2};1/k{st2};#Delta z (z{st1}-z{st2})",100,-5,5,100,-300,300)
    h.SetDirectory(0)
    m1={}
    for muid,k2,z2 in zip(data["mu_id"][st2],data["stub_k"][st2],data["stub_z"][st2]):
        if muid not in m1:
            m1[muid]=(k2,z2)
    #print("m1: ",m1)
    for muid,z1 in zip(data["mu_id"][st1],data["stub_z"][st1]):
        if muid not in m1:
            continue
        k2,z2=m1[muid]
        #min value to avoiud blowing up 1/k
        if abs(k2)<0.001:
            continue
        h.Fill(1.0/k2,z1-z2)
    c1=ROOT.TCanvas(f"c_dz_1_over_k{st2}_MB{st1}_MB{st2}","",900,700)
    h.SetStats(0)
    h.Draw("COLZ")
    c1.Update()
    c1.SaveAs(f"dz_vs_1_over_k{st2}_MB{st2}_MB{st1}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"<z{st1}-z{st2}> vs stub 1/k{st2};1/k{st2};<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprofile_1_over_k{st2}_MB{st1}_MB{st2}","",900,700)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"dz_profile_vs_1_over_k{st2}_MB{st2}_MB{st1}.png")
    store_plots["canvas"][f"dz_vs_1_over_k{st2}_MB{st1}_MB{st2}"]=c1
    store_plots["canvas"][f"dzprof_vs_1_over_k{st2}_MB{st1}_MB{st2}"]=c2
    store_plots["histos"][f"dz_vs_1_over_k{st2}_MB{st1}_MB{st2}"]=h
    store_plots["profiles"][f"dzprof_vs_1_over_k{st2}_{st1}_{st2}"]=p
    return c1,c2,h,p

def plot_st1_to_vtx_vs_k1(data,conv_z=True, conv_k=True):
    if conv_k and conv_z:
        h=ROOT.TH2F(f"h_dz_invk_st1_vtx",f"z_vtx-z1 vs k1;k1;#Delta z (z_vtx-z1)",100,-1.75,1.75,100,-250,250)
    else:
        h=ROOT.TH2F(f"h_dz_invk_st1_vtx",f"z_vtx-z1 vs k1;k1;#Delta z (z_vtx-z1)",100,-50000,50000,100,-40000,40000)
    h.SetDirectory(0)
    m1={}
    for muid,k2,z2 in zip(data["mu_id"][1],data["stub_k"][1],data["stub_z"][1]):
        if muid not in m1:
            m1[muid]=(k2,z2)
    for muid,z1 in zip(data["mu_id"][1],data["gen_vz"][1]):
        if muid not in m1:
            continue
        k2,z2=m1[muid]
        h.Fill(k2,z1-z2)
    c1=ROOT.TCanvas(f"c_dz_vs_k1_MB1_vtx","",900,700)
    h.SetStats(0)
    h.Draw("COLZ")
    c1.Update()
    c1.SaveAs(f"dz_vs_k1_MB1_vtx.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"<z_vtx-z1> vs k1;k1;<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprofile_k1_MB1_vtx","",900,700)
    p.Draw("E")
    x_min=p.GetXaxis().GetXmin()
    x_max=p.GetXaxis().GetXmax()
    f_lin=ROOT.TF1(f"f_lin_MB1_vtx", "pol1", x_min, x_max)
    f_lin.SetLineColor(ROOT.kRed)
    f_lin.SetLineWidth(2)
    fitres=p.Fit(f_lin, "SQR")
    p0=f_lin.GetParameter(0)
    p1=f_lin.GetParameter(1)
    p0e=f_lin.GetParError(0)
    p1e=f_lin.GetParError(1)
    p.Draw("E")
    f_lin.Draw("same")
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextSize(0.035)
    if conv_k and conv_z:
        txt.DrawLatex(0.15,0.80,f"slope={p1:.3g} cm #pm {p1e:.2g} cm")
        txt.DrawLatex(0.15,0.75,f"int={p0:.3g} cm #pm {p0e:.2g} cm")
    else:
        txt.DrawLatex(0.15,0.80,f"slope={p1:.3g} #pm {p1e:.2g}")
        txt.DrawLatex(0.15,0.75,f"int={p0:.3g} #pm {p0e:.2g}")

    c2.Update()
    c2.SaveAs(f"dz_profile_vs_k1_MB1_vtx.png")
    store_plots["canvas"][f"dz_vs_k1_MB1_vtx"]=c1
    store_plots["canvas"][f"dzprof_vs_k1_MB1_vtx"]=c2
    store_plots["histos"][f"dz_vs_k1_MB1_vtx"]=h
    store_plots["profiles"][f"dzprof_vs_k1_MB1_vtx"]=p
    return c1,c2,h,p
 
    
if __name__=="__main__":
    data=event_loop()
