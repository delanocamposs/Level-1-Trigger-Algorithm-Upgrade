import ROOT 
import numpy as np 
from dummy_z import event_loop
from trigger_helpers import *
import os 

store_plots={"canvas":{}, "histos":{}, "fits":{}, "profiles":{}}

def plot_slice_and_gaussian(h2, xbin, st, xvar,xunit,yvar,yunit,y_min=-0.5, y_max=0.5):
    direc=make_plot_dir(f"gaussianfit_slices_{yvar}_vs_{xvar}")
    name = h2.GetName()
    h_int=ROOT.gDirectory.Get(f"{name}_0")
    h_mean=ROOT.gDirectory.Get(f"{name}_1")
    h_sigma=ROOT.gDirectory.Get(f"{name}_2")
    A=h_int.GetBinContent(xbin)
    mu=h_mean.GetBinContent(xbin)
    si=h_sigma.GetBinContent(xbin)
    slice_name= f"{name}_slice_xbin{xbin}"
    h_slice =h2.ProjectionY(slice_name, xbin, xbin)
    g_name =f"{name}_gaus_xbin{xbin}"
    g = ROOT.TF1(g_name, "gaus", y_min, y_max)
    g.SetParameters(A, mu, si)
    c = ROOT.TCanvas(f"c_{name}_xbin{xbin}", "", 800, 600)
    xaxis= h2.GetXaxis()
    x_low= xaxis.GetBinLowEdge(xbin)
    x_high= xaxis.GetBinUpEdge(xbin)
    h_slice.SetTitle(f"Station {st}: {x_low:.1f}<{xvar}<{x_high:.1f} [{xunit}];{yvar} [{yunit}];Entries")
    h_slice.SetLineColor(ROOT.kBlack)
    h_slice.SetMarkerStyle(20)
    h_slice.SetMarkerSize(1.0)
    h_slice.SetStats(0)
    h_slice.Draw()
    g.SetLineColor(ROOT.kRed)
    g.SetLineWidth(2)
    g.Draw("SAME")
    leg = ROOT.TLegend(0.65, 0.75, 0.89, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(0, f"#sigma = {si:.3f} {yunit}", "")
    leg.AddEntry(0, f"mean = {mu:.3f} {yunit}", "")
    leg.SetTextSize(0.035)
    leg.Draw()
    c.leg=leg
    c.SaveAs(f"{direc}/slice{yvar}_vs_{xvar}_bin{xbin}_station_{st}.png")
    c.Close()
    return c, h_slice, g

def plot_rho(data,station):
    direc=make_plot_dir("rho")
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
    c.SaveAs(f"{direc}/rho_{station}.png")
    return c,h,fit,(mean,mean_err,sigma,sigma_err)

def plot_delta_z_vs_pT(data,station):
    direc=make_plot_dir("deltaz_vs_pT")
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
    c.SaveAs(f"{direc}/delta_z_vs_pT_{station}.png")
    return c,h

def plot_sigma_vs_var(data,station,var1,var2,key):
    direc=make_plot_dir(f"sigma{var1}_vs_{var2}")
    if key not in store_plots["histos"]:
        raise RuntimeError(f"missing histogram {key}. run function that makes histo first.")
    h=store_plots["histos"][key]
    h.SetDirectory(0)
    ymean=h.GetMean(2)
    yrms=h.GetRMS(2)
    h.RebinX(2)
    c=ROOT.TCanvas(f"c_sigma{var1}_vs_{var2}_{station}","",800,600)
    fit_gaus=ROOT.TF1(f"fit_gauss_{var2}_{station}","gaus",ymean-3*yrms,ymean+3*yrms)
    h.FitSlicesY(fit_gaus,1,-1,0,"QNR")
    h_name=f"{h.GetName()}_2"
    h_sigma=ROOT.gDirectory.Get(h_name)
    if not h_sigma:
        raise RuntimeError(f"missing FitSlicesY output {h_name}")
    h_sigma=h_sigma.Clone(f"sigma{var1}_vs_{var2}_{station}")
    h_sigma.SetDirectory(0)
    h_sigma.SetTitle(f"Station {station};{h.GetXaxis().GetTitle()};#sigma")
    h_sigma.SetStats(0)
    h_sigma.SetMarkerStyle(8)
    h_sigma.SetMarkerSize(1.2)
    h_sigma.SetLineWidth(2)
    c.SetLeftMargin(0.15)
    h_sigma.Draw("PE")
    store_plots["canvas"][f"sigma{var1}_vs_{var2}_{station}"]=c
    store_plots["histos"][f"sigma{var1}_vs_{var2}_{station}"]=h_sigma
    c.Update()
    c.SaveAs(f"{direc}/sigma{var1}_vs_{var2}_{station}.png")
    return

def plot_yslices(data,station,key,xvar,xunit,yvar,yunit):
    ROOT.gROOT.SetBatch(True)
    if key not in store_plots["histos"]:
        raise RuntimeError(f"missing histogram {key}. run function that makes histo first.")
    h=store_plots["histos"][key]
    ymean=h.GetMean(2)
    yrms=h.GetRMS(2)
    h.RebinX(2)
    ROOT.gROOT.cd()
    fit_gaus=ROOT.TF1(f"fit_gauss_slices_{station}","gaus",ymean-3*yrms,ymean+3*yrms)
    h.FitSlicesY(fit_gaus,1,-1,0,"QNR")
    name=h.GetName()
    h_int=ROOT.gDirectory.Get(f"{name}_0")
    h_mean=ROOT.gDirectory.Get(f"{name}_1")
    h_sigma=ROOT.gDirectory.Get(f"{name}_2")
    if not h_int or not hasattr(h_int,"GetBinContent") or not h_mean or not hasattr(h_mean,"GetBinContent") or not h_sigma or not hasattr(h_sigma,"GetBinContent"):
        raise RuntimeError(f"FitSlicesY outputs missing for {name}: {type(h_int)}, {type(h_mean)}, {type(h_sigma)}")
    for i in range(1,h.GetNbinsX()+1):
        plot_slice_and_gaussian(h,i,station,xvar,xunit,yvar,yunit,y_min=ymean-3*yrms,y_max=ymean+3*yrms)
    print("all y slice plots saved.")
    return

def plot_deltak_and_profile_vs_eta(data,st1,st2):
    direc=make_plot_dir("deltak_vs_eta")
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
    c1.SaveAs(f"{direc}/dk_vs_eta1_MB{st2}_MB{st1}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"<k_{st2}-k_{st1}> vs stub #eta_{st1};stub #eta;<#Delta k>")
    c2=ROOT.TCanvas(f"c_dkprof_eta_{st1}_{st2}","",900,700)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"{direc}/dk_profile_vs_eta1_MB{st2}_MB{st1}.png")
    store_plots["canvas"][f"dk_vs_eta1_{st1}_{st2}"]=c1
    store_plots["canvas"][f"dkprof_vs_eta1_{st1}_{st2}"]=c2
    store_plots["histos"][f"dk_vs_eta1_{st1}_{st2}"]=h
    store_plots["profiles"][f"dkprof_vs_eta1_{st1}_{st2}"]=p
    return c1,c2,h,p

def plot_dz_vs_1_over_k(data,st1,st2):
    direc=make_plot_dir("dz_vs_1_over_k")
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
    c1.SaveAs(f"{direc}/dz_vs_1_over_k{st2}_MB{st2}_MB{st1}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"<z{st1}-z{st2}> vs stub 1/k{st2};1/k{st2};<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprofile_1_over_k{st2}_MB{st1}_MB{st2}","",900,700)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"{direc}/dz_profile_vs_1_over_k{st2}_MB{st2}_MB{st1}.png")
    store_plots["canvas"][f"dz_vs_1_over_k{st2}_MB{st1}_MB{st2}"]=c1
    store_plots["canvas"][f"dzprof_vs_1_over_k{st2}_MB{st1}_MB{st2}"]=c2
    store_plots["histos"][f"dz_vs_1_over_k{st2}_MB{st1}_MB{st2}"]=h
    store_plots["profiles"][f"dzprof_vs_1_over_k{st2}_{st1}_{st2}"]=p
    return c1,c2,h,p

def plot_st1_to_vtx_vs_k1(data,conv_z=True, conv_k=True):
    direc=make_plot_dir("dz_vs_k")
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
    c1=ROOT.TCanvas(f"{direc}/c_dz_vs_k1_MB1_vtx","",900,700)
    h.SetStats(0)
    h.Draw("COLZ")
    c1.Update()
    c1.SaveAs(f"{direc}/dz_vs_k1_MB1_vtx.png")
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

def plot_deltak_vs_curv(data,st):
    direc=make_plot_dir("deltak_vs_curv")
    c=ROOT.TCanvas(f"c_dk_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dk_curv_{st}",f"(k_{{pred,{st}}}-k_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta k (k_{{pred,{st}}}-k_{{meas,{st}}})",100,-10000,10000,100,-15000,15000)
    h.SetDirectory(0)
    m={}
    if st==3:
        raise RuntimeError("ERROR: cannot propagate to st3. (no st4 stub info)")
    for muid,k0,curv in zip(data["mu_id"][st+1],data["stub_k"][st+1], data["gen_curv"][st+1]):
        if muid not in m:
            m[muid]=(k0,curv)
    for muid,k_meas in zip(data["mu_id"][st],data["stub_k"][st]):
        if muid not in m:
            continue
        k0,curv=m[muid]
        k_pred=1*k0
        h.Fill(curv, k_pred-k_meas)
    h.SetStats(0)
    h.Draw("COLZ")
    store_plots["canvas"][f"deltak_vs_curv_{st}"]=c
    store_plots["histos"][f"deltak_vs_curv_{st}"]=h
    print(f"histo name stored: deltak_vs_curv_{st}")
    c.Update()
    c.SaveAs(f"{direc}/deltak_vs_curv_{st}.png")
    return c,h 



def plot_deltaz_vs_curv(data,st,show=True):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltaz_vs_curv")
    c=ROOT.TCanvas(f"c_dz_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dz_curv_{st}",f"(z_{{pred,{st}}}-z_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta z (z_{{pred,{st}}}-z_{{meas,{st}}})",100,-10000,10000,100,-15000,15000)
    h.SetDirectory(0)
    m={}
    if st==1:
        dR=-.101
    elif st==2:
        dR=-.133
    else:
        print("cannot propagate to station 3 (no stubs at st 4).")
    for muid,z0,curv,slope in zip(data["mu_id"][st+1],data["stub_z"][st+1],data["gen_curv"][st+1], data["stub_k"][st+1]):
        if muid not in m:
            m[muid]=(z0,curv, slope)
    for muid,z_meas in zip(data["mu_id"][st],data["stub_z"][st]):
        if muid not in m:
            continue
        z0,curv,slope=m[muid]
        z_pred=z0+slope*dR
        h.Fill(curv, z_pred-z_meas)
    h.SetStats(0)
    c.SetLeftMargin(0.15)
    if show:
        h.Draw("COLZ")
    store_plots["canvas"][f"deltaz_vs_curv_{st}"]=c
    store_plots["histos"][f"deltaz_vs_curv_{st}"]=h
    c.Update()
    c.SaveAs(f"{direc}/deltaz_vs_curv_{st}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"station: {st} <z_pred-z_meas> vs curvature;gen q/pT;<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprof_curv_{st}","",900,700)
    c2.SetLeftMargin(0.15)  
    if show:
        p.Draw("E")
    c2.Update()
    c2.SaveAs(f"{direc}/dz_profile_vs_curv_MB{st}.png")
    store_plots["canvas"][f"dzprof_vs_curv_MB{st}"]=c2
    store_plots["profiles"][f"dzprof_vs_curv_MB{st}"]=p
    print(f"histo name stored: deltaz_vs_curv_{st}")
    return c,c2,h,p

    
if __name__=="__main__":
    data=event_loop(1000, False, False)
