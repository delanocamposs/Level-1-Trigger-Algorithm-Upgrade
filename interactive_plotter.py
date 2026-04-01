import ROOT
import glob
from array import array 
import numpy as np 
from main import event_loop
from trigger_helpers import *
import os 

store_plots={"canvas":{}, "histos":{}, "fits":{}, "profiles":{}}


def clear_old_slice_pngs(station, xvar, yvar):
    direc=make_plot_dir(f"gaussianfit_slices_{yvar}_vs_{xvar}")
    pattern=os.path.join(direc, f"slice{yvar}_vs_{xvar}_bin*_station_{station}.png")
    for path in glob.glob(pattern):
        os.remove(path)

def plot_slice_and_gaussian(h2,xbin,st,xvar,xunit,yvar,yunit,min_entries=10,nsig=2,slice_xrange=None,slice_yrange=None):
    direc=make_plot_dir(f"gaussianfit_slices_{yvar}_vs_{xvar}")
    name=h2.GetName()
    xaxis=h2.GetXaxis()
    x_low=xaxis.GetBinLowEdge(xbin)
    x_high=xaxis.GetBinUpEdge(xbin)
    slice_name=f"{name}_slice_xbin{xbin}"
    h_slice=h2.ProjectionY(slice_name,xbin,xbin)
    h_slice.SetDirectory(0)
    c=ROOT.TCanvas(f"c_{name}_xbin{xbin}","",800,600)
    h_slice.SetTitle(f"Station {st}: {x_low:.1f}<{xvar}<{x_high:.1f} {xunit};{yvar} {yunit};Entries")
    h_slice.SetLineColor(ROOT.kBlack)
    h_slice.SetMarkerStyle(20)
    h_slice.SetMarkerSize(1.0)
    h_slice.SetStats(0)
    if h_slice.GetEntries()<min_entries:
        if slice_xrange is not None:
            h_slice.GetXaxis().SetRangeUser(slice_xrange[0],slice_xrange[1])
        if slice_yrange is not None:
            h_slice.GetYaxis().SetRangeUser(slice_yrange[0],slice_yrange[1])
        h_slice.Draw()
        c.SaveAs(f"{direc}/slice{yvar}_vs_{xvar}_bin{xbin}_station_{st}.png")
        c.Close()
        return c,h_slice,None,(None,None,None,None,1)
    mu0=h_slice.GetMean()
    rms0=h_slice.GetRMS()
    if rms0<=0:
        if slice_xrange is not None:
            h_slice.GetXaxis().SetRangeUser(slice_xrange[0],slice_xrange[1])
        if slice_yrange is not None:
            h_slice.GetYaxis().SetRangeUser(slice_yrange[0],slice_yrange[1])
        h_slice.Draw()
        c.SaveAs(f"{direc}/slice{yvar}_vs_{xvar}_bin{xbin}_station_{st}.png")
        c.Close()
        return c,h_slice,None,(None,None,None,None,2)
    y_min=mu0-nsig*rms0
    y_max=mu0+nsig*rms0
    g=ROOT.TF1(f"{name}_gaus_xbin{xbin}","gaus",y_min,y_max)
    A0=h_slice.GetMaximum()
    g.SetParameters(A0,mu0,max(rms0,1e-6))
    res=h_slice.Fit(g,"RQS")
    if slice_xrange is not None:
        h_slice.GetXaxis().SetRangeUser(slice_xrange[0],slice_xrange[1])
    if slice_yrange is not None:
        h_slice.GetYaxis().SetRangeUser(slice_yrange[0],slice_yrange[1])
    h_slice.Draw()
    g.SetLineColor(ROOT.kRed)
    g.SetLineWidth(2)
    g.Draw("SAME")
    mu=g.GetParameter(1)
    si=g.GetParameter(2)
    mu_e=g.GetParError(1)
    si_e=g.GetParError(2)
    leg=ROOT.TLegend(0.62,0.75,0.89,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(0,f"#sigma = {si:.3f} #pm {si_e:.2g} {yunit}","")
    leg.AddEntry(0,f"mean = {mu:.3f} #pm {mu_e:.2g} {yunit}","")
    leg.SetTextSize(0.035)
    leg.Draw()
    c.leg=leg
    c.SaveAs(f"{direc}/slice{yvar}_vs_{xvar}_bin{xbin}_station_{st}.png")
    c.Close()
    status=res.Status() if res else 999
    return c,h_slice,g,(mu,mu_e,si,si_e,status)

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

def plot_yslices(data,station,key,xvar,xunit,yvar,yunit,rebinx=20,min_entries=10,nsig=2,sigma_xrange=None,sigma_yrange=None,slice_xrange=None,slice_yrange=None):
    ROOT.gROOT.SetBatch(True)
    if key not in store_plots["histos"]:
        raise RuntimeError(f"missing histogram {key}. run function that makes histo first.")
    h0=store_plots["histos"][key]
    h=h0.Clone(f"{h0.GetName()}_tmp_yslices_st{station}")
    h.SetDirectory(0)
    if rebinx and rebinx>1:
        h.RebinX(rebinx)
    direc=make_plot_dir(f"yslices_{yvar}_vs_{xvar}")
    if station=="vtx":
        prop_info = "1 -> vtx"
    else:
        prop_info = f"{station+1} -> {station}"
    hx=ROOT.TH1F("h_sigma",f"{prop_info}: Gaussian fit #sigma vs curv ;{xvar} [{xunit}];#sigma({yvar}) {yunit}",h.GetNbinsX(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
    hx.SetDirectory(0)
    hx.GetXaxis().SetTitle(h.GetXaxis().GetTitle() if h.GetXaxis().GetTitle() else f"{xvar} {xunit}")
    clear_old_slice_pngs(station, xvar, yvar)
    for i in range(1,h.GetNbinsX()+1):
        c,h_slice,g,pars=plot_slice_and_gaussian(h,i,station,xvar,xunit,yvar,yunit,min_entries=min_entries,nsig=nsig,slice_xrange=slice_xrange,slice_yrange=slice_yrange)
        mu,mu_e,si,si_e,status=pars
        if si is None:
            continue
        if status!=0:
            continue
        if si<=0:
            continue
        hx.SetBinContent(i,si)
        hx.SetBinError(i,si_e if si_e is not None else 0.0)
    c_sig=ROOT.TCanvas("c_sigma","",900,700)
    c_sig.SetLeftMargin(0.15)
    hx.SetStats(0)
    hx.SetMarkerStyle(8)
    hx.SetMarkerSize(1.2)
    hx.SetLineWidth(2)
    hx.Draw("PE")
    if sigma_xrange is not None:
        hx.GetXaxis().SetRangeUser(sigma_xrange[0],sigma_xrange[1])
    if sigma_yrange is not None:
        hx.GetYaxis().SetRangeUser(sigma_yrange[0],sigma_yrange[1])
    xlo=hx.GetXaxis().GetXmin()
    xhi=hx.GetXaxis().GetXmax()
    if sigma_xrange is not None:
        xlo=sigma_xrange[0]
        xhi=sigma_xrange[1]
    f=ROOT.TF1(f"f_ms_{key}_st{station}","TMath::Sqrt([0]*[0]+([1]*[1])*(x*x))",xlo,xhi)
    a0=0.0
    for ib in range(1,hx.GetNbinsX()+1):
        if hx.GetBinContent(ib)>0:
            a0=hx.GetBinContent(ib)
            break
    b0=0.0
    for ib in range(hx.GetNbinsX(),0,-1):
        y=hx.GetBinContent(ib)
        x=hx.GetXaxis().GetBinCenter(ib)
        if y>0 and abs(x)>0:
            b0=y/abs(x)
            break
    f.SetParameters(a0,b0)
    f.SetParNames("a","b")
    f.SetLineColor(ROOT.kRed)
    f.SetLineWidth(2)
    fitres=hx.Fit(f,"RQS")
    f.Draw("SAME")
    a=f.GetParameter(0)
    b=f.GetParameter(1)
    ae=f.GetParError(0)
    be=f.GetParError(1)
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextSize(0.035)
    txt.DrawLatex(0.18,0.72,"#sigma(k)=#sqrt{a^{2}+b^{2}k^{2}}")
    txt.DrawLatex(0.18,0.82,f"a={a:.3g} #pm {ae:.2g}")
    txt.DrawLatex(0.18,0.77,f"b={b:.3g} #pm {be:.2g}")
    store_plots["fits"][f"yslices_sigma_fit_{key}_st{station}"]=f
    store_plots["fits"][f"yslices_sigma_fitres_{key}_st{station}"]=fitres
    c_sig.Update()
    c_sig.SaveAs(f"{direc}/sigma_{yvar}_vs_{xvar}_{station}.png")
    store_plots["canvas"][f"yslices_sigma_{yvar}_vs_{xvar}_{station}"]=c_sig
    store_plots["histos"][f"yslices_sigma_{yvar}_vs_{xvar}_{station}"]=hx
    print("all y slice plots saved.")
    return c_sig,hx

def plot_deltak_and_profile_vs_eta(data,st1,st2,xbins=60,ybins=80):
    direc=make_plot_dir("deltak_vs_eta")
    h=ROOT.TH2F(f"h_dk_eta_{st1}_{st2}",f"k_{st2}-k_{st1} vs stub #eta_{st1};stub #eta;#Delta k",xbins,-1.3,1.3,ybins,-1.0,1.0)
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

def plot_deltaz_vs_k(data,st,conv_z=False, conv_k=False, xrange=(-50000,50000), yrange=(-6000,6000),show=False,xbins=100,ybins=100):
    #st refers to the station we are propagating to, not from.
    #for ex: want to plot deltaz vs k for 2->1, we use st=1
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("dz_vs_k")
    h=ROOT.TH2F(f"h_dz_vs_k_MB{st+1}_to_MB{st}",f"z{st}-z{st+1} vs k{st+1};k{st+1};#Delta z (z{st}-z{st+1})",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    m1={}
    for muid,k2,z2 in zip(data["mu_id"][st+1],data["stub_k"][st+1],data["stub_z"][st+1]):
        if muid not in m1:
            m1[muid]=(k2,z2)
    for muid,z1 in zip(data["mu_id"][st],data["stub_z"][st]):
        if muid not in m1:
            continue
        k2,z2=m1[muid]
        h.Fill(k2,z1-z2)
    c1=ROOT.TCanvas(f"{direc}/c_dz_vs_k{st+1}_MB{st+1}_to_MB{st}","",900,700)
    c1.SetLeftMargin(0.15)
    h.SetStats(0)
    h.Draw("COLZ")
    c1.Update()
    c1.SaveAs(f"{direc}/dz_vs_k{st+1}_MB{st+1}_to_MB{st}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"<z{st}-z{st+1}> vs k{st+1};k{st+1};<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprofile_k{st+1}_MB{st+1}_to_MB{st}","",900,700)
    c2.SetLeftMargin(0.15)
    p.Draw("E")
    x_min=p.GetXaxis().GetXmin()
    x_max=p.GetXaxis().GetXmax()
    f_lin=ROOT.TF1(f"f_lin_MB{st+1}_to_MB{st}", "pol1", x_min, x_max)
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
    c2.SaveAs(f"{direc}/dz_profile_vs_k{st+1}_MB{st+1}_to_MB{st}.png")
    store_plots["canvas"][f"dz_vs_k{st+1}"]=c1
    store_plots["canvas"][f"dzprofile_vs_k{st+1}_MB{st+1}_to_MB{st}"]=c2
    store_plots["histos"][f"dz_vs_k{st+1}_MB{st+1}_to_MB{st}"]=h
    store_plots["profiles"][f"dzprof_vs_k{st+1}_MB{st+1}_to_MB{st}"]=p
    return c1,c2,h,p

def plot_deltaz_vs_k1_to_vtx(data,conv_z=False, conv_k=False, xrange=(-50000,50000), yrange=(-30000, 30000), show=False,xbins=100,ybins=100):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("dz_vs_k")
    h=ROOT.TH2F(f"h_dz_invk_st1_vtx",f"z_vtx-z1 vs k1;k1;#Delta z (z_vtx-z1)",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    m1={}
    for muid,k1,z1 in zip(data["mu_id"][1],data["stub_k"][1],data["stub_z"][1]):
        if muid not in m1:
            m1[muid]=(k1,z1)
    for muid,z_vtx in zip(data["mu_id"][1],data["gen_vz"][1]):
        if muid not in m1:
            continue
        k1,z1=m1[muid]
        if conv_z!=conv_k:
            continue
        if conv_z==False:
            z_vtx=z_vtx*(65536.0/1500.0)
        else:
            z_vtx=z_vtx
        h.Fill(k1,z_vtx-z1)
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
    c2.SaveAs(f"{direc}/dz_profile_vs_k1_MB1_vtx.png")
    store_plots["canvas"][f"dz_vs_k1_MB1_vtx"]=c1
    store_plots["canvas"][f"dzprof_vs_k1_MB1_vtx"]=c2
    store_plots["histos"][f"dz_vs_k1_MB1_vtx"]=h
    store_plots["profiles"][f"dzprof_vs_k1_MB1_vtx"]=p
    return c1,c2,h,p

def plot_deltak_vs_curv(data,st,show=False,xrange=(-7000,7000),yrange=(-20000,20000),xbins=100,ybins=100):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltak_vs_curv")
    c=ROOT.TCanvas(f"c_dk_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dk_curv_{st}",f"{st+1} -> {st}: (k_{{pred,{st}}}-k_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta k (k_{{pred,{st}}}-k_{{meas,{st}}})",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
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
    c.SetLeftMargin(0.18)
    h.Draw("COLZ")
    store_plots["canvas"][f"deltak_vs_curv_{st}"]=c
    store_plots["histos"][f"deltak_vs_curv_{st}"]=h
    print(f"histo name stored: deltak_vs_curv_{st}")
    c.Update()
    c.SaveAs(f"{direc}/deltak_vs_curv_{st}.png")
    f=ROOT.TFile(f"deltak_vs_curv_{st}.root","UPDATE")
    f.cd()
    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    c.Write(c.GetName(),ROOT.TObject.kOverwrite)
    f.Close()

    return c,h 

def plot_deltak_vs_curv_to_vtx(data,show=False,xrange=(-7000,7000),yrange=(-10000,10000), conv_k=False, conv_z=False,xbins=100,ybins=100):
    st="vtx"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltak_vs_curv")
    c=ROOT.TCanvas(f"c_dk_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dk_curv_{st}",f"1 -> vtx: (k_{{pred,{st}}}-k_{{gen,{st}}}) vs curvature;gen q/p_{{T}};#Delta k (k_{{pred,{st}}}-k_{{gen,{st}}})",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    m={}
    if conv_k!=conv_z:
        raise ValueError("conv_k and conv_z must be the same. units need to be consistent")
    for muid,curv1,slope1 in zip(data["mu_id"][1],data["gen_curv"][1],data["stub_k"][1]):
        if muid not in m:
            m[muid]=(curv1,slope1)
    for muid,gen_eta in zip(data["mu_id"][1],data["gen_eta"][1]):
        if muid not in m:
            continue
        curv1,slope1=m[muid]
        k_pred=slope1
        k_vtx=np.sinh(gen_eta)
        if conv_k==False:
            #if using digital units, convert back to digitized units since gen eta returns physical units from data dict
            k_vtx=k_vtx*(65536.0/2.0)
        else:
            k_vtx=k_vtx
        h.Fill(curv1,k_pred-k_vtx)
    h.SetStats(0)
    c.SetLeftMargin(0.18)
    py=h.ProjectionY(f"py_dk_{st}")
    mu=py.GetMean()
    rms=py.GetRMS()
    if rms>0:
        h.GetYaxis().SetRangeUser(mu-6*rms,mu+6*rms)
    h.Draw("COLZ")
    store_plots["canvas"][f"deltak_vs_curv_{st}"]=c
    store_plots["histos"][f"deltak_vs_curv_{st}"]=h
    c.Update()
    c.SaveAs(f"{direc}/deltak_vs_curv_{st}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"vtx <k_pred-k_gen> vs curvature (1 -> vtx);gen q/pT;<#Delta k>")
    c2=ROOT.TCanvas(f"c_dkprof_curv_{st}","",900,700)
    c2.SetLeftMargin(0.15)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"{direc}/dk_profile_vs_curv_{st}.png")
    store_plots["canvas"][f"dkprof_vs_curv_{st}"]=c2
    store_plots["profiles"][f"dkprof_vs_curv_{st}"]=p
    print(f"histo name stored: deltak_vs_curv_{st}")
    f=ROOT.TFile(f"deltak_vs_curv_{st}.root","UPDATE")
    f.cd()
    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    c.Write(c.GetName(),ROOT.TObject.kOverwrite)
    f.Close()
    return c,c2,h,p



def plot_deltaz_vs_curv(data,st,conv_k=False, conv_z=False,show=False,xrange=(-7000,7000),yrange=(-10000,10000),xbins=100,ybins=100):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltaz_vs_curv")
    c=ROOT.TCanvas(f"c_dz_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dz_curv_{st}",f"{st+1} -> {st}: (z_{{pred,{st}}}-z_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta z (z_{{pred,{st}}}-z_{{meas,{st}}})",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    m={}
    if conv_k != conv_z:
        raise ValueError("conv_z and conv_k must be the same for consistent units. ")
    if st==1:
        if conv_z==False:
            dR=-.109
        else:
            dR=-80.3
    elif st==2: 
        if conv_z==False:
            dR=-.142
        else:
            dR=-106
    else:
        print("cannot propagate to station 3 (no stubs at st 4).")
    for muid,z0,curv,slope in zip(data["mu_id"][st+1],data["stub_z"][st+1],data["gen_curv"][st+1],data["stub_k"][st+1]):
        if muid not in m:
            m[muid]=(z0,curv,slope)
    for muid,z_meas in zip(data["mu_id"][st],data["stub_z"][st]):
        if muid not in m:
            continue
        z0,curv,slope=m[muid]
        z_pred=z0+slope*dR
        h.Fill(curv,z_pred-z_meas)
    h.SetStats(0)
    c.SetLeftMargin(0.15)
    py=h.ProjectionY(f"py_dz_{st}")
    mu=py.GetMean()
    rms=py.GetRMS()
    if rms>0:
        h.GetYaxis().SetRangeUser(mu-6*rms,mu+6*rms)
    h.Draw("COLZ")
    store_plots["canvas"][f"deltaz_vs_curv_{st}"]=c
    store_plots["histos"][f"deltaz_vs_curv_{st}"]=h
    c.Update()
    c.SaveAs(f"{direc}/deltaz_vs_curv_{st}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"{st+1} -> {st}: <z_pred-z_meas> vs curvature;gen q/pT;<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprof_curv_{st}","",900,700)
    c2.SetLeftMargin(0.15)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"{direc}/dz_profile_vs_curv_MB{st}.png")
    store_plots["canvas"][f"dzprof_vs_curv_MB{st}"]=c2
    store_plots["profiles"][f"dzprof_vs_curv_MB{st}"]=p
    print(f"histo name stored: deltaz_vs_curv_{st}")
    f=ROOT.TFile(f"deltaz_vs_curv_{st}.root","UPDATE")
    f.cd()
    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    c.Write(c.GetName(),ROOT.TObject.kOverwrite)
    f.Close()
    return c,c2,h,p

def plot_deltaz_vs_curv_to_vtx(data,show=False,xrange=(-7000,7000),yrange=(-10000,10000), conv_k=False, conv_z=False,xbins=100,ybins=100):
    #currently i dont have the option to convert k and z to real coords. 
    #will have to obtain dR for the converted case (dR=-0.579 is for non-converted)
    st="vtx"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltaz_vs_curv")
    c=ROOT.TCanvas(f"c_dz_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dz_curv_{st}",f"1 -> vtx: (z_{{pred,{st}}}-z_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta z (z_{{pred,{st}}}-z_{{meas,{st}}})",xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    m={}
    if conv_k!=conv_z:
        raise ValueError("conv_k and conv_z must be the same. units need to be consistent")
    if conv_z==False:
        dR=-.601
    else:
        dR=-438
    for muid,z1,curv1,slope1 in zip(data["mu_id"][1],data["stub_z"][1],data["gen_curv"][1],data["stub_k"][1]):
        if muid not in m:
            m[muid]=(z1,curv1,slope1)
    for muid,z_vtx in zip(data["mu_id"][1],data["gen_vz"][1]):
        if muid not in m:
            continue
        z1,curv1,slope1=m[muid]
        z_pred=z1+slope1*dR
        if conv_z==False:
            #if using digital, raw units, i must change the z_vtx to raw units. from the data dictionary it is returned in physical units as cm from genparticles collection
            z_vtx=z_vtx*(65536.0/1500.0)
        else:
            z_vtx=z_vtx
        h.Fill(curv1,z_pred-z_vtx)
    h.SetStats(0)
    c.SetLeftMargin(0.15)
    py=h.ProjectionY(f"py_dz_{st}")
    mu=py.GetMean()
    rms=py.GetRMS()
    if rms>0:
        h.GetYaxis().SetRangeUser(mu-6*rms,mu+6*rms)
    h.Draw("COLZ")
    store_plots["canvas"][f"deltaz_vs_curv_{st}"]=c
    store_plots["histos"][f"deltaz_vs_curv_{st}"]=h
    c.Update()
    c.SaveAs(f"{direc}/deltaz_vs_curv_{st}.png")
    p=h.ProfileX()
    p.SetStats(0)
    p.SetLineWidth(2)
    p.SetTitle(f"{st+1} -> {st}: <z_pred-z_meas> vs curvature;gen q/pT;<#Delta z>")
    c2=ROOT.TCanvas(f"c_dzprof_curv_{st}","",900,700)
    c2.SetLeftMargin(0.15)
    p.Draw("E")
    c2.Update()
    c2.SaveAs(f"{direc}/dz_profile_vs_curv_MB{st}.png")
    store_plots["canvas"][f"dzprof_vs_curv_MB{st}"]=c2
    store_plots["profiles"][f"dzprof_vs_curv_MB{st}"]=p
    print(f"histo name stored: deltaz_vs_curv_{st}")
    f=ROOT.TFile(f"deltaz_vs_curv_{st}.root","UPDATE")
    f.cd()
    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    c.Write(c.GetName(),ROOT.TObject.kOverwrite)
    f.Close()
    return c,c2,h,p

def plot_eff_vs_pT(data,eff_data,show=False,pt_bins=(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70, 75,80,85,90,95,100,105,110,115,120), title="KMTF efficiency vs gen p_{T};gen p_{T} [GeV];Efficiency (prompt)"):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("kmtf_eff_vs_pt")
    bins=array('d',pt_bins)
    hDen=ROOT.TH1D("hDen",";gen p_{T} [GeV];Efficiency (matched to L1 Muon)",len(pt_bins)-1,bins)
    hNum=ROOT.TH1D("hNum",";gen p_{T} [GeV];Efficiency (matched to L1 Muon)",len(pt_bins)-1,bins)
    hDen.SetDirectory(0);hNum.SetDirectory(0)
    for pt in data["gen_pt_unmatched"]:
        hDen.Fill(float(pt))
    for pt in eff_data:
        hNum.Fill(float(pt))
    eff=ROOT.TEfficiency(hNum,hDen)
    c=ROOT.TCanvas("c_kmtf_eff","",800,600)
    eff.SetTitle(title)
    eff.Draw("AP")
    c.SaveAs(f"{direc}/kmtf_eff_vs_pt.png")
    f=ROOT.TFile("kmtf_eff_vs_pt.root","RECREATE")
    hDen.Write();hNum.Write();eff.Write();c.Write()
    f.Close()
    store_plots["canvas"]["kmtf_eff_vs_pt"]=c
    store_plots["histos"]["hDen_kmtf"]=hDen
    store_plots["histos"]["hNum_kmtf"]=hNum
    store_plots["histos"]["eff_kmtf"]=eff
    return c,eff

def plot_eff_vs_pT_prompt_displaced(data,prompt_data,displaced_data,show=False,pt_bins=(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,45,50,55,60,65,70,75,80,85,90,120,150,200),title="KMTF efficiency vs gen p_{T};gen p_{T} [GeV];Efficiency",color_prompt=ROOT.kRed,color_displaced=ROOT.kBlue):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("kmtf_eff_vs_pt")
    bins=array('d',pt_bins)
    hDen=ROOT.TH1D("hDen","; ; ",len(pt_bins)-1,bins)
    hNumP=ROOT.TH1D("hNumP","; ; ",len(pt_bins)-1,bins)
    hNumD=ROOT.TH1D("hNumD","; ; ",len(pt_bins)-1,bins)
    hDen.SetDirectory(0);hNumP.SetDirectory(0);hNumD.SetDirectory(0)
    for pt in data["gen_pt_unmatched"]:
        hDen.Fill(float(pt))
    for pt in prompt_data:
        hNumP.Fill(float(pt))
    for pt in displaced_data:
        hNumD.Fill(float(pt))
    effP=ROOT.TEfficiency(hNumP,hDen)
    effD=ROOT.TEfficiency(hNumD,hDen)
    effP.SetName("eff_prompt");effD.SetName("eff_displaced")
    effP.SetTitle(title)
    effP.SetLineColor(color_prompt);effP.SetMarkerColor(color_prompt);effP.SetMarkerStyle(20)
    effD.SetLineColor(color_displaced);effD.SetMarkerColor(color_displaced);effD.SetMarkerStyle(21)
    c=ROOT.TCanvas("c_kmtf_eff_pd","",800,600)
    effP.Draw("AP")
    effD.Draw("P SAME")
    leg=ROOT.TLegend(0.60,0.20,0.88,0.35)
    leg.SetBorderSize(0);leg.SetFillStyle(0)
    leg.AddEntry(effP,"prompt","lp")
    leg.AddEntry(effD,"displaced","lp")
    leg.Draw()
    c.SaveAs(f"{direc}/kmtf_eff_vs_pt_prompt_displaced.png")
    f=ROOT.TFile("kmtf_eff_vs_pt_prompt_displaced.root","RECREATE")
    hDen.Write();hNumP.Write();hNumD.Write();effP.Write();effD.Write();c.Write();f.Close()
    store_plots["canvas"]["kmtf_eff_vs_pt_prompt_displaced"]=c
    store_plots["histos"]["hDen_kmtf"]=hDen
    store_plots["histos"]["hNumP_kmtf"]=hNumP
    store_plots["histos"]["hNumD_kmtf"]=hNumD
    store_plots["histos"]["eff_prompt_kmtf"]=effP
    store_plots["histos"]["eff_displaced_kmtf"]=effD
    return c,effP,effD

def plot_eff_vs_eta(data, eff_data, show=False, n_bins=None, eta_bins=(-1.5,-1.4,-1.3,-1.2,-1.1, -1, -.9, -.8, -.7, -.6, -.5,-.4,-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.1, 1.2, 1.3, 1.4, 1.5), title="KMTF efficiency vs gen #eta;gen #eta [GeV];Efficiency (prompt)", y_range=(0, 1.1)):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("kmtf_eff_vs_eta")
    if n_bins is not None:
        import numpy as np
        eta_bins = tuple(np.linspace(eta_bins[0], eta_bins[-1], n_bins + 1))
    bins=array('d',eta_bins)
    hDen=ROOT.TH1D("hDen",";gen p_{T} [GeV];Efficiency (matched to L1 Muon)",len(eta_bins)-1,bins)
    hNum=ROOT.TH1D("hNum",";gen p_{T} [GeV];Efficiency (matched to L1 Muon)",len(eta_bins)-1,bins)
    hDen.SetDirectory(0);hNum.SetDirectory(0)
    for eta in data["gen_eta_unmatched"]:
        hDen.Fill(float(eta))
    for eta in eff_data:
        hNum.Fill(float(eta))
    eff=ROOT.TEfficiency(hNum,hDen)
    c=ROOT.TCanvas("c_kmtf_eff","",800,600)
    eff.SetTitle(title)
    eff.Draw("AP")
    ROOT.gPad.Update()
    eff.GetPaintedGraph().GetYaxis().SetRangeUser(y_range[0], y_range[1])
    c.SaveAs(f"{direc}/kmtf_eff_vs_eta.png")
    f=ROOT.TFile("kmtf_eff_vs_eta.root","RECREATE")
    hDen.Write();hNum.Write();eff.Write();c.Write()
    f.Close()
    store_plots["canvas"]["kmtf_eff_vs_eta"]=c
    store_plots["histos"]["hDen_kmtf_eta"]=hDen
    store_plots["histos"]["hNum_kmtf_eta"]=hNum
    store_plots["histos"]["eff_kmtf_eta"]=eff
    return c,eff

def plot_eff_vs_eta_prompt_displaced(data,prompt_data,displaced_data,show=False,n_bins=None,eta_bins=(-1.5,-1.4,-1.3,-1.2,-1.1, -1, -.9, -.8, -.7, -.6, -.5,-.4,-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.1, 1.2, 1.3, 1.4, 1.5),title="KMTF efficiency vs gen #eta;gen #eta [GeV];Efficiency",color_prompt=ROOT.kRed,color_displaced=ROOT.kBlue,y_range=(0,1.1)):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("kmtf_eff_vs_eta")
    if n_bins is not None:
        import numpy as np
        eta_bins = tuple(np.linspace(eta_bins[0], eta_bins[-1], n_bins + 1))
    bins=array('d',eta_bins)
    hDen=ROOT.TH1D("hDen","; ; ",len(eta_bins)-1,bins)
    hNumP=ROOT.TH1D("hNumP","; ; ",len(eta_bins)-1,bins)
    hNumD=ROOT.TH1D("hNumD","; ; ",len(eta_bins)-1,bins)
    hDen.SetDirectory(0);hNumP.SetDirectory(0);hNumD.SetDirectory(0)
    for eta in data["gen_eta_unmatched"]:
        hDen.Fill(float(eta))
    for eta in prompt_data:
        hNumP.Fill(float(eta))
    for eta in displaced_data:
        hNumD.Fill(float(eta))
    effP=ROOT.TEfficiency(hNumP,hDen)
    effD=ROOT.TEfficiency(hNumD,hDen)
    effP.SetName("eff_prompt");effD.SetName("eff_displaced")
    effP.SetTitle(title)
    effP.SetLineColor(color_prompt);effP.SetMarkerColor(color_prompt);effP.SetMarkerStyle(20)
    effD.SetLineColor(color_displaced);effD.SetMarkerColor(color_displaced);effD.SetMarkerStyle(21)
    c=ROOT.TCanvas("c_kmtf_eff_pd","",800,600)
    effP.Draw("AP")
    ROOT.gPad.Update()
    effP.GetPaintedGraph().GetYaxis().SetRangeUser(y_range[0], y_range[1])
    effD.Draw("P SAME")
    leg=ROOT.TLegend(0.60,0.20,0.88,0.35)
    leg.SetBorderSize(0);leg.SetFillStyle(0)
    leg.AddEntry(effP,"prompt","lp")
    leg.AddEntry(effD,"displaced","lp")
    leg.Draw()
    c.SaveAs(f"{direc}/kmtf_eff_vs_eta_prompt_displaced.png")
    f=ROOT.TFile("kmtf_eff_vs_eta_prompt_displaced.root","RECREATE")
    hDen.Write();hNumP.Write();hNumD.Write();effP.Write();effD.Write();c.Write();f.Close()
    store_plots["canvas"]["kmtf_eff_vs_eta_prompt_displaced"]=c
    store_plots["histos"]["hDen_kmtf"]=hDen
    store_plots["histos"]["hNumP_kmtf"]=hNumP
    store_plots["histos"]["hNumD_kmtf"]=hNumD
    store_plots["histos"]["eff_prompt_kmtf"]=effP
    store_plots["histos"]["eff_displaced_kmtf"]=effD
    return c,effP,effD

_zvtx_call_count = 0
def plot_overlay_prompt_displaced(prompt_data, displaced_data1,displaced_data2=None,displaced_data3=None, show=False, n_bins=50, xrange=(-1500, 1500),title="Normalized KMTF Track z_{vtx};z_{vtx} [cm];a.u.",color_prompt=ROOT.kRed, color_displaced=ROOT.kBlue):
    global _zvtx_call_count
    _zvtx_call_count += 1
    suf = _zvtx_call_count
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc = make_plot_dir("kmtf_zvtx_distribution")
    z_bins = tuple(np.linspace(xrange[0], xrange[1], n_bins + 1))
    bins = array('d', z_bins)
    hPrompt = ROOT.TH1D(f"hPrompt_{suf}", title, n_bins, bins)
    hDisp = ROOT.TH1D(f"hDisp_{suf}",   title, n_bins, bins)
    hPrompt.SetDirectory(0); hDisp.SetDirectory(0)
    for z in prompt_data["kmtf_zvtx"]:
        hPrompt.Fill(float(z))
    for z in displaced_data1["kmtf_zvtx"]:
        hDisp.Fill(float(z))
    if displaced_data2 is not None:
        for z in displaced_data2["kmtf_zvtx"]:
            hDisp.Fill(float(z))
    if displaced_data3 is not None:
        for z in displaced_data3["kmtf_zvtx"]:
            hDisp.Fill(float(z))
    hPrompt.SetLineColor(color_prompt);  hPrompt.SetMarkerColor(color_prompt);  hPrompt.SetMarkerStyle(20)
    hDisp.SetLineColor(color_displaced); hDisp.SetMarkerColor(color_displaced); hDisp.SetMarkerStyle(21)
    if hPrompt.Integral() > 0: hPrompt.Scale(1.0 / hPrompt.Integral())
    if hDisp.Integral()> 0: hDisp.Scale(1.0 / hDisp.Integral())
    hPrompt.SetStats(0)
    hDisp.SetStats(0)
    gP = ROOT.TF1(f"gP_{suf}", "gaus", xrange[0], xrange[1])
    gD = ROOT.TF1(f"gD_{suf}", "gaus", xrange[0], xrange[1])
    gP.SetLineColor(color_prompt);  gP.SetLineWidth(2); gP.SetLineStyle(2)
    gD.SetLineColor(color_displaced); gD.SetLineWidth(2); gD.SetLineStyle(2)
    hPrompt.Fit(gP, "RQ")
    hDisp.Fit(gD, "RQ")
    c = ROOT.TCanvas(f"c_kmtf_zvtx_{suf}", "", 800, 600)
    y_max = max(hPrompt.GetMaximum(), hDisp.GetMaximum()) * 1.4
    hPrompt.GetYaxis().SetRangeUser(0, y_max)
    hPrompt.Draw("HIST")
    hDisp.Draw("HIST SAME")
    #gP.Draw("SAME")
    #gD.Draw("SAME")
    leg = ROOT.TLegend(0.58, 0.62, 0.89, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextSize(0.030)
    leg.AddEntry(hPrompt, "prompt DY", "l")
    #leg.AddEntry(gP, f"#sigma={gP.GetParameter(2):.1f} cm", "l")
    leg.AddEntry(hDisp,"displaced", "l")
    #leg.AddEntry(gD, f"#sigma={gD.GetParameter(2):.1f} cm", "l")
    leg.Draw()
    print("sigma", gP.GetParameter(2))
    c.SaveAs(f"{direc}/kmtf_zvtx_prompt_displaced.png")
    f = ROOT.TFile(f"{direc}/kmtf_zvtx_prompt_displaced.root", "RECREATE")
    hPrompt.Write(); hDisp.Write(); gP.Write(); gD.Write(); c.Write(); f.Close()
    store_plots["canvas"][f"kmtf_zvtx_prompt_displaced_{suf}"] = c
    store_plots["histos"][f"hPrompt_zvtx_{suf}"] = hPrompt
    store_plots["histos"][f"hDisp_zvtx_{suf}"] = hDisp
    store_plots["histos"][f"leg_zvtx_{suf}"] = leg
    store_plots["fits"][f"gP_zvtx_{suf}"] = gP
    store_plots["fits"][f"gD_zvtx_{suf}"]= gD
    return c, hPrompt, hDisp, gP, gD

def plot_ztrack_minus_zgen(data,show=False,nbins=100,xrange=(-300,300),title="z_{track}-z_{gen};z_{track}-z_{gen} (cm);Entries",out_name="ztrack_minus_zgen_1d"):
    ztrack_vals=data["kmtf_zvtx_KMTFTracks_matched"]
    zgen_vals=data["gen_vz_KMTFTracks_matched"]
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(ztrack_vals)!=len(zgen_vals):
        raise ValueError("ztrack_vals and zgen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("ztrack_minus_zgen")
    c=ROOT.TCanvas(f"c_{out_name}","",800,600)
    h=ROOT.TH1F(f"h_{out_name}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    dz=np.array(ztrack_vals,dtype=float)-np.array(zgen_vals,dtype=float)
    for val in dz:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(2)
    h.Draw("HIST")
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c,h,dz
