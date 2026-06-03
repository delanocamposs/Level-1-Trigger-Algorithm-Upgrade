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
            dR=.109
        else:
            dR=80.3
    elif st==2: 
        if conv_z==False:
            dR=.142
        else:
            dR=106
    else:
        print("cannot propagate to station 3 (no stubs at st 4).")
    for muid,z0,curv,slope in zip(data["mu_id"][st+1],data["stub_z"][st+1],data["gen_curv"][st+1],data["stub_k"][st+1]):
        if muid not in m:
            m[muid]=(z0,curv,slope)
    for muid,z_meas in zip(data["mu_id"][st],data["stub_z"][st]):
        if muid not in m:
            continue
        z0,curv,slope=m[muid]
        z_pred=z0-slope*dR
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
def plot_overlay_prompt_displaced(prompt_data, displaced_data1, displaced_data2=None, displaced_data3=None,minbias_data=None, show=False, n_bins=50, xrange=(-1500, 1500),title="Normalized KMTF Track z_{vtx};z_{vtx} [cm];a.u.", color_prompt=ROOT.kRed, color_displaced=ROOT.kBlue, color_minbias=ROOT.kGreen+2):
    global _zvtx_call_count
    _zvtx_call_count += 1
    suf = _zvtx_call_count
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc = make_plot_dir("kmtf_zvtx_distribution")
    z_bins = tuple(np.linspace(xrange[0], xrange[1], n_bins + 1))
    bins = array('d', z_bins)

    hPrompt = ROOT.TH1D(f"hPrompt_{suf}", title, n_bins, bins)
    hDisp   = ROOT.TH1D(f"hDisp_{suf}",   title, n_bins, bins)
    hMB     = ROOT.TH1D(f"hMB_{suf}",     title, n_bins, bins)
    hPrompt.SetDirectory(0); hDisp.SetDirectory(0); hMB.SetDirectory(0)

    for z in prompt_data["kmtf_zvtx"]:
        hPrompt.Fill((1500/65536)*float(z))
    for z in displaced_data1["kmtf_zvtx"]:
        hDisp.Fill((1500/65536)*float(z))
    if displaced_data2 is not None:
        for z in displaced_data2["kmtf_zvtx"]:
            hDisp.Fill((1500/65536)*float(z))
    if displaced_data3 is not None:
        for z in displaced_data3["kmtf_zvtx"]:
            hDisp.Fill((1500/65536)*float(z))
    if minbias_data is not None:
        for z in minbias_data["kmtf_zvtx"]:
            hMB.Fill((1500/65536)*float(z))

    hPrompt.SetLineColor(color_prompt);   hPrompt.SetMarkerColor(color_prompt);   hPrompt.SetMarkerStyle(20); hPrompt.SetLineWidth(2)
    hDisp.SetLineColor(color_displaced);  hDisp.SetMarkerColor(color_displaced);  hDisp.SetMarkerStyle(21); hDisp.SetLineWidth(2)
    hMB.SetLineColor(color_minbias);      hMB.SetMarkerColor(color_minbias);      hMB.SetMarkerStyle(22); hMB.SetLineWidth(2)

    if hPrompt.Integral() > 0: hPrompt.Scale(1.0 / hPrompt.Integral())
    if hDisp.Integral()   > 0: hDisp.Scale(1.0 / hDisp.Integral())
    if hMB.Integral()     > 0: hMB.Scale(1.0 / hMB.Integral())

    hPrompt.SetStats(0); hDisp.SetStats(0); hMB.SetStats(0)

    gP = ROOT.TF1(f"gP_{suf}", "gaus", xrange[0], xrange[1])
    gD = ROOT.TF1(f"gD_{suf}", "gaus", xrange[0], xrange[1])
    gP.SetLineColor(color_prompt);    gP.SetLineWidth(2); gP.SetLineStyle(2)
    gD.SetLineColor(color_displaced); gD.SetLineWidth(2); gD.SetLineStyle(2)
    hPrompt.Fit(gP, "RQ")
    hDisp.Fit(gD, "RQ")

    c = ROOT.TCanvas(f"c_kmtf_zvtx_{suf}", "", 800, 600)

    histos_to_draw = [hPrompt, hDisp]
    if minbias_data is not None:
        histos_to_draw.append(hMB)
    y_max = max(h.GetMaximum() for h in histos_to_draw) * 1.4
    hPrompt.GetYaxis().SetRangeUser(0, y_max)

    hPrompt.Draw("HIST")
    hDisp.Draw("HIST SAME")
    if minbias_data is not None:
        hMB.Draw("HIST SAME")

    leg = ROOT.TLegend(0.58, 0.62, 0.89, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(hPrompt, "DY + Jets",             "l")
    leg.AddEntry(hDisp,   "Displaced Muon Gun",    "l")
    if minbias_data is not None:
        leg.AddEntry(hMB, "MinBias",               "l")
    leg.Draw()

    print("sigma prompt", gP.GetParameter(2))
    print("sigma displaced", gD.GetParameter(2))

    c.SaveAs(f"{direc}/kmtf_zvtx_prompt_displaced.png")
    f = ROOT.TFile(f"{direc}/kmtf_zvtx_prompt_displaced.root", "RECREATE")
    hPrompt.Write(); hDisp.Write(); hMB.Write(); gP.Write(); gD.Write(); c.Write(); f.Close()

    store_plots["canvas"][f"kmtf_zvtx_prompt_displaced_{suf}"] = c
    store_plots["histos"][f"hPrompt_zvtx_{suf}"] = hPrompt
    store_plots["histos"][f"hDisp_zvtx_{suf}"]   = hDisp
    store_plots["histos"][f"hMB_zvtx_{suf}"]     = hMB
    store_plots["histos"][f"leg_zvtx_{suf}"]      = leg
    store_plots["fits"][f"gP_zvtx_{suf}"]         = gP
    store_plots["fits"][f"gD_zvtx_{suf}"]         = gD

    return c, hPrompt, hDisp, hMB, gP, gD

_simple_plot_call_count = 0
def plot_ztrack_minus_zgen(data,show=False,showFit=True,fitColor=ROOT.kRed,nbins=100,xrange=(-300,300),title="z_{track}-z_{gen};z_{track}-z_{gen} (cm);Entries",out_name="ztrack_minus_zgen_1d"):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    ztrack_vals=(1500/65536)*np.array(data["kmtf_zvtx_KMTFTracks_matched"], dtype=float)
    zgen_vals=np.array(data["gen_vz_KMTFTracks_matched"], dtype=float)
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(ztrack_vals)!=len(zgen_vals):
        raise ValueError("ztrack_vals and zgen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("ztrack_minus_zgen")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    dz=ztrack_vals-zgen_vals
    for val in dz:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(2)
    h.Draw("HIST")
    fit=None
    leg=None
    if showFit:
        fit=ROOT.TF1(f"fit_{uniq}","gaus",xrange[0],xrange[1])
        fit.SetLineColor(fitColor)
        fit.SetLineWidth(2)
        h.Fit(fit,"RQS")
        fit.Draw("SAME")
        sigma=fit.GetParameter(2)
        sigma_err=fit.GetParError(2)
        leg=ROOT.TLegend(0.6,0.80,0.86,0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"Gaussian fit", "L")
        leg.AddEntry(0,f"#sigma = {sigma:.3g} #pm {sigma_err:.2g} cm","")
        leg.Draw()

    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    if fit is not None:
        store_plots["fits"][out_name]=fit
    if leg is not None:
        store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c,h,dz

def plot_etatrack_minus_etagen(data,show=False,showFit=True,fitColor=ROOT.kRed,nbins=120,xrange=(-1.2,1.2),title="#eta_{track}-#eta_{gen};#eta_{track}-#eta_{gen};Entries",out_name="etatrack_minus_etagen_1d",overlay_out_name=None,overlay_title="#eta_{track} and #eta_{gen} distributions;#eta;Entries"):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    etatrack_vals=data["kmtf_eta_KMTFTracks_matched"]
    etagen_vals=data["gen_eta_KMTFTracks_allmatched"]
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(etatrack_vals)!=len(etagen_vals):
        raise ValueError("etatrack_vals and etagen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("etatrack_minus_etagen")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    deta=np.array(etatrack_vals,dtype=float)-np.array(etagen_vals,dtype=float)
    for val in deta:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(2)
    h.Draw("HIST")
    fit=None
    leg=None
    if showFit:
        fit=ROOT.TF1(f"fit_{uniq}","gaus",xrange[0],xrange[1])
        fit.SetLineColor(fitColor)
        fit.SetLineWidth(2)
        h.Fit(fit,"RQS")
        fit.Draw("SAME")
        sigma=fit.GetParameter(2)
        sigma_err=fit.GetParError(2)
        leg=ROOT.TLegend(0.6,0.80,0.86,0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"Gaussian fit", "L")
        leg.AddEntry(0,f"#sigma = {sigma:.3g} #pm {sigma_err:.2g}","")
        leg.Draw()
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    if fit is not None:
        store_plots["fits"][out_name]=fit
    if leg is not None:
        store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")

    if overlay_out_name is None:
        overlay_out_name=f"{out_name}_track_vs_gen"
    uniq_overlay=f"{overlay_out_name}_{_simple_plot_call_count}"
    c_overlay=ROOT.TCanvas(f"c_{uniq_overlay}","",800,600)
    h_etatrack=ROOT.TH1F(f"h_etatrack_{uniq_overlay}",overlay_title,nbins,xrange[0],xrange[1])
    h_etagen=ROOT.TH1F(f"h_etagen_{uniq_overlay}",overlay_title,nbins,xrange[0],xrange[1])
    h_etatrack.SetDirectory(0)
    h_etagen.SetDirectory(0)
    for val in etatrack_vals:
        h_etatrack.Fill(float(val))
    for val in etagen_vals:
        h_etagen.Fill(float(val))
    h_etatrack.SetStats(0)
    h_etagen.SetStats(0)
    h_etatrack.SetLineWidth(2)
    h_etagen.SetLineWidth(2)
    h_etatrack.SetLineColor(ROOT.kBlue+1)
    h_etagen.SetLineColor(ROOT.kOrange+7)
    max_y=max(h_etatrack.GetMaximum(),h_etagen.GetMaximum())
    h_etatrack.SetMaximum(1.15*max_y if max_y>0 else 1.0)
    h_etatrack.Draw("HIST")
    h_etagen.Draw("HIST SAME")
    leg_overlay=ROOT.TLegend(0.62,0.78,0.88,0.88)
    leg_overlay.SetBorderSize(0)
    leg_overlay.SetFillStyle(0)
    leg_overlay.AddEntry(h_etatrack,"#eta_{track}","l")
    leg_overlay.AddEntry(h_etagen,"#eta_{gen}","l")
    leg_overlay.Draw()
    store_plots["canvas"][overlay_out_name]=c_overlay
    store_plots["histos"][f"{overlay_out_name}_etatrack"]=h_etatrack
    store_plots["histos"][f"{overlay_out_name}_etagen"]=h_etagen
    store_plots["histos"][f"{overlay_out_name}_legend"]=leg_overlay
    c_overlay.Update()
    c_overlay.SaveAs(f"{direc}/{overlay_out_name}.png")
    return c,h,deta

def plot_etagen_vs_kslope(data, k_phys, eta_phys,show=False,xbins=120,ybins=120,xrange=(-65536,65536),yrange=(-1600,1600),title="#eta_{gen} vs kSlope_{track};kSlope_{track};#eta_{gen}",profile_title="profile: #eta_{gen} vs kSlope_{track};kSlope_{track};< #eta_{gen} >",out_name="etagen_vs_kslope_2d", ytitle_offset_2d=1.05, ytitle_offset_profile=1.35):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    kslope_vals=data["kmtf_kslope_KMTFTracks_matched"]
    etagen_vals=data["gen_eta_KMTFTracks_allmatched"]
    if len(kslope_vals)!=len(etagen_vals):
        raise ValueError("kslope_vals and etagen_vals must have the same length for event-by-event plotting.")
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("etagen_vs_kslope")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",850,700)
    c.SetLeftMargin(0.14)
    h=ROOT.TH2F(f"h_{uniq}",title,xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    h.GetYaxis().SetTitleOffset(ytitle_offset_2d)
    for ks,eta in zip(kslope_vals,etagen_vals):
        if k_phys:
            k_LSB = 2.0/65536.0
            ks = ks*k_LSB
        if not eta_phys:
            eta_LSB = (2*np.pi)/(2**13)
            eta =eta/eta_LSB
        h.Fill(float(ks),float(eta))
    h.SetStats(0)
    h.Draw("COLZ")
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    store_plots["canvas"][uniq]=c
    store_plots["histos"][uniq]=h
    c.Update()
    c.SaveAs(f"{direc}/{uniq}.png")

    p=h.ProfileX(f"p_{uniq}")
    p.SetDirectory(0)
    p.SetTitle(profile_title)
    p.SetStats(0)
    p.SetLineWidth(2)
    c_prof=ROOT.TCanvas(f"cprof_{uniq}","",850,700)
    c_prof.SetLeftMargin(0.14)
    p.GetYaxis().SetTitleOffset(ytitle_offset_profile)
    p.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    p.Draw("E")

    f_lin=ROOT.TF1(f"f_lin_{uniq}", "pol1", xrange[0], xrange[1])
    f_lin.SetLineColor(ROOT.kRed+1)
    f_lin.SetLineWidth(2)
    fitres=p.Fit(f_lin, "RQS")
    f_lin.Draw("SAME")
    p0=f_lin.GetParameter(0)
    p1=f_lin.GetParameter(1)
    p0e=f_lin.GetParError(0)
    p1e=f_lin.GetParError(1)

    leg=ROOT.TLegend(0.18,0.74,0.64,0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.AddEntry(f_lin,"< #eta_{gen} > = p0 + p1 #upoint kSlope_{track}","l")
    leg.AddEntry(0,f"p0 = {p0:.4g} #pm {p0e:.2g}","")
    leg.AddEntry(0,f"p1 = {p1:.4g} #pm {p1e:.2g}","")
    leg.Draw()

    c_prof.Update()
    c_prof.SaveAs(f"{direc}/{uniq}_profile.png")
    store_plots["canvas"][f"{out_name}_profile"]=c_prof
    store_plots["profiles"][f"{out_name}_profile"]=p
    store_plots["fits"][f"{out_name}_profile_fit"]=f_lin
    store_plots["fits"][f"{out_name}_profile_fitres"]=fitres
    store_plots["histos"][f"{out_name}_profile_legend"]=leg
    store_plots["canvas"][f"{uniq}_profile"]=c_prof
    store_plots["profiles"][f"{uniq}_profile"]=p
    store_plots["fits"][f"{uniq}_profile_fit"]=f_lin
    store_plots["fits"][f"{uniq}_profile_fitres"]=fitres
    store_plots["histos"][f"{uniq}_profile_legend"]=leg
    return c,h,c_prof,p,f_lin,fitres

def plot_dz_overlay(data1, data2,show=False,showFit=True, fitColor1=ROOT.kRed-9, fitColor2 = ROOT.kBlue - 10, nbins=100,xrange=(-200,200),title="z_{track}-z_{gen};z_{track}-z_{gen} (cm);Entries",out_name="ztrack_minus_zgen_1d", h1color = ROOT.kRed - 3, h2color = ROOT.kBlue - 7):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    ztrack_vals=(1500/65536)*np.array(data1["kmtf_zvtx_KMTFTracks_matched"], dtype=float)
    zgen_vals=np.array(data1["gen_vz_KMTFTracks_matched"], dtype=float)
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(ztrack_vals)!=len(zgen_vals):
        raise ValueError("ztrack_vals and zgen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("ztrack_minus_zgen")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    dz=ztrack_vals-zgen_vals
    for val in dz:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(1)
    h.SetLineColor(h1color)
    h.Draw("HIST")
    fit=None
    leg=ROOT.TLegend(0.6,0.70,0.86,0.88)
    leg.AddEntry(h, f"Offline KMTF", "L")
    if showFit:
        fit=ROOT.TF1(f"fit_{uniq}","gaus",xrange[0],xrange[1])
        fit.SetLineColor(fitColor1)
        fit.SetLineWidth(2)
        fit.SetNpx(1000)   
        h.Fit(fit,"RQS")
        fit.Draw("SAME")
        sigma=fit.GetParameter(2)
        sigma_err=fit.GetParError(2)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"fit #sigma = {sigma:.3g} #pm {sigma_err:.2g} cm", "")
    ztrack_vals1=(1500/65536)*np.array(data2["kmtf_zvtx_KMTFTracks_matched"], dtype=float)
    zgen_vals1=np.array(data2["gen_vz_KMTFTracks_matched"], dtype=float)
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(ztrack_vals1)!=len(zgen_vals1):
        raise ValueError("ztrack_vals and zgen_vals must have the same length for event-by-event subtraction.")
    h2=ROOT.TH1F(f"h_{uniq}_1",title,nbins,xrange[0],xrange[1])
    h2.SetDirectory(0)
    dz1=ztrack_vals1-zgen_vals1
    for val in dz1:
        h2.Fill(float(val))
    h2.SetStats(0)
    h2.SetLineWidth(1)
    h2.SetLineColor(h2color)
    leg.AddEntry(h2, f"LUT-based KMTF", "L")
    h2.Draw("HIST, SAME")
    if showFit:
        fit1=ROOT.TF1(f"fit_{uniq}_1","gaus",xrange[0],xrange[1])
        fit1.SetLineColor(fitColor2)
        fit1.SetLineWidth(2)
        fit1.SetNpx(1000)   
        h2.Fit(fit1,"RQS")
        fit1.Draw("SAME")
        sigma1=fit1.GetParameter(2)
        sigma_err1=fit1.GetParError(2)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"fit #sigma = {sigma1:.3g} #pm {sigma_err1:.2g} cm", "")
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    if fit is not None:
        store_plots["fits"][out_name]=fit
    if leg is not None:
        store_plots["histos"][f"{out_name}_legend"]=leg
    leg.Draw()
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c


def plot_etatrack_minus_etagen_overlay(data1, data2, show=False, showFit=True, fitColor1=ROOT.kRed-9, fitColor2=ROOT.kBlue-10, nbins=120, xrange=(-1.2,1.2), title="#eta_{track}-#eta_{gen};#eta_{track}-#eta_{gen};Entries", out_name="etatrack_minus_etagen_1d_overlay", h1color = ROOT.kRed - 3, h2color = ROOT.kBlue -7):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    etatrack_vals=np.array(data1["kmtf_eta_KMTFTracks_matched"], dtype=float)
    etagen_vals=np.array(data1["gen_eta_KMTFTracks_allmatched"], dtype=float)
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(etatrack_vals)!=len(etagen_vals):
        raise ValueError("etatrack_vals and etagen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("etatrack_minus_etagen")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    deta=etatrack_vals-etagen_vals
    for val in deta:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(1)
    h.SetLineColor(h1color)
    h.Draw("HIST")
    leg=ROOT.TLegend(0.55,0.70,0.86,0.88)
    leg.AddEntry(h, f"Offline KMTF", "L")
    if showFit:
        fit=ROOT.TF1(f"fit_{uniq}","gaus",xrange[0],xrange[1])
        fit.SetLineColor(fitColor1)
        fit.SetLineWidth(2)
        fit.SetNpx(1000)
        h.Fit(fit,"RQS")
        fit.Draw("SAME")
        sigma=fit.GetParameter(2)
        sigma_err=fit.GetParError(2)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"fit #sigma = {sigma:.3g} #pm {sigma_err:.2g}", "")
    etatrack_vals1=np.array(data2["kmtf_eta_KMTFTracks_matched"], dtype=float)
    etagen_vals1=np.array(data2["gen_eta_KMTFTracks_allmatched"], dtype=float)
    if len(etatrack_vals1)!=len(etagen_vals1):
        raise ValueError("etatrack_vals and etagen_vals must have the same length for event-by-event subtraction.")
    h2=ROOT.TH1F(f"h_{uniq}_1",title,nbins,xrange[0],xrange[1])
    h2.SetDirectory(0)
    deta1=etatrack_vals1-etagen_vals1
    for val in deta1:
        h2.Fill(float(val))
    h2.SetStats(0)
    h2.SetLineWidth(1)
    h2.SetLineColor(h2color)
    max_y=max(h.GetMaximum(), h2.GetMaximum())
    h.SetMaximum(1.15*max_y if max_y>0 else 1.0)
    leg.AddEntry(h2, f"LUT-based KMTF", "L")
    h2.Draw("HIST SAME")
    if showFit:
        fit1=ROOT.TF1(f"fit_{uniq}_1","gaus",xrange[0],xrange[1])
        fit1.SetLineColor(fitColor2)
        fit1.SetLineWidth(2)
        fit1.SetNpx(1000)
        h2.Fit(fit1,"RQS")
        fit1.Draw("SAME")
        sigma1=fit1.GetParameter(2)
        sigma_err1=fit1.GetParError(2)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"fit #sigma = {sigma1:.3g} #pm {sigma_err1:.2g}", "")
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    store_plots["histos"][f"{out_name}_1"]=h2
    if fit and fit1 is not False:
        store_plots["fits"][out_name]=fit
        store_plots["fits"][f"{out_name}_1"]=fit
        store_plots["fits"][f"{out_name}_2"]=fit1
    if leg is not None:
        store_plots["histos"][f"{out_name}_legend"]=leg
        leg.Draw()
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c

def plot_zSAMuon_minus_zgen(data,vertex="prompt",show=False,showFit=True,fitColor=ROOT.kRed,nbins=100,xrange=(-300,300),title=None,out_name=None,overlay_out_name=None,overlay_title=None,overlay_xrange=None,overlay_nbins=None):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    #phZ0 is already in cm, so (unlike the KMTF track zPosition) no hw->cm conversion is needed.
    zSAMuon_vals=np.array(data[f"samuon_z_{vertex}_matched"], dtype=float)
    zgen_vals=np.array(data[f"gen_vz_SAMuons_{vertex}_matched"], dtype=float)
    if title is None:
        title=f"z_{{SAMuon}}-z_{{gen}} ({vertex});z_{{SAMuon}}-z_{{gen}} (cm);Entries"
    if out_name is None:
        out_name=f"zSAMuon_minus_zgen_{vertex}_1d"
    if overlay_title is None:
        overlay_title=f"z_{{SAMuon}} and z_{{gen}} distributions ({vertex});z (cm);Entries"
    if overlay_xrange is None:
        overlay_xrange=xrange
    if overlay_nbins is None:
        overlay_nbins=nbins
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(zSAMuon_vals)!=len(zgen_vals):
        raise ValueError("zSAMuon_vals and zgen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("zSAMuon_minus_zgen")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    dz=4*zSAMuon_vals-zgen_vals
    for val in dz:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(2)
    h.Draw("HIST")
    fit=None
    leg=None
    if showFit:
        fit=ROOT.TF1(f"fit_{uniq}","gaus",xrange[0],xrange[1])
        fit.SetLineColor(fitColor)
        fit.SetLineWidth(2)
        h.Fit(fit,"RQS")
        fit.Draw("SAME")
        sigma=fit.GetParameter(2)
        sigma_err=fit.GetParError(2)
        leg=ROOT.TLegend(0.6,0.80,0.86,0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"Gaussian fit", "L")
        leg.AddEntry(0,f"#sigma = {sigma:.3g} #pm {sigma_err:.2g} cm","")
        leg.Draw()

    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    if fit is not None:
        store_plots["fits"][out_name]=fit
    if leg is not None:
        store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")

    if overlay_out_name is None:
        overlay_out_name=f"{out_name}_SAMuon_vs_gen"
    uniq_overlay=f"{overlay_out_name}_{_simple_plot_call_count}"
    c_overlay=ROOT.TCanvas(f"c_{uniq_overlay}","",800,600)
    h_zSAMuon=ROOT.TH1F(f"h_zSAMuon_{uniq_overlay}",overlay_title,overlay_nbins,overlay_xrange[0],overlay_xrange[1])
    h_zgen=ROOT.TH1F(f"h_zgen_{uniq_overlay}",overlay_title,overlay_nbins,overlay_xrange[0],overlay_xrange[1])
    h_zSAMuon.SetDirectory(0)
    h_zgen.SetDirectory(0)
    for val in zSAMuon_vals:
        h_zSAMuon.Fill(float(16*val))
    for val in zgen_vals:
        h_zgen.Fill(float(val))
    h_zSAMuon.SetStats(0)
    h_zgen.SetStats(0)
    h_zSAMuon.SetLineWidth(2)
    h_zgen.SetLineWidth(2)
    h_zSAMuon.SetLineColor(ROOT.kBlue+1)
    h_zgen.SetLineColor(ROOT.kOrange+7)
    max_y=max(h_zSAMuon.GetMaximum(),h_zgen.GetMaximum())
    h_zSAMuon.SetMaximum(1.15*max_y if max_y>0 else 1.0)
    h_zSAMuon.Draw("HIST")
    h_zgen.Draw("HIST SAME")
    leg_overlay=ROOT.TLegend(0.62,0.78,0.88,0.88)
    leg_overlay.SetBorderSize(0)
    leg_overlay.SetFillStyle(0)
    leg_overlay.AddEntry(h_zSAMuon,"z_{SAMuon}","l")
    leg_overlay.AddEntry(h_zgen,"z_{gen}","l")
    leg_overlay.Draw()
    store_plots["canvas"][overlay_out_name]=c_overlay
    store_plots["histos"][f"{overlay_out_name}_zSAMuon"]=h_zSAMuon
    store_plots["histos"][f"{overlay_out_name}_zgen"]=h_zgen
    store_plots["histos"][f"{overlay_out_name}_legend"]=leg_overlay
    c_overlay.Update()
    c_overlay.SaveAs(f"{direc}/{overlay_out_name}.png")
    return c,h,dz


def plot_etaSAMuon_minus_etagen(data,vertex="prompt",show=False,showFit=True,fitColor=ROOT.kRed,nbins=120,xrange=(-1.2,1.2),title=None,out_name=None,overlay_out_name=None,overlay_title=None):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    etaSAMuon_vals=data[f"samuon_eta_{vertex}_matched"]
    etagen_vals=data[f"gen_eta_SAMuons_{vertex}_allmatched"]
    if title is None:
        title=f"#eta_{{SAMuon}}-#eta_{{gen}} ({vertex});#eta_{{SAMuon}}-#eta_{{gen}};Entries"
    if out_name is None:
        out_name=f"etaSAMuon_minus_etagen_{vertex}_1d"
    if overlay_title is None:
        overlay_title=f"#eta_{{SAMuon}} and #eta_{{gen}} distributions ({vertex});#eta;Entries"
    if not show:
        ROOT.gROOT.SetBatch(True)
    if len(etaSAMuon_vals)!=len(etagen_vals):
        raise ValueError("etaSAMuon_vals and etagen_vals must have the same length for event-by-event subtraction.")
    direc=make_plot_dir("etaSAMuon_minus_etagen")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    deta=np.array(etaSAMuon_vals,dtype=float)-np.array(etagen_vals,dtype=float)
    for val in deta:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(2)
    h.Draw("HIST")
    fit=None
    leg=None
    if showFit:
        fit=ROOT.TF1(f"fit_{uniq}","gaus",xrange[0],xrange[1])
        fit.SetLineColor(fitColor)
        fit.SetLineWidth(2)
        h.Fit(fit,"RQS")
        fit.Draw("SAME")
        sigma=fit.GetParameter(2)
        sigma_err=fit.GetParError(2)
        leg=ROOT.TLegend(0.6,0.80,0.86,0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(fit, f"Gaussian fit", "L")
        leg.AddEntry(0,f"#sigma = {sigma:.3g} #pm {sigma_err:.2g}","")
        leg.Draw()
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    if fit is not None:
        store_plots["fits"][out_name]=fit
    if leg is not None:
        store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")

    if overlay_out_name is None:
        overlay_out_name=f"{out_name}_SAMuon_vs_gen"
    uniq_overlay=f"{overlay_out_name}_{_simple_plot_call_count}"
    c_overlay=ROOT.TCanvas(f"c_{uniq_overlay}","",800,600)
    h_etaSAMuon=ROOT.TH1F(f"h_etaSAMuon_{uniq_overlay}",overlay_title,nbins,xrange[0],xrange[1])
    h_etagen=ROOT.TH1F(f"h_etagen_{uniq_overlay}",overlay_title,nbins,xrange[0],xrange[1])
    h_etaSAMuon.SetDirectory(0)
    h_etagen.SetDirectory(0)
    for val in etaSAMuon_vals:
        h_etaSAMuon.Fill(float(val))
    for val in etagen_vals:
        h_etagen.Fill(float(val))
    h_etaSAMuon.SetStats(0)
    h_etagen.SetStats(0)
    h_etaSAMuon.SetLineWidth(2)
    h_etagen.SetLineWidth(2)
    h_etaSAMuon.SetLineColor(ROOT.kBlue+1)
    h_etagen.SetLineColor(ROOT.kOrange+7)
    max_y=max(h_etaSAMuon.GetMaximum(),h_etagen.GetMaximum())
    h_etaSAMuon.SetMaximum(1.15*max_y if max_y>0 else 1.0)
    h_etaSAMuon.Draw("HIST")
    h_etagen.Draw("HIST SAME")
    leg_overlay=ROOT.TLegend(0.62,0.78,0.88,0.88)
    leg_overlay.SetBorderSize(0)
    leg_overlay.SetFillStyle(0)
    leg_overlay.AddEntry(h_etaSAMuon,"#eta_{SAMuon}","l")
    leg_overlay.AddEntry(h_etagen,"#eta_{gen}","l")
    leg_overlay.Draw()
    store_plots["canvas"][overlay_out_name]=c_overlay
    store_plots["histos"][f"{overlay_out_name}_etaSAMuon"]=h_etaSAMuon
    store_plots["histos"][f"{overlay_out_name}_etagen"]=h_etagen
    store_plots["histos"][f"{overlay_out_name}_legend"]=leg_overlay
    c_overlay.Update()
    c_overlay.SaveAs(f"{direc}/{overlay_out_name}.png")
    return c,h,deta

def plot_railZ0_eta(data,vertex="prompt",rails=(-16,15),show=False,nbins=120,xrange=(-1.2,1.2),title=None,out_name=None):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    eta_vals=np.array(data[f"samuon_eta_{vertex}_matched"],dtype=float)
    hwZ0_vals=np.array(data[f"samuon_hwZ0_{vertex}_matched"],dtype=float)
    if len(eta_vals)!=len(hwZ0_vals):
        raise ValueError("samuon_eta and samuon_hwZ0 arrays must be the same length.")
    rail_set=set(int(r) for r in rails)
    mask=np.array([int(round(z)) in rail_set for z in hwZ0_vals],dtype=bool)
    rail_eta=eta_vals[mask]
    n_total=len(eta_vals)
    n_rail=len(rail_eta)
    if title is None:
        rails_str="/".join(str(int(r)) for r in rails)
        title=f"#eta of SAMuons on z0 rail (hwZ0={rails_str}, {vertex});#eta_{{SAMuon}};Entries"
    if out_name is None:
        out_name=f"railZ0_eta_{vertex}_1d"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("railZ0_eta")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    h=ROOT.TH1F(f"h_{uniq}",title,nbins,xrange[0],xrange[1])
    h.SetDirectory(0)
    for val in rail_eta:
        h.Fill(float(val))
    h.SetStats(0)
    h.SetLineWidth(2)
    h.SetLineColor(ROOT.kRed+1)
    h.Draw("HIST")
    frac=(n_rail/n_total) if n_total>0 else 0.0
    leg=ROOT.TLegend(0.5,0.78,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.032)
    leg.AddEntry(0,f"rail muons: {n_rail}/{n_total} ({100*frac:.1f}%)","")
    leg.Draw()
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c,h,rail_eta

def count_z0_rail(data,vertex,rails=(-16,15)):
    hwZ0_vals=np.array(data[f"samuon_hwZ0_{vertex}_matched"],dtype=float)
    rail_set=set(int(r) for r in rails)
    mask=np.array([int(round(z)) in rail_set for z in hwZ0_vals],dtype=bool)
    return mask,int(mask.sum()),len(hwZ0_vals)

def plot_overlay_prompt_displaced_SAMuon(prompt_data, displaced_data1, displaced_data2=None, displaced_data3=None, minbias_data=None, vertex="prompt", show=False, n_bins=50, xrange=(-30, 30), z_scale=1.0, title=None, color_prompt=ROOT.kRed, color_displaced=ROOT.kBlue, color_minbias=ROOT.kGreen+2, out_name=None):
    global _zvtx_call_count
    _zvtx_call_count += 1
    suf = _zvtx_call_count
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc = make_plot_dir("samuon_z0_distribution")
    if title is None:
        title = f"Normalized SAMuon z_{{0}} ({vertex});z_{{0}} [cm];a.u."
    if out_name is None:
        out_name = f"samuon_z0_prompt_displaced_{vertex}"
    key = f"samuon_z_all_{vertex}"
    z_bins = tuple(np.linspace(xrange[0], xrange[1], n_bins + 1))
    bins = array('d', z_bins)

    hPrompt = ROOT.TH1D(f"hPrompt_SA_{suf}", title, n_bins, bins)
    hDisp   = ROOT.TH1D(f"hDisp_SA_{suf}",   title, n_bins, bins)
    hMB     = ROOT.TH1D(f"hMB_SA_{suf}",     title, n_bins, bins)
    hPrompt.SetDirectory(0); hDisp.SetDirectory(0); hMB.SetDirectory(0)

    for z in prompt_data[key]:
        hPrompt.Fill(z_scale*float(z))
    for z in displaced_data1[key]:
        hDisp.Fill(z_scale*float(z))
    if displaced_data2 is not None:
        for z in displaced_data2[key]:
            hDisp.Fill(z_scale*float(z))
    if displaced_data3 is not None:
        for z in displaced_data3[key]:
            hDisp.Fill(z_scale*float(z))
    if minbias_data is not None:
        for z in minbias_data[key]:
            hMB.Fill(z_scale*float(z))

    hPrompt.SetLineColor(color_prompt);   hPrompt.SetMarkerColor(color_prompt);   hPrompt.SetMarkerStyle(20); hPrompt.SetLineWidth(2)
    hDisp.SetLineColor(color_displaced);  hDisp.SetMarkerColor(color_displaced);  hDisp.SetMarkerStyle(21); hDisp.SetLineWidth(2)
    hMB.SetLineColor(color_minbias);      hMB.SetMarkerColor(color_minbias);      hMB.SetMarkerStyle(22); hMB.SetLineWidth(2)

    if hPrompt.Integral() > 0: hPrompt.Scale(1.0 / hPrompt.Integral())
    if hDisp.Integral()   > 0: hDisp.Scale(1.0 / hDisp.Integral())
    if hMB.Integral()     > 0: hMB.Scale(1.0 / hMB.Integral())

    hPrompt.SetStats(0); hDisp.SetStats(0); hMB.SetStats(0)

    gP = ROOT.TF1(f"gP_SA_{suf}", "gaus", xrange[0], xrange[1])
    gD = ROOT.TF1(f"gD_SA_{suf}", "gaus", xrange[0], xrange[1])
    gP.SetLineColor(color_prompt);    gP.SetLineWidth(2); gP.SetLineStyle(2)
    gD.SetLineColor(color_displaced); gD.SetLineWidth(2); gD.SetLineStyle(2)
    hPrompt.Fit(gP, "RQ")
    hDisp.Fit(gD, "RQ")

    c = ROOT.TCanvas(f"c_samuon_z0_{suf}", "", 800, 600)

    histos_to_draw = [hPrompt, hDisp]
    if minbias_data is not None:
        histos_to_draw.append(hMB)
    y_max = max(h.GetMaximum() for h in histos_to_draw) * 1.4
    hPrompt.GetYaxis().SetRangeUser(0, y_max)

    hPrompt.Draw("HIST")
    hDisp.Draw("HIST SAME")
    if minbias_data is not None:
        hMB.Draw("HIST SAME")

    leg = ROOT.TLegend(0.58, 0.62, 0.89, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(hPrompt, "DY + Jets",          "l")
    leg.AddEntry(hDisp,   "Displaced Muon Gun", "l")
    if minbias_data is not None:
        leg.AddEntry(hMB, "MinBias",            "l")
    leg.Draw()

    print("sigma prompt", gP.GetParameter(2))
    print("sigma displaced", gD.GetParameter(2))

    c.SaveAs(f"{direc}/{out_name}.png")
    f = ROOT.TFile(f"{direc}/{out_name}.root", "RECREATE")
    hPrompt.Write(); hDisp.Write(); hMB.Write(); gP.Write(); gD.Write(); c.Write(); f.Close()

    store_plots["canvas"][f"samuon_z0_prompt_displaced_{vertex}_{suf}"] = c
    store_plots["histos"][f"hPrompt_SAz0_{vertex}_{suf}"] = hPrompt
    store_plots["histos"][f"hDisp_SAz0_{vertex}_{suf}"]   = hDisp
    store_plots["histos"][f"hMB_SAz0_{vertex}_{suf}"]     = hMB
    store_plots["histos"][f"leg_SAz0_{vertex}_{suf}"]     = leg
    store_plots["fits"][f"gP_SAz0_{vertex}_{suf}"]        = gP
    store_plots["fits"][f"gD_SAz0_{vertex}_{suf}"]        = gD

    return c, hPrompt, hDisp, hMB, gP, gD


def _fill_hwD0_vs_hwZ0_TH2(datasets,legend_label,out_name,title,abs_val,show=False,xrange=None,yrange=None,xbins=100,ybins=100):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    z_list=[]
    d_list=[]
    for data in datasets:
        z=np.array(data["samuon_hwZ0_displaced"],dtype=float)
        d=np.array(data["samuon_hwD0_displaced"],dtype=float)
        if len(z)!=len(d):
            raise ValueError("displaced hwZ0 and hwD0 arrays must be the same length.")
        z_list.append(z)
        d_list.append(d)
    hwZ0=np.concatenate(z_list) if z_list else np.array([],dtype=float)
    hwD0=np.concatenate(d_list) if d_list else np.array([],dtype=float)
    if not show:
        ROOT.gROOT.SetBatch(True)
    z_for_range=hwZ0 if len(hwZ0)>0 else np.array([0.0])
    d_for_range=hwD0 if len(hwD0)>0 else np.array([0.0])
    if xrange is None:
        zmin,zmax=float(z_for_range.min()),float(z_for_range.max())
        zpad=0.05*(zmax-zmin) if zmax>zmin else 1.0
        xrange=(zmin-zpad,zmax+zpad)
    if yrange is None:
        dmin,dmax=float(d_for_range.min()),float(d_for_range.max())
        dpad=0.05*(dmax-dmin) if dmax>dmin else 1.0
        yrange=(dmin-dpad,dmax+dpad)
    direc=make_plot_dir("hwD0_vs_hwZ0_SAMuon")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.15)
    h=ROOT.TH2D(f"h_{uniq}",title,xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    h.SetStats(0)
    for z,d in zip(hwZ0,hwD0):
        if abs_val == True:
            h.Fill(float(abs(z)),float(abs(d)))
        else:
            h.Fill(float(z),float(d))
    h.Draw("COLZ")
    leg=ROOT.TLegend(0.50,0.80,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.032)
    leg.AddEntry(0,legend_label,"")
    #leg.Draw()
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c,h


def plot_hwD0_vs_hwZ0_SAMuon(prompt_data=None,disp_data1=None,disp_data2=None,disp_data3=None,MB_data=None,show=False,xrange=None,yrange=None,xbins=100,ybins=100, abs_val=False):
    results={}
    if prompt_data is not None:
        results["prompt"]=_fill_hwD0_vs_hwZ0_TH2(
            [prompt_data],"Prompt DY + Jets","hwD0_vs_hwZ0_SAMuon_prompt_2d",
            "Prompt SAMuon hwD0 vs hwZ0;hwZ0;hwD0",abs_val=abs_val,
            show=show,xrange=xrange,yrange=yrange,xbins=xbins,ybins=ybins)
    disp_datasets=[d for d in (disp_data1,disp_data2,disp_data3) if d is not None]
    if disp_datasets:
        results["displaced"]=_fill_hwD0_vs_hwZ0_TH2(
            disp_datasets,"Displaced Muon Gun sample","hwD0_vs_hwZ0_SAMuon_displaced_2d",
            "Displaced SAMuon hwD0 vs hwZ0;hwZ0;hwD0", abs_val=abs_val,
            show=show,xrange=xrange,yrange=yrange,xbins=xbins,ybins=ybins)
    if MB_data is not None:
        results["minbias"]=_fill_hwD0_vs_hwZ0_TH2(
            [MB_data],"MinBias","hwD0_vs_hwZ0_SAMuon_minbias_2d",
            "MinBias SAMuon hwD0 vs hwZ0;hwZ0;hwD0",abs_val=abs_val,
            show=show,xrange=xrange,yrange=yrange,xbins=xbins,ybins=ybins)
    return results

def _fill_dxy_vs_z_KMTF_TH2(datasets,legend_label,out_name,title,show=False,xrange=None,yrange=None,xbins=100,ybins=100):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    z_list=[]
    d_list=[]
    for data in datasets:
        z=np.array(data["kmtf_zvtx"],dtype=float)
        d=np.array(data["kmtf_dxy"],dtype=float)
        if len(z)!=len(d):
            raise ValueError("kmtf_zvtx and kmtf_dxy arrays must be the same length.")
        z_list.append(z)
        d_list.append(d)
    z_vals=np.concatenate(z_list) if z_list else np.array([],dtype=float)
    dxy_vals=np.concatenate(d_list) if d_list else np.array([],dtype=float)
    if not show:
        ROOT.gROOT.SetBatch(True)
    z_for_range=z_vals if len(z_vals)>0 else np.array([0.0])
    d_for_range=dxy_vals if len(dxy_vals)>0 else np.array([0.0])
    if xrange is None:
        zmin,zmax=float(z_for_range.min()),float(z_for_range.max())
        zpad=0.05*(zmax-zmin) if zmax>zmin else 1.0
        xrange=(zmin-zpad,zmax+zpad)
    if yrange is None:
        dmin,dmax=float(d_for_range.min()),float(d_for_range.max())
        dpad=0.05*(dmax-dmin) if dmax>dmin else 1.0
        yrange=(dmin-dpad,dmax+dpad)
    direc=make_plot_dir("dxy_vs_z_KMTFTrack")
    uniq=f"{out_name}_{_simple_plot_call_count}"
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.15)
    h=ROOT.TH2D(f"h_{uniq}",title,xbins,xrange[0],xrange[1],ybins,yrange[0],yrange[1])
    h.SetDirectory(0)
    h.SetStats(0)
    for z,d in zip(z_vals,dxy_vals):
        h.Fill(float(z),float(d))
    h.Draw("COLZ")
    leg=ROOT.TLegend(0.50,0.80,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.032)
    leg.AddEntry(0,legend_label,"")
    leg.Draw()
    store_plots["canvas"][out_name]=c
    store_plots["histos"][out_name]=h
    store_plots["histos"][f"{out_name}_legend"]=leg
    c.Update()
    c.SaveAs(f"{direc}/{out_name}.png")
    return c,h

def plot_dxy_vs_z_KMTFTrack(prompt_data=None,disp_data1=None,disp_data2=None,disp_data3=None,MB_data=None, show=False,xrange=None,yrange=None,xbins=100,ybins=100):
    results={}
    if prompt_data is not None:
        results["prompt"]=_fill_dxy_vs_z_KMTF_TH2(
            [prompt_data],"DY + Jets","dxy_vs_z_KMTFTrack_prompt_2d",
            "Prompt KMTFTrack dxy vs z;z;dxy",
            show=show,xrange=xrange,yrange=yrange,xbins=xbins,ybins=ybins)
    disp_datasets=[d for d in (disp_data1,disp_data2,disp_data3) if d is not None]
    if disp_datasets:
        results["displaced"]=_fill_dxy_vs_z_KMTF_TH2(
            disp_datasets,"Displaced Muon Gun sample","dxy_vs_z_KMTFTrack_displaced_2d",
            "Displaced KMTFTrack dxy vs z;z;dxy",
            show=show,xrange=xrange,yrange=yrange,xbins=xbins,ybins=ybins)
    if MB_data is not None:
        results["minbias"]=_fill_dxy_vs_z_KMTF_TH2(
            [MB_data],"MinBias","dxy_vs_z_KMTFTrack_minbias_2d",
            "MinBias KMTFTrack dxy vs z;z;dxy",
            show=show,xrange=xrange,yrange=yrange,xbins=xbins,ybins=ybins)
    return results

def plot_eff_vs_pT_SAMuon_d0z0(data1=None,data2=None,data3=None,d0_cut=10,z0_cut=5,eta_cut=0.83,pt_cut=5,show=False,pt_bins=(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120),title=None,legend_text="displaced SAMuon",color=ROOT.kBlue):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    uniq=f"samuon_eff_d0z0_{_simple_plot_call_count}"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("samuon_eff_vs_pt")
    if title is None:
        title=f"SAMuon efficiency vs gen p_{{T}} (p_{{T}}>{pt_cut} & |#eta|<{eta_cut} & |hwD0|>{d0_cut} OR |hwZ0|>{z0_cut});gen p_{{T}} [GeV];Efficiency"
    datasets=[d for d in (data1,data2,data3) if d is not None]
    if not datasets:
        raise ValueError("at least one of data1, data2, data3 must be provided.")
    bins=array('d',pt_bins)
    hDen=ROOT.TH1D(f"hDen_{uniq}",";gen p_{T} [GeV];Efficiency",len(pt_bins)-1,bins)
    hNum=ROOT.TH1D(f"hNum_{uniq}",";gen p_{T} [GeV];Efficiency",len(pt_bins)-1,bins)
    hDen.SetDirectory(0);hNum.SetDirectory(0)
    n_pass=0;n_tot=0
    for data in datasets:
        gen_pt=np.array(data["gen_pt_SAMuons_displaced_drmatched"],dtype=float)
        gen_eta=np.array(data["gen_eta_SAMuons_displaced_allmatched"],dtype=float)
        sa_pt=np.array(data["samuon_phPt_displaced_matched"],dtype=float)
        hwD0=np.array(data["samuon_hwD0_displaced_matched"],dtype=float)
        hwZ0=np.array(data["samuon_hwZ0_displaced_matched"],dtype=float)
        if not (len(gen_pt)==len(gen_eta)==len(sa_pt)==len(hwD0)==len(hwZ0)):
            raise ValueError("matched gen_pt, gen_eta, sa_pt, hwD0, hwZ0 arrays must be the same length.")
        mask=(sa_pt>pt_cut)&(np.abs(gen_eta)<eta_cut)&((np.abs(hwD0)>d0_cut)|(np.abs(hwZ0)>z0_cut))
        num_pt=gen_pt[mask]
        n_pass+=int(mask.sum());n_tot+=len(mask)
        den_pt=np.array(data["gen_pt_unmatched"],dtype=float)
        den_eta=np.array(data["gen_eta_unmatched_full"],dtype=float)
        if len(den_pt)!=len(den_eta):
            raise ValueError("gen_pt_unmatched and gen_eta_unmatched_full must be the same length.")
        for pt in den_pt[np.abs(den_eta)<eta_cut]:
            hDen.Fill(float(pt))
        for pt in num_pt:
            hNum.Fill(float(pt))
    eff=ROOT.TEfficiency(hNum,hDen)
    eff.SetName(f"eff_{uniq}")
    eff.SetTitle(title)
    eff.SetLineColor(color);eff.SetMarkerColor(color);eff.SetMarkerStyle(20)
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    eff.Draw("AP")
    leg=ROOT.TLegend(0.40,0.20,0.88,0.32)
    leg.SetBorderSize(0.1);leg.SetFillStyle(0);leg.SetTextSize(0.030)
    leg.AddEntry(eff,legend_text,"lp")
    #leg.AddEntry(0,f"SAMuon p_{{T}}>{pt_cut}, |#eta|<{eta_cut}, |hwD0|>{d0_cut} OR |hwZ0|>{z0_cut} ({n_pass}/{n_tot})","")
    leg.Draw()
    c.SaveAs(f"{direc}/samuon_eff_vs_pt_d0z0.png")
    f=ROOT.TFile("samuon_eff_vs_pt_d0z0.root","RECREATE")
    hDen.Write();hNum.Write();eff.Write();c.Write();f.Close()
    store_plots["canvas"]["samuon_eff_vs_pt_d0z0"]=c
    store_plots["histos"]["hDen_samuon_d0z0"]=hDen
    store_plots["histos"]["hNum_samuon_d0z0"]=hNum
    store_plots["histos"]["eff_samuon_d0z0"]=eff
    store_plots["histos"]["eff_samuon_d0z0_legend"]=leg
    return c,eff

def plot_eff_vs_d3d_SAMuon(data1=None,data2=None,data3=None,d3d_cut=10,eta_cut=0.83,gen_pt_cut=20,sa_pt_cut=5,show=False,d3d_bins=None,n_bins=30,title=None,legend_text="displaced SAMuon",color=ROOT.kBlue,color_nocut=ROOT.kRed):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    uniq=f"samuon_eff_d3d_{_simple_plot_call_count}"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("samuon_eff_vs_d3d")
    if title is None:
        title=f"SAMuon efficiency vs gen d3d (gen p_{{T}}>{gen_pt_cut}, SAMuon p_{{T}}>{sa_pt_cut}, |#eta|<{eta_cut});gen d3d [cm];Efficiency"
    datasets=[d for d in (data1,data2,data3) if d is not None]
    if not datasets:
        raise ValueError("at least one of data1, data2, data3 must be provided.")
    #first pass: collect the gen d3d that pass the cuts across all datasets, for the denominator
    #and for TWO numerators - one WITH the SAMuon d3d cut and one WITHOUT it. we bin only after
    #seeing the data so the default range covers everything (no silent overflow).
    num_with_chunks=[];num_nocut_chunks=[];den_chunks=[]
    n_with=0;n_nocut=0;n_tot=0
    for data in datasets:
        #numerator base: gen muons dr-matched to a displaced SAMuon, with gen pT>gen_pt_cut,
        #the matched SAMuon's reco pT>sa_pt_cut (numerator-only), and |gen eta|<eta_cut. the
        #"with" numerator additionally requires reco d3d=sqrt(hwD0^2+hwZ0^2)>d3d_cut. all arrays
        #are filled together per matched muon upstream, so they are index-aligned. plotted vs gen d3d.
        gen_d3d=np.array(data["gen_d3d_SAMuons_displaced_drmatched"],dtype=float)
        gen_pt=np.array(data["gen_pt_SAMuons_displaced_drmatched"],dtype=float)
        gen_eta=np.array(data["gen_eta_SAMuons_displaced_allmatched"],dtype=float)
        sa_pt=np.array(data["samuon_phPt_displaced_matched"],dtype=float)
        hwD0=np.array(data["samuon_hwD0_displaced_matched"],dtype=float)
        hwZ0=np.array(data["samuon_hwZ0_displaced_matched"],dtype=float)
        if not (len(gen_d3d)==len(gen_pt)==len(gen_eta)==len(sa_pt)==len(hwD0)==len(hwZ0)):
            raise ValueError("matched gen_d3d, gen_pt, gen_eta, sa_pt, hwD0, hwZ0 arrays must be the same length.")
        sa_d3d=np.sqrt(hwD0**2+hwZ0**2)
        base=(gen_pt>gen_pt_cut)&(sa_pt>sa_pt_cut)&(np.abs(gen_eta)<eta_cut)
        mask_nocut=base
        mask_with=base&(sa_d3d>d3d_cut)
        num_nocut_chunks.append(gen_d3d[mask_nocut])
        num_with_chunks.append(gen_d3d[mask_with])
        n_nocut+=int(mask_nocut.sum());n_with+=int(mask_with.sum());n_tot+=len(base)
        #denominator: all gen muons with gen pT>gen_pt_cut and |gen eta|<eta_cut (NO pT-reco, NO
        #d3d cut), binned by gen d3d. gen_d3d_unmatched_full / gen_pt_unmatched /
        #gen_eta_unmatched_full are all aligned 1:1 (same all-gen selection).
        den_d3d=np.array(data["gen_d3d_unmatched_full"],dtype=float)
        den_pt=np.array(data["gen_pt_unmatched"],dtype=float)
        den_eta=np.array(data["gen_eta_unmatched_full"],dtype=float)
        if not (len(den_d3d)==len(den_pt)==len(den_eta)):
            raise ValueError("gen_d3d_unmatched_full, gen_pt_unmatched, gen_eta_unmatched_full must be the same length.")
        den_mask=(den_pt>gen_pt_cut)&(np.abs(den_eta)<eta_cut)
        den_chunks.append(den_d3d[den_mask])
    num_with_vals=np.concatenate(num_with_chunks) if num_with_chunks else np.array([],dtype=float)
    num_nocut_vals=np.concatenate(num_nocut_chunks) if num_nocut_chunks else np.array([],dtype=float)
    den_vals=np.concatenate(den_chunks) if den_chunks else np.array([],dtype=float)
    #binning: if d3d_bins not given, auto-range from 0 to the max gen d3d in the denominator
    #(the full population) so no entry falls in the overflow bin and gets dropped by TEfficiency.
    if d3d_bins is None:
        dmax=float(den_vals.max()) if len(den_vals)>0 else 1.0
        if not dmax>0:
            dmax=1.0
        edges=np.linspace(0.0,dmax*1.02,n_bins+1)
        bins=array('d',edges);nb=n_bins
    else:
        bins=array('d',d3d_bins);nb=len(d3d_bins)-1
    #count anything beyond the last edge so silent overflow loss is visible (only possible with
    #user-supplied d3d_bins; the auto range above already covers the data).
    last_edge=float(bins[len(bins)-1])
    n_over_den=int(np.sum(den_vals>last_edge));n_over_num=int(np.sum(num_nocut_vals>last_edge))
    hDen=ROOT.TH1D(f"hDen_{uniq}",";gen d3d [cm];Efficiency",nb,bins)
    hNumWith=ROOT.TH1D(f"hNumWith_{uniq}",";gen d3d [cm];Efficiency",nb,bins)
    hNumNoCut=ROOT.TH1D(f"hNumNoCut_{uniq}",";gen d3d [cm];Efficiency",nb,bins)
    hDen.SetDirectory(0);hNumWith.SetDirectory(0);hNumNoCut.SetDirectory(0)
    for d in den_vals:
        hDen.Fill(float(d))
    for d in num_with_vals:
        hNumWith.Fill(float(d))
    for d in num_nocut_vals:
        hNumNoCut.Fill(float(d))
    if n_over_den>0 or n_over_num>0:
        print(f"[plot_eff_vs_d3d_SAMuon] WARNING: {n_over_den} denom / {n_over_num} num entries are above the last d3d bin edge ({last_edge:g} cm) and were NOT counted (overflow). Widen d3d_bins to keep them.")
    effNoCut=ROOT.TEfficiency(hNumNoCut,hDen)
    effWith=ROOT.TEfficiency(hNumWith,hDen)
    effNoCut.SetName(f"eff_nocut_{uniq}");effWith.SetName(f"eff_with_{uniq}")
    effNoCut.SetTitle(title)
    effNoCut.SetLineColor(color_nocut);effNoCut.SetMarkerColor(color_nocut);effNoCut.SetMarkerStyle(21)
    effWith.SetLineColor(color);effWith.SetMarkerColor(color);effWith.SetMarkerStyle(20)
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    effNoCut.Draw("AP")
    effWith.Draw("P SAME")
    leg=ROOT.TLegend(0.40,0.78,0.88,0.88)
    leg.SetBorderSize(1);leg.SetFillStyle(0);leg.SetTextSize(0.028)
    leg.AddEntry(0,f"gen p_{{T}}>{gen_pt_cut} && SAMuon p_{{T}}>{sa_pt_cut} && |#eta|<{eta_cut}","")
    leg.AddEntry(effNoCut,f"{legend_text}, no d3d cut","lp")
    leg.AddEntry(effWith,f"{legend_text}, d3d>{d3d_cut}","lp")
    leg.Draw()
    c.SaveAs(f"{direc}/samuon_eff_vs_d3d.png")
    f=ROOT.TFile("samuon_eff_vs_d3d.root","RECREATE")
    hDen.Write();hNumWith.Write();hNumNoCut.Write();effWith.Write();effNoCut.Write();c.Write();f.Close()
    store_plots["canvas"]["samuon_eff_vs_d3d"]=c
    store_plots["histos"]["hDen_samuon_d3d"]=hDen
    store_plots["histos"]["hNumWith_samuon_d3d"]=hNumWith
    store_plots["histos"]["hNumNoCut_samuon_d3d"]=hNumNoCut
    store_plots["histos"]["eff_with_samuon_d3d"]=effWith
    store_plots["histos"]["eff_nocut_samuon_d3d"]=effNoCut
    store_plots["histos"]["eff_samuon_d3d_legend"]=leg
    return c,effWith,effNoCut

def plot_rate_SAMuon(data,vertex="displaced",bunchfactor=40000*2760.0/3564.0,n_bins=20,pt_min=0,pt_max=100,show=False,logy=True,title=None,legend_text=None,color=ROOT.kBlue, draw_opt="HIST"):
    #  1) per event take the leading (max pT) SAMuon pT (collected upstream in event_loop),
    #  2) fill a pT histogram once per event,
    #  3) reverse-cumulative (GetCumulative(False)) -> # events with leading pT >= threshold,
    #  4) Scale(bunchfactor/n_events) -> rate in kHz.
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    uniq=f"samuon_rate_{vertex}_{_simple_plot_call_count}"
    key=f"samuon_maxpt_event_{vertex}"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("samuon_rate")
    if key not in data:
        raise KeyError(f"'{key}' not in data - rerun event_loop so the per-event leading SAMuon pT is collected.")
    maxpt=np.array(data[key],dtype=float)
    n_events=int(data.get("n_events",0))
    if n_events<=0:
        raise ValueError("data['n_events'] must be >0 to scale the rate (rerun event_loop).")
    if title is None:
        title=f"{vertex} SAMuon trigger rate;L1 p_{{T}} threshold [GeV];Rate [kHz]"
    if legend_text is None:
        legend_text=f"{vertex} SAMuon"
    hRate=ROOT.TH1D(f"hRate_{uniq}",title,n_bins,pt_min,pt_max)
    hRate.SetDirectory(0)
    hRate.SetStats(0)
    for pt in maxpt:
        hRate.Fill(float(pt))
    rate=hRate.GetCumulative(False)
    rate.SetDirectory(0)
    rate.Scale(float(bunchfactor)/float(n_events))
    rate.SetName(f"rate_{uniq}")
    rate.SetTitle(title)
    rate.SetLineColor(color);rate.SetMarkerColor(color);rate.SetMarkerStyle(20);rate.SetLineWidth(2)
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    if logy:
        c.SetLogy()
    rate.Draw(draw_opt)
    leg=ROOT.TLegend(0.6,0.74,0.88,0.88)
    leg.SetBorderSize(1);leg.SetFillStyle(0);leg.SetTextSize(0.030)
    leg.AddEntry(rate,legend_text,"l")
    leg.Draw()
    
    c.SaveAs(f"{direc}/samuon_rate_{vertex}.png")
    f=ROOT.TFile(f"samuon_rate_{vertex}.root","RECREATE")
    hRate.Write();rate.Write();c.Write();f.Close()
    store_plots["canvas"][f"samuon_rate_{vertex}"]=c
    store_plots["histos"][f"hRate_samuon_{vertex}"]=hRate
    store_plots["histos"][f"rate_samuon_{vertex}"]=rate
    store_plots["histos"][f"rate_samuon_{vertex}_legend"]=leg
    return c,rate

def plot_rate_SAMuon_overlay(data,bunchfactor=40000*2760.0/3564.0,n_bins=100,pt_min=0,pt_max=100,show=False,logy=True,draw_opt="HIST",title=None,color_prompt=ROOT.kRed,color_displaced=ROOT.kBlue):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    uniq=f"samuon_rate_overlay_{_simple_plot_call_count}"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("samuon_rate")
    n_events=int(data.get("n_events",0))
    if n_events<=0:
        raise ValueError("data['n_events'] must be >0 to scale the rate (rerun event_loop).")
    if title is None:
        title="SAMuon trigger rate (prompt vs displaced);L1 p_{T} threshold [GeV];Rate [kHz]"
    def _make_rate(vertex):
        key=f"samuon_maxpt_event_{vertex}"
        if key not in data:
            raise KeyError(f"'{key}' not in data - rerun event_loop so the per-event leading SAMuon pT is collected.")
        maxpt=np.array(data[key],dtype=float)
        h=ROOT.TH1D(f"hRate_{vertex}_{uniq}",title,n_bins,pt_min,pt_max)
        h.SetDirectory(0)
        for pt in maxpt:
            h.Fill(float(pt))
        nb_x=h.GetNbinsX()
        h.SetBinContent(nb_x,h.GetBinContent(nb_x)+h.GetBinContent(nb_x+1))
        h.SetBinContent(nb_x+1,0.0)
        r=h.GetCumulative(False)
        r.SetDirectory(0)
        r.Scale(float(bunchfactor)/float(n_events))
        r.SetName(f"rate_{vertex}_{uniq}")
        r.SetTitle(title)
        h.SetStats(0)
        r.SetStats(0)
        return h,r,int(len(maxpt))
    hP,rateP,nP=_make_rate("prompt")
    hD,rateD,nD=_make_rate("displaced")
    rateP.SetLineColor(color_prompt);rateP.SetMarkerColor(color_prompt);rateP.SetMarkerStyle(20);rateP.SetLineWidth(2)
    rateD.SetLineColor(color_displaced);rateD.SetMarkerColor(color_displaced);rateD.SetMarkerStyle(21);rateD.SetLineWidth(2)
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    if logy:
        c.SetLogy()
    rateD.Draw(draw_opt)
    rateP.Draw(draw_opt+" SAME")
    leg=ROOT.TLegend(0.6,0.74,0.88,0.88)
    leg.SetBorderSize(1);leg.SetFillStyle(0);leg.SetTextSize(0.030)
    leg.AddEntry(rateD,"displaced SAMuon","l")
    leg.AddEntry(rateP,"prompt SAMuon","l")
    leg.AddEntry(0,f"{n_events} events","")
    leg.Draw()
    c.SaveAs(f"{direc}/samuon_rate_overlay.png")
    f=ROOT.TFile("samuon_rate_overlay.root","RECREATE")
    hP.Write();hD.Write();rateP.Write();rateD.Write();c.Write();f.Close()
    store_plots["canvas"]["samuon_rate_overlay"]=c
    store_plots["histos"]["hRate_prompt_overlay"]=hP
    store_plots["histos"]["hRate_displaced_overlay"]=hD
    store_plots["histos"]["rate_prompt_overlay"]=rateP
    store_plots["histos"]["rate_displaced_overlay"]=rateD
    store_plots["histos"]["rate_overlay_legend"]=leg
    return c,rateP,rateD

def plot_leading_pt_spectrum_SAMuon(data,n_bins=100,pt_min=0,pt_max=100,show=False,logy=True,draw_opt="HIST",title=None,color_prompt=ROOT.kRed,color_displaced=ROOT.kBlue):
    global _simple_plot_call_count
    _simple_plot_call_count += 1
    uniq=f"samuon_leadpt_{_simple_plot_call_count}"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("samuon_leading_pt")
    if title is None:
        title="Leading SAMuon p_{T} per event (prompt vs displaced);leading L1 p_{T} [GeV];events / bin"
    def _make_spectrum(vertex):
        key=f"samuon_maxpt_event_{vertex}"
        if key not in data:
            raise KeyError(f"'{key}' not in data - rerun event_loop so the per-event leading SAMuon pT is collected.")
        maxpt=np.array(data[key],dtype=float)
        h=ROOT.TH1D(f"hLeadPt_{vertex}_{uniq}",title,n_bins,pt_min,pt_max)
        h.SetDirectory(0);h.SetStats(0)
        for pt in maxpt:
            h.Fill(float(pt))
        return h,int(len(maxpt))
    hP,nP=_make_spectrum("prompt")
    hD,nD=_make_spectrum("displaced")
    hP.SetLineColor(color_prompt);hP.SetMarkerColor(color_prompt);hP.SetLineWidth(2)
    hD.SetLineColor(color_displaced);hD.SetMarkerColor(color_displaced);hD.SetLineWidth(2)
    c=ROOT.TCanvas(f"c_{uniq}","",800,600)
    if logy:
        c.SetLogy()
    ymax=1.4*max(hP.GetMaximum(),hD.GetMaximum(),1.0)
    hD.SetMaximum(ymax)
    if logy:
        hD.SetMinimum(0.5)
    hD.Draw(draw_opt)
    hP.Draw(draw_opt+" SAME")
    leg=ROOT.TLegend(0.45,0.72,0.88,0.88)
    leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(0.030)
    leg.AddEntry(hD,f"displaced SAMuon ({nD})","l")
    leg.AddEntry(hP,f"prompt SAMuon ({nP})","l")
    leg.Draw()
    c.SaveAs(f"{direc}/samuon_leading_pt_spectrum.png")
    f=ROOT.TFile("samuon_leading_pt_spectrum.root","RECREATE")
    hP.Write();hD.Write();c.Write();f.Close()
    store_plots["canvas"]["samuon_leading_pt_spectrum"]=c
    store_plots["histos"]["hLeadPt_prompt"]=hP
    store_plots["histos"]["hLeadPt_displaced"]=hD
    store_plots["histos"]["hLeadPt_legend"]=leg
    return c,hP,hD
