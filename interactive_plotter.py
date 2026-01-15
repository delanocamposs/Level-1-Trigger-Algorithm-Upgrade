import ROOT 
import numpy as np 
from dummy_z import event_loop
from trigger_helpers import *
import os 

store_plots={"canvas":{}, "histos":{}, "fits":{}, "profiles":{}}

def plot_slice_and_gaussian(h2,xbin,st,xvar,xunit,yvar,yunit,min_entries=20,nsig=2.5,slice_xrange=None,slice_yrange=None):
    direc=make_plot_dir(f"gaussianfit_slices_{yvar}_vs_{xvar}")
    name=h2.GetName()
    xaxis=h2.GetXaxis()
    x_low=xaxis.GetBinLowEdge(xbin)
    x_high=xaxis.GetBinUpEdge(xbin)
    slice_name=f"{name}_slice_xbin{xbin}"
    h_slice=h2.ProjectionY(slice_name,xbin,xbin)
    h_slice.SetDirectory(0)
    c=ROOT.TCanvas(f"c_{name}_xbin{xbin}","",800,600)
    h_slice.SetTitle(f"Station {st}: {x_low:.1f}<{xvar}<{x_high:.1f} [{xunit}];{yvar} [{yunit}];Entries")
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

def plot_yslices(data,station,key,xvar,xunit,yvar,yunit,rebinx=2,min_entries=20,nsig=2.5,sigma_xrange=None,sigma_yrange=None,slice_xrange=None,slice_yrange=None):
    ROOT.gROOT.SetBatch(True)
    if key not in store_plots["histos"]:
        raise RuntimeError(f"missing histogram {key}. run function that makes histo first.")
    h0=store_plots["histos"][key]
    h=h0.Clone(f"{h0.GetName()}_tmp_yslices_st{station}")
    h.SetDirectory(0)
    if rebinx and rebinx>1:
        h.RebinX(rebinx)
    direc=make_plot_dir(f"yslices_{yvar}_vs_{xvar}")
    hx=ROOT.TH1F("h_sigma",f"Station {station};{xvar} [{xunit}];#sigma({yvar}) [{yunit}]",h.GetNbinsX(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
    hx.SetDirectory(0)
    hx.GetXaxis().SetTitle(h.GetXaxis().GetTitle() if h.GetXaxis().GetTitle() else f"{xvar} [{xunit}]")
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
    c_sig.SaveAs(f"{direc}/sigma_{yvar}_vs_{xvar}.png")
    store_plots["canvas"]["yslices_sigma"]=c_sig
    store_plots["histos"]["yslices_sigma"]=hx
    print("all y slice plots saved.")
    return c_sig,hx

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

def plot_deltak_vs_curv(data,st,show=True,xrange=(-7000,7000),yrange=(-20000,20000)):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltak_vs_curv")
    c=ROOT.TCanvas(f"c_dk_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dk_curv_{st}",f"(k_{{pred,{st}}}-k_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta k (k_{{pred,{st}}}-k_{{meas,{st}}})",100,xrange[0],xrange[1],100,yrange[0],yrange[1])
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
    f=ROOT.TFile(f"deltak_vs_curv_{st}.root","UPDATE")
    f.cd()
    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    c.Write(c.GetName(),ROOT.TObject.kOverwrite)
    f.Close()

    return c,h 

def plot_deltaz_vs_curv(data,st,conv_k,show=True,xrange=(-7000,7000),yrange=(-10000,10000)):
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltaz_vs_curv")
    c=ROOT.TCanvas(f"c_dz_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dz_curv_{st}",f"(z_{{pred,{st}}}-z_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta z (z_{{pred,{st}}}-z_{{meas,{st}}})",100,xrange[0],xrange[1],100,yrange[0],yrange[1])
    h.SetDirectory(0)
    m={}
    if st==1:
        if not conv_k:
            dR=-.101
        else:
            dR=-76.7
    elif st==2: 
        if not conv_k:
            dR=-.133
        else:
            dR=-101
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
    p.SetTitle(f"station: {st} <z_pred-z_meas> vs curvature;gen q/pT;<#Delta z>")
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

def plot_deltaz_vs_curv_to_vtx(data,show=True,xrange=(-7000,7000),yrange=(-10000,10000)):
    #currently i dont have the option to convert k and z to real coords. 
    #will have to obtain dR for the converted case (dR=-0.579 is for non-converted)
    st="vtx"
    if not show:
        ROOT.gROOT.SetBatch(True)
    direc=make_plot_dir("deltaz_vs_curv")
    c=ROOT.TCanvas(f"c_dz_curv_{st}","",800,600)
    h=ROOT.TH2F(f"h_dz_curv_{st}",f"(z_{{pred,{st}}}-z_{{meas,{st}}}) vs curvature;gen q/p_{{T}};#Delta z (z_{{pred,{st}}}-z_{{meas,{st}}})",100,xrange[0],xrange[1],100,yrange[0],yrange[1])
    h.SetDirectory(0)
    m={}
    dR=-.579
    for muid,z1,curv1,slope1 in zip(data["mu_id"][1],data["stub_z"][1],data["gen_curv"][1],data["stub_k"][1]):
        if muid not in m:
            m[muid]=(z1,curv1,slope1)
    for muid,z_vtx in zip(data["mu_id"][1],data["gen_vz"][1]):
        if muid not in m:
            continue
        z1,curv1,slope1=m[muid]
        z_pred=z1+slope1*dR
        z_vtx=z_vtx*(65536.0/1500.0)
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
    p.SetTitle(f"station: {st} <z_pred-z_meas> vs curvature;gen q/pT;<#Delta z>")
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
#def plot_deltak_vs_curv_to_vtx():
    
if __name__=="__main__":
    data=event_loop(1000, False, False)
