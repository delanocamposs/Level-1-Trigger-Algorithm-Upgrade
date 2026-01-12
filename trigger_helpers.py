import math
import os
import numpy as np
import ROOT
import matplotlib.pyplot as plt
from DataFormats.FWLite import Events, Handle
from trigger_helpers import *
from itertools import permutations

ROOT.FWLiteEnabler.enable()
events=Events("DY_Phase2_200_merged.root")
thetahandle=Handle("L1Phase2MuDTThContainer")
genhandle=Handle("vector<reco::GenParticle>")


def eta_from_z(z_cm,R_cm):
    hyp=math.hypot(z_cm,R_cm)
    if hyp == 0:
        return None
    theta=math.atan2(R_cm,z_cm)
    return -math.log(math.tan(theta/2))

def theta_from_z(z_cm,R_cm):
    hyp=math.hypot(z_cm,R_cm)
    if hyp == 0:
        return None
    theta=math.atan2(R_cm,z_cm)
    return theta

def multiple_scattering(pT, eta, st):
    rad_lengths={"ECAL":[0.8903],
                "HCAL":[1.49],
                "SOLENOID":[8.879],
                "YOKE":[1.757]}

    thickness={"ECAL":[181.1-129.0],
                "HCAL":[286.4-181.1],
                "SOLENOID":[380.0-295.0],
                "YOKE1":[490.5-402.0-38.0],
                "YOKE2":[597.5-490.5-38.0],
                "YOKE3":[597.5-490.5-38.0]}

    ecal=thickness["ECAL"][0]/rad_lengths["ECAL"][0]
    hcal=thickness["HCAL"][0]/rad_lengths["HCAL"][0]
    solenoid=thickness["SOLENOID"][0]/rad_lengths["SOLENOID"][0]
    yoke1=thickness["YOKE1"][0]/rad_lengths["YOKE"][0]
    yoke2=thickness["YOKE2"][0]/rad_lengths["YOKE"][0]
    yoke3=thickness["YOKE3"][0]/rad_lengths["YOKE"][0]

    base=ecal+hcal+solenoid

    if st==1:
        x_over_x0=base
    if st==2:
        x_over_x0=base+yoke1
    if st==3:
        x_over_x0=base+yoke1+yoke2

    #p=pT/(np.cosh(eta))
    theta=(0.0136/pT)*(np.sqrt(x_over_x0))*(1+0.088*np.log(x_over_x0))
    return theta

def match_indices_local(listA, listB, max_diff=float("inf")):
    pairs = []
    for i, a in enumerate(listA):
        for j, b in enumerate(listB):
            diff = abs(a - b)
            if diff <= max_diff:
                pairs.append((diff, i, j))
    pairs.sort(key=lambda x: x[0])
    match_idx = [None] * len(listA)
    used_muons = set()
    used_stubs = set()
    for diff, i, j in pairs:
        if i in used_muons or j in used_stubs:
            continue
        match_idx[i] = j
        used_muons.add(i)
        used_stubs.add(j)
    return match_idx

def match_indices_global(listA,listB,max_diff=0.3):
    nA=len(listA)
    nB=len(listB)
    match_idx=[None]*nA
    if nA==0 or nB==0:
        return match_idx
    pairs=[]
    for i,a in enumerate(listA):
        for j,b in enumerate(listB):
            d=abs(a-b)
            if d<=max_diff:
                pairs.append((d,i,j))
    pairs.sort(key=lambda x:x[0])
    usedA=set()
    usedB=set()
    for d,i,j in pairs:
        if i in usedA or j in usedB:
            continue
        match_idx[i]=j
        usedA.add(i)
        usedB.add(j)
    return match_idx
    

def get_gen_muons_eta(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.eta()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_z(event,st,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    if st==1:
        r=402
    if st==2:
        r=490.5
    if st==3:
        r=597.5
    denom=np.tan(2*np.arctan(np.exp(-1*float(g.eta))))
    gen_z=(r)/(denom)
    val=[float(gen_z) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_pt(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.pt()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_theta(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.p4().theta()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_vz(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.vz()) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_vx(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.vx()) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_vy(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.vy()) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_curv(event,pt_min=0,pt_max=1000):
    event.getByLabel("genParticles", genhandle)
    val=[float((g.charge())/(g.pt())) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<1.3 and pt_min<g.pt()<pt_max]
    return val

def make_plot_dir(name):
    outdir = os.path.join("plot_images", name)
    os.makedirs(outdir, exist_ok=True)
    return outdir

def fit(infile,hname,outfile,fitlo,fithi,exlo,exhi,title,xaxis,yaxis):
    f=ROOT.TFile.Open(infile)
    h=f.Get(hname)
    if not h:
        raise RuntimeError("missing histogram")
    h=h.Clone("h")
    h.SetDirectory(0)
    f.Close()
    ax=h.GetXaxis()
    g=ROOT.TGraphErrors()
    n=0
    for i in range(1,h.GetNbinsX()+1):
        x=ax.GetBinCenter(i)
        if x<fitlo or x>fithi:
            continue
        if exlo<=x<=exhi:
            continue
        y=h.GetBinContent(i)
        e=h.GetBinError(i)
        if e<=0:
            continue
        g.SetPoint(n,x,y)
        g.SetPointError(n,0,e)
        n+=1
    g.SetName("included")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.0)
    fun=ROOT.TF1("fit","[0]+[1]*x+[2]*x*x",fitlo,fithi)
    g.Fit(fun,"RQ")
    c=ROOT.TCanvas("c","",800,600)
    g.GetXaxis().SetLimits(fitlo,fithi)
    g.SetName("included")
    g.SetTitle(f"{title};{xaxis};{yaxis}")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.0)    
    g.Draw("AP")
    fun.Draw("SAME")
    lat=ROOT.TLatex()
    lat.SetNDC(True)
    lat.SetTextSize(0.035)
    lat.DrawLatex(0.15,0.86,f"ax^{{2}} + bx + c")
    lat.DrawLatex(0.15,0.81,f"a = {fun.GetParameter(2):.2e}")
    lat.DrawLatex(0.15,0.76,f"b = {fun.GetParameter(1):.2e}")
    lat.DrawLatex(0.15,0.71,f"c = {fun.GetParameter(0):.2f}")
    o=ROOT.TFile.Open(outfile,"RECREATE")
    g.Write()
    fun.Write("fit")
    c.Write()
    o.Close()

def get_dR(st1,st2,)
