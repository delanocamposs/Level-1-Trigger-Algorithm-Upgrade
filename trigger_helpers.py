import math
import os
import numpy as np
import ROOT
import matplotlib.pyplot as plt
from DataFormats.FWLite import Events, Handle
from itertools import permutations

ROOT.FWLiteEnabler.enable()
genhandle=Handle("vector<reco::GenParticle>")
KMTFhandle=Handle("vector<l1t::SAMuon>")
SAMuonshandle=Handle("vector<l1t::SAMuon>")
Trackshandle=Handle("vector<l1t::KMTFTrack>")

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

def match_indices_global_tracks(gen_eta, gen_phi, track_eta, track_phi, max_dr=0.6):
    n_gen=len(gen_phi)
    n_trk=len(track_phi)
    match_idx=[None]*n_gen
    if n_gen==0 or n_trk==0:
        return match_idx
    if len(gen_eta)!=n_gen:
        raise ValueError("len(gen_eta) must equal len(gen_phi).")
    if len(track_eta)!=n_trk:
        raise ValueError("len(track_eta) must equal len(track_phi).")
    pairs=[]
    for i,(eta_gen,phi_gen) in enumerate(zip(gen_eta,gen_phi)):
        for j,(eta_trk,phi_trk) in enumerate(zip(track_eta,track_phi)):
            dphi=abs(np.arctan2(np.sin(phi_gen-phi_trk), np.cos(phi_gen-phi_trk)))
            deta=abs(eta_gen-eta_trk)
            dr=np.hypot(deta,dphi)
            if dr>max_dr:
                continue
            pairs.append((dr,i,j))
    pairs.sort(key=lambda x:x[0])
    used_gen=set()
    used_trk=set()
    for d,i,j in pairs:
        if i in used_gen or j in used_trk:
            continue
        match_idx[i]=j
        used_gen.add(i)
        used_trk.add(j)
    return match_idx

def get_gen_muons_eta(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.eta()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_z(event,st,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    if st==1:
        r=402
    if st==2:
        r=490.5
    if st==3:
        r=597.5
    denom=np.tan(2*np.arctan(np.exp(-1*float(g.eta))))
    gen_z=(r)/(denom)
    val=[float(gen_z) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_pt(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.pt()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_theta(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.p4().theta()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_vz(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.vz()) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_vx(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.vx()) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_vy(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.vy()) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_curv(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float((g.charge())/(g.pt())) for g in genhandle.product() if abs(g.pdgId())==13 and g.status()==1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def get_gen_muons_phi(event,pt_min=0,pt_max=1000,eta_max=1.3):
    event.getByLabel("genParticles", genhandle)
    val=[float(g.phi()) for g in genhandle.product() if abs(g.pdgId()) == 13 and g.status() == 1 and abs(g.eta())<eta_max and pt_min<g.pt()<pt_max]
    return val

def make_plot_dir(name):
    outdir = os.path.join("plot_images", name)
    os.makedirs(outdir, exist_ok=True)
    return outdir

def get_KMTF_muons_phPt(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", KMTFhandle)
    val=[float(g.phPt()) for g in KMTFhandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_KMTF_muons_phEta(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", KMTFhandle)
    val=[float(g.phEta()) for g in KMTFhandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_phPt(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.phPt()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_phEta(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.phEta()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_phZ0(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.phZ0()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_phPhi(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.phPhi()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_hwPt(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.hwPt()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_hwEta(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.hwEta()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_hwZ0(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.hwZ0()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_hwD0(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.hwD0()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_hwPhi(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    val=[float(g.hwPhi()) for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return val

def get_SAMuons_rate_tuples(event, vertex, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons",vertex,"L1P2GT2", SAMuonshandle)
    out=[(float(g.phPt()),float(g.phEta()),float(g.hwD0()),float(g.hwZ0()))
         for g in SAMuonshandle.product() if abs(g.phEta())<eta_max and pt_min<g.phPt()<pt_max]
    return out

def get_KMTFTrack_dxy(event, thetaDigi=True, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons", "kmtfTracks", "L1P2GT2", Trackshandle)
    tracks = Trackshandle.product()
    if thetaDigi:
        tracks = [g for g in tracks if any(s.etaQuality() > 0 for s in g.stubs()) and abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    else:
        tracks = [g for g in tracks if abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    val = [float(g.dxy()) for g in tracks]
    return val

def get_KMTFTrack_zVtx(event, thetaDigi=True, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons", "kmtfTracks", "L1P2GT2", Trackshandle)
    tracks = Trackshandle.product()
    if thetaDigi:
        tracks = [g for g in tracks if any(s.etaQuality() > 0 for s in g.stubs()) and abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    else:
        tracks = [g for g in tracks if abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    val = [float(g.zPosition()) for g in tracks]
    return val

def get_KMTFTrack_eta(event, thetaDigi=True, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons", "kmtfTracks", "L1P2GT2", Trackshandle)
    tracks = Trackshandle.product() 
    if thetaDigi:
        tracks = [g for g in tracks if any(s.etaQuality() > 0 for s in g.stubs()) and abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    else:
        tracks = [g for g in tracks if abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    val = [float(g.eta()) for g in tracks]
    return val

def get_KMTFTrack_phi(event, thetaDigi=True, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons", "kmtfTracks", "L1P2GT2", Trackshandle)
    tracks = Trackshandle.product() 
    if thetaDigi:
        tracks = [g for g in tracks if any(s.etaQuality() > 0 for s in g.stubs()) and abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    else:
        tracks = [g for g in tracks if abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    val = [float(g.phi()) for g in tracks]
    return val

def get_KMTFTrack_kSlope(event, thetaDigi=True, pt_min=0, pt_max=1000, eta_max=1.3):
    event.getByLabel("gmtKMTFMuons", "kmtfTracks", "L1P2GT2", Trackshandle)
    tracks = Trackshandle.product()
    if thetaDigi:
        tracks = [g for g in tracks if any(s.etaQuality() > 0 for s in g.stubs()) and abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    else:
        tracks = [g for g in tracks if abs(g.eta()) < eta_max and pt_min < g.pt() < pt_max]
    val = [float(g.kSlope()) for g in tracks]
    return val

