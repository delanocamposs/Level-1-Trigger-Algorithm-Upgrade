import math
import numpy as np
import ROOT
import matplotlib.pyplot as plt
from DataFormats.FWLite import Events, Handle
from trigger_helpers import *

ROOT.FWLiteEnabler.enable()
events=Events("DY_Phase2_200_merged.root")
thetahandle=Handle("L1Phase2MuDTThContainer")
genhandle=Handle("vector<reco::GenParticle>")

#the radii are obtained from interative_plotter.py plot_rho function per stattion. obtained from data.
R_MB_CM={1:445.0,2:526,3:635}
stations=(1,2,3)
ZRES_CONV=65536.0/1500.0
KRES_CONV=65536.0/2
CURV_CONV=(1<<15)/1.25

def event_loop(event_num=100, conv_z=False, conv_k=False):
    gen_eta_glob, gen_z_glob, gen_pt_glob, gen_vz_glob, gen_vr_glob, gen_curv_glob={st: [] for st in stations}, {st: [] for st in stations}, {st: [] for st in stations}, {st:[] for st in stations}, {st:[] for st in stations}, {st:[] for st in stations}
    stub_eta_glob, stub_z_glob, stub_k_glob={st: [] for st in stations}, {st: [] for st in stations}, {st:[] for st in stations}
    delta_z_glob={st:[] for st in stations}
    mu_id_glob={st:[] for st in stations}
    for i, event in enumerate(events):
        if i==event_num+1:
            break
        event.getByLabel("dtTriggerPhase2PrimitiveDigis", "", "L1P2GT", thetahandle)
        container=thetahandle.product().getContainer()
    
        #grab vector of gen muon pT, eta, vz,vy,vx, etc per event
        #also save a gen_z_event which is station-dependent. filled later based on matching (genParticles has no z information - must build the z).
        gen_pt_event=get_gen_muons_pt(event,pt_min=0,pt_max=1000)
        gen_eta_event=get_gen_muons_eta(event,pt_min=0,pt_max=1000)
        gen_vz_event=get_gen_muons_vz(event,pt_min=0,pt_max=1000)
        gen_vy_event=get_gen_muons_vy(event,pt_min=0,pt_max=1000)
        gen_vx_event=get_gen_muons_vx(event,pt_min=0,pt_max=1000)
        gen_vr_event=np.sqrt((np.array(gen_vx_event))**2+(np.array(gen_vy_event))**2)
        gen_curv_event=get_gen_muons_curv(event, pt_min=0, pt_max=1000)
    
        if len(gen_pt_event)!=len(gen_eta_event):
            print("ERROR: size mismatch")
    
        #storing per event information about stub hits organized by station. 
        stub_eta_event={st: [] for st in stations} #filled with stub eta values from: station and z value.
        stub_indices_event={} #filled with indices of matched stubs to gen muons.
        stub_z_event={st: [] for st in stations} #filled with stub z value from knowing base z value from theta primitive.  
        stub_k_event={st: [] for st in stations} 

        gen_pt_matched_by_sta={st: [] for st in stations} #need this to keep stations and pT organized after matching
        gen_eta_matched_by_sta={st: [] for st in stations}
        gen_vz_matched_by_sta={st: [] for st in stations}
        gen_z_matched_by_sta={st: [] for st in stations}
        gen_vr_matched_by_sta={st: [] for st in stations}
        gen_curv_matched_by_sta={st: [] for st in stations}

        stub_z_matched_by_sta={st:[] for st in stations} 
        stub_k_matched_by_sta={st:[] for st in stations} 
        stub_eta_matched_by_sta={st:[] for st in stations} 
        
        delta_z_event={st: [] for st in stations}
        mu_id_by_st={st:[] for st in stations}
    
        #loop through stub hits in each event to obtain:
        #stub eta, stub theta and store in above arrays
        for stub in container:
            st=int(stub.stNum())
            z_raw=stub.z()
            k_raw=stub.k()
            z_phys=z_raw/ZRES_CONV
            k_phys=k_raw/KRES_CONV 

            eta_stub=eta_from_z(z_phys, R_MB_CM[st])
    
            stub_z_event[st].append(z_phys if conv_z else z_raw)
            stub_k_event[st].append(k_phys if conv_k else k_raw)
            stub_eta_event[st].append(eta_stub)
    
        #this loop begins the matching of gen to stub by closest eta. loop through each station and match stub to gen.
        for st in stations:
            #crucial: match_indices_global returns an array the same size as the gen_eta_event array. if the muon doesnt have a stub to match, it reutrns index=None at the same index of the muon.
            #example: gen_eta_event=[0.5,0.75,1.1], stub_eta_event[st]=[0.77,0.45]. the matching function will return: [1, 0, None]. muon 0 matched with stub 1. muon 1 matches with stub 0. muon 2 didnt match.
            match_idx=match_indices_global(gen_eta_event, stub_eta_event[st]) 
            stub_indices_event[st]=match_idx 
        
            #this is the loop where i ensure correct matching between gen muons and stub hits. indices are tracked:
            for mu_idx, stub_idx in enumerate(match_idx):
                if stub_idx is not None:
                    mu_eta=gen_eta_event[mu_idx]            
                    mu_vz=gen_vz_event[mu_idx]
                    mu_vy=gen_vy_event[mu_idx]
                    mu_vx=gen_vx_event[mu_idx]
                    mu_vr=gen_vr_event[mu_idx]
                    mu_pt=gen_pt_event[mu_idx]
                    mu_curv=CURV_CONV*gen_curv_event[mu_idx]
                    gen_z_val=mu_vz+R_MB_CM[st]/(np.tan(2*np.arctan(np.exp(-1*mu_eta))))
    
                    gen_z_matched_by_sta[st].append(gen_z_val)
                    gen_vz_matched_by_sta[st].append(mu_vz)
                    gen_vr_matched_by_sta[st].append(mu_vr)
                    gen_pt_matched_by_sta[st].append(mu_pt)
                    gen_eta_matched_by_sta[st].append(mu_eta)
                    gen_curv_matched_by_sta[st].append(mu_curv)
    
                    stub_eta_matched_by_sta[st].append(stub_eta_event[st][stub_idx])
                    stub_z_matched_by_sta[st].append(stub_z_event[st][stub_idx])
                    stub_k_matched_by_sta[st].append(stub_k_event[st][stub_idx])

                    #this is so each muon can have a unique identifier that depends on the event #, mu_idx in the event.    
                    #of course only true if num muons per event <= 2**16 which obviously will be true
                    mu_id=(i<<16)|mu_idx
                    mu_id_by_st[st].append(mu_id)
    
        for st in stations:  
            gen_eta_glob[st].extend(np.array(gen_eta_matched_by_sta[st]))
            gen_vz_glob[st].extend(np.array(gen_vz_matched_by_sta[st]))
            gen_vr_glob[st].extend(np.array(gen_vr_matched_by_sta[st]))
            gen_z_glob[st].extend(np.array(gen_z_matched_by_sta[st]))
            gen_pt_glob[st].extend(np.array(gen_pt_matched_by_sta[st]))
            gen_curv_glob[st].extend(np.array(gen_curv_matched_by_sta[st]))
    
            stub_z_glob[st].extend(np.array(stub_z_matched_by_sta[st]))
            stub_k_glob[st].extend(np.array(stub_k_matched_by_sta[st]))
            delta_z_glob[st].extend(np.array(gen_z_matched_by_sta[st])-np.array(stub_z_matched_by_sta[st]))
            stub_eta_glob[st].extend(np.array(stub_eta_matched_by_sta[st]))
            mu_id_glob[st].extend(np.array(mu_id_by_st[st])) 
    print(f"successful event loop. events: {event_num}")
    return_dict={"gen_eta":gen_eta_glob, "gen_pt":gen_pt_glob, "gen_z":gen_z_glob, "stub_z":stub_z_glob, "delta_z":delta_z_glob, "gen_vz":gen_vz_glob, "gen_vr":gen_vr_glob, "stub_k":stub_k_glob, "stub_eta":stub_eta_glob, "mu_id":mu_id_glob, "gen_curv":gen_curv_glob}
    return return_dict
