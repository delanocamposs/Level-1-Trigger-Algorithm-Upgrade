import math
import re, os
import numpy as np
import ROOT
import matplotlib.pyplot as plt
from DataFormats.FWLite import Events, Handle
from trigger_helpers import *

ROOT.FWLiteEnabler.enable()

#the radii are obtained from interative_plotter.py plot_rho function per stattion. obtained from simulation.
R_MB_CM={1:445.0,2:526,3:635}
stations=(1,2,3)
ZRES_CONV=65536.0/1500.0
KRES_CONV=65536.0/2
CURV_CONV=(1<<15)/1.25

def event_loop(event_num, thetaDigi, eta_max, local_dataset=None,eos_dataset=None, conv_z=False, conv_k=False):
    if local_dataset is None and eos_dataset is None:
        print("you have to choose a dataset to process. set local_dataset or eos_dataset")
    if eos_dataset is not None:
        dataset_path={
        "prompt":"/eos/uscms/store/user/dacampos/L1KMTF/PromptDY_PU200/DYToLL_M-50_TuneCP5_14TeV-pythia8/PHASEII_PromptDY_PU200/260324_080705/0000/",
        "displaced2-10":"/eos/uscms/store/user/dacampos/L1KMTF/DisplacedMu_Pt2To10_PU140/DisplacedMuons_Pt-2To10_Dxy-0To3000-gun/PHASEII_DisplacedMu_Pt2To10_PU140/260324_080742/0000/",
        "displaced10-30":"/eos/uscms/store/user/dacampos/L1KMTF/DisplacedMu_Pt10To30_PU140/DisplacedMuons_Pt-10To30_Dxy-0To3000-gun/PHASEII_DisplacedMu_Pt10To30_PU140/260324_080802/0000/",
        "displaced30-100":"/eos/uscms/store/user/dacampos/L1KMTF/DisplacedMu_Pt30To100_PU140/DisplacedMuons_Pt-30To100_Dxy-0To3000-gun/PHASEII_DisplacedMu_Pt30To100_PU140/260324_080823/0000/",
        "minbias":"/eos/uscms/store/user/dacampos/L1KMTF/MinBias_PU200/MinBias_TuneCP5_14TeV-pythia8/PHASEII_MinBias_PU200/260324_081247/"}
        file_format={
        "prompt":r'^output_Phase2_L1T_(\d+)\.root$',
        "displaced2-10":r'^output_Phase2_L1T_(\d+)\.root$',
        "displaced10-30":r'^output_Phase2_L1T_(\d+)\.root$',
        "displaced30-100":r'^output_Phase2_L1T_(\d+)\.root$',
        "minbias":r'^output_Phase2_L1T_MinBias_(\d+)\.root$'}
        EOS_DIR=dataset_path[eos_dataset]
        pattern = re.compile(file_format[eos_dataset])
        def is_good(filepath):
            try:
                f = ROOT.TFile.Open(filepath)
                if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered) or f.GetNkeys() == 0:
                    return False
                f.Close()
                return True
            except Exception:
                return False
        all_files = sorted([os.path.join(EOS_DIR, f) for f in os.listdir(EOS_DIR) if pattern.match(f)],key=lambda f: int(pattern.match(os.path.basename(f)).group(1)))
        good_files = [f for f in all_files if is_good(f)]
        print(f"Using {len(good_files)} files")
        events=Events(good_files)
    else:
        events=Events(local_dataset)
    thetahandle=Handle("L1Phase2MuDTThContainer")
    genhandle=Handle("vector<reco::GenParticle>")
    KMTFhandle=Handle("vector<l1t::SAMuon>")
    SAMuonshandle=Handle("vector<l1t::SAMuon>")

    gen_eta_glob, gen_z_glob, gen_pt_glob, gen_vz_glob, gen_vr_glob, gen_curv_glob={st: [] for st in stations}, {st: [] for st in stations}, {st: [] for st in stations}, {st:[] for st in stations}, {st:[] for st in stations}, {st:[] for st in stations}
    stub_eta_glob, stub_z_glob, stub_k_glob={st: [] for st in stations}, {st: [] for st in stations}, {st:[] for st in stations}
    delta_z_glob={st:[] for st in stations}
    mu_id_glob={st:[] for st in stations}

    #global arrays outside event loop for KMTF muon information.
    #adding unmatched pt dictionary. this means gen muon pt straight from gen particles without any matching algo to stubs. no station split because genparticles doesnt know about stations.
    gen_pt_unmatched_glob=[]
    gen_pt_KMTFTracks_matched_glob=[]
    gen_pt_KMTF_matched_displaced_glob=[]
    gen_pt_KMTF_matched_prompt_glob=[]
    gen_pt_SAMuons_matched_displaced_glob=[]
    gen_pt_SAMuons_matched_prompt_glob=[]

    gen_eta_unmatched_glob=[]
    gen_eta_KMTFTracks_matched_glob=[]
    gen_eta_KMTF_matched_displaced_glob=[]
    gen_eta_KMTF_matched_prompt_glob=[]
    gen_eta_SAMuons_matched_displaced_glob=[]
    gen_eta_SAMuons_matched_prompt_glob=[]
    kmtf_zvtx_glob = []
    kmtf_dxy_glob = []
    samuon_hwD0_displaced_matched_glob=[]
    gen_pt_SAMuons_displaced_drmatched_glob=[]

    samuon_z_all_prompt_glob = []
    samuon_z_all_displaced_glob = []
    kmtf_zvtx_KMTFTracks_matched_glob=[]
    gen_vz_KMTFTracks_matched_glob=[]
    kmtf_eta_KMTFTracks_matched_glob=[]
    gen_eta_KMTFTracks_allmatched_glob=[]
    kmtf_kslope_KMTFTracks_matched_glob=[]

    samuon_eta_prompt_matched_glob=[]
    gen_eta_SAMuons_prompt_allmatched_glob=[]
    samuon_z_prompt_matched_glob=[]
    gen_vz_SAMuons_prompt_matched_glob=[]
    samuon_eta_displaced_matched_glob=[]
    gen_eta_SAMuons_displaced_allmatched_glob=[]
    samuon_z_displaced_matched_glob=[]
    gen_vz_SAMuons_displaced_matched_glob=[]
    samuon_hwZ0_prompt_matched_glob=[]
    samuon_hwZ0_displaced_matched_glob=[]
    samuon_hwZ0_prompt_glob=[]
    samuon_hwZ0_displaced_glob=[]
    samuon_hwD0_prompt_glob=[]
    samuon_hwD0_displaced_glob=[]

    for i, event in enumerate(events):
        if i>=event_num:
            break
        event.getByLabel("dtTriggerPhase2PrimitiveDigis", "", "L1P2GT", thetahandle)
        
        #this container is needed for stub information because the C++ class which defines the dataformat is not iterable as thetahandle.product(). need to get container then iterate later in loop
        container=thetahandle.product().getContainer()
    
        #grab vector of gen muon pT, eta, vz,vy,vx, etc per event
        #also save a gen_z_event which is station-dependent. filled later based on matching (genParticles has no z information - must build the z).
        gen_pt_event=get_gen_muons_pt(event,pt_min=3,pt_max=1000,eta_max=eta_max)
        gen_eta_event=get_gen_muons_eta(event,pt_min=3,pt_max=1000,eta_max=eta_max)
        gen_phi_event=get_gen_muons_phi(event,pt_min=3,pt_max=1000,eta_max=eta_max)
        gen_vz_event=get_gen_muons_vz(event,pt_min=3,pt_max=1000,eta_max=eta_max)
        gen_vy_event=get_gen_muons_vy(event,pt_min=3,pt_max=1000,eta_max=eta_max)
        gen_vx_event=get_gen_muons_vx(event,pt_min=3,pt_max=1000,eta_max=eta_max)
        gen_vr_event=np.sqrt((np.array(gen_vx_event))**2+(np.array(gen_vy_event))**2)
        gen_curv_event=get_gen_muons_curv(event, pt_min=3, pt_max=1000,eta_max=eta_max)
    
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
            if st not in R_MB_CM:
                continue
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


        kmtf_zvtx_event = get_KMTFTrack_zVtx(event, thetaDigi, pt_min=3,eta_max=eta_max)
        kmtf_zvtx_glob.extend(kmtf_zvtx_event)
        kmtf_dxy_event = get_KMTFTrack_dxy(event, thetaDigi, pt_min=3,eta_max=eta_max)
        kmtf_dxy_glob.extend(kmtf_dxy_event)
        #this loop will attempt to match genmuons to KMTF muons if eta<eta_max and pT>20GeV. used for efficiency plot. 
        #this piece of the code is used to return information about efficiency vs pT in "===================".
        #putting it at the end since it uses different collection (l1tKMTFGmt vs L1Phase2MuDTThContainer). both studies use genParticles however.
#=====================================================================================================================
        gen_pt_unmatched_glob.extend(np.array(gen_pt_event)) #add unmatched gen muon pt to global array WITHOUT 20 GEV PT CUT! denom of efficiency curve. 

        KMTFTrack_eta_event=get_KMTFTrack_eta(event, thetaDigi, pt_min=3, eta_max=eta_max)
        KMTFTrack_phi_event=get_KMTFTrack_phi(event, thetaDigi,pt_min=3, eta_max=eta_max)
        KMTFTrack_zvtx_event=get_KMTFTrack_zVtx(event, thetaDigi, pt_min=3,eta_max=eta_max)
        KMTFTrack_kslope_event=get_KMTFTrack_kSlope(event, thetaDigi, pt_min=3,eta_max=eta_max)
        gen_pt_KMTFTracks_matched_event=[]
        KMTF_phEta_displaced_event=get_KMTF_muons_phEta(event, "displaced",pt_min=3,eta_max=eta_max)
        KMTF_phEta_prompt_event=get_KMTF_muons_phEta(event, "prompt",pt_min=3,eta_max=eta_max)
        gen_pt_KMTF_matched_displaced_event=[]
        gen_pt_KMTF_matched_prompt_event=[]

        SAMuons_phEta_displaced_event=get_SAMuons_phEta(event, "displaced",pt_min=3,eta_max=eta_max)
        SAMuons_phEta_prompt_event=get_SAMuons_phEta(event, "prompt",pt_min=3,eta_max=eta_max)
        SAMuons_phPhi_displaced_event=get_SAMuons_phPhi(event, "displaced",pt_min=3,eta_max=eta_max)
        SAMuons_phPhi_prompt_event=get_SAMuons_phPhi(event, "prompt",pt_min=3,eta_max=eta_max)
        SAMuons_phZ0_displaced_event=get_SAMuons_phZ0(event, "displaced",pt_min=3,eta_max=eta_max)
        SAMuons_phZ0_prompt_event=get_SAMuons_phZ0(event, "prompt",pt_min=3,eta_max=eta_max)
        samuon_z_all_prompt_glob.extend(SAMuons_phZ0_prompt_event)
        samuon_z_all_displaced_glob.extend(SAMuons_phZ0_displaced_event)
        SAMuons_hwZ0_displaced_event=get_SAMuons_hwZ0(event, "displaced",pt_min=3,eta_max=eta_max)
        SAMuons_hwZ0_prompt_event=get_SAMuons_hwZ0(event, "prompt",pt_min=3,eta_max=eta_max)
        SAMuons_hwD0_displaced_event=get_SAMuons_hwD0(event, "displaced",pt_min=3,eta_max=eta_max)
        SAMuons_hwD0_prompt_event=get_SAMuons_hwD0(event, "prompt",pt_min=3,eta_max=eta_max)
        samuon_hwZ0_prompt_glob.extend(SAMuons_hwZ0_prompt_event)
        samuon_hwZ0_displaced_glob.extend(SAMuons_hwZ0_displaced_event)
        samuon_hwD0_prompt_glob.extend(SAMuons_hwD0_prompt_event)
        samuon_hwD0_displaced_glob.extend(SAMuons_hwD0_displaced_event)
        gen_pt_SAMuons_matched_displaced_event=[]
        gen_pt_SAMuons_matched_prompt_event=[]

        #match all gen muons in acceptance to KMTF tracks by global deltaR.
        #IMPORTANT: this is additive to (not a replacement for) the legacy
        #KMTF/SAMuons eta-based matching blocks below.
        match_idx_KMTFTracks=match_indices_global_tracks(gen_eta_event, gen_phi_event, KMTFTrack_eta_event, KMTFTrack_phi_event, max_dr=0.6)
        match_idx_KMTF_displaced=match_indices_global(gen_eta_event, KMTF_phEta_displaced_event)
        match_idx_KMTF_prompt=match_indices_global(gen_eta_event, KMTF_phEta_prompt_event)

        match_idx_SAMuons_displaced=match_indices_global(gen_eta_event, SAMuons_phEta_displaced_event)
        match_idx_SAMuons_prompt=match_indices_global(gen_eta_event, SAMuons_phEta_prompt_event)
        match_idx_SAMuons_prompt_dr=match_indices_global_tracks(gen_eta_event, gen_phi_event, SAMuons_phEta_prompt_event, SAMuons_phPhi_prompt_event, max_dr=0.6)
        match_idx_SAMuons_displaced_dr=match_indices_global_tracks(gen_eta_event, gen_phi_event, SAMuons_phEta_displaced_event, SAMuons_phPhi_displaced_event, max_dr=0.6)

        for mu_idx, KMTFTrack_idx in enumerate(match_idx_KMTFTracks):
            if KMTFTrack_idx is not None:
                gen_pt_KMTFTracks_matched_event.append(gen_pt_event[mu_idx])
                kmtf_zvtx_KMTFTracks_matched_glob.append(KMTFTrack_zvtx_event[KMTFTrack_idx])
                gen_vz_KMTFTracks_matched_glob.append(gen_vz_event[mu_idx])
                kmtf_eta_KMTFTracks_matched_glob.append(KMTFTrack_eta_event[KMTFTrack_idx])
                kmtf_kslope_KMTFTracks_matched_glob.append(KMTFTrack_kslope_event[KMTFTrack_idx])
                gen_eta_KMTFTracks_allmatched_glob.append(gen_eta_event[mu_idx])
        gen_pt_KMTFTracks_matched_glob.extend(gen_pt_KMTFTracks_matched_event)

        for mu_idx, KMTF_idx in enumerate(match_idx_KMTF_displaced):
            if KMTF_idx is not None:
                gen_pt_KMTF_matched_displaced_event.append(gen_pt_event[mu_idx])
        gen_pt_KMTF_matched_displaced_glob.extend(gen_pt_KMTF_matched_displaced_event)

        for mu_idx, KMTF_idx in enumerate(match_idx_KMTF_prompt):
            if KMTF_idx is not None:
                gen_pt_KMTF_matched_prompt_event.append(gen_pt_event[mu_idx])
        gen_pt_KMTF_matched_prompt_glob.extend(gen_pt_KMTF_matched_prompt_event)

        for mu_idx, SAMuons_idx in enumerate(match_idx_SAMuons_displaced):
            if SAMuons_idx is not None:
                gen_pt_SAMuons_matched_displaced_event.append(gen_pt_event[mu_idx])
        gen_pt_SAMuons_matched_displaced_glob.extend(gen_pt_SAMuons_matched_displaced_event)

        for mu_idx, SAMuons_idx in enumerate(match_idx_SAMuons_prompt):
            if SAMuons_idx is not None:
                gen_pt_SAMuons_matched_prompt_event.append(gen_pt_event[mu_idx])
        gen_pt_SAMuons_matched_prompt_glob.extend(gen_pt_SAMuons_matched_prompt_event)

        for mu_idx, SAMuons_idx in enumerate(match_idx_SAMuons_prompt_dr):
            if SAMuons_idx is not None:
                samuon_eta_prompt_matched_glob.append(SAMuons_phEta_prompt_event[SAMuons_idx])
                gen_eta_SAMuons_prompt_allmatched_glob.append(gen_eta_event[mu_idx])
                samuon_z_prompt_matched_glob.append(SAMuons_phZ0_prompt_event[SAMuons_idx])
                gen_vz_SAMuons_prompt_matched_glob.append(gen_vz_event[mu_idx])
                samuon_hwZ0_prompt_matched_glob.append(SAMuons_hwZ0_prompt_event[SAMuons_idx])

        for mu_idx, SAMuons_idx in enumerate(match_idx_SAMuons_displaced_dr):
            if SAMuons_idx is not None:
                samuon_eta_displaced_matched_glob.append(SAMuons_phEta_displaced_event[SAMuons_idx])
                gen_eta_SAMuons_displaced_allmatched_glob.append(gen_eta_event[mu_idx])
                samuon_z_displaced_matched_glob.append(SAMuons_phZ0_displaced_event[SAMuons_idx])
                gen_vz_SAMuons_displaced_matched_glob.append(gen_vz_event[mu_idx])
                samuon_hwZ0_displaced_matched_glob.append(SAMuons_hwZ0_displaced_event[SAMuons_idx])
                samuon_hwD0_displaced_matched_glob.append(SAMuons_hwD0_displaced_event[SAMuons_idx])
                gen_pt_SAMuons_displaced_drmatched_glob.append(gen_pt_event[mu_idx])


        for mu_idx in range(len(gen_pt_event)):
            if gen_pt_event[mu_idx]>25:
                gen_eta_unmatched_glob.append(gen_eta_event[mu_idx])

        for mu_idx, KMTFTrack_idx in enumerate(match_idx_KMTFTracks):
            if KMTFTrack_idx is not None and gen_pt_event[mu_idx]>25:
                gen_eta_KMTFTracks_matched_glob.append(gen_eta_event[mu_idx])

        for mu_idx, KMTF_idx in enumerate(match_idx_KMTF_displaced):
            if KMTF_idx is not None and gen_pt_event[mu_idx]>25:
                gen_eta_KMTF_matched_displaced_glob.append(gen_eta_event[mu_idx])

        for mu_idx, KMTF_idx in enumerate(match_idx_KMTF_prompt):
            if KMTF_idx is not None and gen_pt_event[mu_idx]>25:
                gen_eta_KMTF_matched_prompt_glob.append(gen_eta_event[mu_idx])

        for mu_idx, SAMuons_idx in enumerate(match_idx_SAMuons_displaced):
            if SAMuons_idx is not None and gen_pt_event[mu_idx]>25:
                gen_eta_SAMuons_matched_displaced_glob.append(gen_eta_event[mu_idx])

        for mu_idx, SAMuons_idx in enumerate(match_idx_SAMuons_prompt):
            if SAMuons_idx is not None and gen_pt_event[mu_idx]>25:
                gen_eta_SAMuons_matched_prompt_glob.append(gen_eta_event[mu_idx])

#=====================================================================================================================
    print(f"successful event loop. events: {event_num}")
    return_dict={
        "gen_eta":gen_eta_glob, 
        "gen_pt":gen_pt_glob, 
        "gen_z":gen_z_glob, 
        "stub_z":stub_z_glob, 
        "delta_z":delta_z_glob, 
        "gen_vz":gen_vz_glob, 
        "gen_vr":gen_vr_glob, 
        "stub_k":stub_k_glob, 
        "stub_eta":stub_eta_glob, 
        "mu_id":mu_id_glob, 
        "gen_curv":gen_curv_glob, 
        "gen_pt_unmatched":gen_pt_unmatched_glob, 
        "gen_pt_KMTFTracks_matched":gen_pt_KMTFTracks_matched_glob,
        "gen_pt_KMTF_displaced_matched":gen_pt_KMTF_matched_displaced_glob, 
        "gen_pt_KMTF_prompt_matched":gen_pt_KMTF_matched_prompt_glob,
        "gen_pt_SAMuons_displaced_matched":gen_pt_SAMuons_matched_displaced_glob,
        "gen_pt_SAMuons_prompt_matched":gen_pt_SAMuons_matched_prompt_glob,
        "gen_eta_unmatched":gen_eta_unmatched_glob,
        "gen_eta_KMTFTracks_matched":gen_eta_KMTFTracks_matched_glob,
        "gen_eta_KMTF_displaced_matched":gen_eta_KMTF_matched_displaced_glob,
        "gen_eta_KMTF_prompt_matched":gen_eta_KMTF_matched_prompt_glob,
        "gen_eta_SAMuons_displaced_matched":gen_eta_SAMuons_matched_displaced_glob,
        "gen_eta_SAMuons_prompt_matched":gen_eta_SAMuons_matched_prompt_glob,
        "kmtf_zvtx": kmtf_zvtx_glob,
        "kmtf_zvtx_KMTFTracks_matched":kmtf_zvtx_KMTFTracks_matched_glob,
        "gen_vz_KMTFTracks_matched":gen_vz_KMTFTracks_matched_glob,
        "kmtf_eta_KMTFTracks_matched":kmtf_eta_KMTFTracks_matched_glob,
        "gen_eta_KMTFTracks_allmatched":gen_eta_KMTFTracks_allmatched_glob,
        "kmtf_kslope_KMTFTracks_matched":kmtf_kslope_KMTFTracks_matched_glob,
        "samuon_eta_prompt_matched":samuon_eta_prompt_matched_glob,
        "gen_eta_SAMuons_prompt_allmatched":gen_eta_SAMuons_prompt_allmatched_glob,
        "samuon_z_prompt_matched":samuon_z_prompt_matched_glob,
        "gen_vz_SAMuons_prompt_matched":gen_vz_SAMuons_prompt_matched_glob,
        "samuon_eta_displaced_matched":samuon_eta_displaced_matched_glob,
        "gen_eta_SAMuons_displaced_allmatched":gen_eta_SAMuons_displaced_allmatched_glob,
        "samuon_z_displaced_matched":samuon_z_displaced_matched_glob,
        "gen_vz_SAMuons_displaced_matched":gen_vz_SAMuons_displaced_matched_glob,
        "samuon_hwZ0_prompt_matched":samuon_hwZ0_prompt_matched_glob,
        "samuon_hwZ0_displaced_matched":samuon_hwZ0_displaced_matched_glob,
        "samuon_z_all_prompt":samuon_z_all_prompt_glob,
        "samuon_z_all_displaced":samuon_z_all_displaced_glob,
        "samuon_hwZ0_prompt":samuon_hwZ0_prompt_glob,
        "samuon_hwZ0_displaced":samuon_hwZ0_displaced_glob,
        "samuon_hwD0_prompt":samuon_hwD0_prompt_glob,
        "samuon_hwD0_displaced":samuon_hwD0_displaced_glob,
        "kmtf_dxy": kmtf_dxy_glob,
        "samuon_hwD0_displaced_matched":samuon_hwD0_displaced_matched_glob,
        "gen_pt_SAMuons_displaced_drmatched":gen_pt_SAMuons_displaced_drmatched_glob,

        }
    return return_dict
