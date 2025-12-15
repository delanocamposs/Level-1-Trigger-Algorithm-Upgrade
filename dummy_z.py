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
R_MB_CM={1:450.0,2:525,3:630,4:700.0}
stations=(1,2,3)
ZRES_CONV=65536.0/1500.0

def event_loop(event_num=100):
    gen_eta_glob, gen_z_glob, gen_pt_glob, gen_vz_glob, gen_vr_glob={st: [] for st in stations}, {st: [] for st in stations}, {st: [] for st in stations}, {st:[] for st in stations}, {st:[] for st in stations}
    stub_eta_glob, stub_z_glob={st: [] for st in stations}, {st: [] for st in stations}
    delta_z_glob={st:[] for st in stations}
    for i, event in enumerate(events):
        if i==event_num+1:
            break
        event.getByLabel("dtTriggerPhase2PrimitiveDigis", "", "L1P2GT", thetahandle)
        container=thetahandle.product().getContainer()
    
        #grab vector of gen muon pT, eta, theta per event
        #also save a gen_z_event which is station-dependent. filled later based on matching (genParticles has no z information - must build the z).
        gen_pt_event=get_gen_muons_pt(event,pt_min=0,pt_max=1000)
        gen_eta_event=get_gen_muons_eta(event,pt_min=0,pt_max=1000)
        gen_vz_event=get_gen_muons_vz(event,pt_min=0,pt_max=1000)
        gen_vy_event=get_gen_muons_vy(event,pt_min=0,pt_max=1000)
        gen_vx_event=get_gen_muons_vx(event,pt_min=0,pt_max=1000)
        gen_vr_event=np.sqrt((np.array(gen_vx_event))**2+(np.array(gen_vy_event))**2)
    
        if len(gen_pt_event)!=len(gen_eta_event)!=len(gen_vz_event):
            print("ERROR: size mismatch")
    
        gen_z_event_by_sta={st: [] for st in stations}
        gen_vz_event_by_sta={st: [] for st in stations}
        gen_vr_event_by_sta={st: [] for st in stations}
        gen_eta_event_by_sta={st: [] for st in stations}
        gen_pt_event_by_sta={st: [] for st in stations} #need this to keep stations and pT organized after matching
    
        #storing per event information about stub hits organized by station. 
        stub_eta_event={st: [] for st in stations} #filled with stub eta values from: station and z value.
        stub_indices_event={} #filled with indices of matched stubs to gen muons.
        stub_z_event={st: [] for st in stations} #filled with stub z value from knowing base z value from theta primitive.  
        stub_eta_matched={st:len(gen_eta_event)*[None] for st in stations} #new stub eta values after matching. IMPORTANT TO KEEP NONE TYPES.
        stub_z_matched={st:len(gen_eta_event)*[None] for st in stations} 
        stub_z_matched_by_sta={st:[] for st in stations} 
        
        delta_z_event={st: [] for st in stations}
    
        #loop through stub hits in each event to obtain:
        #stub eta, stub theta and store in above arrays
        for stub in container:
            st=int(stub.stNum())
            z_stub_cm=float(stub.z()/ZRES_CONV)
    
            stub_z_event[st].append(z_stub_cm)
            eta_stub=eta_from_z(float(stub.z()/ZRES_CONV), R_MB_CM[st])
            stub_eta_event[st].append(eta_stub)
    
        #this loop begins the matching of gen to stub by closest eta. loop through each station and match stub to gen.
        for st in stations:
            match_idx=match_indices_local(gen_eta_event, stub_eta_event[st]) #returns an array the size of gen_eta_event with the best indices of stub eta values that match.
            stub_indices_event[st]=match_idx #throw the best indices per station to their respecitve location.
        
            for mu_idx, stub_idx in enumerate(match_idx): #mu_idx is index of gen muon. stub_idx is best stub index.
                mu_eta=gen_eta_event[mu_idx]            
                mu_vz=gen_vz_event[mu_idx]
                mu_vy=gen_vy_event[mu_idx]
                mu_vx=gen_vx_event[mu_idx]
                mu_vr=gen_vr_event[mu_idx]
                gen_z_val=mu_vz+R_MB_CM[st]/(np.tan(2*np.arctan(np.exp(-1*mu_eta))))
    
                #this must be here to preserve size of gen_event arrays. None types should stay in the case of less stub hits than gen muons in a given station!
                if stub_idx==None:
                    continue
    
                gen_z_event_by_sta[st].append(gen_z_val)
                gen_vz_event_by_sta[st].append(mu_vz)
                gen_vr_event_by_sta[st].append(mu_vr)
                gen_pt_event_by_sta[st].append(gen_pt_event[mu_idx])
                gen_eta_event_by_sta[st].append(gen_eta_event[mu_idx])
    
                stub_eta_matched[st][mu_idx]=stub_eta_event[st][stub_idx]
                stub_z_matched[st][mu_idx]=stub_z_event[st][stub_idx]
                stub_z_matched_by_sta[st].append(stub_z_matched[st][mu_idx])
    
        for st in stations:  
            gen_eta_glob[st].extend(np.array(gen_eta_event_by_sta[st]))
            gen_vz_glob[st].extend(np.array(gen_vz_event_by_sta[st]))
            gen_vr_glob[st].extend(np.array(gen_vr_event_by_sta[st]))
            gen_z_glob[st].extend(np.array(gen_z_event_by_sta[st]))
            gen_pt_glob[st].extend(np.array(gen_pt_event_by_sta[st]))
    
            stub_z_glob[st].extend(np.array(stub_z_matched_by_sta[st]))

            delta_z_glob[st].extend(np.array(gen_z_event_by_sta[st])-np.array(stub_z_matched_by_sta[st]))
    
    print(f"successful event loop. events: {event_num}")
    return_dict={"gen_eta":gen_eta_glob, "gen_pt":gen_pt_glob, "gen_z":gen_z_glob, "stub_z":stub_z_glob, "delta_z":delta_z_glob, "gen_vz":gen_vz_glob, "gen_vr":gen_vr_glob}
    return return_dict
               
#sigma_hists = {}
#
#for st in range(1,4):
#    c_st = ROOT.TCanvas(f"c_{st}", "", 800, 600)
#    h_st = ROOT.TH2F(f"h_{st}", f"Station {st};pT (GeV);#Delta z (cm)",50, 0, 100, 50, -100, 100)
#
#    for n1, n2 in zip(gen_pt_glob[st], delta_z_glob[st]):
#        h_st.Fill(n1, n2)
#
#    h_st.SetStats(0)
#    h_st.Draw("COLZ")
#    c_st.SaveAs(f"delta_z_vs_pT_{st}.png")
#    c_st.Close()
#
#    c_st_2 = ROOT.TCanvas(f"c_{st}_2", "", 800, 600)
#    fit_gaus_st = ROOT.TF1(f"fit_gauss_{st}", "gaus", -100, 100)
#    h_st.FitSlicesY(fit_gaus_st, 1, -1, 0, "QNR")
#    h_sigma_st = ROOT.gDirectory.Get(f"h_{st}_2")
#    h_sigma_st = h_sigma_st.Clone(f"h_sigma_st{st}")
#    h_sigma_st.SetDirectory(0)
#    sigma_hists[st] = h_sigma_st
#    h_sigma_st.SetTitle(f"Station {st};pT (GeV);#sigma [cm]")
#    h_sigma_st.SetStats(0)
#    h_sigma_st.SetMarkerStyle(8)
#    h_sigma_st.SetMarkerSize(1.2)
#    h_sigma_st.SetLineWidth(2)
#    h_sigma_st.Draw("PE")
#    c_st_2.SaveAs(f"sigma_vs_pT_station_{st}.png")
#    c_st_2.Close()
#    
#    for i in range(1, h_st.GetNbinsX()+1):
#        c, h_slice, g = plot_slice_and_gaussian(h_st, i, st, -100, 100)  
#   
#c_all = ROOT.TCanvas("c_all", "Sigma vs pT (all stations)", 800, 600)
#leg = ROOT.TLegend(0.6, 0.7, 0.85, 0.88)
#colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
#markers = [20, 21, 22]
#xmin=10.0
#xmax=50.0
#ymin=0.0
#y_max = max(h.GetMaximum() for h in sigma_hists.values()) * 1.1
#first = True
#
#for i, st in enumerate([1, 2, 3]):
#    h = sigma_hists[st]
#    h.SetMarkerStyle(markers[i])
#    h.SetMarkerSize(1)
#    h.SetMarkerColor(colors[i])
#    h.SetLineColor(colors[i])
#    h.SetLineWidth(2)
#
#    if first:
#        h.SetTitle(";pT (GeV);#sigma [cm]")
#        h.GetXaxis().SetRangeUser(xmin, xmax)
#        h.GetYaxis().SetRangeUser(ymin, y_max)    
#        h.Draw("P")
#        first = False
#    else:
#        h.Draw("P SAME")
#
#    leg.AddEntry(h, f"Station {st}", "PE")
#
#leg.Draw()
#c_all.SaveAs("sigma_vs_pT_all_stations.png") 















#c_ms = ROOT.TCanvas("c_ms", "theta_rms vs pT", 800, 600)
#colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
#graphs = {}
#leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
#
#for i, st in enumerate([1, 2, 3]):
#    if len(gen_pt[st]) == 0:
#        continue
#
#    x = np.array(gen_pt[st], dtype="float64")
#    y = np.array(theta_rms[st], dtype="float64")
#
#    g = ROOT.TGraph(len(x), x, y)
#    graphs[st] = g
#
#    g.SetMarkerStyle(20 + i)
#    g.SetMarkerSize(0.8)
#    g.SetMarkerColor(colors[i])
#    g.SetLineColor(colors[i])
#    g.GetYaxis().SetLimits(0, 0.5)
#    g.GetXaxis().SetLimits(0, 15)
#
#    if st == 1:
#        g.SetTitle(r"$\theta_{\text{RMS}}$ vs $p_T$; pT [GeV]; $\theta_{\text{RMS}}$ [rad]")
#        g.Draw("AP")
#    else:
#        g.Draw("P SAME")
#
#    leg.AddEntry(g, f"Station {st}", "p")
#
#leg.Draw()
#c_ms.SaveAs("theta_rms_vs_pt_all_stations.png")
