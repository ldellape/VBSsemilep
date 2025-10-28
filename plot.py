import awkward as ak
import numpy as np
import os
import glob
import ROOT
from sample import samples
path = "/eos/user/l/ldellape/VBS/parquet"
norm_luminosity=False
target_lumi = 140
norm_unity = True
polarized=True


tag="WW"
#tag="ZZ"


#category="baseline"
#category="signal_AK8"
category="signal_AK4"
#category="side"
#category="CR_L_Wjets"
#category="CR_R_Wjets"



parquet_patterns, labels = samples(category, tag, polarized)

def load_awkward_parquet(pattern):
    files = glob.glob(pattern)
    arrays = [ak.from_parquet(f) for f in files]
    if len(arrays) == 0:
        print(f"Warning: no files found for pattern {pattern}")
        return None
    return ak.concatenate(arrays)


def sum_lumi(tag):
    if tag == "ssWW_LL":
        log_path = "log_ssWW_LL_mg5_madspin.txt"
    elif tag=="ssWW_TT":
        log_path = "log_ssWW_TT_mg5_madspin.txt"
    elif tag == "TTtoLNu2Q":
        log_path="log_TTtoLNu2Q.txt"
    elif tag == "WLtoLNu_2Jets":
        log_path ="log_WtoLNu-XJets.txt"
    elif tag == "TTto2L2Nu-2Jets":
        log_path = "log_TTto2L2Nu.txt"
    elif tag == "DYto2L-4Jets":
        log_path = "log_DY.txt"
    elif tag == "ZZ_LL":
        log_path = "log_ZZLL_mg5_madspin.txt"
    elif tag == "ZZ_TT":
        log_path = "log_ZZTT_mg5_madspin.txt"
    elif tag == "ZZZ":
        log_path = "log_ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix.txt"
    elif tag == "WZZ":
        log_path = "log_WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix.txt"
    else: 
        return 0.0
    total_lumi = 0.0

    if not os.path.exists(log_path):
        print(f"No log file found for dataset {log_path}")
        return 0.0

    with open(log_path, "r") as f:
        for line in f:
            parts = line.strip().split(',')
            for part in parts:
                if part.strip().startswith("lumi="):
                    lumi_value = float(part.split('=')[1])
                    total_lumi += lumi_value
    print(total_lumi)
    return total_lumi*1e-03

datasets = [load_awkward_parquet(p) for p in parquet_patterns]



lumis = [sum_lumi(label) for label in labels]

lumis = np.array(lumis, dtype=float)
print(lumis)
scaling_factors = target_lumi / lumis
for label, lumi, scale in zip(labels, lumis, scaling_factors):
    print(f"{label:<10} Lumi: {lumi:.3f}  Scale factor: {scale:.6f}")

c = ROOT.TCanvas("c", "Comparison", 800, 600)
ROOT.gStyle.SetOptStat(0)

field_name = ["lumi", "CleanFatJet_pt", "CleanJet_pt", "CleanFatJet_eta", "CleanFatJet_msoftdrop", "CleanFatJet_mass", "events_zepp_lep", "CleanSubJet_pair_zg", 
              "CleanFatJet_tau21", "events_nCleanJets", "events_nElectronGood", "events_nMuonGood", "events_nFatJet", "events_nCleanFatJets", "events_nBJetGood", "ElectronGood_pt", "MuonGood_pt", "V_dijet_candidate_mass", "V_dijet_candidate_pt",
              "VBS_dijet_system_pt", "VBS_dijet_system_mass", "VBS_dijet_system_deltaEta", "events_MT_ele_miss", "events_MT_mu_miss", "events_MT_lep_miss"]

number_events = []
for (data,label) in zip(datasets, labels):
    number_events.append(len(data))  
    print(f"{len(data)} events in {label}")  

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleFontSize(0.04)
ROOT.gStyle.SetLabelSize(0.04, "XY")
ROOT.gStyle.SetTitleSize(0.045, "XY")
ROOT.gStyle.SetLegendBorderSize(0)

for var in field_name:
    hists = {}

    global_min = None
    global_max = None
    is_discrete = var in [
        "events_nCleanJet",
        "events_nCleanFatJets",
        "events_nFatJet",
        "events_nJet",
        "events_nBJetGood",
    ]

    # --- find global min/max for continuous variables ---
    for data in datasets:
        if var in data.fields:
            vals = ak.flatten(data[var], axis=None)
            if len(vals) > 0 and not is_discrete:
                vmin = float(min(vals))
                vmax = float(max(vals))
                global_min = vmin if global_min is None else min(global_min, vmin)
                global_max = vmax if global_max is None else max(global_max, vmax)

    # --- build histograms ---
    for i, (data, label) in enumerate(zip(datasets, labels)):
        if var not in data.fields:
            continue

        vals = ak.flatten(data[var], axis=None)
        if len(vals) == 0:
            print(f"[WARNING] No entries for variable '{var}' in dataset '{label}'")
            continue

        if not is_discrete:
            if "miss" in var:
                h = ROOT.TH1F(f"h_{var}_{i}", label, 50, 0, global_max)
            else: 
                h = ROOT.TH1F(f"h_{var}_{i}", label, 50, global_min, global_max)
        else:
            h = ROOT.TH1F(f"h_{var}_{i}", label, 8, 0, 8)

        for v in vals:
            h.Fill(v)

        if norm_luminosity:
            h.Scale(scaling_factors[i])

        if label in hists:
            hists[label].Add(h)
        else:
            h.SetLineWidth(2)
            hists[label] = h
        if norm_unity:
            sum_content = sum(h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1))
            h.Scale(1./sum_content)
        print(var)
       

    if not hists:
        continue
    

    # --- Overlay plot ---
    c_overlay = ROOT.TCanvas(f"c_{var}_overlay", f"Overlay: {var}", 800, 600)
    legend_overlay = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    sorted_hists = sorted(hists.items(), key=lambda kv: kv[1].GetMaximum(), reverse=True)

    for idx, (label, h) in enumerate(sorted_hists):
        h.SetLineColor(idx + 1)
        draw_option = "hist" if idx == 0 else "hist same"
        h.Draw(draw_option)
        if idx == 0:
            h.GetXaxis().SetTitle(var)
            h.GetYaxis().SetTitle("Entries")
            h.SetTitle(f"{var}   {category}")
        if var== "VBS_dijet_system_mass":
            h.GetXaxis().SetTitle("M_{jj} [vbs]")
        elif var == "VBS_dijet_system_deltaEta":
            h.GetXaxis().SetTitle("|#Delta #eta| [vbs]")
            c_overlay.SetLogy()
        elif var == "events_zepp_lep":
            h.GetXaxis().SetTitle("zepp")
            c_overlay.SetLogy()

        elif var == "MuonGood_eta":
            h.GetXaxis().SetTitle("lepton #eta")
            c_overlay.SetLogy()
        legend_overlay.AddEntry(h, label, "f")

    legend_overlay.Draw()
    os.makedirs(f"./plots/{tag}/{category}", exist_ok=True)
    c_overlay.SaveAs(f"./plots/{tag}/{category}/{var}_{category}_overlay.pdf")

    # --- Stack plot ---
    c_stack = ROOT.TCanvas(f"c_{var}_stack", f"Stack: {var}", 800, 600)
    stack = ROOT.THStack(f"hs_{var}", f"{var}   {category}")
    legend_stack = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend_stack.SetFillStyle(0)

    for idx, (label, h) in enumerate(hists.items()):
        h.SetFillColorAlpha(idx + 1, 0.4)
        stack.Add(h)
        legend_stack.AddEntry(h, label, "f")

    stack.Draw("hist")
    stack.GetXaxis().SetTitle(var)
    stack.GetYaxis().SetTitle("Entries")
    stack.SetMinimum(0)
    stack.SetMaximum(stack.GetMaximum() * 1.2)
    legend_stack.Draw()

    c_stack.SaveAs(f"./plots/{tag}/{category}/{var}_{category}_stack.pdf")
