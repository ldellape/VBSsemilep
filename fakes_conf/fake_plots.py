import pandas as pd
import numpy as np
import awkward as ak
import glob
import os
import ROOT



# --------------------------------------------------------
# User configuration
# --------------------------------------------------------

# WITHOUT cuts
loose_path = [
    "/eos/user/l/ldellape/VBS/parquet_forfake_MC/Muon1_2023_EraC/Fakes_loose_muon/*.parquet"
]
good_path = [
    "/eos/user/l/ldellape/VBS/parquet_forfake_MC/Muon1_2023_EraC/Fakes_tight_muon/*.parquet"
]
loose_path_MC = [
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/Fakes_loose_muon/*.parquet"
]
good_path_MC = [
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/Fakes_tight_muon/*.parquet"]


vars_to_plot = ["pt", "eta", "phi"]
outdir = "fake_plots_root"
os.makedirs(outdir, exist_ok=True)

MET_THRESHOLD = 30
MT_THRESHOLD  = 20

# --------------------------------------------------------
# Helper: load all parquet files matching a pattern
# --------------------------------------------------------
def load_parquet_collection(patterns, prefix, max_files=10):
    if isinstance(patterns, str):
        patterns = [patterns]

    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))

    files = sorted(files)

    if max_files is not None:
        files = files[:max_files]

    if not files:
        print(f"⚠ No parquet files found in: {patterns}")
        return pd.DataFrame()

    print(f"📂 Loading {len(files)} parquet files for prefix {prefix}")
    jets_field = ["JetForFakes_tight_eta", "JetForFakes_tight_phi"]
    muon_fiels = ["MuonGood_eta", "MuonGood_phi"]
    dfs = []
    for f in files:
        df = pd.read_parquet(f)

        # get muon columns
        cols = [c for c in df.columns if c.startswith(prefix)]

        # also load MET and MT
        if "MET_pt" in df.columns:
            cols.append("MET_pt")
            cols.append("MET_phi")
        if "events_MT_lep_miss" in df.columns:
            cols.append("events_MT_lep_miss")
            cols.append("events_MT_lepLoose_miss")
        if "JetForFakes_tight_eta" in df.columns:
            cols.append("JetForFakes_loose_pt")
            cols.append("JetForFakes_tight_pt")
            cols.append("events_deltaR_jetMuon")

        dfs.append(df[cols])

    return pd.concat(dfs, ignore_index=True)


# --------------------------------------------------------
# Load Loose and Good collections
# --------------------------------------------------------
df_loose = load_parquet_collection(loose_path, "MuonLoose_")
df_good  = load_parquet_collection(good_path,  "MuonGood_")
df_looseMC = load_parquet_collection(loose_path_MC, "MuonLoose_")
df_goodMC = load_parquet_collection(good_path_MC, "MuonGood_")


print(df_looseMC["events_deltaR_jetMuon"])
print(print(df_loose))


df_loose = df_loose[
    (df_loose["MET_pt"] < MET_THRESHOLD) &
    (df_loose["MuonLoose_pt"] > 30) &
    (df_loose["events_MT_lepLoose_miss"] < 30)
]
df_good = df_good[
    (df_good["MET_pt"] < MET_THRESHOLD) &
    (df_good["MuonGood_pt"] > 30)  &
    (df_good["events_MT_lep_miss"] < 30)
]
df_looseMC = df_looseMC[
    (df_looseMC["MET_pt"] < MET_THRESHOLD) &
    (df_looseMC["MuonLoose_pt"] > 30) &
    (df_looseMC["events_MT_lepLoose_miss"] < 30)
]
df_goodMC = df_goodMC[
    (df_goodMC["MET_pt"] < MET_THRESHOLD) &
    (df_goodMC["MuonGood_pt"] > 30)  &
    (df_goodMC["events_MT_lep_miss"] < 30)
]

if df_loose.empty or df_good.empty:
    print("❌ Cannot continue: missing Loose or Good data.")
    exit(0)

# --------------------------------------------------------
# Ratio helper
# --------------------------------------------------------
def compute_ratio(h_good, h_loose):
    h_ratio = h_good.Clone()
    h_ratio.SetName(h_good.GetName() + "_ratio")
    h_ratio.Divide(h_good, h_loose, 1.0, 1.0, "B")
    return h_ratio

# --------------------------------------------------------
# Loop over variables
# --------------------------------------------------------
for var in vars_to_plot:

    loose_field = f"MuonLoose_{var}"
    good_field  = f"MuonGood_{var}"

    if loose_field not in df_loose.columns or good_field not in df_good.columns:
        print(f"⏩ Skipping {var}: missing field")
        continue

    try:
        loose_flat = np.concatenate(df_loose[loose_field].dropna().to_numpy())
        good_flat  = np.concatenate(df_good[good_field].dropna().to_numpy())
        loose_flatMC = np.concatenate(df_looseMC[loose_field].dropna().to_numpy())
        good_flatMC = np.concatenate(df_goodMC[good_field].dropna().to_numpy())
    except Exception as e:
        print(f"⚠ Skipping {var}: {e}")
        continue

    if var == "eta":
        loose_flat = np.abs(loose_flat)
        loose_flatMC = np.abs(loose_flatMC)
        good_flatMC = np.abs(good_flatMC)
        good_flat  = np.abs(good_flat)


    print(f"✔ {var}: Loose={len(loose_flat)}, Good={len(good_flat)}")

    # ----------------------------------------------
    # Histograms
    # ----------------------------------------------
    nbins = 100

    if var == "pt":
        h_loose = ROOT.TH1F(f"h_loose_{var}", f"Loose {var}", nbins, 15, 100)
        h_good  = ROOT.TH1F(f"h_good_{var}",  f"Good {var}",  nbins, 15, 100)
        h_goodMC = ROOT.TH1F(f"h_good_MC_{var}", f"Good_{var}_MC", nbins, 15, 100)
        h_looseMC = ROOT.TH1F(f"h_loose_MC_{var}", f"Good_{var}_MC", nbins, 15, 100)
    else:
        xmin = min(loose_flat.min(), good_flat.min())
        xmax = max(loose_flat.max(), good_flat.max())
        h_loose = ROOT.TH1F(f"h_loose_{var}", f"Loose {var}", nbins, xmin, xmax)
        h_good  = ROOT.TH1F(f"h_good_{var}",  f"Good {var}",  nbins, xmin, xmax)
        h_goodMC = ROOT.TH1F(f"h_good_MC_{var}", f"Good_{var}_MC", nbins, xmin, xmax)
        h_looseMC = ROOT.TH1F(f"h_loose_MC_{var}", f"Good_{var}_MC", nbins, xmin, xmax)

    for x in loose_flat:
        h_loose.Fill(x)
    for x in good_flat:
        h_good.Fill(x)
    for x in loose_flatMC: 
        h_looseMC.Fill(x)
    for x in good_flatMC:
        h_goodMC.Fill(x)

    h_loose.SetLineColor(ROOT.kBlue)
    h_good.SetLineColor(ROOT.kRed)
    h_loose.SetStats(0)
    h_good.SetStats(0)

    # ----------------------------------------------
    # Plot
    # ----------------------------------------------
    c1 = ROOT.TCanvas(f"c_{var}", f"Muon {var}", 700, 600)
    h_loose.Draw("HIST")
    h_good.Draw("HIST SAME")

    legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    legend.AddEntry(h_loose, "Loose", "l")
    legend.AddEntry(h_good, "Good", "l")
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.12, 0.92, "CMS Preliminary")

    c1.SaveAs(f"{outdir}/Muon_{var}.png")
    c1.Close()

    # ----------------------------------------------
    # Ratio
    # ----------------------------------------------
    h_ratioData = compute_ratio(h_good, h_loose)
    h_ratioMC = compute_ratio(h_goodMC, h_looseMC)
    h_ratio = h_ratioData.Clone()
    #h_ratio.Add(h_ratioMC,-1)
    h_ratio.GetYaxis().SetTitle("Fake rate Muon (Run2023C)")
    h_ratio.SetMarkerStyle(20)

    c2 = ROOT.TCanvas(f"c2_{var}", f"Muon {var} Ratio", 700, 500)
    h_ratio.Draw("E1")

    latex.DrawLatex(0.12, 0.92, "CMS Preliminary")
    c2.SaveAs(f"{outdir}/Muon_{var}_ratio.pdf")
    c2.Close()

    # ----------------------------------------------
    # 2D FAKE MAP (pt vs |eta|)
    # ----------------------------------------------
    if var != "pt":
        continue

    loose_eta = np.abs(np.concatenate(df_loose["MuonLoose_eta"].dropna().to_numpy()))
    good_eta  = np.abs(np.concatenate(df_good["MuonGood_eta"].dropna().to_numpy()))


    pt_bins  = np.linspace(30, 100, 10)
    eta_bins = np.linspace(0, 2.5, 10)

    h2_fake = ROOT.TH2F(
        "h2_fake",
        "Fake rate; p_{T} [GeV]; |#eta|",
        len(pt_bins)-1, pt_bins,
        len(eta_bins)-1, eta_bins
    )

    counts_loose = np.zeros((len(pt_bins)-1, len(eta_bins)-1))
    counts_good  = np.zeros_like(counts_loose)

    for pt, eta in zip(loose_flat, loose_eta):
        i = np.searchsorted(pt_bins, pt) - 1
        j = np.searchsorted(eta_bins, eta) - 1
        if 0 <= i < counts_loose.shape[0] and 0 <= j < counts_loose.shape[1]:
            counts_loose[i, j] += 1

    for pt, eta in zip(good_flat, good_eta):
        i = np.searchsorted(pt_bins, pt) - 1
        j = np.searchsorted(eta_bins, eta) - 1
        if 0 <= i < counts_good.shape[0] and 0 <= j < counts_good.shape[1]:
            counts_good[i, j] += 1

    for i in range(counts_good.shape[0]):
        for j in range(counts_good.shape[1]):
            denom = counts_good[i, j] + counts_loose[i, j]
            if denom > 0:
                h2_fake.SetBinContent(i+1, j+1, counts_good[i, j] / denom)

    c3 = ROOT.TCanvas("c3", "Fake rate map", 1500, 1000)
    h2_fake.SetStats(0)
    h2_fake.Draw("COLZ TEXT")

    latex.DrawLatex(0.12, 0.92, "CMS Preliminary, Run2023C")
    c3.SaveAs(f"{outdir}/Muon_pt_eta_fake2D.pdf")
    c3.SaveAs(f"{outdir}/Muon_pt_eta_fake2D.png")
    c3.Close()

print("\n✅ Done — entry-limited ROOT histograms and maps saved.")
