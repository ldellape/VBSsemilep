import pandas as pd
import numpy as np
import glob
import os
import ROOT

# --------------------------------------------------------
# User configuration
# --------------------------------------------------------
# met < 30 e MT < 20
#loose_path = ["/eos/user/l/ldellape/VBS/parquet_forfake/JetMET0_2023_EraC/Fakes_loose_muon/*.parquet", "/eos/user/l/ldellape/VBS/parquet_forfake/JetMET1_2023_EraC/Fakes_loose_muon/*.parquet"]
#good_path  = ["/eos/user/l/ldellape/VBS/parquet_forfake/JetMET1_2023_EraC/Fakes_tight_muon/*.parquet", "/eos/user/l/ldellape/VBS/parquet_forfake/JetMET0_2023_EraC/Fakes_tight_muon/*.parquet"]

#without cuts 
loose_path = ["/eos/user/l/ldellape/VBS/parquet_forfake_nocuts/JetMET0_2023_EraC/Fakes_loose_muon/*.parquet", "/eos/user/l/ldellape/VBS/parquet_forfake_nocuts/JetMET1_2023_EraC/Fakes_loose_muon/*.parquet" ]
good_path = ["/eos/user/l/ldellape/VBS/parquet_forfake_nocuts/JetMET0_2023_EraC/Fakes_tight_muon/*.parquet", "/eos/user/l/ldellape/VBS/parquet_forfake_nocuts/JetMET1_2023_EraC/Fakes_tight_muon/*.parquet" ]

# MET < 30
#loose_path = ["/eos/user/l/ldellape/VBS/parquet_forfake_METonly/JetMET0_2023_EraC/Fakes_loose_muon/*.parquet", "/eos/user/l/ldellape/VBS/parquet_forfake_METonly/JetMET1_2023_EraC/Fakes_loose_muon/*.parquet" ]
#good_path = ["/eos/user/l/ldellape/VBS/parquet_forfake_METonly/JetMET0_2023_EraC/Fakes_tight_muon/*.parquet", "/eos/user/l/ldellape/VBS/parquet_forfake_METonly/JetMET1_2023_EraC/Fakes_tight_muon/*.parquet" ]

# met and mt variables 


Wjets_loose_path = ""
tbarWplus_loose_path = ""
tbarWplus_tight_path = ""
Wjets_tight_path = ""



vars_to_plot = ["pt", "eta", "phi"]
outdir = "fake_plots_root"
os.makedirs(outdir, exist_ok=True)

# --------------------------------------------------------
# Helper: load all parquet files matching a pattern
# --------------------------------------------------------
def load_parquet_collection(patterns, prefix):
    if isinstance(patterns, str):
        patterns = [patterns]

    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))

    files = sorted(files)

    if not files:
        print(f"⚠ No parquet files found in: {patterns}")
        return pd.DataFrame()

    dfs = []
    for f in files:
        df = pd.read_parquet(f)
        # get muon columns
        cols = [c for c in df.columns if c.startswith(prefix)]
        # ALSO load MET_pt
        if "MET_pt" in df.columns:
            cols.append("MET_pt")
        cols.append("events_MT_lep_miss")
        dfs.append(df[cols])

    return pd.concat(dfs, ignore_index=True)

# --------------------------------------------------------
# Load Loose and Good collections
# --------------------------------------------------------
df_loose = load_parquet_collection(loose_path, "MuonLoose_")
df_good  = load_parquet_collection(good_path,  "MuonGood_")
MET_THRESHOLD = 20
MT_THRESHOLD = 20000000000


# Filter events based on MET threshold
df_loose = df_loose[
    (df_loose["MET_pt"] < MET_THRESHOLD) &
    (df_loose["events_MT_lep_miss"] < MT_THRESHOLD)
]
df_good = df_good[
    (df_good["MET_pt"] < MET_THRESHOLD) &
    (df_good["events_MT_lep_miss"] < MT_THRESHOLD)
]


if df_loose.empty or df_good.empty:
    print("❌ Cannot continue: missing Loose or Good data.")
    exit(0)

# --------------------------------------------------------
# Ratio helper using ROOT histogram divide
# --------------------------------------------------------
def compute_ratio(h_good, h_loose):
    h_ratio = h_good.Clone()
    h_ratio.SetName(h_good.GetName() + "_ratio")
    h_ratio.Divide(h_good, h_loose, 1.0, 1.0, "B")  # binomial errors
    return h_ratio

# --------------------------------------------------------
# Loop over variables (1D plots + ratio)
# --------------------------------------------------------
for var in vars_to_plot:

    loose_field = f"MuonLoose_{var}"
    good_field  = f"MuonGood_{var}"

    if loose_field not in df_loose.columns or good_field not in df_good.columns:
        print(f"⏩ Skipping {var}: missing field")
        continue

    # Flatten ragged arrays
    try:
        loose_flat = np.concatenate(df_loose[loose_field].dropna().to_numpy())
        good_flat  = np.concatenate(df_good[good_field].dropna().to_numpy())
    except Exception as e:
        print(f"⚠ Skipping {var}: error flattening → {e}")
        continue

    if var == "eta":
        loose_flat = np.abs(loose_flat)
        good_flat  = np.abs(good_flat)

    print(f"✔ {var}: Loose={len(loose_flat)}, Good={len(good_flat)}")

    # ----------------------------------------------
    # Create ROOT histograms (1D)
    # ----------------------------------------------
    nbins = 5
    xmin = min(loose_flat.min(), good_flat.min())
    xmax = max(loose_flat.max(), good_flat.max())
    if var=="pt":
        h_loose = ROOT.TH1F(f"h_loose_{var}", f"Loose {var}", nbins, 15, 100)
        h_good  = ROOT.TH1F(f"h_good_{var}",  f"Good {var}",  nbins, 15, 100)
    else: 
        h_loose = ROOT.TH1F(f"h_loose_{var}", f"Loose {var}", nbins, xmin, xmax)
        h_good  = ROOT.TH1F(f"h_good_{var}",  f"Good {var}",  nbins, xmin, xmax)
    for x in loose_flat: h_loose.Fill(x)
    for x in good_flat:  h_good.Fill(x)

    h_loose.SetLineColor(ROOT.kBlue)
    h_good.SetLineColor(ROOT.kRed)
    h_loose.SetStats(0)
    h_good.SetStats(0)

    # ----------------------------------------------
    # Plot Loose vs Good (1D)
    # ----------------------------------------------
    c1 = ROOT.TCanvas(f"c_{var}", f"Muon {var}", 700, 600)
    h_loose.Draw("HIST")
    h_good.Draw("HIST SAME")

    legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    legend.AddEntry(h_loose, "Loose", "l")
    legend.AddEntry(h_good, "Good", "l")
    legend.Draw()

    # CMS Preliminary label
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.12, 0.92, "CMS Preliminary")

    c1.SaveAs(f"{outdir}/Muon_{var}.png")
    c1.Close()

    # ----------------------------------------------
    # Ratio histogram (Good / (Good + Loose))
    # ----------------------------------------------
    h_ratio = compute_ratio(h_good, h_loose + h_good)
    h_ratio.SetTitle(f"Muon {var} Ratio")
    h_ratio.GetYaxis().SetTitle("Fake rate Muon (Run2023C)")
    h_ratio.SetMarkerStyle(20)

    c2 = ROOT.TCanvas(f"c2_{var}", f"Muon {var} Ratio", 700, 500)
    h_ratio.Draw("E1")

    # CMS label
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.12, 0.92, "CMS Preliminary")

    c2.SaveAs(f"{outdir}/Muon_{var}_ratio.pdf")
    c2.Close()

    # --------------------------------------------------------
    # 2D FAKE RATE MAP (ONLY FOR var == "pt")
    # --------------------------------------------------------
    if var != "pt":
        continue

    print("👉 Building 2D fake-rate map: pt vs |eta|")

    # Also load eta fields
    loose_eta = np.abs(np.concatenate(df_loose["MuonLoose_eta"].dropna().to_numpy()))
    good_eta  = np.abs(np.concatenate(df_good["MuonGood_eta"].dropna().to_numpy()))

    loose_pt = loose_flat
    good_pt  = good_flat

    # Binning: pt up to 100 GeV, |eta| up to 2.5
    pt_bins  = np.linspace(15, 100, 6)    # 20 bins
    eta_bins = np.linspace(0, 2.5, 5)    # 25 bins

    # Create TH2F
    h2_fake = ROOT.TH2F("h2_fake",
                        "Fake rate, tight/(tight+loose); p_{T} [GeV]; |#eta|",
                        len(pt_bins)-1, pt_bins,
                        len(eta_bins)-1, eta_bins)

    # 2D counters
    counts_loose = np.zeros((len(pt_bins)-1, len(eta_bins)-1))
    counts_good  = np.zeros((len(pt_bins)-1, len(eta_bins)-1))

    # Fill loose
    for pt_val, eta_val in zip(loose_pt, loose_eta):
        i_pt  = np.searchsorted(pt_bins, pt_val) - 1
        i_eta = np.searchsorted(eta_bins, eta_val) - 1
        if 0 <= i_pt < len(pt_bins)-1 and 0 <= i_eta < len(eta_bins)-1:
            counts_loose[i_pt, i_eta] += 1

    # Fill good
    for pt_val, eta_val in zip(good_pt, good_eta):
        i_pt  = np.searchsorted(pt_bins, pt_val) - 1
        i_eta = np.searchsorted(eta_bins, eta_val) - 1
        if 0 <= i_pt < len(pt_bins)-1 and 0 <= i_eta < len(eta_bins)-1:
            counts_good[i_pt, i_eta] += 1

    # Fill TH2F with good/(good+loose)
    for i in range(len(pt_bins)-1):
        for j in range(len(eta_bins)-1):
            denom = counts_good[i,j] + counts_loose[i,j]
            if denom > 0:
                fake = counts_good[i,j] / denom
                h2_fake.SetBinContent(i+1, j+1, fake)

    # Draw 2D map
    c3 = ROOT.TCanvas("c3_pt_eta", "Muon pt vs eta Fake Rate", 1500, 1000)
    h2_fake.SetStats(0)
    h2_fake.Draw("COLZ TEXT")

    # CMS label
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.12, 0.92, "CMS Preliminary, Run2023C")

    c3.SaveAs(f"{outdir}/Muon_pt_eta_fake2D.pdf")
    c3.Close()

print("\n✅ Done — ROOT histograms and 2D maps saved.")
