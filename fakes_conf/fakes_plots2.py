import pandas as pd
import numpy as np
import awkward as ak
import glob
import ROOT
import os

ROOT.gROOT.SetBatch(True)
MAX_FILES=10

# --------------------------------------------------------
# Input
# --------------------------------------------------------
paths = [
    "/eos/user/l/ldellape/VBS/parquet_forfake_MC/Muon1_2023_EraC/baseline/*.parquet"
]
paths_mc = [
    "/eos/user/l/ldellape/VBS/parquet_forfake_MC/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/baseline/*.parquet",
"/eos/user/l/ldellape/VBS/parquet_forfake_MC/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/baseline/*.parquet"]

outdir = "fake_rate_maps"
os.makedirs(outdir, exist_ok=True)

# --------------------------------------------------------
# Load parquets
# --------------------------------------------------------
# --------------------------------------------------------
# Load parquets (with file limit)
# --------------------------------------------------------
files = []
files_mc = []
for p in paths:
    files.extend(glob.glob(p))
for p in paths_mc:
    files_mc.extend(glob.glob(p))

files = sorted(files)
files_mc = sorted(files_mc)

if MAX_FILES is not None:
    files = files[:MAX_FILES]
    files_mc = files[:MAX_FILES]

if not files:
    raise RuntimeError("No parquet files found")

print(f"📂 Loading {len(files)} parquet files")

df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
df_mc = pd.concat([pd.read_parquet(f) for f in files_mc], ignore_index=True)

print(df["MuonGood_pt"][:2])
# --------------------------------------------------------
# Event selection
# --------------------------------------------------------
event_mask = (
    (df["MET_pt"] < 30) &
    (df["events_MT_lepLoose_miss"] < 30) 
)
event_mask_mc = (
    (df["MET_pt"] <30) & 
    (df["events_MT_lepLoose_miss"] < 30)
)

df = df[event_mask]
df_mc = df_mc[event_mask_mc]
'''
df["selected_muons_pt"] = df.apply(
    lambda row: [pt for pt, dr_list in zip(row["MuonGood_pt"], row["events_deltaR_jetMuon"]) 
                 if any(dr > 1 for dr in dr_list)],
    axis=1
)
'''
# --------------------------------------------------------
# Flatten muons
# --------------------------------------------------------
mu_loose_pt  = ak.flatten(df["MuonLoose_pt"])
mu_loose_eta = ak.flatten(df["MuonLoose_eta"])

mu_good_pt   = ak.flatten(df["MuonGood_pt"])
mu_good_eta  = ak.flatten(df["MuonGood_eta"])

mu_loose_pt  = ak.to_numpy(mu_loose_pt)
mu_loose_eta = np.abs(ak.to_numpy(mu_loose_eta))

mu_good_pt   = ak.to_numpy(mu_good_pt)
mu_good_eta  = np.abs(ak.to_numpy(mu_good_eta))

mu_loose_pt_mc  = ak.flatten(df_mc["MuonLoose_pt"])
mu_loose_eta_mc = ak.flatten(df_mc["MuonLoose_eta"])

mu_good_pt_mc   = ak.flatten(df_mc["MuonGood_pt"])
mu_good_eta_mc  = ak.flatten(df_mc["MuonGood_eta"])

mu_loose_pt_mc  = ak.to_numpy(mu_loose_pt_mc)
mu_loose_eta_mc = np.abs(ak.to_numpy(mu_loose_eta_mc))

mu_good_pt_mc   = ak.to_numpy(mu_good_pt_mc)
mu_good_eta_mc  = np.abs(ak.to_numpy(mu_good_eta_mc))

print(f"✔ Loose muons: {len(mu_loose_pt)}")
print(f"✔ Good  muons: {len(mu_good_pt)}")

# --------------------------------------------------------
# Binning
# --------------------------------------------------------
pt_bins  = np.linspace(30, 100, 20)
eta_bins = np.linspace(0, 2.5, 10)

# --------------------------------------------------------
# Count loose and good
# --------------------------------------------------------
counts_loose, _, _ = np.histogram2d(
    mu_loose_pt, mu_loose_eta,
    bins=[pt_bins, eta_bins]
)

counts_good, _, _ = np.histogram2d(
    mu_good_pt, mu_good_eta,
    bins=[pt_bins, eta_bins]
)


counts_loose_mc, _, _ = np.histogram2d(
    mu_loose_pt_mc, mu_loose_eta_mc,
    bins=[pt_bins, eta_bins]
)

counts_good_mc, _, _ = np.histogram2d(
    mu_good_pt_mc, mu_good_eta_mc,
    bins=[pt_bins, eta_bins]
)

# --------------------------------------------------------
# Fake rate = good / loose
# --------------------------------------------------------
fake_rate = np.zeros_like(counts_loose)
fake_rate_mc = np.zeros_like(counts_loose_mc)

mask = counts_loose > 0
mask_mc = counts_loose_mc > 0
fake_rate[mask] = counts_good[mask] / counts_loose[mask]
fake_rate_mc[mask_mc] = counts_good_mc[mask_mc]/counts_loose_mc[mask_mc]

# --------------------------------------------------------
# ROOT 2D histogram
# --------------------------------------------------------
h2 = ROOT.TH2F(
    "h2_fake_rate",
    "Muon fake rate; p_{T} [GeV]; |#eta|",
    len(pt_bins)-1, pt_bins,
    len(eta_bins)-1, eta_bins
)

for i in range(fake_rate.shape[0]):
    for j in range(fake_rate.shape[1]):
        h2.SetBinContent(i+1, j+1, fake_rate[i,j])

# --------------------------------------------------------
# Plot
# --------------------------------------------------------
c = ROOT.TCanvas("c", "", 800, 600)
h2.SetStats(0)
h2.SetMinimum(0)
h2.SetMaximum(1)
h2.Draw("COLZ TEXT")

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.DrawLatex(0.12, 0.92, "CMS Preliminary")

c.SaveAs(f"{outdir}/Muon_fakeRate_pt_eta.png")
c.SaveAs(f"{outdir}/Muon_fakeRate_pt_eta.pdf")
c.Close()

print("✅ Fake-rate map produced.")
