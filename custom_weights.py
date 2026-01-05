from pocket_coffea.lib.weights.weights import WeightLambda, WeightLambda
from fakes_conf.fakes_rates import Fakes_bkg
import numpy as np
import awkward as ak



######################################################################################
# temporary for fakes 
######################################################################################
Muon_samples = ["Muon1_forFake"]
fakes_mu_weight = WeightLambda.wrap_func(
    name="fakesMU",
    function=lambda params, metadata, events, size, shape_variations:
        apply_mu_fakes(events, metadata["year"], metadata["sample"]),
    has_variations=True,
    isMC_only=False
    )
def apply_mu_fakes(events, year, sample):
    print("weightsss")
    if sample in Muon_samples:
        rates = Fakes_bkg["Muon"][year]
        pt_bins = rates["pt"]
        pt_min  = np.array([b["min"]   for b in pt_bins])
        pt_max  = np.array([b["max"]   for b in pt_bins])
        pt_val  = np.array([b["value"] for b in pt_bins])
        print(f"pt_min: {pt_min}")
        print(f"pt_max: {pt_max}")
        eta_bins = rates["eta"]
        eta_min  = np.array([b["min"]   for b in eta_bins])
        eta_max  = np.array([b["max"]   for b in eta_bins])
        eta_val  = np.array([b["value"] for b in eta_bins])
        mu_pt  = events.MuonLoose.pt
        mu_eta = abs(events.MuonLoose.eta)
        pt_edges = np.array([15,20,30,45,70,100])
        eta_edges = np.array([0.0,0.8,1.479,2.0,2.4])

        pt_idx  = ak.digitize(mu_pt, pt_edges) - 1 
        eta_idx = ak.digitize(mu_eta, eta_edges) - 1

        pt_idx = ak.where((pt_idx>=0) & (pt_idx<len(pt_val)), pt_idx, 0)
        eta_idx = ak.where((eta_idx>=0) & (eta_idx<len(eta_val)), eta_idx, 0)

        fr_pt  = pt_val[pt_idx]
        fr_eta = eta_val[eta_idx]

        fr = fr_pt * fr_eta
        fr = ak.where(fr > 0, fr, 0)
        fake_mu = fr / (1.0 - fr)
        fake_mu_up   = fake_mu * 1.30
        fake_mu_down = fake_mu * 0.70
        n_mu = ak.num(fake_mu)

        w_nom = ak.where(
            n_mu > 0,
            ak.prod(fake_mu, axis=1),
            1.0
        )
        w_up = ak.where(
            n_mu > 0,
            ak.prod(fake_mu_up, axis=1),
            1.0
        )
        w_down = ak.where(
            n_mu > 0,
            ak.prod(fake_mu_down, axis=1),
            1.0
        )
        print(f"weight: {w_nom}")
        return w_nom, w_up, w_down
    else: 
        return np.ones(len(events), dtype=np.float64)
######################################################################################
######################################################################################



######################################################################################
#https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
######################################################################################
SF_pt_top = WeightLambda.wrap_func(
    name="sf_top_pt",
    function=lambda params, metadata, events, size, shape_variations: 
        get_sf_top_pt(events,metadata["sample"]),
    has_variations=False
    )

tt_samples=["TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", 
            "TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8", 
            "TTtoLNu2Q_HT-500_NJet-9_Hdamp-418_TuneCP5_13p6TeV_powheg-pythia8",
            "TTtoLNu2Q_HT-500_NJet-9_TuneCP5_13p6TeV_powheg-pythia8",
            "TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8"
            ]
def get_sf_top_pt(events, sample):
    if sample in tt_samples:
        mask_top = (
            events.GenPart.statusFlags == 10497 & 
            abs(events.GenPart.pdgId) == 6 &
            events.GenPart.status == 62
        )
        top_part = events.GenPart[mask_top]
        top_part = top_part[ak.argsort(top_part.pdgId, ascending=False)]
        
        # primo è il top, second tbar
        weight_t = 0.103*np.exp(-0.0118*top_part.pt[:,0])  -0.000134*top_part.pt[:,0] + 0.973
        weight_tbar = 0.103*np.exp(-0.0118*top_part.pt[:,1]) - 0.000134*top_part.pt[:,1] + 0.973
        weight = np.sqrt(ak.prod([weight_t, weight_tbar], axis=0))
        return weight 
    else: 
        return np.ones(len(events), dtype=np.float64)
        
        