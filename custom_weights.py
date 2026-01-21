from pocket_coffea.lib.weights.weights import WeightLambda, WeightLambda
from parameters.fakes_rates import Fakes_bkg
import numpy as np
import awkward as ak

######################################################################################
# temporary for fakes                                                                #                                       
######################################################################################
Muon_samples = ["Muon1_forFake", "SingleMuon_forFake"]
Electron_samples = ["EGamma_forFake", "SingleElectron_forFake"]

fakesMuon = WeightLambda.wrap_func(
    name="fakesMU",
    function=lambda params, metadata, events, size, shape_variations:
        apply_fakes(events,metadata["year"], metadata["sample"], "Muon"),
    isMC_only=False,
    has_variations=True
    )
fakesEle = WeightLambda.wrap_func(
    name="fakesEle",
    function=lambda params, metadata, events, size, shape_variations:
        apply_fakes(events,metadata["year"], metadata["sample"], "Electron"),
        isMC_only=False,
        has_variations=True
    )
def apply_fakes(events, year, sample, particle):
    if particle == "Muon":
        if sample not in Muon_samples:
            return np.ones(len(events)), np.zeros(len(events)), np.zeros(len(events))
    elif particle == "Electron":
        if sample not in Electron_samples: 
            return np.ones(len(events)), np.zeros(len(events)), np.zeros(len(events))

    if particle == "Muon":
        lepton = "Muon"
        lep_pt  = events.MuonLoose.pt
        lep_eta = abs(events.MuonLoose.eta)
    else:
        lepton = "Electron"
        lep_pt = events.ElectronLoose.pt
        lep_eta = abs(events.ElectronLoose.eta) 
        

    fr_pt  = ak.zeros_like(lep_pt, dtype=np.float64)
    fr_eta = ak.zeros_like(lep_eta, dtype=np.float64)

    for b in Fakes_bkg[lepton][year]["pt"]:
        lo, hi, val = b["min"], b["max"], b["value"]
        fr_pt = ak.where((lep_pt >= lo) & (lep_pt < hi), val, fr_pt)

    for b in Fakes_bkg[lepton][year]["eta"]:
        lo, hi, val = b["min"], b["max"], b["value"]
        fr_eta = ak.where((lep_eta >= lo) & (lep_eta < hi), val, fr_eta)

    fr = fr_pt * fr_eta
    fr = ak.where(fr > 0, fr, 0)
    print(fr)
    fake_lep = fr / (1.0 - fr)
    fake_lep_up   = fake_lep* 1.30
    fake_lep_down = fake_lep* 0.70

    n_lep = ak.num(fake_lep)

    w_nom = ak.where(n_lep > 0, ak.prod(fake_lep, axis=1), 1.0)
    w_up  = ak.where(n_lep > 0, ak.prod(fake_lep_up, axis=1), 1.0)
    w_dn  = ak.where(n_lep > 0, ak.prod(fake_lep_down, axis=1), 1.0)
    print(n_lep)
    print(f"w nom: {w_nom}")
    return w_nom, w_up, w_dn
######################################################################################
######################################################################################



######################################################################################
#https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting                           #
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
            ]
ST_samples = ["TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8", 
              "TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
              "TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8"
            ]
top_samples = tt_samples + ST_samples
def get_sf_top_pt(events, sample):
    if sample in top_samples:
        print(f"statusFlags: {events.GenPart.statusFlags}")
        is_lastcopy = (events.GenPart.statusFlags & (1 << 13)) != 0

        mask_top = (
            is_lastcopy & 
            (abs(events.GenPart.pdgId) == 6) &
            (events.GenPart.status == 62)
        )
        top_part = events.GenPart[mask_top]
        
        # primo è il top, second tbar
        if sample in tt_samples:
            top_part = top_part[ak.argsort(top_part.pdgId, ascending=False)]
            weight_t = 0.103*np.exp(-0.0118*top_part.pt[:,0])  -0.000134*top_part.pt[:,0] + 0.973
            weight_tbar = 0.103*np.exp(-0.0118*top_part.pt[:,1]) - 0.000134*top_part.pt[:,1] + 0.973
            weight = np.sqrt(ak.prod([weight_t, weight_tbar], axis=0))
            print(weight)
            return weight
        elif sample in ST_samples: 
            weight_t = 0.103*np.exp(-0.0118*top_part.pt)-0.000134*top_part.pt + 0.973
            weight = ak.prod(weight_t, axis=1)
            print(weight)
            return weight
    else: 
        return np.ones(len(events), dtype=np.float64)
        
        