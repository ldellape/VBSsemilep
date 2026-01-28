import awkward as ak
import math
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq, apply_golden_json,get_JetVetoMap
from pocket_coffea.lib.objects import get_dilepton





#############################################################
# cut functions for fakes                                   #
#############################################################
def closure_test_non_prompt_deepMETresponse(events, params, year, sample, **kwargs):
    events.CleanJet_forNonPrompt = events.CleanJet[events.CleanJet.pt > 30]
    print(f"closure: {events.CleanJet_forNonPrompt.pt}")
    events.nCleanJet_forNonPrompt = ak.num(events.CleanJet_forNonPrompt)
    print(f"cleanjets: {events.CleanJet}")
    print(f"puppiMET: {events.PuppiMET.pt}")
    print(f"muon pt: {events.MuonGood.pt}")
    mask=( 
          (events.nCleanJet_forNonPrompt >= params["nJet"]) &
          (events.MET.pt < params["met"]) & 
          (ak.firsts(events.MuonGood.pt) > params["pt"])
        )
    print(f"mask closure: {mask}")
    return ak.where(ak.is_none(mask), False, mask)
closure_test_deepMETresponse = Cut(
    name="closure_test_1",
    params={
        "nJet" : 4,
        "met": 30,
        "pt" : 30,
    },
    function=closure_test_non_prompt_deepMETresponse,
)
def closure_test_non_prompt_deepMETresolution(events, params, year, sample, **kwargs):
    events.CleanJet_forNonPrompt = events.CleanJet[events.CleanJet.pt > 30]
    print(f"closure: {events.CleanJet_forNonPrompt.pt}")
    events.nCleanJet_forNonPrompt = ak.num(events.CleanJet_forNonPrompt)
    print(f"cleanjets: {events.CleanJet}")
    print(f"puppiMET: {events.PuppiMET.pt}")
    print(f"muon pt: {events.MuonGood.pt}")
    mask=( 
          (events.nCleanJet_forNonPrompt >= params["nJet"]) &
          (events.MET.pt < params["met"]) & 
          (ak.firsts(events.MuonGood.pt) > params["pt"])
        )
    print(f"mask closure: {mask}")
    return ak.where(ak.is_none(mask), False, mask)
closure_test_deepMETres = Cut(
    name="closure_test_1",
    params={
        "nJet" : 4,
        "met": 30,
        "pt" : 30,
    },
    function=closure_test_non_prompt_deepMETresolution,
)

def closure_test_non_prompt_puppiMET(events, params, year, sample, **kwargs):
    events.CleanJet_forNonPrompt = events.CleanJet[events.CleanJet.pt > 30]
    print(f"closure: {events.CleanJet_forNonPrompt.pt}")
    events.nCleanJet_forNonPrompt = ak.num(events.CleanJet_forNonPrompt)
    print(f"cleanjets: {events.CleanJet}")
    print(f"puppiMET: {events.PuppiMET.pt}")
    print(f"muon pt: {events.MuonGood.pt}")
    mask=( 
          (events.nCleanJet_forNonPrompt >= params["nJet"]) &
          (events.MET.pt < params["met"]) & 
          (ak.firsts(events.MuonGood.pt) > params["pt"])
        )
    print(f"mask closure: {mask}")
    return ak.where(ak.is_none(mask), False, mask)
closure_test_1 = Cut(
    name="closure_test_1",
    params={
        "nJet" : 4,
        "met": 30,
        "pt" : 30,
    },
    function=closure_test_non_prompt_puppiMET,
)




def met_for_fakes_loose(events, params, year, sample, **kwargs):
    mask = (
        (events.PuppiMET.pt < params["met"]) & 
        (events.MT_lep_miss < params["mt"])
    )
    return ak.where(ak.is_none(mask), False, mask)
FakeLoose = Cut(
    name="FakeLoose",
    params={
        "met" : 30,
        "mt" : 1000000000000,
        },
    function=met_for_fakes_loose,
)
def jet_for_fakes_loose(events, params, year, sample, **kwargs):
    single_jet = events.nJetForFakes_loose == 1
    mask = single_jet
    return ak.where(ak.is_none(mask), False, mask)
JetForFakes_loose = Cut(
    name="JetForFakes_loose",
    params={},
    function=jet_for_fakes_loose,
)
def jet_for_fakes_tight(events, params, year, sample, **kwargs):
    single_jet = events.nJetForFakes_tight == 1
    mask = single_jet
    print(mask)
    return ak.where(ak.is_none(mask), False, mask)
JetForFakes_tight = Cut(
    name="JetForFakes_tight",
    params={},
    function=jet_for_fakes_tight,
)
def met_for_fakes_tight(events, params, year, sample, **kwargs):
    single_lepton = events.nLeptonGood == 1
    mask = (
        single_lepton &
        (events.MET.pt < params["met"]) & 
        (events.MT_lep_miss < params["mt"])
    )
    return ak.where(ak.is_none(mask), False, mask)
FakeTight = Cut(
    name="FakeLoose",
    params={
        "met" : 20,
        "mt" : 20,
        },
    function=met_for_fakes_tight,
)
#############################################################
#############################################################




#############################################################
# single good lepton (electron or muon)                     #
#############################################################
def single_good_lepton(events, params, year, sample, **kwargs):
    single_lepton = events.nLeptonGood == 1
    mask = single_lepton
    return ak.where(ak.is_none(mask), False, mask)   
SingleLepton = Cut(
    name="SingleGoodLeptons",
    params={},
    function=single_good_lepton,
)
#############################################################
#############################################################




#############################################################
# VBS TOPOLOGY                                              #
#############################################################
def VBS_topology(events, params, year, sample, **kwargs):
    mask = (
        (events.VBS_dijet_system.mass > params["mass"])
        & 
        (events.VBS_dijet_system.deltaEta > params["deltaEta"])
        &
        (events.VBS_dijet_system.pt1 > params["pt_leadingJet"])
        &
        (
            (
                (events.nCleanFatJets == 1) &
                (events.nCleanJet >= params["nJet_with_FatJet"])
            )
            |
            (
                (events.nCleanJet >= params["nJet"]) &
                (events.nCleanFatJets == 0)
            )
        )
    )
    return ak.where(ak.is_none(mask), False, mask)
VBS_jets_presel = Cut(
    name="VBS_jets_presel",
    params = {
        "mass" : 500,
        "deltaEta" : 2.5,
        "pt_leadingJet" : 50,
        "nJet" : 4,
        "nJet_with_FatJet" : 2,
    },
    function=VBS_topology,
)
#############################################################
#############################################################



#############################################################
# semileptonic, requirements on leptons and/or MET          #
#############################################################
def semileptonic_deepMETresponse(events, params, year, sample, **kwargs):
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    single_lepton = events.nLeptonGood == 1
    double_electron = events.nElectronGood == 2
    double_muon = events.nMuonGood == 2
    print(f" lepton good: {events.LeptonGood.pt}")
    print(f" electron good: {events.ElectronGood.pt}")
    if params["W"] is True:
        mask = (
                ( single_electron
                  & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_electron"])
                  & (events.DeepMETResponseTune.pt > params["met_electron"]) 
                )
                | 
                (  single_muon
                   & (ak.firsts(events.MuonGood.pt) > params["pt_leading_muon"])
                   & (events.DeepMETResponseTune.pt > params["met_muon"])
                )
        )
    elif params["Z"] is True:
        muon_pt = ak.pad_none(events.MuonGood.pt, 2)
        muon_charge = ak.pad_none(events.MuonGood.charge, 2)
        electron_pt = ak.pad_none(events.ElectronGood.pt, 2)
        electron_charge = ak.pad_none(events.ElectronGood.charge, 2)
        mu_mask = (
                    (muon_pt[:, 0] > params["pt_leading_muon"]) &
                    (muon_pt[:, 1] > params["pt_subleading_muon"]) &
                    (muon_charge[:,0] != muon_charge[:, 1])
        )
        mu_mask = ak.fill_none(mu_mask, False)
        el_mask = (
            (electron_pt[:, 0] > params["pt_leading_electron"]) &
            (electron_pt[:, 1] > params["pt_subleading_electron"]) &
            (electron_charge[:,0] != electron_charge[:,1])
        )
        el_mask = ak.fill_none(el_mask, False)
        mask = mu_mask | el_mask
    return ak.where(ak.is_none(mask), False, mask)
def semileptonic_deepMETresolution(events, params, year, sample, **kwargs):
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    single_lepton = events.nLeptonGood == 1
    double_electron = events.nElectronGood == 2
    double_muon = events.nMuonGood == 2
    print(f" lepton good: {events.LeptonGood.pt}")
    print(f" electron good: {events.ElectronGood.pt}")
    if params["W"] is True:
        mask = (
                ( single_electron
                  & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_electron"])
                  & (events.DeepMETResolutionTune.pt > params["met_electron"]) 
                )
                | 
                (  single_muon
                   & (ak.firsts(events.MuonGood.pt) > params["pt_leading_muon"])
                   & (events.DeepMETResolutionTune.pt > params["met_muon"])
                )
        )
    elif params["Z"] is True:
        muon_pt = ak.pad_none(events.MuonGood.pt, 2)
        muon_charge = ak.pad_none(events.MuonGood.charge, 2)
        electron_pt = ak.pad_none(events.ElectronGood.pt, 2)
        electron_charge = ak.pad_none(events.ElectronGood.charge, 2)
        mu_mask = (
                    (muon_pt[:, 0] > params["pt_leading_muon"]) &
                    (muon_pt[:, 1] > params["pt_subleading_muon"]) &
                    (muon_charge[:,0] != muon_charge[:, 1])
        )
        mu_mask = ak.fill_none(mu_mask, False)
        el_mask = (
            (electron_pt[:, 0] > params["pt_leading_electron"]) &
            (electron_pt[:, 1] > params["pt_subleading_electron"]) &
            (electron_charge[:,0] != electron_charge[:,1])
        )
        el_mask = ak.fill_none(el_mask, False)
        mask = mu_mask | el_mask
    return ak.where(ak.is_none(mask), False, mask)
def semileptonic(events, params, year, sample, **kwargs):
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    single_lepton = events.nLeptonGood == 1
    double_electron = events.nElectronGood == 2
    double_muon = events.nMuonGood == 2
    print(f" lepton good: {events.LeptonGood.pt}")
    print(f" electron good: {events.ElectronGood.pt}")
    if params["W"] is True:
        mask = (
                ( single_electron
                  & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_electron"])
                  & (events.MET.pt > params["met_electron"]) 
                )
                | 
                (  single_muon
                   & (ak.firsts(events.MuonLoose.pt) > params["pt_leading_muon"])
                   & (events.MET.pt > params["met_muon"])
                )
        )
    elif params["Z"] is True:
        muon_pt = ak.pad_none(events.MuonGood.pt, 2)
        muon_charge = ak.pad_none(events.MuonGood.charge, 2)
        electron_pt = ak.pad_none(events.ElectronGood.pt, 2)
        electron_charge = ak.pad_none(events.ElectronGood.charge, 2)
        mu_mask = (
                    (muon_pt[:, 0] > params["pt_leading_muon"]) &
                    (muon_pt[:, 1] > params["pt_subleading_muon"]) &
                    (muon_charge[:,0] != muon_charge[:, 1])
        )
        mu_mask = ak.fill_none(mu_mask, False)
        el_mask = (
            (electron_pt[:, 0] > params["pt_leading_electron"]) &
            (electron_pt[:, 1] > params["pt_subleading_electron"]) &
            (electron_charge[:,0] != electron_charge[:,1])
        )
        el_mask = ak.fill_none(el_mask, False)
        mask = mu_mask | el_mask
    return ak.where(ak.is_none(mask), False, mask)
semileptonic_preselW_puppiMET = Cut(
    name="semileptonic_preselW", 
    params = {
        "W" : True,
        "pt_leading_electron" : 30,
        "pt_leading_muon" : 30,
        "eta_max_lep" : 2.5, 
        "nJet" : 4, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 2,
        "met_electron" : 30,
        "met_muon" : 30,
        "max_mt":185,
    },
    function = semileptonic,
) 
semileptonic_preselW_deepMETresolution= Cut(
    name="semileptonic_preselW_deepMETresolution", 
    params = {
        "W" : True,
        "pt_leading_electron" : 30,
        "pt_leading_muon" : 30,
        "eta_max_lep" : 2.5, 
        "nJet" : 4, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 2,
        "met_electron" : 30,
        "met_muon" : 30,
        "max_mt":185,
    },
    function = semileptonic_deepMETresolution,
) 
semileptonic_preselW_deepMETresponse = Cut(
    name="semileptonic_preselW_deepMETresponse", 
    params = {
        "W" : True,
        "pt_leading_electron" : 30,
        "pt_leading_muon" : 30,
        "eta_max_lep" : 2.5, 
        "nJet" : 4, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 2,
        "met_electron" : 30,
        "met_muon" : 30,
        "max_mt":185,
    },
    function = semileptonic_deepMETresponse,    
)
semileptonic_preselW = Cut(
    name="semileptonic_preselW", 
    params = {
        "W" : True,
        "pt_leading_electron" : 30,
        "pt_leading_muon" : 35,
        "eta_max_lep" : 2.5, 
        "nJet" : 4, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 2,
        "met_electron" : 30,
        "met_muon" : 30,
        "max_mt": 185,
    },
    function = semileptonic,
) 
semileptonic_preselW_3 = Cut(
    name="semileptonic_preselW", 
    params = {
        "W" : True,
        "pt_leading_electron" : 0,
        "pt_leading_muon" : 0,
        "eta_max_lep" : 2.5, 
        "nJet" : 4, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 2,
        "met_electron" : 30,
        "met_muon" : 0,
        "max_mt":185,
    },
    function = semileptonic,
) 


semileptonic_preselZ = Cut(
    name="semileptonic_preselZ", 
    params = {
        "W" : False,
        "Z" : True,
        "pt_leading_electron" : 20,
        "pt_subleading_electron" : 20,
        "pt_leading_muon" : 20,
        "pt_subleading_muon" : 20,
        "nJet" : 1, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 1,
        "met" : 20,
    },
    function = semileptonic,
) 
semileptonic_preselDY = Cut(
    name="semileptonic_preselDY",
    params = {
        "pt_leading_electron" : 20,
        "pt_subleading_electron" : 20,
        "pt_leading_muon" : 20,
        "pt_subleading_muon" : 20,
        "nJet" : 1, 
        "nFatJets" : 1,
        "nJet_with_FatJet" : 1,
        "met" : 20,        
    },
    function = semileptonic,
)
#############################################################
#############################################################



#############################################################
# V jet mass from AK8 or 2AK4 (pair wuth the mass closer to W/Z mass)
#############################################################
def Vjet_mass_AK8(events, params, year, sample, **kwargs):
    fj_mask = (
        (events.CleanFatJet.msoftdrop > params["mass_min"]) & 
        (events.CleanFatJet.msoftdrop < params["mass_max"]) & 
        (events.CleanFatJet.pt > params["pt_min"])
    )
    mask = ak.any(fj_mask, axis=1)
    return ak.where(ak.is_none(mask), False, mask)
def Vjet_mass_AK4(events, params, year, sample, **kwargs):
    if "V_dijet_candidate" not in events.fields:
        return ak.zeros_like(events.event, dtype=bool)
    fj_mask = (
        (events.nCleanFatJets == 0) &
        (events.V_dijet_candidate.mass > params["mass_min"]) &
        (events.V_dijet_candidate.mass < params["mass_max"]) &
        (events.V_dijet_candidate.pt < params["pt_max"])
    )
    return ak.where(ak.is_none(fj_mask), False, fj_mask)

Vjet_massZ_boosted = Cut(
    name="Vjet_massZ",
    params={
        "VV": True,
        "mass_min": 70,
        "mass_max": 110,
        "nJet_min" : 4,
        "tau21":0.45,
    },
    function=Vjet_mass_AK8,
)
Vjet_massZ_resolved = Cut(
    name="Vjet_massZ",
    params={
        "VV": True,
        "mass_min": 70,
        "mass_max": 110,
        "nJet_min" : 4,
    },
    function=Vjet_mass_AK8,
)
Vjet_massW_resolved = Cut(
    name="Vjet_massW",
    params={
        "VV": True,
        "mass_min": 65,
        "mass_max": 105,
        "pt_max" : 200,
    },
    function=Vjet_mass_AK4,
)
Vjet_massW_boosted = Cut(
    name="Vjet_massW",
    params={
        "VV": True,
        "mass_min": 70,
        "mass_max": 115,
        "pt_min" : 200,
    },
    function=Vjet_mass_AK8,
)
#############################################################
#############################################################



#############################################################
# number of b-tagged (central) jets (for ttbar etc.)        #
#############################################################
def Bjets_presel(events, params, year, sample, **kwargs):
    b_mask = (events["nBJetGood"] >= params["nBjets"])
    return ak.where(ak.is_none(b_mask), False, b_mask)
Bjets_presel = Cut(
    name="Bjets_presel",
    params={
        "nBjets" : 1,
    },
    function=Bjets_presel,
)
#############################################################
#############################################################



#############################################################
# jet mass-side W/Z
#############################################################

def Vjet_massSide_resolved(events, params, year, sample, **kwargs):
    if params["L"] is True: 
        mask = (
            (events.V_dijet_candidate.mass < params["mass_min"])
        )
    else:
        mask = (
            (events.V_dijet_candidate.mass > params["mass_max"])
        )
    return ak.where(ak.is_none(mask), False, mask)


def Vjet_massSide_boosted_ALL(events, params, year, sample, **kwargs):
    mask = (
        events.CleanFatJet.msoftdrop < params["mass_min"] & 
        events.CleanFatJet.msoftdrop > params["mass_max"]
    )
    return ak.where(ak.is_none(mask), False, mask)
Wjet_side_ALL = Cut(
    name="Wjet_side_ALL",
    params={
        "mass_min" : 65,
        "mass_max" : 105,
    },
    function=Vjet_massSide_boosted_ALL,
)
    



def Vjet_massSide_boosted(events, params, year, sample, **kwargs):
    if params["L"] is True:
        mask1 = (
            (events.CleanFatJet.msoftdrop < params["mass_min"])
        )
    else: 
        mask1 = (
            (events.CleanFatJet.msoftdrop > params["mass_max"])
        )
    mask = ak.any(mask1, axis=1)
    return ak.where(ak.is_none(mask), False, mask) 
Wjet_sideL_resolved = Cut(
    name="Vjet_sideL_resolved",
    params={
        "L" : True,
        "mass_min" : 60,
    },
    function=Vjet_massSide_resolved,    
)
Wjet_sideR_resolved = Cut(
    name="Vjet_sideL_resolved",
    params={
        "L" : False,
        "mass_max" : 100,
        "mass_min" : 60,
    },
    function=Vjet_massSide_resolved,    
)
Wjet_sideL_boosted = Cut(
    name="Vjet_sideL_boosted",
    params={
        "L" : True,
        "mass_min" : 60,
        "mass_max" : 100,
    },
    function=Vjet_massSide_resolved,    
)
Wjet_sideR_boosted = Cut(
    name="Vjet_sideL_boosted",
    params={
        "L" : False,
        "mass_max" : 100,
        "mass_min" : 60,
    },
    function=Vjet_massSide_boosted,    
)
#############################################################
#############################################################



#############################################################
# transverse mass requirements form W->lv                   #
#############################################################
def Wtransverse_mass(events, params, year, sample, **kwargs):
    mask = (
        (events.MT_lep_miss < params["transverse_max"]) &
        (events.MT_lep_miss > 0)
    )
    return ak.where(ak.is_none(mask), False, mask)
Wtransverse_mass_presel = Cut(
    name="Wtransverse_mass_presel",
    params={
        "transverse_max" : 185,
    },
    function=Wtransverse_mass,
)
Wtransverse_forFakes = Cut(
    name="Wtransverse_forFakes",
    params={"transverse_max" : 30},
    function=Wtransverse_mass,
)
#############################################################
#############################################################





#############################################################
# Z->ll dilepton mass                                       #
#############################################################
def dilepton_mass(events, params, year, sample, **kwargs):
    double_muon = events.nMuonGood == 2
    double_electron = events.nElectronGood == 2
    if params["VV"] is True:
        mask = (
            ( 
                (double_muon) &
                (events.dimuon_candidate.mass > params["mass_min"]) & 
                (events.dimuon_candidate.mass < params["mass_max"])
            )
            |
            ( 
                (double_electron) &
                (events.dielectron_candidate.mass > params["mass_min"]) &
                (events.dielectron_candidate.mass < params["mass_max"])
            )
        )
    elif params["DY"] is True:
        mask = (
            (events.dilepton_candidate.mass > params["mass_max"] | events.dilepton_candidate.mass < params["mass_min"])
        )
    return ak.where(ak.is_none(mask), False, mask)

dilepton_massZ = Cut(
    name="dilepton_massZ",
    params = {
        "VV" : True,
        "mass_min" : 70,
        "mass_max" : 110,
    },
    function= dilepton_mass,
)
dilepton_massDY = Cut(
    name="dilepton_massDY",
    params = {
        "DY" : True,
        "mass_min" : 70,
        "mass_max" : 110,
    },
    function = dilepton_mass,
)   
#############################################################
#############################################################





#############################################################
# leptons flavour (Z->ll vs tt)                             #
#############################################################
def check_flavour(events, params, year, sample, **kwargs):
    double_leptons = events.nLeptonGood == 2
    single_electron = events.nElectronGood == 1
    single_muon = events.nMuonGood == 1
    if params["VV"] is True:
        mask=(
            double_leptons
            &
            (
                (single_electron & ~single_muon)
                |
                (single_muon & ~single_electron)
            )
        )
    else: 
        mask=(
            double_leptons
            &
            single_electron
            & 
            single_muon
        )
    return ak.where(ak.is_none(mask), False, mask)
check_flavour_SF = Cut(
    name="check_flavour_SF",
    params={ "VV" : True, },
    function=check_flavour,
)
check_flavour_OF =Cut(
    name="check_flavour_OF",
    params={ "VV":False,},
    function=check_flavour,
)
#############################################################
#############################################################















def cut_function(events, params, year, sample, **kwargs):
    masks = []
    for i, c in enumerate(params["coll"]):
        objs = getattr(events, c)
        ptmin = params["minpt"][i]
        nmin = params["nMin"][i]

        if ptmin is not None:
            objs = objs[objs.pt > ptmin]

        mask = ak.num(objs) >= nmin
        masks.append(mask)
    combined_mask = masks[0]
    for m in masks[1:]:
        combined_mask = combined_mask | m
    return ak.where(ak.is_none(combined_mask), False, combined_mask)

def cut_function_2(events, params, year, sample, **kwargs):
    masks = []
    for i, c in enumerate(params["coll"]):
        objs = getattr(events, c)
        nmin = params["nMin"][i]
        mask = ak.num(objs) == nmin
        masks.append(mask)
    combined_mask = masks[0]
    for m in masks[1:]:
        combined_mask = combined_mask | m

    return ak.where(ak.is_none(combined_mask), False, combined_mask)


def skim_single_leptons(events, params, year, sample, **kwargs):
    mask = (
        (ak.num(events.Muon) >=1)
        |
        (ak.num(events.Electron) >=1)
    )
    return ak.where(ak.is_none(mask), False, mask)
skim_sigle_lepton= Cut(
    name="skim_single_leptons",
    params={},
    function=skim_single_leptons
)

def skim_double_leptons(events, params, year, sample, **kwargs):
    mask = (
        (ak.num(events.Muon) >= 2)
        |
        (ak.num(events.Electron) >= 2)
    )
    return ak.where(ak.is_none(mask), False, mask)

skim_double_lepton = Cut(
    name="skim_double_leptons",
    params={},
    function=skim_double_leptons,
)




def get_nObj_min_or(nMin, minpt,coll, name=None):
    if not (len(coll) == len(nMin) and (minpt is None or len(minpt) == len(coll))):
        raise ValueError("coll, nMin, and minpt (if provided) must have same length")

    if name is None:
        name = "cut_" + "_or_".join([
            f"{c}_min{nMin[i]}" + (f"_pt{minpt[i]}" if minpt else "")
            for i, c in enumerate(coll)
        ])
    return Cut(
        name=name,
        params={
            "coll": coll,
            "nMin": nMin,
            "minpt": minpt if minpt else [None] * len(coll),
        },
        function=cut_function,
    )

def single_good_electron(events, params, year, sample, **kwargs):
    mask = (
        ak.num(events.ElectronGood) >= 2
        )
    return ak.where(ak.is_none(mask), False, mask)  

def single_good_muon(events, params, year, sample, **kwargs):
    mask = (
        ak.num(events.MuonGood) >= 2
    )  
    return ak.where(ak.is_none(mask), False, mask)





SingleEle = Cut(
    name="SingleGoodEle",
    params = {},
    function=single_good_electron,
)
SingleMuon = Cut(
    name="SingleGoodMuon", 
    params = {},
    function=single_good_muon,
)
SingleLepton = Cut(
    name="SingleGoodLeptons",
    params={},
    function=single_good_lepton,
)
def get_nObj_eq_or(nMin, coll, name=None):
    if name is None:
        name = "cut_" + "_or_".join([
            f"{c}_min{nMin[i]}"
            for i, c in enumerate(coll)
        ])
    return Cut(
        name=name,
        params={
            "coll": coll,
            "nMin": nMin,
        },
        function=cut_function_2,
    )