from custom_cut_functions import *
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq, count_objects_eq
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
import workflowVBS
from workflowVBS import VBS_WV_Processor
from skim_dictionary import skim_dict
from pocket_coffea.lib.weights.common import common_weights
import cloudpickle
from pocket_coffea.lib.columns_manager import ColOut
import os
import custom_cut_functions 
cloudpickle.register_pickle_by_value(workflowVBS)
cloudpickle.register_pickle_by_value(custom_cut_functions)
localdir = os.path.dirname(os.path.abspath(__file__))
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/parameters")


parameters = defaults.merge_parameters_from_files(default_parameters,
                                    f"{localdir}/parameters/object_presel.yaml",
                                    f"{localdir}/parameters/btagging.yaml",
                                    update=True                                    
                                    )
cfg = Configurator(
    parameters=parameters,
    datasets = {
        "tag" : "VBS_ssWW",
        "jsons" : [#f"{localdir}/datasets/VBS_mc_132X_Summer23wmLHEGS.json",
                   f"{localdir}/datasets/TTbar/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8.json",
                   f"{localdir}/datasets/TTbar/TTtoLNu2Q_HT-500_NJet-9_Hdamp-418_TuneCP5_13p6TeV_powheg-pythia8.json",
                   f"{localdir}/datasets/TTbar/TTtoLNu2Q_HT-500_NJet-9_TuneCP5_13p6TeV_powheg-pythia8.json",
                   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.json",
                   f"{localdir}/datasets/TbarWplus/TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8.json",
                   #f"{localdir}/datasets/VVV/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8.json",
                   #f"{localdir}/datasets/VVV/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8.json",
                   ],
        "filter" : {
            "samples" : [
               "ssWWLL",
               "ssWWTT",
              # "ssWWTL",
              # "ssWW_unpolarized",
                "TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8", 
                "TTtoLNu2Q_HT-500_NJet-9_TuneCP5_13p6TeV_powheg-pythia8",
                "WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", 
                "TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8", 
                "WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8", 
                "ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8",
                
            ],
            "year" : ["2023_preBPix"]
        }
        }, 
    workflow=VBS_WV_Processor,
    skim = [
            get_nPVgood(1), 
            get_nObj_min_or(
                [skim_dict["Wlep_V"]["Muon"]["nMin"],skim_dict["Wlep_V"]["Electron"]["nMin"]], 
                [skim_dict["Wlep_V"]["Muon"]["ptMin"], skim_dict["Wlep_V"]["Electron"]["ptMin"]],
                ["Muon", "Electron"]),            
            ],
    preselections=[SingleLepton, VBS_jets_presel, semileptonic_preselW],
    categories= {
        "baseline" : [passthrough],
        
        "SingleEle_AK8" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_boosted],
        "SingleEle_AK4" : [get_nElectron(1, coll="ElectronGood"),  get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_resolved],
        "SingleMuon_AK8" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"),  get_nObj_eq(0, coll="BJetGood"), Vjet_massW_boosted],
        "SingleMuon_AK4" : [get_nMuon(1, coll="MuonGood"),  get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_resolved],
        "SingleLepton_AK8" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet") , get_nObj_eq(0, coll="BJetGood"),  Vjet_massW_boosted],
        "SingleLepton_AK4" : [get_nObj_eq(1, coll="LeptonGood"),get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_resolved],

        # ttbar/tbarWplus
        "SingleEle_AK8_bjets_ttbar" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1 , coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        "SingleEle_AK4_bjets_ttbar" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(0 , coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        "SingleMuon_AK8_bjets_ttbar" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1 , coll="CleanFatJet"), get_nObj_eq(1, coll="BJetGood")],
        "SingleMuon_AK4_bjets_ttbar" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        "SingleLepton_AK8_bjets_ttbar" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        "SingleLepton_AK4_bjets_ttbar" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        
        #WtoLNu-XJets (check contamination in the jet mass window of SR)
        "SingleEle_AK8_sideL_Wjets" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"),  Wjet_sideL_boosted],
        "SingleEle_AK8_sideR_Wjets" : [get_nElectron(1, coll="ElectronGood"),  get_nObj_eq(1, coll="CleanFatJet"),   Wjet_sideR_boosted],
        "SingleMuon_AK8_sideL_Wjets" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideL_boosted],
        "SingleMuon_AK8_sideR_Wjets" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideR_boosted],
        "SingleLepton_AK8_sideL_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideL_boosted],
        "SingleLepton_AK8_sideR_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideR_boosted],

        "SingleLepton_AK4_sideL_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(0, coll="CleanFatJet"), Wjet_sideL_resolved],
        "SingleLepton_AK4_sideR_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(0, coll="CleanFatJet"), Wjet_sideR_resolved],
    },    
    weights_classes = common_weights,
    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
                          ],
            "bycategory" : {
            }
        },
        "bysample": {
        }
    },
    variations = {
        "weights": {
            "common": {
                "inclusive": [ 
                              ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
    },
    variables = {},
    workflow_options = {
        "dump_columns_as_arrays_per_chunk": "root://eosuser.cern.ch//eos/user/l/ldellape/VBS/parquet/"
    },
    columns = {
        "common" : {
            "inclusive" : [ColOut("ElectronGood", ["pt","eta", "phi"], flatten=False),
                           ColOut("MuonGood", ["pt", "eta", "phi"], flatten=False),
                           ColOut("LeptonGood", ["pt", "eta", "phi"], flatten=False),
                           ColOut("CleanFatJet", ["pt", "eta", "phi", "tau1", "tau2", "tau21", "msoftdrop", "mass"], flatten=False),
                           ColOut("CleanJet", ["pt", "eta", "phi"], flatten=False),
                           ColOut("CleanJet_noVBS", ["pt", "eta", "phi"], flatten=False),
                           ColOut("VBS_dijet_system", ["mass", "pt", "deltaEta", "pt1", "pt2"], flatten=False),    
                           ColOut("V_dijet_candidate", ["mass", "pt", "deltaEta"], flatten=False),  
                           ColOut("events", ["zepp_ele", "zepp_muon", "zepp_lep"], flatten=False),
                           ColOut("events", ["MT_mu_miss", "MT_ele_miss", "MT_lep_miss"], flatten=False),
                           ColOut("events", ["nCleanJets", "nCleanFatJets", "nBJetGood", "nJet", "nFatJet"], flatten=False),
                           ColOut("CleanSubJet_pair", ["zg"], flatten=False),
                           ],
            "bycategory" : {},
            },            
    },   
)