from custom_cut_functions import *
from custom_weights import *
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel, get_HLTsel_custom, get_nPVgood, goldenJson, eventFlags, get_nElectron, get_nMuon, get_nObj_eq, count_objects_eq, apply_golden_json, get_JetVetoMap
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
import custom_weights
cloudpickle.register_pickle_by_value(workflowVBS)
cloudpickle.register_pickle_by_value(custom_cut_functions)
cloudpickle.register_pickle_by_value(custom_weights)
localdir = os.path.dirname(os.path.abspath(__file__))
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/parameters")

from custom_weights import SF_pt_top, fakesMuon

parameters = defaults.merge_parameters_from_files(default_parameters,
                                    f"{localdir}/parameters/object_presel.yaml",
                                    f"{localdir}/parameters/btagging.yaml",
                                  #  f"{localdir}/parameters/triggers_forfake.yaml",
                                    update=True                                    
                                    )
cfg = Configurator(
    parameters=parameters,
    datasets = {
        "tag" : "VBS_ssWW",
        "jsons" : [
                   
                   #DY
                 #  f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-4to10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-50_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-50_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-2Jets_MLL-50_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-4Jets_MLL-50to120_HT-2500_TuneCP5_13p6TeV_madgraphMLM-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-4Jets_MLL-50to120_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-4Jets_MLL-50to120_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8_redirector.json",
                 #  f"{localdir}/datasets/DY/DYto2L-4Jets_MLL-50to120_HT-800to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8_redirector.json",
                   
                   
                   # W->lv, 2 jets
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",     
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                #   f"{localdir}/datasets/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
         #          # ttbar 
                 #  f"{localdir}/datasets/TTbar/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",
                 #  f"{localdir}/datasets/TTbar/TTtoLNu2Q_HT-500_NJet-9_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",
               #    f"{localdir}/datasets/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",

                   
                  # f"{localdir}/datasets/tW/TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8_redirector.json",
                 #  f"{localdir}/datasets/tW/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",
                 #  f"{localdir}/datasets/tW/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",
          #         f"{localdir}/datasets/TTto4Q_Hdamp-418_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",
           #        f"{localdir}/datasets/TTto4Q_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_redirector.json",

            #       f"{localdir}/datasets/JetMET0.json",
             #      f"{localdir}/datasets/JetMET1.json",
                #   f"{localdir}/datasets/DATA/MUON/2023/Muon1_Run2023D_v2_redirector.json",
                #   f"{localdir}/datasets/DATA/MUON/2023/Muon1_Run2023C_2_v12_redirector.json",
              #      f"{localdir}/datasets/DATA/Muon1_Run2023D_v1_redirector.json",

                    f"{localdir}/datasets/DATA/MUON/2023/Muon1_2023C_v1.json",
                    f"{localdir}/datasets/DATA/MUON/2023/Muon1_2023C_v2.json",
                    f"{localdir}/datasets/DATA/MUON/2023/Muon1_2023C_v3.json",
                    f"{localdir}/datasets/DATA/MUON/2023/Muon1_2023C_v4.json",

                

               #     f"{localdir}/datasets/TTZ/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_redirector.json",
                    
                #   f"{localdir}/datasets/VVV/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_redirector.json",
                #   f"{localdir}/datasets/VVV/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_redirector.json",
                #   f"{localdir}/datasets/VVV/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8_redirector.json",
                #   f"{localdir}/datasets/VVV/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_redirector.json",
                    
               #    f"{localdir}/datasets/DATA/EGamma1.json",
                   ],
        "filter" : {
            "samples" : [
         #       "WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #       "WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #        "WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #        "WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #        "WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #        "WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #        "WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        #        "WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", 
       #        "WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",   
        
       #       "TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
       #       "TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8",
       #       "TTtoLNu2Q_HT-500_NJet-9_Hdamp-418_TuneCP5_13p6TeV_powheg-pythia8",
      #         "TTtoLNu2Q_HT-500_NJet-9_TuneCP5_13p6TeV_powheg-pythia8",
        
     #   "TTto4Q_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8",
         #   "Muon0",
           "Muon1",
          #  "Muon1_forFake",
            #  "EGamma1",
            # "ssWWLL",
            #  "ssWWTT",
            #  "ssWWLL",
            #  "ssWWTL",
            #  "JetMET0",
            #  "JetMET1",
            #  "ssWWTT",
            #  "ssWWTL",
            #  "ssWW_unpolarized",
            #  "TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8", 
              
            #  "WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", 
           #   "TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8", 
            #
       #    "TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        #    "TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
              
      #          "WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8", 
      #          "ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8",
      #          "WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8",
       #         "WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8",
             
            ],
            "year" : ["2023_preBPix", "2023_postBPix"],
        }
        }, 
    workflow=VBS_WV_Processor,
    skim = [
            get_nPVgood(1), 
            eventFlags,
            goldenJson,
            get_HLTsel_custom(["IsoMu27"]),
            get_JetVetoMap("JetVetoMaps"),
            get_nObj_min_or(
                [skim_dict["Wlep_V"]["Muon"]["nMin"],skim_dict["Wlep_V"]["Electron"]["nMin"]], 
                [skim_dict["Wlep_V"]["Muon"]["ptMin"], skim_dict["Wlep_V"]["Electron"]["ptMin"]],
                ["Muon", "Electron"]),  
            ],
    #preselections=[SingleLepton, VBS_jets_presel, semileptonic_preselW],
    # for fakes 
    preselections=[SingleLepton, VBS_jets_presel, Wtransverse_mass_presel],
    categories= {
        "baseline" : [passthrough],
        #"SingleEle_AK8" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="JetGood"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_boosted],
        #"SingleEle_AK4" : [get_nElectron(1, coll="ElectronGood"),  get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="JetGood"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_resolved],
        #"SingleMuon_AK8" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"),  get_nObj_eq(0, coll="BJetGood"), Vjet_massW_boosted],
        #"SingleMuon_AK4" : [get_nMuon(1, coll="MuonGood"),  get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_resolved],
      #  "SingleLepton_AK8" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet") , get_nObj_eq(0, coll="BJetGood"),  Vjet_massW_boosted],
      #  "SingleLepton_AK4" : [get_nObj_eq(1, coll="LeptonGood"),get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Vjet_massW_resolved],
      #  "Fakes_loose_ele" : [FakeLoose,  get_nObj_eq(1, coll="ElectronLoose")],
      #  "Fakes_loose_muon" : [get_nObj_eq(1, coll="MuonLoose")],
        
     #   "Fakes_tight_ele" : [FakeTight, get_nObj_eq(1, coll="ElectronGood")],
     #   "Fakes_tight_muon" : [get_nObj_eq(1, coll="MuonGood")],


        # ttbar/tbarWplus
        
      #  "Fakes_tight_ele" : [FakeTight, get_nObj_eq(1, coll="MuonGood"), get_nObj_min(1, coll="JetForFakes_tight")],
       # "Fakes_tight_muon" : [FakeTight, get_nObj_eq(1, coll="ElectronGood"), get_nObj_min(1, coll="JetForFakes_tight")],


        # ttbar/tbarWplus
        #"SingleEle_AK8_bjets_ttbar" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1 , coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        #"SingleEle_AK4_bjets_ttbar" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(0 , coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        #"SingleMuon_AK8_bjets_ttbar" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1 , coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        #"SingleMuon_AK4_bjets_ttbar" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(1, coll="BJetGood")],
        #"SingleMuon_bjets_ttbar_inclusive" : [get_nMuon(1, coll="MuonGood"), get_nObj_min(1, coll="BJetGood")], 
 #       "SingleLepton_AK8_bjets_mediumWP_ttbar" : [get_nObj_eq(1, coll="MuonLoose"), get_nObj_eq(0, coll="MuonGood"), semileptonic_preselW, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGood")],
 #       "SingleLepton_AK8_bjets_mediumWP_ttbar" : [get_nObj_eq(1, coll="MuonLoose"), get_nObj_eq(0, coll="MuonGood"), semileptonic_preselW, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGood")],


        "SingleLepton_AK8_bjets_mediumWP_ttbar" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGood")],
     #   "SingleLepton_AK4_bjets_mediumWP_ttbar" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW, get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), Vjet_massW_resolved, get_nObj_min(1, coll="BJetGood")],
        "SingleLepton_AK8_looseWP_bjets_ttbar" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGoodLoose")],
     #   "SingleLepton_AK4_looseWP_bjets_ttbar" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW, get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), Vjet_massW_resolved, get_nObj_min(1, coll="BJetGoodLoose")],
        "SingleLepton_AK8_bjets_mediumWP_ttbar_preselW_2" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_2, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGood")],
     #   "SingleLepton_AK4_bjets_mediumWP_ttbar_preselW_2" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_2, get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), Vjet_massW_resolved, get_nObj_min(1, coll="BJetGood")],
        "SingleLepton_AK8_looseWP_bjets_ttbar_preselW_2" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_2, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGoodLoose")],
     #   "SingleLepton_AK4_looseWP_bjets_ttbar_preselW_2" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_2, get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), Vjet_massW_resolved, get_nObj_min(1, coll="BJetGoodLoose")],
        "SingleLepton_AK8_bjets_mediumWP_ttbar_preselW_3" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_3, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGood")],
     #   "SingleLepton_AK4_bjets_mediumWP_ttbar_preselW_3" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_3, get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), Vjet_massW_resolved, get_nObj_min(1, coll="BJetGood")],
        "SingleLepton_AK8_looseWP_bjets_ttbar_preselW_3" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_3, get_nObj_eq(1, coll="CleanFatJet"), get_nObj_min(2, coll="CleanJet"), Vjet_massW_boosted, get_nObj_min(1, coll="BJetGoodLoose")],
      #  "SingleLepton_AK4_looseWP_bjets_ttbar_preselW_3" : [get_nObj_eq(1, coll="MuonGood"), semileptonic_preselW_3, get_nObj_eq(0, coll="CleanFatJet"), get_nObj_min(4, coll="CleanJet"), Vjet_massW_resolved, get_nObj_min(1, coll="BJetGoodLoose")],
       # "SingleLepton_AK8_sideBand" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Wjet_side_boosted],
       # "SingleLepton_AK4_sideBand" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(0, coll="CleanFatJet"), get_nObj_eq(0, coll="BJetGood"), Wjet_side_resolved],
       # "SingleLepton_bjets_ttbar_inclusive" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(1, coll="BJetGood")],
        
        #WtoLNu-XJets (check contamination in the jet mass window of SR)
        #"SingleEle_AK8_sideL_Wjets" : [get_nElectron(1, coll="ElectronGood"), get_nObj_eq(1, coll="CleanFatJet"),  Wjet_sideL_boosted],
        #"SingleEle_AK8_sideR_Wjets" : [get_nElectron(1, coll="ElectronGood"),  get_nObj_eq(1, coll="CleanFatJet"),   Wjet_sideR_boosted],
        #"SingleMuon_AK8_sideL_Wjets" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideL_boosted],
        #"SingleMuon_AK8_sideR_Wjets" : [get_nMuon(1, coll="MuonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideR_boosted],
        #"SingleLepton_AK8_sideL_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideL_boosted],
        #"SingleLepton_AK8_sideR_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_eq(1, coll="CleanFatJet"), Wjet_sideR_boosted],

        #"SingleLepton_AK4_sideL_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(0, coll="CleanFatJet"), Wjet_sideL_resolved],
        #"SingleLepton_AK4_sideR_Wjets" : [get_nObj_eq(1, coll="LeptonGood"), get_nObj_min(0, coll="CleanFatJet"), Wjet_sideR_resolved],
    },    
    weights_classes = common_weights + [SF_pt_top, fakesMuon, fakesEle],
    weights = {
        "common": { "inclusive": ["genWeight","lumi","XS", "pileup","sf_mu_id", "sf_mu_iso", "sf_top_pt", "sf_btag", "fakesMU", "fakesEle"],
            "bycategory" : {
            }
        },
    },
    variations={"weights": {"common": {"inclusive": ["pileup","sf_mu_id","sf_mu_iso" , "sf_btag", "fakesMU", "fakesEle"]}}},

    workflow_options = {
        "dump_columns_as_arrays_per_chunk": "root://eosuser.cern.ch//eos/user/l/ldellape/VBS/parquet_forfake/"
    },
    columns = {},

    variables = {
        "ElectronGood_pt"  : HistConf([Axis(coll="ElectronGood", field="pt",  bins=20, start=0, stop=500, label="Electron pT")]),
        "ElectronGood_eta" : HistConf([Axis(coll="ElectronGood", field="eta", bins=20, start=-5, stop=5, label="Electron η")]),
        "ElectronGood_phi" : HistConf([Axis(coll="ElectronGood", field="phi", bins=20, start=-3.2, stop=3.2, label="Electron φ")]),
        "MuonGood_pt"  : HistConf([Axis(coll="MuonGood", field="pt",  bins=30, start=0, stop=220, label="Muon pT")]),
        "MuonGood_eta" : HistConf([Axis(coll="MuonGood", field="eta", bins=30, start=-2.5, stop=2.5, label="Muon η")]),
        "MuonGood_phi" : HistConf([Axis(coll="MuonGood", field="phi", bins=30, start=-3.2, stop=3.2, label="Muon φ")]),
        "ElectronLoose_pt"  : HistConf([Axis(coll="ElectronLoose", field="pt",  bins=30, start=0, stop=500, label="Muon pT")]),
        "ElectronLoose_eta" : HistConf([Axis(coll="ElectronLoose", field="eta", bins=30, start=-5, stop=5, label="Muon η")]),
        "ElectronLoose_phi" : HistConf([Axis(coll="ElectronLoose", field="phi", bins=30, start=-3.2, stop=3.2, label="Muon φ")]),
        "MuonLoose_pt"  : HistConf([Axis(coll="MuonLoose", field="pt",  bins=30, start=0, stop=500, label="Muon pT")]),
        "MuonLoose_eta" : HistConf([Axis(coll="MuonLoose", field="eta", bins=30, start=-5, stop=5, label="Muon η")]),
        "MuonLoose_phi" : HistConf([Axis(coll="MuonLoose", field="phi", bins=10, start=-3.2, stop=3.2, label="Muon φ")]),
        "LeptonGood_pt"  : HistConf([Axis(coll="LeptonGood", field="pt",  bins=30, start=0, stop=500, label="Lepton pT")]),
        "LeptonGood_eta" : HistConf([Axis(coll="LeptonGood", field="eta", bins=30, start=-5, stop=5, label="Lepton η")]),
        "LeptonGood_phi" : HistConf([Axis(coll="LeptonGood", field="phi", bins=30, start=-3.2, stop=3.2, label="Lepton φ")]),
        "CleanFatJet_pt"        : HistConf([Axis(coll="CleanFatJet", field="pt", bins=30, start=170, stop=800, label="FatJet pT")]),
        "CleanFatJet_eta"       : HistConf([Axis(coll="CleanFatJet", field="eta", bins=30, start=-2.5, stop=2.5, label="FatJet η")]),
        "CleanFatJet_phi"       : HistConf([Axis(coll="CleanFatJet", field="phi", bins=30, start=-3.2, stop=3.2, label="FatJet φ")]),
        "CleanFatJet_tau1"      : HistConf([Axis(coll="CleanFatJet", field="tau1", bins=30, start=0, stop=1, label="τ₁")]),
        "CleanFatJet_tau2"      : HistConf([Axis(coll="CleanFatJet", field="tau2", bins=30, start=0, stop=1, label="τ₂")]),
        "CleanFatJet_tau21"     : HistConf([Axis(coll="CleanFatJet", field="tau21", bins=30, start=0, stop=1, label="τ₂₁")]),
        "CleanFatJet_msoftdrop" : HistConf([Axis(coll="CleanFatJet", field="msoftdrop", bins=30, start=50, stop=120, label="SoftDrop Mass")]),
        "CleanFatJet_mass"      : HistConf([Axis(coll="CleanFatJet", field="mass", bins=30, start=0, stop=300, label="FatJet Mass")]),
        "MET_pt"      : HistConf([Axis(coll="MET", field="pt", bins=30, start=0, stop=300, label="MET  pt")]),
        "JetGood_pt"  : HistConf([Axis(coll="JetGood", field="pt", bins=30, start=0, stop=500, label="JetGood pT")]),
        "JetGood_eta" : HistConf([Axis(coll="JetGood", field="eta", bins=30, start=-5, stop=5, label="JetGood η")]),
        "JetGood_phi" : HistConf([Axis(coll="JetGood", field="phi", bins=30, start=-3.2, stop=3.2, label="JetGood φ")]),     
        "VBSjet_pt" : HistConf([Axis(coll="VBS_dijet_system", field="pt", bins=30, start=0, stop=800, label="VBS jet pT")]),
        "VBSjet_deltaEta" : HistConf([Axis(coll="VBS_dijet_system", field="pt", bins=30, start=0, stop=800, label="VBS jet delta|η|")]),
        "VBSjet_mass" : HistConf([Axis(coll="VBS_dijet_system", field="mass", bins=30, start=400, stop=1000, label="VBS jet m_{jj}")]),
    },
)
