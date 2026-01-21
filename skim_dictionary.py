skim_dict = {
    "Wlep_V" : {
         "Muon" : {
             "ptMin" : 15,
             "nMin" : 1,
         },
         "Electron" : {
             "ptMin" : 1,
             "nMin" : 1,
         },
         "Jet" : {
             "ptMin" : 1,
             "nMin" : 2,
         },
         "triggers_forfakes" : [
             "PFJet40",
             "PFJet60",
             "PFJet80",
             "PFJet140",
             "PFJet200",
             "PFJet200",
             "PFJet260",
             "PFHT125",
             "PFHT200",
             "PFHT250",
             "PFHT300",
             "PFHT350",
            ],
         "triggers_muon" : [
             "IsoMu24",
             "IsoMu24_eta2p1",
             "IsoMu24_TwoProngs35",
             "IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1",
             "IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
             "IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1"
             
         ]
         
    },
    "ZZ" : {
         "Muon" : {
             "ptMin" : 1,
             "nMin" : 2,
         },
         "Electron" : {
             "ptMin" : 1,
             "nMin" : 2,
         },
         "Jet" : {
             "ptMin" : 1,
             "nMin" : 2,
         }
    },
    "DY" : {
        "Muon" : {
            "ptMin" : 1,
            "nMIn" : 2,
        },
        "Electron" : {
            "ptMin" : 1,
            "nMin" : 2, 
        },
        "Jet" : {
            "ptMin" : 1, 
            "nMin": 2,
        },
    },
    "ttsemilep" : {
        "Muon" : {
            "ptMin" : 1,
            "nMin" : 1,
        },
        "Electron" : {
            "ptMin" : 1,
            "nMin" : 1, 
        },
        "Jet" : {
            "ptMin" : 1, 
            "nMin": 1,
        },
    },
    "Wjets" : {
        "Muon" : {
            "ptMin" : 1,
            "nMin" : 1,
        },
        "Electron" : {
            "ptMin" : 1,
            "nMin" : 1, 
        },
        "Jet" : {
            "ptMin" : 1, 
            "nMin": 1,
        },
    }
    
}