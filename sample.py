def samples(category, sample, polarized, path="/eos/user/l/ldellape/VBS/parquet"): 
    parquet_patterns = []
    labels = []
    if category == "baseline":
        if sample=="WW":
            parquet_patterns = [
                f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/{category}/*.parquet",
                f"{path}/ssWW_TT_mg5_madspin/{category}/*.parquet",
                f"{path}/ssWW_LL_mg5_madspin/{category}/*.parquet",
            #    f"{path}/TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8_2023_preBPix/{category}/*.parquet"
                f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/{category}/*.parquet",
                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/{category}/*.parquet",
                f"{path}/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix",
            ]
            if polarized is True:
                labels = ["TTtoLNu2Q", "ssWW_TT", "ssWW_LL", "WLtoLNu_2Jets", "ZZZ", "WZZ"]
            else: 
                labels = ["TTtoLNu2Q", "ssWW", "ssWW", "WLtoLNu_2Jets", "WZZ+ZZZ", "WZZ+ZZZ"]
        elif sample=="ZZ":
            parquet_patterns = [
                f"{path}/ZZLL_mg5_madspin/{category}/*.parquet",
                f"{path}/ZZTT_mg5_madspin/{category}/*.parquet",
                f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/{category}/*.parquet",
                f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/{category}/*.parquet",
            ]
            if polarized is True:
                labels = ["ZZ_LL", "ZZ_TT", "TTto2L2Nu-2Jets", "DYto2L-4Jets"]
            else: 
                labels = ["ZZ", "ZZ", "TTto2L2Nu-2Jets", "WLtoLNu_2Jets" , "DYto2L-4Jets"]
                
    elif category=="signal_AK8":
        if sample == "WW":
            parquet_patterns = [
                f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK8/*.parquet",
                f"{path}/ssWW_LL_mg5_madspin/SingleLepton_AK8/*.parquet",
                f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",
                f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",
                                                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",
                        f"{path}/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",
                                f"{path}/TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK8/*.parquet",


            


            ]
            if polarized is True:
                labels = ["ssWW_TT", "ssWW_LL", "TTtoLNu2Q", "WLtoLNu_2Jets", "ZZZ", "WZZ", "tWtoLNu-2Q"]
            else:
                labels = ["ssWW", "ssWW", "TTtoLNu2Q", "WLtoLNu_2Jets"]
        elif sample == "ZZ":
            parquet_patterns = [
                f"{path}/ZZLL_mg5_madspin/DoubleMuon_AK8/*.parquet",
                f"{path}/ZZLL_mg5_madspin/DoubleEle_AK8/*.parquet",
                f"{path}/ZZTT_mg5_madspin/DoubleMuon_AK8/*.parquet",
                f"{path}/ZZTT_mg5_madspin/DoubleEle_AK8/*.parquet",
                f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/DoubleEle_AK8/*.parquet",
                f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/DoubleMuon_AK8/*.parquet",
                f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/DoubleEle_AK8/*.parquet",
                f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/DoubleMuon_AK8/*.parquet",
            ]
            if polarized is True:
                labels = ["ZZ_LL", "ZZ_LL", "ZZ_TT", "ZZ_TT",  "TTto2L2Nu-2Jets", "TTto2L2Nu-2Jets", "WLtoLNu_2Jets", "DYto2L-4Jets", "DYto2L-4Jets"]
            else: 
                labels = ["ZZ", "ZZ", "ZZ", "ZZ", "TTto2L2Nu-2Jets", "TTto2L2Nu-2Jets" , "WLtoLNu_2Jets", "DYto2L-4Jets", "DYto2L-4Jets"]
    elif category=="side":
        if sample== "WW":
            parquet_patterns = [
                f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK8_sideR_Wjets/*.parquet",
                f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK8_sideL_Wjets/*.parquet",
                f"{path}/ssWW_LL_mg5_madspin/SingleLepton_AK8_sideL_Wjets/*.parquet",
                f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK8_sideR_Wjets/*.parquet",
                f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8_sideR_Wjets/*.parquet",
                f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8_sideL_Wjets/*.parquet",
             #   f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8_sideR_Wjets/*.parquet",
              #  f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8_sideL_Wjets/*.parquet",
                                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_sideL_Wjets/*.parquet",
             #   f"{path}/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_sideR_Wjets/*.parquet",
                                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_sideL_Wjets/*.parquet",
              #  f"{path}/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/SingleLepton_AK8_sideR_Wjets/*.parquet",
            ]
            if polarized is True:
                labels = ["ssWW_TT", "ssWW_TT", "ssWW_TT", "ssWW_TT", "ssWW_LL", "ssWW_LL", "ssWW_LL", "ssWW_LL","TTtoLNu2Q", "TTtoLNu2Q", "TTtoLNu2Q" , "TTtoLNu2Q", "WLtoLNu_2Jets","WLtoLNu_2Jets", "WLtoLNu_2Jets", "WtoLNu_2Jets" , 
                          "tWtoLNu2Q", "tWtoLNu2Q", "tWtoLNu2Q", "tWtoLNu2Q", "ZZZ", "ZZZ", "ZZZ", "ZZZ", "WZZ", "WZZ", "WZZ", "WZZ"]
            else: 
                labels = ["ssWW", "ssWW" ,"ssWW","ssWW", "WtoLNu-2Jets",   "WtoLNu-2Jets", "WtoLNu-2Jets", "WLtoLNu_2Jets","TTbar","TTbar" , "WZZ", "WZZ", "ZZZ", "ZZZ"] 
            
    elif category=="signal_AK4":
        if sample == "WW":
            parquet_patterns = [
                      f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK4/*.parquet",
                f"{path}/ssWW_LL_mg5_madspin/SingleLepton_AK4/*.parquet",
                f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
                f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
                                                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
                        f"{path}/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
                                f"{path}/TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK4/*.parquet",
            ]
            if polarized is True:
                labels = ["ssWW_TT", "ssWW_LL", "TTtoLNu2Q", "WLtoLNu_2Jets", "ZZZ", "WZZ", "tWtoLNu-2Q"]
            else:
                labels = ["ssWW", "ssWW", "TTtoLNu2Q", "WLtoLNu_2Jets"]     
        elif sample == "ZZ":
            parquet_patterns = [
                    f"{path}/ZZLL_mg5_madspin/DoubleMuon_AK4/*.parquet",
                    f"{path}/ZZLL_mg5_madspin/DoubleEle_AK4/*.parquet",
                    f"{path}/ZZTT_mg5_madspin/DoubleMuon_AK4/*.parquet",
                    f"{path}/ZZTT_mg5_madspin/DoubleEle_AK4/*.parquet",
                    f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/DoubleEle_AK4/*.parquet",
                    f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/DoubleMuon_AK4/*.parquet",
                    f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/DoubleEle_AK4/*.parquet",
                    f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/DoubleMuon_AK4/*.parquet",
            ]
            if polarized is True:
                labels = ["ZZ_LL", "ZZ_LL", "ZZ_TT", "ZZ_TT", "TTto2L2Nu-2Jets", "TTto2L2Nu-2Jets" , "DYto2L-4Jets", "DYto2L-4Jets"]
            else: 
                labels = ["ZZ", "ZZ", "ZZ", "ZZ",  "TTto2L2Nu-2Jets", "TTto2L2Nu-2Jets" ,"WLtoLNu_2Jets", "DYto2L-4Jets", "DYto2L-4Jets" ]
    elif category=="side":
            parquet_pattern
    elif category=="CR_TTbar":
        if sample == "WW":
            parquet_patterns = [
                f"{path}/ssWW_TT_mg5_madspin/SingleLepton_AK8_bjets_ttbar/*.parquet",
                f"{path}/ssWW_LL_mg5_madspin/SingleLepton_AK4_bjets_ttbar/*.parquet",
                f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/*.parquet",
           #     f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleMuon_AK4_bjets_ttbar/*.parquet",
           #     f"{path}/TTtoLNu2Q_HT-500_NJet-9_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8_2023_preBPix/SingleMuon_AK4_bjets_ttbar/*.parquet",
                f"{path}/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/*.parquet",
                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/*.parquet",
                f"{path}/TbarWplustoLNu2Q_TuneCP5Down_13p6TeV_powheg-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/*.parquet",
                f"{path}/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/*.parquet",
                f"{path}/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_2023_preBPix/SingleLepton_AK8_bjets_ttbar/*.parquet",

            ]
            if polarized is True:
                labels = ["ssWW_TT", "ssWW_TT", "ssWW_TT", "ssWW_TT", "ssWW_LL", "ssWW_LL", "ssWW_LL", "ssWW_LL","TTtoLNu2Q", "TTtoLNu2Q", "TTtoLNu2Q" , "TTtoLNu2Q", "WLtoLNu_2Jets","WLtoLNu_2Jets", "WLtoLNu_2Jets", "WtoLNu_2Jets" , 
                          "tWtoLNu2Q", "tWtoLNu2Q", "tWtoLNu2Q", "tWtoLNu2Q", "ZZZ", "ZZZ", "ZZZ", "ZZZ", "WZZ", "WZZ", "WZZ", "WZZ"]
            else: 
                labels = ["ssWW", "ssWW" ,"TTtoLNu2Q",  "WLtoLNu_2Jets","tWtoLNu2Q", "ZZZ", "WZZ"] 
        elif sample == "ZZ":
            parquet_patterns = [
                f"{path}/ZZLL_mg5_madspin/AK8_OF/*.parquet",
                f"{path}/ZZLL_mg5_madspin/AK4_OF/*.parquet",
                f"{path}/ZZTT_mg5_madspin/AK8_OF/*.parquet",
                f"{path}/ZZTT_mg5_madspin/AK4_OF/*.parquet",         
                f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/AK8_OF/*.parquet",
                f"{path}/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_2023_preBPix/AK4_OF/*.parquet",
                f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/AK4_OF/*.parquet",
                f"{path}/DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_2023_preBPix/AK8_OF/*.parquet",

            ]
            if polarized is True:
                labels = ["ZZ_LL",  "ZZ_LL","ZZ_TT", "ZZ_TT", "TTto2L2Nu-2Jets","TTto2L2Nu-2Jets",  "DYto2L-4Jets", "DYto2L-4Jets"  ]
            else:
                labels = ["ZZ", "ZZ", "ZZ", "ZZ","TTto2L2Nu-2Jets","TTto2L2Nu-2Jets", "DYto2L-4Jets", "DYto2L-4Jets" ]

    return parquet_patterns, labels
    