void zg_plots() {
    // Open ROOT files
    TFile *f_LL = TFile::Open("output_diffMassAK4_AK8 ssWW_LL_mg5+ms.root");
    TFile *f_TT = TFile::Open("output_diffMassAK4_AK8 ssWW_TT_mg5+ms.root");

    // --- Get histograms (AK8) ---
    TH1F *h_zg_LL_ak8 = (TH1F*) f_LL->Get("h_subjet_zg_ak8_ssWW_LL_mg5+ms");
    TH1F *h_zg_TT_ak8 = (TH1F*) f_TT->Get("h_subjet_zg_ak8_ssWW_TT_mg5+ms");
    TH1F *h_ratio_LL_ak8 = (TH1F*) f_LL->Get("h_subjet_ratio_ak8_ssWW_LL_mg5+ms");
    TH1F *h_ratio_TT_ak8 = (TH1F*) f_TT->Get("h_subjet_ratio_ak8_ssWW_TT_mg5+ms");

    // --- Get histograms (AK4) ---
    TH1F *h_zg_LL_ak4 = (TH1F*) f_LL->Get("h_subjet_zg_ak4_ssWW_LL_mg5+ms");
    TH1F *h_zg_TT_ak4 = (TH1F*) f_TT->Get("h_subjet_zg_ak4_ssWW_TT_mg5+ms");
    TH1F *h_ratio_LL_ak4 = (TH1F*) f_LL->Get("h_subjet_ratio_ak4_ssWW_LL_mg5+ms");
    TH1F *h_ratio_TT_ak4 = (TH1F*) f_TT->Get("h_subjet_ratio_ak4_ssWW_TT_mg5+ms");

    // --- Rebin all histograms once ---
    h_zg_LL_ak8->Rebin(3);
    h_zg_TT_ak8->Rebin(3);
    h_ratio_LL_ak8->Rebin(3);
    h_ratio_TT_ak8->Rebin(3);

    h_zg_LL_ak4->Rebin(3);
    h_zg_TT_ak4->Rebin(3);
    h_ratio_LL_ak4->Rebin(3);
    h_ratio_TT_ak4->Rebin(3);

    h_zg_LL_ak8->Scale(1./h_zg_LL_ak8->Integral());
    h_zg_LL_ak4->Scale(1./h_zg_LL_ak4->Integral());
    h_zg_TT_ak8->Scale(1./h_zg_TT_ak8->Integral());
    h_zg_TT_ak4->Scale(1./h_zg_TT_ak4->Integral());
    h_ratio_LL_ak8->Scale(1./h_ratio_LL_ak8->Integral());
    h_ratio_TT_ak8->Scale(1./h_ratio_TT_ak8->Integral());
    h_ratio_LL_ak4->Scale(1./h_ratio_LL_ak4->Integral());
    h_ratio_TT_ak4->Scale(1./h_ratio_TT_ak4->Integral());


    // --- Style settings ---
    h_zg_LL_ak8->SetLineColor(kRed);
    h_zg_TT_ak8->SetLineColor(kBlue);
    h_zg_LL_ak4->SetLineColor(kRed+2);
    h_zg_TT_ak4->SetLineColor(kBlue+2);

    h_ratio_LL_ak8->SetLineColor(kRed);
    h_ratio_TT_ak8->SetLineColor(kBlue);
    h_ratio_LL_ak4->SetLineColor(kRed+2);
    h_ratio_TT_ak4->SetLineColor(kBlue+2);

    // --- Create canvases ---
    TCanvas *c1 = new TCanvas("c1", "zg comparison", 1000, 800);
    c1->Divide(1,2);

    // zg plots (left: AK8, right: AK4)
    c1->cd(1);
    h_zg_LL_ak8->SetTitle("zg (AK8);z_{g};Events");
    h_zg_LL_ak8->Draw("HIST");
    h_zg_LL_ak8->SetStats(0);
    h_zg_TT_ak8->Draw("HIST SAME");
    TLegend *leg1 = new TLegend(0.2,0.7,0.3,0.8);
    leg1->AddEntry(h_zg_LL_ak8, "LL (AK8)", "l");
    leg1->AddEntry(h_zg_TT_ak8, "TT (AK8)", "l");
    leg1->Draw();

    c1->cd(2);
    h_zg_LL_ak4->SetTitle("zg (AK4);z_{g};Events");
    h_zg_LL_ak4->Draw("HIST");
    h_zg_LL_ak4->SetStats(0);
    h_zg_TT_ak4->Draw("HIST SAME");
    TLegend *leg2 = new TLegend(0.2,0.7,0.3,0.8);
    leg2->AddEntry(h_zg_LL_ak4, "LL (AK4)", "l");
    leg2->AddEntry(h_zg_TT_ak4, "TT (AK4)", "l");
    leg2->Draw();

    // --- Ratio plots ---
    TCanvas *c2 = new TCanvas("c2", "Subjet ratio comparison", 1000, 800);
    c2->Divide(1,2);

    c2->cd(1);
    h_ratio_LL_ak8->SetTitle("AK8 Subjet ;p^{min}_{T}/p^{max}_{T};Events");
    h_ratio_LL_ak8->Draw("HIST");
        h_ratio_LL_ak8->SetStats(0);

    h_ratio_TT_ak8->Draw("HIST SAME");
    TLegend *leg3 = new TLegend(0.2,0.7,0.3,0.8);
    leg3->AddEntry(h_ratio_LL_ak8, "ssWW LL (AK8)", "l");
    leg3->AddEntry(h_ratio_TT_ak8, "ssWW TT (AK8)", "l");
    leg3->Draw();

    c2->cd(2);
    h_ratio_LL_ak4->SetTitle(" AK4 jets ;p^{min}_{T}/p^{max}_{T};Events");
    h_ratio_LL_ak4->Draw("HIST");
    h_ratio_LL_ak4->SetStats(0);
    h_ratio_TT_ak4->Draw("HIST SAME");
    TLegend *leg4 = new TLegend(0.2,0.7,0.3,0.8);
    leg4->AddEntry(h_ratio_LL_ak4, "ssWW LL (AK4)", "l");
    leg4->AddEntry(h_ratio_TT_ak4, "ssWW TT (AK4)", "l");
    leg4->Draw();

    c1->Update();
    c2->Update();
    c1->SaveAs("./polarization_plots/zg_LL.pdf");
    c2->SaveAs("./polarization_plots/zg_TT.pdf");
}
