// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"


void hist_setting(TH1D *h, TLegend *l, Color_t color, TString label){
    h->SetLineColor(color);
    h->SetLineWidth(1.5);
    l->AddEntry(h, label.Data(), "l");
}

void analyze(TString path, Int_t focus_pdg_code){
    Config& conf = Config::getInstance();
    
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x"); // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y"); // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);
    
    // +-----------+
    // | load file |
    // +-----------+
    auto *f = new TFile( path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << path << std::endl;
        return;
    }
    TTreeReader reader("g4hyptpc", f);
    TTreeReaderValue<Int_t> generator(reader, "generator");
    TTreeReaderValue<Int_t> decay_particle_code(reader, "decay_particle_code"); 
    TTreeReaderValue<std::vector<TParticle>> HTOF(reader, "HTOF");

    // +---------------------------------------+
    // | prepare output pdf file and root file |
    // +---------------------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];

    // -- pdf file -----
    TString pdf_name = Form("%s/img/htof_%s_%d.pdf", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);
    if (std::ifstream(pdf_name.Data())) std::remove(pdf_name.Data());

    // -- root file -----
    TString output_path = Form("%s/root/htof_%s_%d.root", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    auto h_hitpat_all      = new TH1D("hitpat_all", "hitpat_all", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_hitpat_proton   = new TH1D("hitpat_proton", "hitpat_proton", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_hitpat_piplus   = new TH1D("hitpat_piplus", "hitpat_piplus", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_hitpat_piminus  = new TH1D("hitpat_piminus", "hitpat_piminus", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_hitpat_kminus   = new TH1D("hitpat_kminus", "hitpat_kminus", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_hitpat_electron = new TH1D("hitpat_electron", "hitpat_electron", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_hitpat_muon     = new TH1D("hitpat_muon", "hitpat_muon", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);

    auto h_edep_all      = new TH1D("edep_all", "edep_all", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_edep_proton   = new TH1D("edep_proton", "edep_proton", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_edep_piplus   = new TH1D("edep_piplus", "edep_piplus", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_edep_piminus  = new TH1D("edep_piminus", "edep_piminus", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_edep_kminus   = new TH1D("edep_kminus", "edep_kminus", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_edep_electron = new TH1D("edep_electron", "edep_electron", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_edep_muon     = new TH1D("edep_muon", "edep_muon", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    
    auto h_hitpat_vs_edep = new TH2D("hitpat_vs_edep", "HTOF hitpat vs edep", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_hitpat_vs_edep_proton = new TH2D("hitpat_vs_edep_proton", "HTOF hitpat vs edep", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi = new TH1D("multiplicity", "Multiplicity", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);

    auto h_multi1_hitpat_all      = new TH1D("multi1_hitpat_all", "multi1_hitpat_all", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_multi1_hitpat_proton   = new TH1D("multi1_hitpat_proton", "multi1_hitpat_proton", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_multi1_hitpat_piplus   = new TH1D("multi1_hitpat_piplus", "multi1_hitpat_piplus", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_multi1_hitpat_piminus  = new TH1D("multi1_hitpat_piminus", "multi1_hitpat_piminus", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_multi1_hitpat_kminus   = new TH1D("multi1_hitpat_kminus", "multi1_hitpat_kminus", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_multi1_hitpat_electron = new TH1D("multi1_hitpat_electron", "multi1_hitpat_electron", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);
    auto h_multi1_hitpat_muon     = new TH1D("multi1_hitpat_muon", "multi1_hitpat_muon", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max);

    auto h_multi1_edep_all      = new TH1D("multi1_edep_all", "multi1_edep_all", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_edep_proton   = new TH1D("multi1_edep_proton", "multi1_edep_proton", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_edep_piplus   = new TH1D("multi1_edep_piplus", "multi1_edep_piplus", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_edep_piminus  = new TH1D("multi1_edep_piminus", "multi1_edep_piminus", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_edep_kminus   = new TH1D("multi1_edep_kminus", "multi1_edep_kminus", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_edep_electron = new TH1D("multi1_edep_electron", "multi1_edep_electron", conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_edep_muon     = new TH1D("multi1_edep_muon", "multi1_edep_muon", conf.edep_bin_num, conf.edep_min, conf.edep_max);
        
    auto h_multi1_hitpat_vs_edep = new TH2D("multi1_hitpat_vs_edep", "HTOF (multi < 2) hitpat vs edep", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.edep_bin_num, conf.edep_min, conf.edep_max);
    auto h_multi1_hitpat_vs_edep_proton = new TH2D("multi1_hitpat_vs_edep_proton", "HTOF (multi < 2) hitpat vs edep", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.edep_bin_num, conf.edep_min, conf.edep_max);
    

    auto h_hitpos_all      = new TH2D("hitpos_all", "hitpos_all", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.hitpos_bin_num, conf.hitpos_min, conf.hitpos_max);
    auto h_hitpos_proton   = new TH2D("hitpos_proton", "hitpos_proton", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.hitpos_bin_num, conf.hitpos_min, conf.hitpos_max);
    auto h_multi1_hitpos_all      = new TH2D("multi1_hitpos_all", "multi1_hitpos_all", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.hitpos_bin_num, conf.hitpos_min, conf.hitpos_max);
    auto h_multi1_hitpos_proton   = new TH2D("multi1_hitpos_proton", "multi1_hitpos_proton", conf.max_htof_ch, conf.htof_seg_min, conf.htof_seg_max, conf.hitpos_bin_num, conf.hitpos_min, conf.hitpos_max);


    // +----------------------+
    // | check and fill event |
    // +----------------------+
    reader.Restart();
    while (reader.Next()){
        if (*generator != conf.beam_generator && (*decay_particle_code == focus_pdg_code || focus_pdg_code == 9999)) {
            std::set<Int_t> htof_seg_unique;
            for (const auto& item : (*HTOF)) {
                if (item.GetWeight() > conf.edep_threshold) {
                    htof_seg_unique.insert(item.GetMother(1));
                    h_hitpat_all->Fill(item.GetMother(1));
                    h_edep_all->Fill(item.GetWeight());
                    h_hitpat_vs_edep->Fill(item.GetMother(1), item.GetWeight());
                    h_hitpos_all->Fill(item.GetMother(1), item.Vy());

                    if (item.GetPdgCode() == 2212) { // proton
                        h_hitpat_proton->Fill(item.GetMother(1));
                        h_edep_proton->Fill(item.GetWeight());
                        h_hitpos_proton->Fill(item.GetMother(1), item.Vy());
                        h_hitpat_vs_edep_proton->Fill(item.GetMother(1), item.GetWeight());
                    }
                    else if (item.GetPdgCode() == 211) { // pi+
                        h_hitpat_piplus->Fill(item.GetMother(1));
                        h_edep_piplus->Fill(item.GetWeight());
                    } 
                    else if (item.GetPdgCode() == -211) { // pi-
                        h_hitpat_piminus->Fill(item.GetMother(1));
                        h_edep_piminus->Fill(item.GetWeight());
                    } 
                    else if (item.GetPdgCode() == -321) { // K-
                        h_hitpat_kminus->Fill(item.GetMother(1));
                        h_edep_kminus->Fill(item.GetWeight());
                    } 
                    else if (item.GetPdgCode() == 11 || item.GetPdgCode() == -11) { // electron
                        h_hitpat_electron->Fill(item.GetMother(1));
                        h_edep_electron->Fill(item.GetWeight());
                    }
                    else if (item.GetPdgCode() == 13 || item.GetPdgCode() == -13) { // muon
                        h_hitpat_muon->Fill(item.GetMother(1));
                        h_edep_muon->Fill(item.GetWeight());
                    }
                }
            }
            Int_t htof_multi = htof_seg_unique.size();
            h_multi->Fill(htof_multi);

            // -- HTOF mp 1 -----
            if (htof_multi < 2) {
                for (const auto& item : (*HTOF)) {
                    if (item.GetWeight() > conf.edep_threshold) {
                        h_multi1_hitpat_all->Fill(item.GetMother(1));
                        h_multi1_edep_all->Fill(item.GetWeight());
                        h_multi1_hitpat_vs_edep->Fill(item.GetMother(1), item.GetWeight());
                        h_multi1_hitpos_all->Fill(item.GetMother(1), item.Vy());

                        if (item.GetPdgCode() == 2212) { // proton
                            h_multi1_hitpat_proton->Fill(item.GetMother(1));
                            h_multi1_edep_proton->Fill(item.GetWeight());
                            h_multi1_hitpos_proton->Fill(item.GetMother(1), item.Vy());
                            h_multi1_hitpat_vs_edep_proton->Fill(item.GetMother(1), item.GetWeight());
                        }
                        else if (item.GetPdgCode() == 211) { // pi+
                            h_multi1_hitpat_piplus->Fill(item.GetMother(1));
                            h_multi1_edep_piplus->Fill(item.GetWeight());
                        } 
                        else if (item.GetPdgCode() == -211) { // pi-
                            h_multi1_hitpat_piminus->Fill(item.GetMother(1));
                            h_multi1_edep_piminus->Fill(item.GetWeight());
                        } 
                        else if (item.GetPdgCode() == -321) { // K-
                            h_multi1_hitpat_kminus->Fill(item.GetMother(1));
                            h_multi1_edep_kminus->Fill(item.GetWeight());
                        } 
                        else if (item.GetPdgCode() == 11 || item.GetPdgCode() == -11) { // electron
                            h_multi1_hitpat_electron->Fill(item.GetMother(1));
                            h_multi1_edep_electron->Fill(item.GetWeight());
                        }
                        else if (item.GetPdgCode() == 13 || item.GetPdgCode() == -13) { // muon
                            h_multi1_hitpat_muon->Fill(item.GetMother(1));
                            h_multi1_edep_muon->Fill(item.GetWeight());
                        }
                    }
                }
            }
        }
    }


    // +----------------+
    // | Draw histogram |
    // +----------------+

    // -- HTOF hitpat -----
    TCanvas *c_hitpat = new TCanvas("", "", 1500, 1200);
    c_hitpat->Print(pdf_name + "[");
    c_hitpat->cd(1);
    TLegend *l_hitpat = new TLegend( 0.7, 0.58, 0.9, 0.88);
    l_hitpat->SetFillStyle(0);
    
    hist_setting(h_hitpat_all, l_hitpat, kBlack, "all");
    h_hitpat_all->Draw();

    hist_setting(h_hitpat_proton, l_hitpat, kRed, "proton");
    h_hitpat_proton->Draw("same");

    hist_setting(h_hitpat_piplus, l_hitpat, kOrange, "pi+");
    h_hitpat_piplus->Draw("same");

    hist_setting(h_hitpat_piminus, l_hitpat, kBlue, "pi-");
    h_hitpat_piminus->Draw("same");

    hist_setting(h_hitpat_kminus, l_hitpat, kViolet, "K-");
    h_hitpat_kminus->Draw("same");

    hist_setting(h_hitpat_electron, l_hitpat, kGreen, "e-, e+");
    h_hitpat_electron->Draw("same");

    hist_setting(h_hitpat_muon, l_hitpat, kCyan, "muon");
    h_hitpat_muon->Draw("same");

    l_hitpat->Draw("same");
    c_hitpat->Print(pdf_name);


    // -- HTOF hitpat -----
    TCanvas *c_multi1_hitpat = new TCanvas("", "", 1500, 1200);
    c_multi1_hitpat->cd(1);
    TLegend *l_multi1_hitpat = new TLegend( 0.7, 0.58, 0.9, 0.88);
    l_multi1_hitpat->SetFillStyle(0);
    
    hist_setting(h_multi1_hitpat_all, l_multi1_hitpat, kBlack, "all");
    h_multi1_hitpat_all->Draw();

    hist_setting(h_multi1_hitpat_proton, l_multi1_hitpat, kRed, "proton");
    h_multi1_hitpat_proton->Draw("same");

    hist_setting(h_multi1_hitpat_piplus, l_multi1_hitpat, kOrange, "pi+");
    h_multi1_hitpat_piplus->Draw("same");

    hist_setting(h_multi1_hitpat_piminus, l_multi1_hitpat, kBlue, "pi-");
    h_multi1_hitpat_piminus->Draw("same");

    hist_setting(h_multi1_hitpat_kminus, l_multi1_hitpat, kViolet, "K-");
    h_multi1_hitpat_kminus->Draw("same");

    hist_setting(h_multi1_hitpat_electron, l_multi1_hitpat, kGreen, "e-, e+");
    h_multi1_hitpat_electron->Draw("same");

    hist_setting(h_multi1_hitpat_muon, l_multi1_hitpat, kCyan, "muon");
    h_multi1_hitpat_muon->Draw("same");

    l_multi1_hitpat->Draw("same");
    c_multi1_hitpat->Print(pdf_name);

    // -- HTOF edep -----
    TCanvas *c_edep = new TCanvas("", "", 1500, 1200);
    c_edep->cd(1);
    TLegend *l_edep = new TLegend( 0.7, 0.58, 0.9, 0.88);
    l_edep->SetFillStyle(0);
    
    hist_setting(h_edep_all, l_edep, kBlack, "all");
    h_edep_all->Draw();

    hist_setting(h_edep_proton, l_edep, kRed, "proton");
    h_edep_proton->Draw("same");

    hist_setting(h_edep_piplus, l_edep, kOrange, "pi+");
    h_edep_piplus->Draw("same");

    hist_setting(h_edep_piminus, l_edep, kBlue, "pi-");
    h_edep_piminus->Draw("same");

    hist_setting(h_edep_kminus, l_edep, kViolet, "K-");
    h_edep_kminus->Draw("same");
    
    hist_setting(h_edep_electron, l_edep, kGreen, "e-, e+");
    h_edep_electron->Draw("same");

    hist_setting(h_edep_muon, l_edep, kCyan, "muon");
    h_edep_muon->Draw("same");

    l_edep->Draw("same");
    c_edep->Print(pdf_name);


    // -- HTOF edep -----
    TCanvas *c_multi1_edep = new TCanvas("", "", 1500, 1200);
    c_multi1_edep->cd(1);
    TLegend *l_multi1_edep = new TLegend( 0.7, 0.58, 0.9, 0.88);
    l_multi1_edep->SetFillStyle(0);
    
    hist_setting(h_multi1_edep_all, l_multi1_edep, kBlack, "all");
    h_multi1_edep_all->Draw();

    hist_setting(h_multi1_edep_proton, l_multi1_edep, kRed, "proton");
    h_multi1_edep_proton->Draw("same");

    hist_setting(h_multi1_edep_piplus, l_multi1_edep, kOrange, "pi+");
    h_multi1_edep_piplus->Draw("same");

    hist_setting(h_multi1_edep_piminus, l_multi1_edep, kBlue, "pi-");
    h_multi1_edep_piminus->Draw("same");

    hist_setting(h_multi1_edep_kminus, l_multi1_edep, kViolet, "K-");
    h_multi1_edep_kminus->Draw("same");
    
    hist_setting(h_multi1_edep_electron, l_multi1_edep, kGreen, "e-, e+");
    h_multi1_edep_electron->Draw("same");

    hist_setting(h_multi1_edep_muon, l_multi1_edep, kCyan, "muon");
    h_multi1_edep_muon->Draw("same");

    l_multi1_edep->Draw("same");
    c_multi1_edep->Print(pdf_name);


    // -- HTOF hitpat vs edep -----
    TCanvas *c_hitpat_vs_edep = new TCanvas("", "", 1500, 1200);
    c_hitpat_vs_edep->cd(1);
    h_hitpat_vs_edep->Draw("colz");
    c_hitpat_vs_edep->Print(pdf_name);


    // -- HTOF hitpat vs edep -----
    TCanvas *c_multi1_hitpat_vs_edep = new TCanvas("", "", 1500, 1200);
    c_multi1_hitpat_vs_edep->cd(1);
    h_multi1_hitpat_vs_edep->Draw("colz");
    c_multi1_hitpat_vs_edep->Print(pdf_name);

    // -- HTOF Multiplicity -----
    TCanvas *c_multi = new TCanvas("", "", 1500, 1200);
    c_multi->cd(1);
    h_multi->Draw();
    c_multi->Print(pdf_name);

    // -- HTOF hitpos -----
    TCanvas *c_hitpos = new TCanvas("", "", 1500, 1200);
    c_hitpos->cd(1);

    // h_hitpos_all->Draw();
    h_hitpos_proton->Draw("colz");

    c_hitpos->Print(pdf_name);

    // -- HTOF hitpos Mp1 -----
    TCanvas *c_multi1_hitpos = new TCanvas("", "", 1500, 1200);
    c_multi1_hitpos->cd(1);

    // h_multi1_hitpos_all->Draw();
    h_multi1_hitpos_proton->Draw("colz");

    c_multi1_hitpos->Print(pdf_name);
    c_multi1_hitpat->Print(pdf_name + "]");


    // -- close output root file -----
    fout.cd(); // 明示的にカレントディレクトリを設定
    h_hitpat_all->Write();
    h_hitpat_proton->Write();
    h_hitpat_piplus->Write();
    h_hitpat_piminus->Write();
    h_hitpat_kminus->Write();
    h_hitpat_electron->Write();
    h_hitpat_muon->Write();

    h_edep_all->Write();
    h_edep_proton->Write();
    h_edep_piplus->Write();
    h_edep_piminus->Write();
    h_edep_kminus->Write();
    h_edep_electron->Write();
    h_edep_muon->Write();

    h_hitpat_vs_edep->Write();
    h_hitpat_vs_edep_proton->Write();
    h_multi->Write();

    h_multi1_hitpat_all->Write();
    h_multi1_hitpat_proton->Write();
    h_multi1_hitpat_piplus->Write();
    h_multi1_hitpat_piminus->Write();
    h_multi1_hitpat_kminus->Write();
    h_multi1_hitpat_electron->Write();
    h_multi1_hitpat_muon->Write();

    h_multi1_edep_all->Write();
    h_multi1_edep_proton->Write();
    h_multi1_edep_piplus->Write();
    h_multi1_edep_piminus->Write();
    h_multi1_edep_kminus->Write();
    h_multi1_edep_electron->Write();
    h_multi1_edep_muon->Write();

    h_multi1_hitpat_vs_edep->Write();
    h_multi1_hitpat_vs_edep_proton->Write();

    h_hitpos_all->Write();
    h_hitpos_proton->Write();
    h_multi1_hitpos_all->Write();
    h_multi1_hitpos_proton->Write();
    

}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file path> [focus_pdg_code]" << std::endl;
        return 1;
    }
    TString path = argv[1];
    
    // -- for beam data -----
    if (path.Contains("beam")) {
        conf.beam_initialize();
        std::cout << " ---" << std::endl;
    }

    Int_t focus_pdg_code = 9999;
    if (argc >= 3) {
        focus_pdg_code = std::atoi(argv[2]);
        std::cout << "focus_pdg_code: " << focus_pdg_code << std::endl;
    } else {
        std::cout << "No optional focus_pdg_code provided." << std::endl;
    }

    analyze(path, focus_pdg_code);
    return 0;

}