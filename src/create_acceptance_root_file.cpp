// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <tuple>

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
#include <TSpline.h>     
#include <TGraphErrors.h>
#include <TKey.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"


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
    TTreeReader reader("g4hyptpc_light", f);
    TTreeReaderValue<Double_t> mom_kaon_lab(reader, "mom_kaon_lab");
    TTreeReaderValue<Double_t> cos_theta(reader, "cos_theta");
    TTreeReaderValue<Int_t> trig_flag(reader, "trig_flag");
    TTreeReaderValue<Int_t> pdg_code(reader, "decay_particle_code");
    Int_t total_entry = reader.GetEntries();


    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/for_cusp/acceptance_%s_%d.root", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "recreate");


    // +-------------------+
    // | prepare histogram |
    // +-------------------+ 
    TH1D *h_cos_theta_raw[conf.n_mom_points];
    TH1D *h_cos_theta_trig[conf.n_mom_points];
    TH1D *h_acceptance[conf.n_mom_points];
    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        h_cos_theta_raw[i]  = new TH1D(
            Form("cos_theta_raw%d", i),
            Form("%f , %f", conf.mom_start+i*conf.mom_step_size, conf.mom_start+(i+1)*conf.mom_step_size), 
            conf.cos_theta_bin_num, -1.0, 1.0
        );
        h_cos_theta_trig[i] = new TH1D(
            Form("cos_theta_trig%d", i),
            Form("%f , %f", conf.mom_start+i*conf.mom_step_size, conf.mom_start+(i+1)*conf.mom_step_size),
            conf.cos_theta_bin_num, -1.0, 1.0
        );
        h_acceptance[i] = new TH1D(
            Form("acceptance%d", i),
            Form("%f , %f", conf.mom_start+i*conf.mom_step_size, conf.mom_start+(i+1)*conf.mom_step_size),
            conf.cos_theta_bin_num, -1.0, 1.0
        );
    }

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t evnum = 0;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum, total_entry);
        if (focus_pdg_code == 9999 || *pdg_code == focus_pdg_code) {
            Int_t index = ana_helper::get_index(*mom_kaon_lab);
            if (index != -1) {
                h_cos_theta_raw[index]->Fill(*cos_theta);
                if (*trig_flag != 0) h_cos_theta_trig[index]->Fill(*cos_theta);
            }
        }
    }

    // +-----------+
    // | Print PDF |
    // +-----------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2;
    Int_t cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_name = Form("%s/img/for_cusp/acceptance_%s_%d.pdf", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);


    auto *c = new TCanvas("acceptance", "", 1500, 1200);
    c->Divide(cols, rows);
    c->Print(pdf_name + "["); // start
    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        h_cos_theta_raw[i]->SetLineColor(kBlack);
        h_cos_theta_raw[i]->GetYaxis()->SetRangeUser(0.0, h_cos_theta_raw[i]->GetMaximum()*1.2);
        h_cos_theta_raw[i]->Draw();
        h_cos_theta_trig[i]->SetLineColor(kRed);
        h_cos_theta_trig[i]->Draw("same");

        c->cd(++nth_pad);
        h_acceptance[i]->Divide( h_cos_theta_trig[i], h_cos_theta_raw[i], 1.0, 1.0 );
        h_acceptance[i]->SetLineColor(kBlack);
        h_acceptance[i]->GetYaxis()->SetRangeUser(0.0, 1.0);
        h_acceptance[i]->Draw();

        auto h_smooth_acceptance = (TH1D*)h_acceptance[i]->Clone();
        h_smooth_acceptance->Smooth(5);
        h_smooth_acceptance->SetLineColor(kRed);
        h_smooth_acceptance->Draw("same");
        
        nth_pad++;
    }
    c->Print(pdf_name);
    c->Print(pdf_name + "]"); // end
    delete c;

    // +------------+
    // | Write data |
    // +------------+
    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        h_cos_theta_raw[i]->Write();
        h_cos_theta_trig[i]->Write();
        h_acceptance[i]->Write();
    }
    fout.Close(); 

}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    // =======================================================
    // Program Arguments:
    // argv[1]: root file path (required)
    // argv[2]: focus PDG code (optional, default: 9999)
    // =======================================================

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path> [focus_pdg_code]" << std::endl;
        return 1;
    }

    // Load parameters
    TString path = argv[1];
    Int_t focus_pdg_code = (argc > 2) ? std::atoi(argv[2]) : 9999; // Default: 9999

    // Output parameters
    std::cout << "\n====================================" << std::endl;
    std::cout << "Path:           " << path << std::endl;
    std::cout << "Focus PDG Code: " << focus_pdg_code << std::endl;
    std::cout << "====================================\n" << std::endl;

    analyze(path, focus_pdg_code);
    return 0;
}