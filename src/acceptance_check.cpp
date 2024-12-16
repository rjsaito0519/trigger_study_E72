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
    TTreeReaderValue<Int_t> decay_particle_code(reader, "decay_particle_code");

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/acceptance_%s_%d.root", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    auto h_cos_theta_raw    = new TH1D("cos_theta_raw", "cos_theta", 50, -1.0, 1.0);
    auto h_cos_theta_trig   = new TH1D("cos_theta_trig", "cos_theta", 50, -1.0, 1.0);
    auto h_cos_theta_mp2    = new TH1D("cos_theta_mp2", "cos_theta", 50, -1.0, 1.0);
    auto h_cos_theta_htofp  = new TH1D("cos_theta_htofp", "cos_theta", 50, -1.0, 1.0);

    auto h_acceptance       = new TH1D("acceptance", "cos_theta", 50, -1.0, 1.0);
    auto h_acceptance_mp2   = new TH1D("acceptance_mp2", "cos_theta", 50, -1.0, 1.0);
    auto h_acceptance_htofp = new TH1D("acceptance_htofp", "cos_theta", 50, -1.0, 1.0);    

    auto h_mom_dist = new TH1D("mom_dist", "mom_dist", 600, 600.0, 900.0);

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    reader.Restart();
    while (reader.Next()){
        if (*decay_particle_code == focus_pdg_code || focus_pdg_code == 9999) {   
            h_mom_dist->Fill(*mom_kaon_lab);
            h_cos_theta_raw->Fill(*cos_theta);
            
            if (*trig_flag != 0) {
                h_cos_theta_trig->Fill(*cos_theta);
                if (*trig_flag == 1) h_cos_theta_mp2->Fill(*cos_theta);
                else if (*trig_flag == 2) h_cos_theta_htofp->Fill(*cos_theta);
            }
        }
    }

    // +-------+
    // | Write |
    // +-------+
    // -- cal acceptance -----
    h_acceptance->Divide( h_cos_theta_trig, h_cos_theta_raw, 1, 1 );
    h_acceptance_mp2->Divide( h_cos_theta_mp2, h_cos_theta_raw, 1, 1 );
    h_acceptance_htofp->Divide( h_cos_theta_htofp, h_cos_theta_raw, 1, 1 );

    // -- write -----
    fout.cd();
    h_cos_theta_raw->Write();
    h_cos_theta_trig->Write();
    h_cos_theta_mp2->Write();
    h_cos_theta_htofp->Write();
    h_acceptance->Write();
    h_acceptance_mp2->Write();
    h_acceptance_htofp->Write();
    h_mom_dist->Write();

    fout.Close(); 
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file path>" << std::endl;
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