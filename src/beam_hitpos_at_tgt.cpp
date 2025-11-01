// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <cctype>

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
#include <TRandom.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

static const TDatabasePDG *pdg_database = new TDatabasePDG();


void analyze(TString path){
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
    
    gRandom->SetSeed(72);

    // +-----------+
    // | load file |
    // +-----------+
    auto *f = new TFile( path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << path << std::endl;
        return;
    }
    TTreeReader reader("g4hyptpc", f);
    TTreeReaderValue<std::vector<TParticle>> BAC(reader,  "BAC");
    TTreeReaderValue<std::vector<TParticle>> TGT(reader,  "TGT");
    Int_t total_entry = reader.GetEntries();


    // +--------------------+
    // | prepare histograms |
    // +--------------------+
    HistPair h_xhit("x_profile", "x_profile", 200, -50.0, 50.0);

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t evnum = 0;
    Int_t n_kaon = 0;
    Int_t n_hit_kaon = 0;
    Int_t n_hit_tgt = 0;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum, total_entry);
        // -- kaon beam -----
        Bool_t is_kaon_at_bac = false;   
        for(const auto& item : (*BAC)) {
            TParticlePDG *particle = pdg_database->GetParticle(item.GetPdgCode());
            if (!particle) continue;
            Double_t mass = particle->Mass()*1000.0; // MeV/c^2
            Double_t mom  = item.P();                // MeV/c^2
            Double_t beta = mom / TMath::Sqrt( mass*mass + mom*mom );
            if (beta < 1.0/conf.refractive_index_bac) is_kaon_at_bac = true;
        }

        Bool_t do_hit_tgt = false;
        for(const auto& item : (*TGT)) if (item.GetPdgCode() == -321) {
            h_xhit.raw->Fill(item.Vx());
            if (is_kaon_at_bac) {
                h_xhit.trig->Fill(item.Vx());
            }
            do_hit_tgt = true;
        }
        if (do_hit_tgt) n_hit_tgt++;
        if (is_kaon_at_bac) {
            n_kaon++;
            if (do_hit_tgt) n_hit_kaon++;
        }
    }

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/bac_pos_optimize/%s.root", OUTPUT_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "RECREATE");
    h_xhit.raw->Write();
    h_xhit.trig->Write();

    TTree output_tree("tree", "");
    output_tree.Branch("n_kaon", &n_kaon, "n_kaon/I");
    output_tree.Branch("n_hit_kaon", &n_hit_kaon, "n_hit_kaon/I");
    output_tree.Branch("n_hit_tgt", &n_hit_tgt, "n_hit_tgt/I");
    output_tree.Fill();   
    output_tree.Write();
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
    analyze(path);
    return 0;
}