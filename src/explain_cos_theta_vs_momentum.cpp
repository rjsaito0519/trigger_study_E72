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
    TTreeReaderValue<Double_t> cos_theta(reader, "cos_theta"); 
    TTreeReaderValue<std::vector<TParticle>> PRM(reader, "PRM");
    
    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/explain_%s.root", OUTPUT_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    auto h_cos_theta_vs_mom_meson = new TH2D("meson", ";cos_theta;mom", 50, -1.0, 1.0, 1000, 0.0, 1000.0);
    auto h_cos_theta_vs_mom_baryon = new TH2D("baryon", ";cos_theta;mom", 50, -1.0, 1.0, 1000, 0.0, 1000.0);

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    reader.Restart();
    while (reader.Next()){
        if (*generator != conf.beam_generator) {
            h_cos_theta_vs_mom_meson->Fill(*cos_theta, (*PRM)[0].P());
            h_cos_theta_vs_mom_baryon->Fill(*cos_theta, (*PRM)[1].P());
        }
    }

    // +-------+
    // | Write |
    // +-------+
    // -- write -----
    fout.cd();
    h_cos_theta_vs_mom_meson->Write();
    h_cos_theta_vs_mom_baryon->Write();
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