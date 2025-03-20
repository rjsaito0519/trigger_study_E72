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
#include "progress_bar.h"

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
    auto *f = new TFile(path.Data());
    TTreeReader reader("g4hyptpc", f);
    TTreeReaderValue<std::vector<TParticle>> BAC(reader, "BAC");
    TTreeReaderValue<std::vector<TParticle>> TGT(reader, "TGT");
    TTreeReaderValue<std::vector<TParticle>> KVC(reader, "KVC");
    Int_t total_entry = reader.GetEntries();


    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    auto h_BAC_profile = new TH2D("BAC_profile", "", 
        300, -300.0, 300.0, 
        300, -300.0, 300.0
    );

    auto h_TGT_profile = new TH2D("TGT_profile", "", 
        300, -300.0, 300.0, 
        300, -300.0, 300.0
    );

    auto h_KVC_profile = new TH2D("KVC_profile", "", 
        200, -400.0, 400.0, 
        250, -125.0, 125.0
    );


    // +------------+
    // | Fill event |
    // +------------+
    Int_t evnum_beam = 0;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum_beam, total_entry);
        for(const auto& item : (*BAC)) if (item.GetPdgCode() == -321) {
            h_BAC_profile->Fill(item.Vx(), item.Vy());
        }

        for(const auto& item : (*TGT)) if (item.GetPdgCode() == -321) {
            h_TGT_profile->Fill(item.Vx(), item.Vy());
        }
        
        for(const auto& item : (*KVC)) if (item.GetPdgCode() == -321) {
            h_KVC_profile->Fill(item.Vx(), item.Vy());
        }
    }


    // +-------------+
    // | save object |
    // +-------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/%s.root", OUTPUT_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    
    h_BAC_profile->Write();
    h_TGT_profile->Write();
    h_KVC_profile->Write();

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