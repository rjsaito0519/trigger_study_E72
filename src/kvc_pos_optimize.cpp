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

void analyze(TString path_beam, TString path_event){
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
    auto *f_beam = new TFile(path_beam.Data());
    TTreeReader reader_beam("g4hyptpc", f_beam);
    TTreeReaderValue<std::vector<TParticle>> KVC_beam(reader_beam, "KVC");
    Int_t total_entry_beam = reader_beam.GetEntries();

    auto *f_event = new TFile(path_event.Data());
    TTreeReader reader_event("g4hyptpc_light", f_event);
    TTreeReaderValue<Int_t> decay_particle_code(reader_event, "decay_particle_code");
    TTreeReaderValue<std::vector<TParticle>> KVC_event(reader_event, "KVC");
    Int_t total_entry_event = reader_event.GetEntries();


    // +------------+
    // | Fill event |
    // +------------+
    std::vector<Double_t> beam_x, beam_y, beam_z, beam_u, beam_v;
    Int_t evnum_beam = 0;
    reader_beam.Restart();
    while (reader_beam.Next()){ displayProgressBar(++evnum_beam, total_entry_beam);
        for(const auto& item : (*KVC_beam)) if (item.GetPdgCode() == -321) {
            beam_x.push_back(item.Vx());
            beam_y.push_back(item.Vy());
            beam_z.push_back(item.Vz());
            Double_t u = item.Px()/item.Pz();
            Double_t v = item.Py()/item.Pz();
            beam_u.push_back(u);
            beam_v.push_back(v);
        }
    }

    std::vector<Double_t> proton_x, proton_y, proton_z, proton_u, proton_v;
    Int_t evnum_event = 0;
    reader_event.Restart();
    while (reader_event.Next()){ displayProgressBar(++evnum_event, total_entry_event);
        if (*decay_particle_code == 2212) {
            for(const auto& item : (*KVC_event)) if (item.GetPdgCode() == 2212) {
                proton_x.push_back(item.Vx());
                proton_y.push_back(item.Vy());
                proton_z.push_back(item.Vz());
                Double_t u = item.Px()/item.Pz();
                Double_t v = item.Py()/item.Pz();
                proton_u.push_back(u);
                proton_v.push_back(v);
            }    
        }
    }

    // +-------------+
    // | save object |
    // +-------------+
    TString output_path = Form("%s/root/kvc_pos_optimize.root", OUTPUT_DIR.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", "");
    output_tree.Branch("beam_x", &beam_x);
    output_tree.Branch("beam_y", &beam_y);
    output_tree.Branch("beam_z", &beam_z);
    output_tree.Branch("beam_u", &beam_u);
    output_tree.Branch("beam_v", &beam_v);

    output_tree.Branch("proton_x", &proton_x);
    output_tree.Branch("proton_y", &proton_y);
    output_tree.Branch("proton_z", &proton_z);
    output_tree.Branch("proton_u", &proton_u);
    output_tree.Branch("proton_v", &proton_v);

    output_tree.Fill();
    output_tree.Write();

    fout.Close();
}


Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // -- check argments -----
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <beam root> <event root>" << std::endl;
        return 1;
    }
    TString path_beam  = argv[1];
    TString path_event = argv[2];
    
    analyze(path_beam, path_event);
    return 0;
}