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

Config& conf = Config::getInstance();

FitResult FillHist(TString path, TH1D *h)
{
    Config& conf = Config::getInstance();
    // +-----------+
    // | load file |
    // +-----------+
    auto *f = new TFile(path.Data());

    TTreeReader reader("g4hyptpc", f);
    TTreeReaderValue<std::vector<TParticle>> KVC(reader, "KVC");
    Int_t total_entry = reader.GetEntries();

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t evnum = 0;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum, total_entry);
        for(const auto& item : (*KVC)) if (item.GetPdgCode() == -321) {
            h->Fill(item.Vx());
        }
    }
    delete f;

    FitResult result;
    result.par.push_back(h->GetMean());
    result.par.push_back(h->GetStdDev());
    result.err.push_back(h->GetMeanError());
    result.err.push_back(h->GetStdDevError());

    std::cout << "finish: " << path << std::endl;
    return result;
}

void analyze(TString dir){
    
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

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    Int_t n_z_pos = 8;
    TH1D* h_mom645[n_z_pos];
    TH1D* h_mom735[n_z_pos];
    TH1D* h_mom805[n_z_pos];
    
    for (Int_t i = 0; i < n_z_pos; i++) {
        Int_t z_pos = 450 + 50*i;
        h_mom645[i] = new TH1D(Form("mom645_%d", z_pos), Form("mom645 z%d", z_pos), 1400, 0., 700.);
        h_mom735[i] = new TH1D(Form("mom735_%d", z_pos), Form("mom735 z%d", z_pos), 1400, 0., 700.);
        h_mom805[i] = new TH1D(Form("mom805_%d", z_pos), Form("mom805 z%d", z_pos), 1400, 0., 700.);
    }

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = Form("%s/root/kvc_profile.root", OUTPUT_DIR.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", "");
    Double_t mom, z_pos, mean_val, mean_err, stdev_val, stdev_err;
    output_tree.Branch("mom", &mom, "mom/D");
    output_tree.Branch("z_pos", &z_pos, "z_pos/D");
    output_tree.Branch("mean_val", &mean_val, "mean_val/D");
    output_tree.Branch("mean_err", &mean_err, "mean_err/D");
    output_tree.Branch("stdev_val", &stdev_val, "stdev_val/D");
    output_tree.Branch("stdev_err", &stdev_err, "stdev_err/D");    
    
    // +----------------+
    // | Fill histogram |
    // +----------------+
    for (Int_t i = 0; i < n_z_pos; i++) {
        z_pos = 450.0 + 50.0*i;
        {
            mom = 645.0;
            FitResult result = FillHist(Form("%s/kvc_pos_mom%.0f_%.0f.root", dir.Data(), mom, z_pos/10), h_mom645[i]);
            mean_val  = result.par[0]; mean_err  = result.err[0];
            stdev_val = result.par[1]; stdev_err = result.err[1];
            output_tree.Fill();
        }
        {
            mom = 735.0;
            FitResult result = FillHist(Form("%s/kvc_pos_mom%.0f_%.0f.root", dir.Data(), mom, z_pos/10), h_mom735[i]);
            mean_val  = result.par[0]; mean_err  = result.err[0];
            stdev_val = result.par[1]; stdev_err = result.err[1];
            output_tree.Fill();
        }
        {
            mom = 805.0;
            FitResult result = FillHist(Form("%s/kvc_pos_mom%.0f_%.0f.root", dir.Data(), mom, z_pos/10), h_mom805[i]);
            mean_val  = result.par[0]; mean_err  = result.err[0];
            stdev_val = result.par[1]; stdev_err = result.err[1];
            output_tree.Fill();
        }        
    }

    // +-------------+
    // | save object |
    // +-------------+
    for (Int_t i = 0; i < n_z_pos; i++) {
        h_mom645[i]->Write();
        h_mom735[i]->Write();
        h_mom805[i]->Write();
    }
    output_tree.Write();
    fout.Close();

}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file dir>" << std::endl;
        return 1;
    }
    TString dir = argv[1];
    
    analyze(dir);
    return 0;
}