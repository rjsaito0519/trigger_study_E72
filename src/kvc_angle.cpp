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

void data_fill(TString path, TH1D *h)
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
        for (const auto& item : (*KVC)) {
            if (item.GetPdgCode() == -321) {
                Double_t mom_x = item.Px();
                Double_t mom_z = item.Pz();
                Double_t angle = TMath::ATan(mom_x/mom_z) * TMath::RadToDeg();
                h->Fill(angle);
            }
            break;
        }
    }
    delete f;

    std::cout << "finish: " << path << std::endl;
}


void analyze(TString dir){
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

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    auto h_angle = new TH1D("angle", "angle", 1800, 0.0, 180.0);

    // +----------------+
    // | Fill histogram |
    // +----------------+
    data_fill(Form("%s/beam_mom685.root", dir.Data()), h_angle);
    data_fill(Form("%s/beam_mom705.root", dir.Data()), h_angle);
    data_fill(Form("%s/beam_mom725.root", dir.Data()), h_angle);
    data_fill(Form("%s/beam_mom745.root", dir.Data()), h_angle);
    data_fill(Form("%s/beam_mom765.root", dir.Data()), h_angle);
    data_fill(Form("%s/beam_mom735.root", dir.Data()), h_angle);

    // +-------------+
    // | save object |
    // +-------------+
    TString output_path = Form("%s/root/angle.root", OUTPUT_DIR.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", "");
    Double_t mean  = h_angle->GetMean();
    Double_t stdev = h_angle->GetStdDev();
    
    output_tree.Branch("mean", &mean, "mean/D");
    output_tree.Branch("stdev", &stdev, "stdev/D");
    output_tree.Fill();    
    output_tree.Write();

    h_angle->Write();
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