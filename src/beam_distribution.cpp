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

static std::vector<Double_t> n_kaon_container;

void data_fill(TString path, TH1D *h, Double_t factor)
{
    Config& conf = Config::getInstance();
    // +-----------+
    // | load file |
    // +-----------+
    auto *f = new TFile(path.Data());

    TTreeReader reader("g4hyptpc", f);
    TTreeReaderValue<Int_t> generator(reader, "generator");
    TTreeReaderValue<std::vector<TParticle>> BH2(reader, "BH2");
    TTreeReaderValue<std::vector<TParticle>> BAC(reader, "BAC");
    TTreeReaderValue<std::vector<TParticle>> TGT(reader, "TGT");

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t n_points = conf.n_mom_points;
    std::vector<Double_t> tmp_n_kaon_container(n_points, 0.0);
    Int_t tot_num = 0;
    reader.Restart();
    while (reader.Next()){
        // -- kaon beam -----
        Bool_t is_kaon_at_bac = false, kaon_beam = false;   
        std::set<Int_t> bh2_seg_unique;
        for(const auto& item : (*BAC)) if (item.GetPdgCode() == -321) is_kaon_at_bac = true;
        for(const auto& item : (*BH2)) if (item.GetWeight() > conf.edep_threshold) bh2_seg_unique.insert(item.GetMother(1));;
        Int_t bh2_multi = bh2_seg_unique.size();
        if ( is_kaon_at_bac && bh2_multi != 0 ) kaon_beam = true;

        if (kaon_beam) {
            if ((*TGT).size() > 0) {
                const auto& item = (*TGT)[0];
                if (item.GetPdgCode() == -321) {
                    Double_t mom = item.P();
                    h->Fill(mom);
                    Int_t index = ana_helper::get_index(mom);
                    if (index != -1) {
                        tmp_n_kaon_container[index]++;
                        tot_num++;
                    }
                }
            }
        }
    }

    for (Int_t i = 0; i < n_points; i++) n_kaon_container[i] += tmp_n_kaon_container[i] * factor / tot_num;
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

    auto kaon_intensity = new TF1("kaon_intensity", "328.860759*x - 202920.253", 600., 900.); // fitting result (unit: /spill)

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    auto h_mom685 = new TH1D("mom685", "mom685", 600, 600., 900.);
    auto h_mom705 = new TH1D("mom705", "mom705", 600, 600., 900.);
    auto h_mom725 = new TH1D("mom725", "mom725", 600, 600., 900.);
    auto h_mom745 = new TH1D("mom745", "mom745", 600, 600., 900.);
    auto h_mom765 = new TH1D("mom765", "mom765", 600, 600., 900.);
    auto h_mom735 = new TH1D("mom735", "mom735", 600, 600., 900.);    
    auto h_momscan = new TH1D("mom_scan", "Kaon momentum", 600, 600., 900.);
    auto h_momall = new TH1D("mom_all", "Kaon momentum", 600, 600., 900.);


    // +----------------+
    // | Fill histogram |
    // +----------------+
    Double_t t_measure_scan = 0.5*24.0*3600.0; // unit: sec
    Double_t t_measure_735  = 5.5*24.0*3600.0; // unit: sec
    
    Double_t factor = kaon_intensity->Eval(685.0) * t_measure_scan / conf.spill_length;
    data_fill(Form("%s/beam_mom685.root", dir.Data()), h_mom685, factor);
    h_mom685->Scale( factor / h_mom685->GetEntries() );
    h_momscan->Add(h_mom685, 1.0);
    h_momall->Add(h_mom685, 1.0);

    factor = kaon_intensity->Eval(705.0) * t_measure_scan / conf.spill_length;
    data_fill(Form("%s/beam_mom705.root", dir.Data()), h_mom705, factor);
    h_mom705->Scale( factor / h_mom705->GetEntries() );
    h_momscan->Add(h_mom705, 1.0);
    h_momall->Add(h_mom705, 1.0);

    factor = kaon_intensity->Eval(725.0) * t_measure_scan / conf.spill_length;
    data_fill(Form("%s/beam_mom725.root", dir.Data()), h_mom725, factor);
    h_mom725->Scale( factor / h_mom725->GetEntries() );
    h_momscan->Add(h_mom725, 1.0);
    h_momall->Add(h_mom725, 1.0);

    factor = kaon_intensity->Eval(745.0) * t_measure_scan / conf.spill_length;
    data_fill(Form("%s/beam_mom745.root", dir.Data()), h_mom745, factor);
    h_mom745->Scale( factor / h_mom745->GetEntries() );
    h_momscan->Add(h_mom745, 1.0);
    h_momall->Add(h_mom745, 1.0);

    factor = kaon_intensity->Eval(765.0) * t_measure_scan / conf.spill_length;
    data_fill(Form("%s/beam_mom765.root", dir.Data()), h_mom765, factor);
    h_mom765->Scale( factor / h_mom765->GetEntries() );
    h_momscan->Add(h_mom765, 1.0);
    h_momall->Add(h_mom765, 1.0);

    factor = kaon_intensity->Eval(735.0) * t_measure_735 / conf.spill_length;
    data_fill(Form("%s/beam_mom735.root", dir.Data()), h_mom735, factor);
    h_mom735->Scale( factor / h_mom735->GetEntries() );
    h_momall->Add(h_mom735, 1.0);


    // +-------------+
    // | save object |
    // +-------------+
    TString output_path = Form("%s/root/n_kaon.root", OUTPUT_DIR.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", "");
    output_tree.Branch("n_kaon", &n_kaon_container);
    output_tree.Fill();    
    output_tree.Write();

    h_mom685->Write();
    h_mom705->Write();
    h_mom725->Write();
    h_mom745->Write();
    h_mom765->Write();
    h_mom735->Write();
    h_momscan->Write();
    h_momall->Write();
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
    
    n_kaon_container.resize(conf.n_mom_points, 0.0);
    analyze(dir);
    return 0;
}