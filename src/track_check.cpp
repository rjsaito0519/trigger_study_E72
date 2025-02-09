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
#include <TH2Poly.h>
#include <TGraph.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

static const TDatabasePDG *pdg_database = new TDatabasePDG();
static Double_t htof_threshold = 3.0;
static std::vector<Int_t> forward_seg;
static Int_t htof_multi_threshold = 2;

void analyze(TString path, Int_t focus_pdg_code, Int_t n_rand){
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
    TTreeReaderValue<Int_t> decay_particle_code(reader, "decay_particle_code");
    TTreeReaderValue<std::vector<TParticle>>  BH2(reader,  "BH2");
    TTreeReaderValue<std::vector<TParticle>>  BAC(reader,  "BAC");
    TTreeReaderValue<std::vector<TParticle>>  TPC(reader,  "TPC");
    TTreeReaderValue<std::vector<TParticle>> HTOF(reader, "HTOF");
    TTreeReaderValue<std::vector<TParticle>>  KVC(reader,  "KVC");
    Int_t total_entry = reader.GetEntries();

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Bool_t beam_flag = false;
    std::vector<TParticle> beam_bac;
    std::vector<TParticle> beam_bh2;
    Int_t n = n_rand % total_entry;
    reader.Restart();
    reader.SetEntry(n-1);
    while (reader.Next()){
        if (*generator == conf.beam_generator){
            beam_bac = *BAC;
            beam_bh2 = *BH2;   
            beam_flag = true;         
        } else {
            if (beam_flag && (*cos_theta > 0.8) && (*decay_particle_code == focus_pdg_code || focus_pdg_code == 9999)) {

                std::vector<Double_t> x;
                std::vector<Double_t> z;

                m_canvas_tpc->Clear();
                auto *m_tpc_adc2d = new TH2Poly("h_tpc_adc2d", "TPC ADC;Z;X", DrawMinZ, DrawMaxZ, DrawMinX, DrawMaxX);
                
                TPC_pad_template(m_tpc_adc2d);
                
                for (Int_t i = 0; i < *nhTpc; i++) {
                    Int_t pad = padHelper::getPadID( (*layerTpc)[i], (*rowTpc)[i] );
                    Bool_t found = std::binary_search(padHelper::noisy_pad.begin(), padHelper::noisy_pad.end(), pad);
                    if (found && removeNoisyPad) continue;
                    m_tpc_adc2d->SetBinContent(pad+1, (*deTpc)[i]);
                }
            

                // -- kaon beam -----
                for(const auto& item : (beam_bac)) {
                    x.push_back(item.Vx());
                    z.push_back(item.Vz());
                }
                
                for(const auto& item : (beam_bh2)) {
                    x.push_back(item.Vx());
                    z.push_back(item.Vz());
                }

                // -- TPC -----
                for(const auto& item : (*TPC)) {
                    x.push_back(item.Vx());
                    z.push_back(item.Vz());
                }

                // -- HTOF -----
                for(const auto& item : (*HTOF)) {
                    x.push_back(item.Vx());
                    z.push_back(item.Vz());
                }

                // -- Cherenkov radiation at KVC -----
                for(const auto& item : (*KVC)) {
                    x.push_back(item.Vx());
                    z.push_back(item.Vz());
                }

                auto *c = new TCanvas("", "", 800, 800);
                c->cd(1);
                auto g = new TGraph(z.size(), z.data(), x.data());
                g->GetXaxis()->SetRangeUser(-750, 500);
                g->GetYaxis()->SetRangeUser(-350, 350);

                g->SetTitle(Form("cos_theta %f;z axis;x axis", *cos_theta)); // タイトルと軸ラベル
                g->SetMarkerStyle(20); // マーカーの形状（丸）
                g->SetMarkerSize(1.0); // マーカーサイズ
                g->SetMarkerColor(kBlue); // マーカー色
                g->SetLineColor(kRed); // 線の色

                g->Draw("AP");
                break;
            }
        }
    }

}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    forward_seg = conf.forward_seg_wide;

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file path> [focus_pdg_code]" << std::endl;
        return 1;
    }
    TString path = argv[1];
    
    Int_t focus_pdg_code = 9999;
    if (argc >= 3) {
        focus_pdg_code = std::atoi(argv[2]);
        std::cout << "focus_pdg_code: " << focus_pdg_code << std::endl;
    } else {
        std::cout << "No optional focus_pdg_code provided." << std::endl;
    }

    TApplication *theApp = new TApplication("App", &argc, argv);    

    std::random_device rd;
    std::mt19937 gen(rd());
    Int_t n_rand = gen();
    while (n_rand < 0) n_rand = gen();
    analyze(path, focus_pdg_code, n_rand);

    theApp->Run();
    
    return 0;

}