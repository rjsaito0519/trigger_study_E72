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
static Double_t htof_threshold = 3.0;
static std::vector<Int_t> forward_seg;
static Int_t htof_multi_threshold = 2;

void print_result(std::unordered_map<Int_t, Int_t>& num){
    std::cout << "\n----------" << std::endl;
    for (const auto& pair : num) {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }
}

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
    TTreeReaderValue<Int_t> decay_particle_code(reader, "decay_particle_code");
    TTreeReaderValue<std::vector<TParticle>>  BH2(reader,  "BH2");
    TTreeReaderValue<std::vector<TParticle>>  BAC(reader,  "BAC");
    TTreeReaderValue<std::vector<TParticle>> HTOF(reader, "HTOF");
    TTreeReaderValue<std::vector<TParticle>>  KVC(reader,  "KVC");
    Int_t total_entry = reader.GetEntries();

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t n_kaon_all = 0, n_trig_all = 0;
    std::unordered_map<Int_t, Int_t> n_kaon;
    std::unordered_map<Int_t, Int_t> n_trig;

    Int_t n_trig_mp2_all = 0, n_trig_htofp_all = 0;
    std::unordered_map<Int_t, Int_t> n_trig_mp2;
    std::unordered_map<Int_t, Int_t> n_trig_htofp;

    Int_t evnum = 0;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum, total_entry);
        // -- kaon beam -----
        Bool_t is_kaon_at_bac = false, kaon_beam = false;   
        std::set<Int_t> bh2_seg_unique;
        for(const auto& item : (*BAC)) {
            // if (item.GetPdgCode() == -321) is_kaon_at_bac = true;
            // -- calc beta -----
            TParticlePDG *particle = pdg_database->GetParticle(item.GetPdgCode());
            if (!particle) continue;
            Double_t mass = particle->Mass()*1000.0; // MeV/c^2
            Double_t mom  = item.P();                // MeV/c^2
            Double_t beta = mom / TMath::Sqrt( mass*mass + mom*mom );
            if (beta < 1.0/conf.refractive_index_bac) is_kaon_at_bac = true;
        }
        for(const auto& item : (*BH2)) 
            if (item.GetWeight() >conf.edep_threshold*conf.counter_thickness.at("bh2")) 
                bh2_seg_unique.insert(item.GetMother(1));
        Int_t bh2_multi = bh2_seg_unique.size();
        if ( is_kaon_at_bac && bh2_multi != 0 ) kaon_beam = true;

        // -- HTOF -----
        std::set<Int_t> htof_seg_unique;
        Bool_t is_proton_forward_htof = false;
        for(const auto& item : (*HTOF))  {
            if (item.GetWeight() > conf.edep_threshold*conf.counter_thickness.at("htof")) htof_seg_unique.insert(item.GetMother(1));
            if (item.GetWeight() > htof_threshold && std::binary_search(forward_seg.begin(), forward_seg.end(), item.GetMother(1))) is_proton_forward_htof =true;
        }
        Int_t htof_multi = htof_seg_unique.size();

        // -- Cherenkov radiation at KVC -----
        Bool_t hit_anyseg_kvc = false;
        std::vector<std::pair<Bool_t, Int_t>> cherenkov_kvc;
        for(const auto& item : (*KVC)) {
            // -- calc beta -----
            TParticlePDG *particle = pdg_database->GetParticle(item.GetPdgCode());
            if (!particle) continue;
            Double_t mass = particle->Mass()*1000.0; // MeV/c^2
            Double_t mom  = item.P();                // MeV/c^2
            Double_t beta = mom / TMath::Sqrt( mass*mass + mom*mom );
            if (beta > 1.0/conf.refractive_index_kvc &&  gRandom->Uniform(0.0, 1.0) < conf.efficiency_kvc) hit_anyseg_kvc = true;
            cherenkov_kvc.emplace_back( beta > 1.0/conf.refractive_index_kvc, item.GetMother(1) );
        }

        // -- check trigger -------
        if (kaon_beam) {
            n_kaon_all++;
            n_kaon[*decay_particle_code]++;
            if ((htof_multi >= htof_multi_threshold || is_proton_forward_htof) && !hit_anyseg_kvc) {
                n_trig_all++;
                n_trig[*decay_particle_code]++;
            }

            // -- mp2 -----
            if (htof_multi >= htof_multi_threshold && !hit_anyseg_kvc) {
                n_trig_mp2_all++;
                n_trig_mp2[*decay_particle_code]++;
            }

            // -- HTOF p -----
            if (is_proton_forward_htof && !hit_anyseg_kvc) {
                n_trig_htofp_all++;
                n_trig_htofp[*decay_particle_code]++;
            }
        }

        // if (n_kaon_all >= 19000) break; // 1 sec (Kaon)
    }

    // +--------------+
    // | Print result |
    // +--------------+
    // auto kaon_intensity = new TF1("kaon_intensity", "328.860759*x - 202920.253", 600., 900.); // fitting result (unit: /spill)
    auto kaon_intensity = new TF1("kaon_intensity", "137.962111*exp(x / 126.840771) - 6979.77790", 600., 1000.); // fitting result (unit: /spill)
    Ssiz_t mom_pos = path.Index("mom");
    Double_t mom_value = 735.0;
    if (mom_pos != kNPOS) {
        mom_pos += 3;
        Ssiz_t len = 0;
        while (mom_pos + len < path.Length() && std::isdigit(path[mom_pos + len])) {
            ++len;
        }
        if (len > 0) {
            TString mom_str = path(mom_pos, len);
            mom_value = mom_str.Atof();
        }
        std::cout << "Momentum value: " << mom_value << std::endl;
    } else {
        std::cerr << "No 'mom' found. Use 735 MeV/c" << std::endl;
    }


    // Double_t Kaon_rate2024_at_735MeV = 19.2; // kHz
    Double_t Kaon_rate2024 = kaon_intensity->Eval(mom_value)/2000.0;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "total entries: " << total_entry << std::endl;
    std::cout << "n_kaon: " << n_kaon_all << std::endl;
    std::cout << "n_trig: " << n_trig_all << std::endl;
    std::cout << "rate: " << Kaon_rate2024 * static_cast<Double_t>(n_trig_all) / static_cast<Double_t>(n_kaon_all) << "\n" << std::endl;
 
    std::cout << "n_trig mp2: " << n_trig_mp2_all << std::endl;
    std::cout << "n_trig f-p: " << n_trig_htofp_all << std::endl;
    std::cout << "------------------------------------------" << std::endl;

    std::vector<Int_t> tmp_n_kaon, tmp_n_trig, tmp_n_trig_mp2, tmp_n_trig_htofp, tmp_pdg_code;
    for (Int_t pdg_code : {0, 13, 111, -211}) {
        tmp_n_kaon.push_back(n_kaon[pdg_code]);
        tmp_n_trig.push_back(n_trig[pdg_code]);
        tmp_n_trig_mp2.push_back(n_trig_mp2[pdg_code]);
        tmp_n_trig_htofp.push_back(n_trig_htofp[pdg_code]);
        tmp_pdg_code.push_back(pdg_code);

        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "decay particle code: " << pdg_code << std::endl;
        std::cout << "n_kaon: " << n_kaon[pdg_code] << std::endl;
        std::cout << "n_trig: " << n_trig[pdg_code] << std::endl;
        std::cout << "rate: " << Kaon_rate2024 * static_cast<Double_t>(n_trig[pdg_code]) / static_cast<Double_t>(n_kaon_all) << "\n" << std::endl;
    
        std::cout << "n_trig mp2: " << n_trig_mp2[pdg_code] << std::endl;
        std::cout << "n_trig f-p: " << n_trig_htofp[pdg_code] << std::endl;
        std::cout << "------------------------------------------" << std::endl;
    }

    // print_result(n_kaon);
    // print_result(n_trig);
    // print_result(n_trig_mp2);
    // print_result(n_trig_htofp);

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/kvc_pos_optimize/%s.root", OUTPUT_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "RECREATE");
    TTree output_tree("tree", "");
    output_tree.Branch("kvc_eff", &conf.efficiency_kvc, "kvc_eff/D");
    output_tree.Branch("kaon_rate2024", &Kaon_rate2024, "Kaon_rate2024/D");
    output_tree.Branch("n_kaon_all", &n_kaon_all, "n_kaon_all/I");
    output_tree.Branch("n_trig_all", &n_trig_all, "n_trig_all/I");
    output_tree.Branch("n_trig_mp2_all", &n_trig_mp2_all, "n_trig_mp2_all/I");
    output_tree.Branch("n_trig_htofp_all", &n_trig_htofp_all, "n_trig_htofp_all/I");

    output_tree.Branch("n_kaon", &tmp_n_kaon);
    output_tree.Branch("n_trig", &tmp_n_trig);
    output_tree.Branch("n_trig_mp2", &tmp_n_trig_mp2);
    output_tree.Branch("n_trig_htofp", &tmp_n_trig_htofp);
    output_tree.Branch("pdg_code", &tmp_pdg_code);
    
    output_tree.Fill();   
    output_tree.Write();
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    forward_seg = conf.forward_seg_wide;

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file path>" << std::endl;
        return 1;
    }
    TString path = argv[1];
    analyze(path);
    return 0;
}