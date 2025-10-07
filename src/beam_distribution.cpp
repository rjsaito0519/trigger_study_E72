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

static std::vector<Double_t> n_kaon_container_scan;
static std::vector<Double_t> n_kaon_container_all;

static const TDatabasePDG *pdg_database = new TDatabasePDG();

void data_fill(TString path, TH1D *h, Double_t factor, Bool_t do_fill_scan = true)
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
    Int_t total_entry = reader.GetEntries();

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t n_points = conf.n_mom_points;
    std::vector<Double_t> tmp_n_kaon_container(n_points, 0.0);
    Int_t tot_num = 0;
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

    for (Int_t i = 0; i < n_points; i++) {
        if (do_fill_scan) n_kaon_container_scan[i] += tmp_n_kaon_container[i] * factor / tot_num;
        n_kaon_container_all[i] += tmp_n_kaon_container[i] * factor / tot_num;
    }
    delete f;

    std::cout << "finish: " << path << std::endl;
}


Double_t kaon_intensity_max() {
    return 137.962111 * TMath::Exp(884.5 / 126.840771) - 6979.77790;
}

Double_t kaon_intensity_func(Double_t* x, Double_t* /*p*/) {
    Double_t val = 137.962111 * TMath::Exp(x[0] / 126.840771) - 6979.77790;
    return TMath::Min(val, kaon_intensity_max());
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

    // auto kaon_intensity = new TF1("kaon_intensity", "328.860759*x - 202920.253", 600., 900.); // fitting result (unit: /spill)
    // auto kaon_intensity = new TF1("kaon_intensity", "137.962111*exp(x / 126.840771) - 6979.77790", 600., 1000.); // fitting result (unit: /spill)
    auto kaon_intensity = new TF1("kaon_intensity", kaon_intensity_func, 600., 1000., 0);

    // +-----------------------+
    // | Histogram preparation |
    // +-----------------------+
    auto make_hist = [](Double_t m) {
        return new TH1D(Form("mom%.0f", m), Form("mom%.0f", m), 1000, 500., 1000.);
    };
    
    auto h_momscan = new TH1D("mom_scan", "Kaon momentum (scan)", 1000, 500., 1000.);
    auto h_momall  = new TH1D("mom_all",  "Kaon momentum (all)" , 1000, 500., 1000.);

    // List of scan momenta (each with 0.5 days measurement time)
    // const std::vector<Double_t> scan_moms = {645.0, 665.0, 685.0, 705.0, 725.0, 745.0, 765.0, 785.0, 805.0, 825.0, 845.0, 865.0, 885.0, 905.0, 925.0};
    const std::vector<Double_t> scan_moms = {645.0, 665.0, 685.0, 705.0, 725.0, 745.0, 765.0, 790.0, 814.0, 842.0, 870.0, 902.0, 933.0};

    // Special case: 735 MeV/c (3.5 days, different data_fill flag, only added to "all")
    const Double_t special_mom = 735.0;

    // Measurement times [sec]
    const Double_t physics_run_days = 8.0;
    const Double_t day    = 24.0 * 3600.0;
    const Double_t hour   = 3600.0; 
    const Double_t minute = 60.0; 

    // const std::unordered_map<Double_t, Double_t> t_measure_scan = {
    //     {645.0, 12.0*hour },
    //     {665.0, 12.0*hour },
    //     {685.0, 12.0*hour },
    //     {705.0, 12.0*hour },
    //     {725.0, 11.0*hour + 20.0*minute },
    //     {745.0,  9.0*hour + 30.0*minute },
    //     {765.0,  7.0*hour + 50.0*minute },
    //     {785.0,  6.0*hour + 40.0*minute },
    //     {805.0,  5.0*hour + 30.0*minute },
    //     {825.0,  4.0*hour + 40.0*minute },
    //     {845.0,  4.0*hour },
    //     {865.0,  3.0*hour + 20.0*minute },
    //     {885.0,  2.0*hour + 50.0*minute },
    //     {905.0,  2.0*hour + 20.0*minute },
    //     {925.0,  2.0*hour },        
    // };
    const std::unordered_map<Double_t, Double_t> t_measure_scan = {
        {645.0, 12.0*hour },
        {665.0, 12.0*hour },
        {685.0, 12.0*hour },
        {705.0, 10.0*hour },                 // 10:00
        {725.0,  8.0*hour + 10.0*minute },   // 08:10
        {745.0,  6.0*hour + 50.0*minute },   // 06:50
        {765.0,  5.0*hour + 40.0*minute },   // 05:40
        {790.0,  4.0*hour + 30.0*minute },   // 04:30
        {814.0,  3.0*hour + 40.0*minute },   // 03:40
        {842.0,  2.0*hour + 50.0*minute },   // 02:50
        {870.0,  2.0*hour + 20.0*minute },   // 02:20
        {902.0,  2.0*hour },                 // 02:00
        {933.0,  2.0*hour },                 // 02:00
    };

    // --- 計算と標準出力 ---
    std::cout << "--- Individual Scan Durations ---" << std::endl;
    // 小数点以下の表示を固定
    std::cout << std::fixed << std::setprecision(4);

    Double_t total_scan_seconds = 0.0;
    for (const auto& pair : t_measure_scan) {
        // pair.first = キー, pair.second = 値（秒）
        Double_t scan_item_days = pair.second / day;
        total_scan_seconds += pair.second;
        
        // 各スキャンが何日分かを出力
        std::cout << "Scan at key " << pair.first << ": " << scan_item_days << " days" << std::endl;
    }

    Double_t total_scan_days = total_scan_seconds / day;
    const Double_t t_measure_735 = physics_run_days - total_scan_days;

    std::cout << "\n--- Summary ---" << std::endl;
    // 合計スキャン日数を出力
    std::cout << "Total scan time: " << total_scan_days << " days" << std::endl;
    // 残りの日数を出力
    std::cout << "Remaining time (t_measure_735): " << t_measure_735 << " days" << std::endl;
    

    // Store produced histograms for later writing
    std::vector<std::unique_ptr<TH1D>> produced_hists;
    produced_hists.reserve(scan_moms.size() + 1);

    // +-------------------+
    // | Fill scan spectra |
    // +-------------------+
    for (Double_t m : scan_moms) {
        auto h = std::unique_ptr<TH1D>(make_hist(m));

        const Double_t factor = kaon_intensity->Eval(m) * t_measure_scan.at(m) / conf.spill_length;

        data_fill(Form("%s/beam_mom%.0f.root", dir.Data(), m), h.get(), factor);

        const Double_t entries = h->GetEntries();
        if (entries > 0.0) h->Scale(factor / entries); // Avoid division by zero

        h_momscan->Add(h.get(), 1.0);
        h_momall ->Add(h.get(), 1.0);

        produced_hists.emplace_back(std::move(h));
    }

    // +-----------------------------+
    // | Fill special 735 MeV/c case |
    // +-----------------------------+
    {
        auto h735 = std::unique_ptr<TH1D>(make_hist(special_mom));

        const Double_t factor = kaon_intensity->Eval(special_mom) * t_measure_735 / conf.spill_length;

        // Only 735 MeV/c: last argument of data_fill set to false
        data_fill(Form("%s/beam_mom%.0f.root", dir.Data(), special_mom), h735.get(), factor, false);

        const Double_t entries = h735->GetEntries();
        if (entries > 0.0) h735->Scale(factor / entries);

        // 735 MeV/c contributes only to "all"
        h_momall->Add(h735.get(), 1.0);

        produced_hists.emplace_back(std::move(h735));
    }

    // +-------------+
    // | Save output |
    // +-------------+
    TString output_path = Form("%s/root/n_kaon.root", OUTPUT_DIR.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());

    TFile fout(output_path.Data(), "RECREATE");
    TTree output_tree("tree", "");
    output_tree.Branch("n_kaon_scan", &n_kaon_container_scan);
    output_tree.Branch("n_kaon_all",  &n_kaon_container_all);
    output_tree.Fill();
    output_tree.Write();

    // Write individual and summary histograms
    for (auto &h : produced_hists) h->Write();
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
    
    n_kaon_container_scan.resize(conf.n_mom_points, 0.0);
    n_kaon_container_all.resize(conf.n_mom_points, 0.0);
    analyze(dir);
    return 0;
}