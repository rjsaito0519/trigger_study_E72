// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <string>
#include <sstream>
#include <cctype>
#include <unordered_map>

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

#include <TString.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cctype>
#include <vector>
#include <unordered_map>
#include <algorithm>

// Helper function to parse a time string like "1h30m" into seconds (Double_t).
Double_t parse_time_string(TString s) {
    Double_t total_seconds = 0.0;
    Double_t current_number = 0.0;
    TString number_buffer; // Buffer to temporarily hold number characters

    for (Int_t i = 0; i < s.Length(); ++i) {
        char c = s[i];
        // If c is a digit or a decimal point, append it to the buffer.
        if (std::isdigit(c) || c == '.') {
            number_buffer += c;
        } else {
            // If c is a unit character, convert the number in the buffer and add to total.
            if (number_buffer.Length() > 0) {
                current_number = number_buffer.Atof(); // Convert TString to Double_t
                number_buffer.Clear();
            }

            if (c == 'h' || c == 'H') total_seconds += current_number * 3600.0;
            else if (c == 'm' || c == 'M') total_seconds += current_number * 60.0;
            else if (c == 's' || c == 'S') total_seconds += current_number;
            current_number = 0.0;
        }
    }
    // Case where the string ends with a number (e.g., "3600" for seconds).
    if (number_buffer.Length() > 0) {
        total_seconds += number_buffer.Atof();
    }
    return total_seconds;
}

// Function to load parameters from the configuration file.
bool load_config(
    const TString& filename,
    Double_t& physics_run_days,
    Double_t& special_mom,
    std::vector<Double_t>& scan_moms,
    std::unordered_map<Double_t, Double_t>& t_measure_scan
) {
    // Get const char* using .Data() for ifstream.
    std::ifstream config_file(filename.Data());
    if (!config_file.is_open()) {
        std::cerr << "Error: Could not open config file: " << filename.Data() << std::endl;
        return false;
    }

    scan_moms.clear();
    t_measure_scan.clear();

    // getline requires std::string, so we receive it here first.
    std::string line_std;
    while (std::getline(config_file, line_std)) {
        TString line = line_std; // Convert from std::string to TString
        
        // Skip empty lines or lines starting with '#'.
        if (line.IsNull() || line.BeginsWith("#")) continue;

        std::stringstream ss(line.Data());
        std::string key_std;
        ss >> key_std;
        TString key = key_std;

        if (key == "physics_run_days") ss >> physics_run_days;
        else if (key == "special_mom") ss >> special_mom;
        else if (key == "scan_point") {
            Double_t mom;
            std::string time_str_std;
            ss >> mom >> time_str_std;
            
            scan_moms.push_back(mom);
            t_measure_scan[mom] = parse_time_string(TString(time_str_std));
        }
    }
    std::sort(scan_moms.begin(), scan_moms.end());
    std::cout << "Successfully loaded config file: " << filename.Data() << std::endl;
    return true;
}

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

void analyze(TString dir, TString config_filename){
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

    // --- Prepare intensity function and histograms ---
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


    // +--------------------------+
    // | Prepare momentum setting |
    // +--------------------------+

    // Measurement times [sec]
    const Double_t day    = 24.0 * 3600.0;
    const Double_t hour   = 3600.0; 
    const Double_t minute = 60.0; 

    // --- Load settings from config file ---
    Double_t physics_run_days, special_mom;
    std::vector<Double_t> scan_moms;
    std::unordered_map<Double_t, Double_t> t_measure_scan;

    if (!load_config(config_filename, physics_run_days, special_mom, scan_moms, t_measure_scan)) {
        std::cerr << "Failed to load config. Aborting." << std::endl;
        return;
    }


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
    const Double_t t_measure_735 = (physics_run_days - total_scan_days) * day;

    std::cout << "\n--- Summary ---" << std::endl;
    // 合計スキャン日数を出力
    std::cout << "Total scan time: " << total_scan_days << " days" << std::endl;
    // 残りの日数を出力
    std::cout << "Remaining time (t_measure_735): " << t_measure_735/day << " days" << std::endl;
    

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

    // --- Check Command-Line Arguments ---
    // Expect 3 arguments: 1. program name, 2. data directory, 3. config file
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <data_directory> <config_file>" << std::endl;
        return 1; // Exit with an error code
    }

    // --- Get Arguments ---
    TString data_dir = argv[1];
    TString config_filename = argv[2];

    // --- Initialize Global Containers ---
    n_kaon_container_scan.resize(conf.n_mom_points, 0.0);
    n_kaon_container_all.resize(conf.n_mom_points, 0.0);

    // --- Run Analysis ---
    analyze(data_dir, config_filename);

    std::cout << "Analysis finished." << std::endl;
    return 0; // Exit successfully
}