// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <tuple>

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
#include <TSpline.h>     
#include <TGraphErrors.h>
#include <TKey.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

static std::vector<Double_t> n_react_container;

void analyze(TString path, Int_t focus_pdg_code){
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
    TTreeReader reader("g4hyptpc_light", f);
    TTreeReaderValue<Double_t> mom_kaon_lab(reader, "mom_kaon_lab");
    TTreeReaderValue<Double_t> cos_theta(reader, "cos_theta");
    TTreeReaderValue<Int_t> trig_flag(reader, "trig_flag");
    TTreeReaderValue<Int_t> pdg_code(reader, "decay_particle_code");
    Int_t total_entry = reader.GetEntries();


    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t sla_index = path.Last('/');
    for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    TString output_path = Form("%s/root/yield_%s_%d_for_cusp.root", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "recreate");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    TH1D *h_cos_theta_raw[conf.n_mom_points];
    TH1D *h_cos_theta_trig[conf.n_mom_points];
    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        h_cos_theta_raw[i]  = new TH1D(
            Form("cos_theta_raw%d", i),
            Form("%f , %f", conf.mom_start+i*conf.mom_step_size, conf.mom_start+(i+1)*conf.mom_step_size), 
            conf.cos_theta_bin_num, -1.0, 1.0
        );
        h_cos_theta_trig[i] = new TH1D(
            Form("cos_theta_trig%d", i),
            Form("%f , %f", conf.mom_start+i*conf.mom_step_size, conf.mom_start+(i+1)*conf.mom_step_size),
            conf.cos_theta_bin_num, -1.0, 1.0
        );
    }
    auto *h_mom_dist_raw = new TH1D("mom_dist_raw", "mom_dist", 1000, 500.0, 1000.0);
    auto *h_mom_dist_trig = new TH1D("mom_dist_trig", "mom_dist", 1000, 500.0, 1000.0);

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t n_points = n_react_container.size();

    std::vector<std::vector<std::tuple<Bool_t, Double_t, Double_t>>> container(n_points, std::vector<std::tuple<Bool_t, Double_t, Double_t>>());
    Int_t evnum = 0;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum, total_entry);
        // -- re-make trig_flag and prepare fill data-----
        Bool_t trig_flag_include_branch = false;
        if (*trig_flag != 0 && (focus_pdg_code == 9999 || *pdg_code == focus_pdg_code) ) trig_flag_include_branch = true;

        Int_t index = ana_helper::get_index( *mom_kaon_lab );
        if (index != -1) container[index].emplace_back(trig_flag_include_branch, *mom_kaon_lab, *cos_theta);        
    }

    // -- check statistics -----
    Bool_t exit_flag = false;
    for (Int_t i = 0; i < n_points; i++) if (container[i].size() < n_react_container[i]*conf.daq_eff) {
        std::cout << "not enough data: " <<  conf.mom_start+i*conf.mom_step_size << " - " << conf.mom_start+(i+1)*conf.mom_step_size << ", " << n_react_container[i]*conf.daq_eff - container[i].size() << std::endl;
        exit_flag = true;
    }
    if (exit_flag) {
        fout.Close(); 
        return;
    }

    // +--------------------+
    // | ramdomize and fill |
    // +--------------------+
    for (Int_t i = 0; i < n_points; i++) {
        unsigned int seed = 72*(i+1);
        std::mt19937 gen(seed);
        std::vector<Int_t> indices(container[i].size());
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);
        for (Int_t j = 0; j < static_cast<Int_t>(n_react_container[i]*conf.daq_eff); j++) {
            Bool_t flag;
            Double_t tmp_mom, tmp_cos_theta;
            std::tie(flag, tmp_mom, tmp_cos_theta) = container[i][indices[j]];
            
            Int_t index = ana_helper::get_index(tmp_mom);
            h_cos_theta_raw[index]->Fill(tmp_cos_theta);
            h_mom_dist_raw->Fill(tmp_mom);
            if (flag) {
                h_cos_theta_trig[index]->Fill(tmp_cos_theta);
                h_mom_dist_trig->Fill(tmp_mom);
            }
        }
    }

    // +-----------+
    // | Print PDF |
    // +-----------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 4;
    Int_t cols = 4;
    Int_t max_pads = rows * cols;
    TString pdf_name = Form("%s/img/yield_%s_%d_for_cusp.pdf", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);

    auto *c = new TCanvas("yield", "", 1500, 1200);
    c->Divide(cols, rows);
    c->Print(pdf_name + "["); // start
    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        h_cos_theta_raw[i]->SetLineColor(kBlack);
        h_cos_theta_raw[i]->GetYaxis()->SetRangeUser(0.0, h_cos_theta_raw[i]->GetMaximum()*1.2);
        h_cos_theta_raw[i]->Draw();
        h_cos_theta_trig[i]->SetLineColor(kRed);
        h_cos_theta_trig[i]->Draw("same");   
        nth_pad++;
    }
    c->Print(pdf_name);
    c->Print(pdf_name + "]"); // end
    delete c;

    // +------------+
    // | Write data |
    // +------------+
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Yield: " << h_mom_dist_trig->GetEntries() << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    

    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        h_cos_theta_raw[i]->Write();
        h_cos_theta_trig[i]->Write();
    }
    h_mom_dist_raw->Write();
    h_mom_dist_trig->Write();
    
    fout.Close(); 

    
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    // =======================================================
    // Program Arguments:
    // argv[1]: root file path (required)
    // argv[2]: focus PDG code (optional, default: 9999)
    // =======================================================

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path> [focus_pdg_code]" << std::endl;
        return 1;
    }

    // Load parameters
    TString path = argv[1];
    Int_t focus_pdg_code = (argc > 2) ? std::atoi(argv[2]) : 9999; // Default: 9999

    // Output parameters
    std::cout << "\n====================================" << std::endl;
    std::cout << "Path:           " << path << std::endl;
    std::cout << "Focus PDG Code: " << focus_pdg_code << std::endl;
    std::cout << "====================================\n" << std::endl;


    // +-------------+
    // | load n_kaon |
    // +-------------+
    auto *f_kaon = new TFile(Form("%s/data/n_kaon.root", WORK_DIR.Data()));
    TTreeReader reader_kaon("tree", f_kaon);
    TTreeReaderValue<std::vector<Double_t>> n_kaon_container(reader_kaon, "n_kaon_all");
    reader_kaon.Restart();
    reader_kaon.Next();


    // +-----------------------+
    // | load spline root file |
    // +-----------------------+
    TString reaction1 = "kp";
    TString reaction2 = "k0n";
    TString reaction3 = "pi+sigma-";
    
    TString spline_root_file_path;
    if (path.Contains(reaction1)) {
        spline_root_file_path = Form("%s/data/spline/Kp_spline_for_cusp.root", WORK_DIR.Data());
    } else if (path.Contains(reaction2)) {
        spline_root_file_path = Form("%s/data/spline/K0n_spline_for_cusp.root", WORK_DIR.Data());
    } else if (path.Contains(reaction3)) {
        spline_root_file_path = Form("%s/data/spline/piPlusSigmaMinus_spline_for_cusp.root", WORK_DIR.Data());
    } else {
        std::cerr << "Error: No matching reaction found in path: " << path << std::endl;
        return 1;
    }
    auto *f_spline = new TFile(spline_root_file_path.Data());
    std::cout << "Loaded spline root file: " << spline_root_file_path << std::endl;

    // -- prepare spline func -----
    std::vector<TSpline3*> splines;
    TIter next(f_spline->GetListOfKeys());
    TKey* spline_key;
    while ((spline_key = (TKey*)next())) {
        TObject* obj = spline_key->ReadObj();
        if (obj->InheritsFrom("TSpline3")) {
            TSpline3* spline = (TSpline3*)obj;
            splines.push_back(spline);
            std::cout << "Loaded spline: " << spline->GetName() << std::endl;
        }
    }

    // -- calc n_react -----
    Double_t mom_min = splines[0]->GetXmin();
    Double_t mom_max = splines[0]->GetXmax(); 
    TString f_legendre_str = "";
    for (Int_t order = 0, n_coeff = splines.size(); order < n_coeff; order++) f_legendre_str += Form(" + [%d]*ROOT::Math::legendre(%d,x)", order, order);
    
    for (Double_t i = 0; i < conf.n_mom_points; i++) {
        Double_t mom = conf.mom_start + (i+0.5)*conf.mom_step_size;
        
        // -- calc cross section use spline coeff -----
        TF1 f_legendre("tmp_legendre", f_legendre_str, -1.0, 1.0);
        for (Int_t order = 0, n_coeff = splines.size(); order < n_coeff; order++) f_legendre.SetParameter(order, splines[order]->Eval(mom));
        Double_t cs = 2.0*TMath::Pi()*f_legendre.Integral(-1.0, 1.0);

        // -- calc n react -----
        Double_t n_kaon = (*n_kaon_container)[ana_helper::get_index(mom)] > 500.0 ? (*n_kaon_container)[ana_helper::get_index(mom)] : 0.0;
        Double_t tmp_n_react = cs*TMath::Power(10.0, -27.0) * n_kaon * conf.density_LH2 * conf.Na * conf.d;
        if (mom_min <= mom && mom <= mom_max) n_react_container.push_back(tmp_n_react);
        else n_react_container.push_back(0.0);
    }

    analyze(path, focus_pdg_code);
    return 0;
}