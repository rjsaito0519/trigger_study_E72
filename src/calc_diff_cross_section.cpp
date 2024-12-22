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
#include <TMinuit.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"


// -- for create hist -----
static Double_t mom_left = 600.0;
static Double_t mom_right = 899.9;
static Int_t n_cluster = 4; // 4であれば2MeV/cごとにまとめる
static Int_t n_rebin = 4; // binの数を変える。4なら100->25

// -- for calc diff cross sec -----
static Double_t branching_ratio = 1.0;

// -- for legendre fit -----
static Int_t n_coeff = 4;
static TString f_legendre_str = "";
static Double_t chi2_value = 0.;
static std::vector<Double_t> fit_cos_theta, fit_cos_theta_err, fit_diff_cs, fit_diff_cs_err;
static Int_t fit_left_index_offset = 0;
static Int_t fit_right_index_offset = 0;

static std::vector<std::vector<Double_t>> old_legendre;

// +--------+
// | Minuit |
// +--------+
// ______________________________________________________________________________________________________
void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{    
    // -- initialize -----
    f = 0;
    // -- calc. chisqare -----
    auto f_legendre = new TF1("legendre", f_legendre_str, -1.0, 1.0);
    for (Int_t order = 0; order < n_coeff; order++) f_legendre->SetParameter(order, par[order]);
    for (Int_t i = 0, n_fit_points = fit_cos_theta.size(); i < n_fit_points; i++) {
        Double_t range_left  = fit_cos_theta[i] - fit_cos_theta_err[i];
        Double_t range_right = fit_cos_theta[i] + fit_cos_theta_err[i];        
        f += TMath::Power( (fit_diff_cs[i]*(range_right - range_left) - f_legendre->Integral( range_left, range_right )) / (fit_diff_cs_err[i]*(range_right - range_left) ), 2.0 );
    }
    delete f_legendre;
    chi2_value = f;
}

FitResult fit_legendre() 
{
    FitResult fit_result;

    // -- minuit setting ----------
    Int_t n_par = n_coeff;
    TMinuit *min = new TMinuit(n_par);
    min->SetPrintLevel(0);  // 0 simple, 1 verbose

    std::vector<Double_t> step, init_par, min_par, max_par;
    step.resize(n_par, 0.0002);
    init_par.resize(n_par, 0.1);
    min_par.resize(n_par, -100.0);
    max_par.resize(n_par,  100.0);

    for(Int_t i =0; i < n_par; i++) min->DefineParameter(i, Form("A%d",i), init_par[i], step[i], min_par[i], max_par[i]);
    min->SetFCN(chi2);

    // -- execute minuit ----------
    Int_t migrad_stats = min->Migrad();
    Int_t ndf = fit_cos_theta.size() - n_par;

    Double_t par, par_err;
    for(Int_t i =0; i < n_par; i++) {
        min->GetParameter(i, par, par_err);
        fit_result.par.push_back(par);
        fit_result.err.push_back(par_err);        
        std::cout << "A" << i << ": " << par << " +/- " << par_err << std::endl;
    }
    fit_result.migrad_stats = migrad_stats;
    fit_result.chi_square = chi2_value;
    fit_result.ndf = ndf;
    std::cout << "Status of Migrad: " << migrad_stats << std::endl;
    std::cout << "chi-square: " << chi2_value << std::endl;
    std::cout << "ndf: " << ndf << std::endl;

    // clean
    delete min;

    return fit_result;
}


// +-------------+
// | Main stream |
// +-------------+
// ______________________________________________________________________________________________________
void analyze(TString path_yield, TString path_acceptance){
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
    auto *f_yield = new TFile(path_yield.Data());
    auto *f_acceptance = new TFile(path_acceptance.Data());
    auto *f_kaon = new TFile(Form("%s/data/n_kaon.root", WORK_DIR.Data()));
    TTreeReader reader_kaon("tree", f_kaon);
    TTreeReaderValue<std::vector<Double_t>> n_kaon_container(reader_kaon, "n_kaon_all");
    reader_kaon.Restart();
    reader_kaon.Next();

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path_yield.Last('.');
    Int_t sla_index = path_yield.Last('/');
    for (Int_t i = sla_index+7; i < dot_index; i++) save_name += path_yield[i];
    TString output_path = Form("%s/root/dcs_%s.root", OUTPUT_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");


    // +-------------------+
    // | prepare histogram |
    // +-------------------+ 
    TH1D *h_cos_theta_raw[conf.n_mom_points];
    TH1D *h_cos_theta_trig[conf.n_mom_points];
    TH1D *h_acceptance[conf.n_mom_points];
    for (Int_t i = 0; i < conf.n_mom_points; i++) {
        h_cos_theta_raw[i]  = (TH1D*)f_yield->Get(Form("cos_theta_raw%d", i));
        h_cos_theta_trig[i] = (TH1D*)f_yield->Get(Form("cos_theta_trig%d", i));
        h_acceptance[i] = (TH1D*)f_acceptance->Get(Form("acceptance%d", i));
        h_acceptance[i]->Smooth(conf.n_smoothing);
    }

    // +-----------------------------------------+
    // | clusterize and calculate diff cross sec |
    // +-----------------------------------------+ 
    // -- prepare container ------
    Int_t n_graph = static_cast<Int_t>( (ana_helper::get_index(mom_right) - ana_helper::get_index(mom_left) + 1) / n_cluster );
    Int_t push_back_index = 0;
    std::vector<std::vector<Double_t>> cos_theta_container(n_graph);
    std::vector<std::vector<Double_t>> cos_theta_err_container(n_graph);
    std::vector<std::vector<Double_t>> diff_cs_container(n_graph);
    std::vector<std::vector<Double_t>> diff_cs_err_container(n_graph);
    std::vector<std::pair<Double_t, Double_t>> mom_container{};
    // -- for consistency check ---
    std::vector<std::vector<Double_t>> calc_n_react_container(n_graph);
    std::vector<std::vector<Double_t>> true_n_react_container(n_graph);
    
    // -- clusterize -----
    for (Int_t mom_index = ana_helper::get_index(mom_left); mom_index+n_cluster <= ana_helper::get_index(mom_right); mom_index+=n_cluster) {
        std::cout << "-----------\n" << conf.mom_start + mom_index*conf.mom_step_size << " - " << conf.mom_start + (mom_index+n_cluster)*conf.mom_step_size << std::endl;

        // -- calc. mom center -----
        Double_t mom_center = 0.0;
        for (Int_t mom_i_offset = 0; mom_i_offset < n_cluster; mom_i_offset++) mom_center += (conf.mom_start + (mom_index+mom_i_offset+0.5)*conf.mom_step_size) / static_cast<Double_t>(n_cluster);
        std::cout << mom_center << ", " << (conf.mom_step_size * n_cluster)/2.0 << std::endl;
        mom_container.emplace_back(mom_center, conf.mom_step_size * n_cluster/2.0);

        // -- rebin histogram and calc. diff cross sec ------
        for (Int_t bin_index = 1, n_bin = h_cos_theta_raw[mom_index]->GetNbinsX(); bin_index <= n_bin; bin_index+=n_rebin) {
            Double_t n_stat = 0.0;
            Double_t n_kaon = 0.0;
            Double_t diff_cs_val = 0.0;
            Double_t cos_theta_center = 0.0;
            Double_t delta_cos_theta = 0.0;
            Double_t n_react = 0.0;
            for (Int_t mom_i_offset = 0; mom_i_offset < n_cluster; mom_i_offset++) {
                n_kaon += (*n_kaon_container)[mom_index+mom_i_offset];
                for (Int_t bin_i_offset = 0; bin_i_offset < n_rebin; bin_i_offset++) {
                    Double_t n_detect_content = h_cos_theta_trig[mom_index+mom_i_offset]->GetBinContent(bin_index+bin_i_offset); 
                    n_stat += n_detect_content;
                    diff_cs_val += n_detect_content / ( branching_ratio * h_acceptance[mom_index+mom_i_offset]->GetBinContent(bin_index+bin_i_offset) );
                    n_react += h_cos_theta_raw[mom_index+mom_i_offset]->GetBinContent(bin_index+bin_i_offset);
                    if (mom_i_offset == 0) {
                        cos_theta_center += h_cos_theta_trig[mom_index+mom_i_offset]->GetBinCenter(bin_index+bin_i_offset) / static_cast<Double_t>(n_rebin);
                        delta_cos_theta  += h_cos_theta_trig[mom_index+mom_i_offset]->GetBinWidth(bin_index+bin_i_offset);
                    }
                }
            }
            // -- error management ---
            if (n_stat == 0 || std::isinf(diff_cs_val)) {
                diff_cs_val = 0.0;
                n_stat = 1.0;
            }
            calc_n_react_container[push_back_index].push_back(diff_cs_val); // at this point, diff_cs_val is corresponding to n_react 
            true_n_react_container[push_back_index].push_back(n_react);

            // unit is mb
            diff_cs_val /= n_kaon * conf.density_LH2 * conf.Na * conf.d * conf.daq_eff * 2*TMath::Pi()*delta_cos_theta * TMath::Power(10.0, -27);
            Double_t diff_cs_err = diff_cs_val / TMath::Sqrt(n_stat);

            diff_cs_container[push_back_index].push_back(diff_cs_val);
            diff_cs_err_container[push_back_index].push_back(diff_cs_err);
            cos_theta_container[push_back_index].push_back(std::round(cos_theta_center * 1000.0) / 1000.0);            
            cos_theta_err_container[push_back_index].push_back(delta_cos_theta/2.0);
        }
        push_back_index++;
    }
    std::cout << "-----------" << std::endl;


    // +-----------------------+
    // | fitting and print PDF |
    // +-----------------------+
    // -- prepare branch -----
    Double_t tmp_mom, tmp_mom_err, chi_square;
    Int_t ndf, migrad_stats;
    std::vector<Double_t> coeff, coeff_err;

    TTree output_tree("tree", "");
    output_tree.Branch("mom", &tmp_mom, "mom/D");
    output_tree.Branch("mom_err", &tmp_mom_err, "mom_err/D");
    output_tree.Branch("cos_theta", &fit_cos_theta);
    output_tree.Branch("cos_theta_err", &fit_cos_theta_err);
    output_tree.Branch("diff_cs", &fit_diff_cs);
    output_tree.Branch("diff_cs_err", &fit_diff_cs_err);

    output_tree.Branch("chi_square", &chi_square, "chi_square/D");
    output_tree.Branch("ndf", &ndf, "ndf/I");
    output_tree.Branch("migrad_stats", &migrad_stats, "migrad_stats/I");
    output_tree.Branch("coeff", &coeff);
    output_tree.Branch("coeff_err", &coeff_err);

    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2;
    Int_t cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_name = Form("%s/img/dcs_%s.pdf", OUTPUT_DIR.Data(), save_name.Data());

    // -- prepare coeff graph -----
    TGraphErrors *g_coeff[n_coeff];
    TGraphErrors *g_coeff_old[n_coeff];
    for (Int_t order = 0; order < n_coeff; order++) {
        g_coeff[order] = new TGraphErrors();
        g_coeff[order]->SetName(Form("A%d", order));
        g_coeff[order]->GetXaxis()->SetTitle("p^{Lab}_{K}");
        g_coeff[order]->GetYaxis()->SetTitle(Form("A%d", order));
        g_coeff[order]->SetLineWidth(2);

        g_coeff_old[order] = new TGraphErrors();
        g_coeff_old[order]->SetName(Form("A%d_old", order));
        g_coeff_old[order]->GetXaxis()->SetTitle("p^{Lab}_{K}");
        g_coeff_old[order]->GetYaxis()->SetTitle(Form("A%d", order));
        g_coeff_old[order]->SetLineColor(kRed);
        g_coeff_old[order]->SetLineWidth(2);

        Int_t nth_point = 0;
        for (const auto &row : old_legendre) {
            g_coeff_old[order]->SetPoint(nth_point, row[0], row[2*order+2]);
            g_coeff_old[order]->SetPointError(nth_point, row[1], row[2*order+3]);
            nth_point++;
        } 
    }

    TGraphErrors *g_cs = new TGraphErrors();
    g_cs->SetName("cross_sec");
    g_cs->GetXaxis()->SetTitle("p^{Lab}_{K}");
    g_cs->GetYaxis()->SetTitle("cross secton [mb]");
    g_cs->SetLineWidth(2); 

    TGraphErrors *g_cs_old = new TGraphErrors();
    g_cs_old->SetName("cross_sec_old");
    g_cs_old->GetXaxis()->SetTitle("p^{Lab}_{K}");
    g_cs_old->GetYaxis()->SetTitle("cross secton [mb]");
    g_cs_old->SetLineColor(kRed);
    g_cs_old->SetLineWidth(2);
    {   
        Int_t nth_point = 0;
        for (const auto &row : old_legendre) {
            auto f_legendre = new TF1("tmp_legendre", f_legendre_str, -1, 1);
            for (Int_t order = 0; order < n_coeff; order++) f_legendre->SetParameter(order, row[2*order+2]);
            Double_t cs = 2.0*TMath::Pi()*f_legendre->Integral(-1.0, 1.0);
            g_cs_old->SetPoint(nth_point, row[0], cs);
            g_cs_old->SetPointError(nth_point, row[1], 0.0);
            nth_point++;
        }

    }

    auto *c = new TCanvas("", "", 1500, 1200);
    c->Divide(cols, rows);
    c->Print(pdf_name + "["); // start
    for (Int_t i = 0; i < n_graph; i++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }

        // -- for consistency check ---
        c->cd(nth_pad);
        Int_t n_points = cos_theta_container[i].size();
        auto *h_true_n_react = new TH1D(Form("h_true_n_react%d", i), "", n_points, -1.0, 1.0);
        for (Int_t j = 0; j < n_points; j++) h_true_n_react->SetBinContent(j+1, true_n_react_container[i][j]);
        h_true_n_react->SetLineColor(kRed);
        h_true_n_react->Draw();

        auto *h_calc_n_react = new TH1D(Form("h_calc_n_react%d", i), "N_{react} #times #epsilon_{daq}", n_points, -1, 1);
        for (Int_t j = 0; j < n_points; j++) h_calc_n_react->SetBinContent(j+1, calc_n_react_container[i][j]);
        h_calc_n_react->SetLineColor(kBlue);
        h_calc_n_react->Draw("same");

        auto *legend = new TLegend(0.4, 0.7, 0.6, 0.9);
        legend->AddEntry(h_calc_n_react, "calc_n_react", "l");
        legend->AddEntry(h_true_n_react, "true_n_react", "l");
        legend->SetBorderSize(0);
        legend->Draw("same");

        // -- diff cross sec -----
        c->cd(++nth_pad);
        auto *g_diff_cs = new TGraphErrors(
            n_points - fit_left_index_offset - fit_right_index_offset,
            &cos_theta_container[i][fit_left_index_offset],
            &diff_cs_container[i][fit_left_index_offset],
            &cos_theta_err_container[i][fit_left_index_offset],
            &diff_cs_err_container[i][fit_left_index_offset]
        );
        g_diff_cs->SetTitle(Form("mom %f;cos#theta;", mom_container[i].first));
        g_diff_cs->Draw("AP");


        // -- legendre fit -----
        // -- initialize -----
        fit_cos_theta.clear(); fit_cos_theta_err.clear(); fit_diff_cs.clear(); fit_diff_cs_err.clear();
        fit_cos_theta.assign(cos_theta_container[i].begin()+fit_left_index_offset, cos_theta_container[i].end()-fit_right_index_offset);
        fit_cos_theta_err.assign(cos_theta_err_container[i].begin()+fit_left_index_offset, cos_theta_err_container[i].end()-fit_right_index_offset);
        fit_diff_cs.assign(diff_cs_container[i].begin()+fit_left_index_offset, diff_cs_container[i].end()-fit_right_index_offset);
        fit_diff_cs_err.assign(diff_cs_err_container[i].begin()+fit_left_index_offset, diff_cs_err_container[i].end()-fit_right_index_offset);

        // -- fitting -----
        FitResult fit_result = fit_legendre();
        // -- fill branch -----
        tmp_mom = mom_container[i].first; tmp_mom_err = mom_container[i].second;
        coeff.clear(); coeff_err.clear();
        for (Int_t order = 0; order < n_coeff; order++) {
            coeff.push_back(fit_result.par[order]);
            coeff_err.push_back(fit_result.err[order]);
            g_coeff[order]->SetPoint(i, tmp_mom, fit_result.par[order]);
            g_coeff[order]->SetPointError(i, tmp_mom_err, fit_result.err[order]);
        }
        migrad_stats = fit_result.migrad_stats;
        ndf = fit_result.ndf;
        chi_square = fit_result.chi_square;
        output_tree.Fill();

        // -- draw -----
        auto f_legendre = new TF1(Form("legendre%d", i), f_legendre_str, -1.0, 1.0);
        for (Int_t order = 0; order < n_coeff; order++) f_legendre->SetParameter(order, fit_result.par[order]);
        f_legendre->SetLineColor(kOrange);
        f_legendre->Draw("same");

        // -- calc total cross section -----
        Double_t cs = 2.0*TMath::Pi()*f_legendre->Integral(-1.0, 1.0);
        g_cs->SetPoint(i, tmp_mom, cs);
        g_cs->SetPointError(i, tmp_mom_err, 0.0);

        nth_pad++;
    }
    c->Print(pdf_name);
    c->Print(pdf_name + "]"); // end
    delete c;

    // +------------+
    // | Write data |
    // +------------+
    for (Int_t order = 0; order < n_coeff; order++) {
        g_coeff[order]->Write();
        g_coeff_old[order]->Write();
    }
    g_cs->Write();
    g_cs_old->Write();
    output_tree.Write();
    fout.Close(); 

}

#include <iostream>
#include <TString.h>

int main(int argc, char** argv) {
    // Check the number of arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <yield_root_file_path> <acceptance_root_file_path>" << std::endl;
        return 1;
    }

    // Retrieve file paths from command-line arguments
    TString path_yield = argv[1];
    TString path_acceptance = argv[2];

    // Print the file paths for verification
    std::cout << "====================================\n" << std::endl;
    std::cout << "Yield ROOT file path: " << path_yield.Data() << std::endl;
    std::cout << "Acceptance ROOT file path: " << path_acceptance.Data() << std::endl;
    std::cout << "====================================\n" << std::endl;


    // +-----------------------+
    // | load spline root file |
    // +-----------------------+
    TString reaction1 = "eta_lambda";
    TString reaction2 = "kp";
    TString reaction3 = "k0n";
    TString reaction4 = "pi+sigma-";
    TString reaction5 = "pi-sigma+";
    TString reaction6 = "pi0sigma0";
    TString reaction7 = "pi0lambda";
    
    TString spline_root_file_path;
    if (path_yield.Contains(reaction1)) {
        // -- for create hist -----
        mom_left = 724.0;
        mom_right = 770.0;
        n_cluster = 4;
        branching_ratio = 0.641;
        // -- for fitting -----
        n_rebin = 2;
        fit_left_index_offset  = 0;
        fit_right_index_offset = 0;
        n_coeff = 4;
        f_legendre_str = "";
        for (Int_t order = 0; order < n_coeff; order++) f_legendre_str += Form(" + [%d]*ROOT::Math::legendre(%d,x)", order, order);
        old_legendre = ana_helper::load_data("data/legendre/etaLambda_crystal_ball_legendre.csv");

    } else if (path_yield.Contains(reaction2)) {
        // -- for create hist -----
        mom_left = 659.0;
        mom_right = 786.0;
        n_cluster = 4;
        branching_ratio = 1.0;
        n_rebin = 2;
        fit_left_index_offset  = 2;
        fit_right_index_offset = 3;
        n_coeff = 6;
        f_legendre_str = "";
        for (Int_t order = 0; order < n_coeff; order++) f_legendre_str += Form(" + [%d]*ROOT::Math::legendre(%d,x)", order, order);
        old_legendre = ana_helper::load_data("data/legendre/Kp_bubble_chamber1970_legendre.csv");

    } else if (path_yield.Contains(reaction3)) {
        // -- for create hist -----
        // mom_left = 659.0;
        mom_left = 659.0+8.0;
        // mom_right = 787.0;
        mom_right = 787.0-8.0;
        n_cluster = 4;
        branching_ratio = 0.692;
        // -- for fitting -----
        n_rebin = 2;
        fit_left_index_offset  = 0;
        fit_right_index_offset = 0;
        n_coeff = 6;
        f_legendre_str = "";
        for (Int_t order = 0; order < n_coeff; order++) f_legendre_str += Form(" + [%d]*ROOT::Math::legendre(%d,x)", order, order);
        old_legendre = ana_helper::load_data("data/legendre/K0n_bubble_chamber1970_legendre.csv");

    } else if (path_yield.Contains(reaction4)) {
        // -- for create hist -----
        // mom_left = 659.0;
        mom_left = 659.0+4.0;
        // mom_right = 787.0;
        mom_right = 787.0-12.0;
        n_cluster = 8;
        branching_ratio = 0.99848;
        // -- for fitting -----
        n_rebin = 2;
        fit_left_index_offset  = 0;
        fit_right_index_offset = 0;
        n_coeff = 5;
        f_legendre_str = "";
        for (Int_t order = 0; order < n_coeff; order++) f_legendre_str += Form(" + [%d]*ROOT::Math::legendre(%d,x)", order, order);
        old_legendre = ana_helper::load_data("data/legendre/pi+Sigma-_bubble_chamber1970_legendre.csv");

    } else if (path_yield.Contains(reaction5)) {
        // spline_root_file_path_yield = Form("%s/data/spline/piMinusSigmaPlus_bubble1970_spline.root", WORK_DIR.Data());
        return 1;
    } else if (path_yield.Contains(reaction6)) {
        // spline_root_file_path_yield = Form("%s/data/spline/pi0Sigma0_bubble1970_spline.root", WORK_DIR.Data());
        return 1;
    } else if (path_yield.Contains(reaction7)) {
        // spline_root_file_path = Form("%s/data/spline/pi0Lambda_bubble1970_spline.root", WORK_DIR.Data());
        return 1;
    }else {
        std::cerr << "Error: No matching reaction found in path: " << path_yield << std::endl;
        return 1;
    }

    analyze(path_yield, path_acceptance);
    return 0;
}