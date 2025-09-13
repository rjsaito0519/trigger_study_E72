// C++
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sys/stat.h>
#include <vector>

// ROOT
#include <Math/IFunction.h>
#include <Math/IntegratorOptions.h>
#include <Math/SpecFuncMathMore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TComplex.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSpline.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>


// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "DiffCrossSection.hh"
#include "fit_functions.h"

// namespace alias
namespace FF = fit_functions;


static std::vector<Double_t> mom_kaons;
static std::unordered_map<Double_t, std::vector<Double_t>> legendre_coeff;

void analyze(TString save_name){
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

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = Form("%s/data/spline/%s_spline_for_cusp.root", WORK_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");

    Int_t n_coeff = legendre_coeff.at(mom_kaons[0]).size();

    // -- prepare graph and spline -----
    TGraphErrors *g_coeff[n_coeff];
    TSpline3 *spline[n_coeff];
    for (Int_t order = 0; order < n_coeff; order++) {
        g_coeff[order] = new TGraphErrors();
        g_coeff[order]->SetName(Form("A%d", order));
        g_coeff[order]->GetXaxis()->SetTitle("p^{Lab}_{K}");
        g_coeff[order]->GetYaxis()->SetTitle(Form("A%d", order));
        g_coeff[order]->SetLineWidth(2); 

        // -- fill -----
        if (order == 0) {
            Double_t x_min = *std::min_element(mom_kaons.begin(), mom_kaons.end());
            Double_t x_max = *std::max_element(mom_kaons.begin(), mom_kaons.end());
            
            Double_t range_min = FF::cal_sqrt_s(x_min);
            Double_t range_max = FF::cal_sqrt_s(x_max);
            
            auto cusp_conv = FF::f_cusp_with_range_conv(5.0, range_min, range_max);
            auto *cusp_func = new TF1("cusp_func", cusp_conv, range_min, range_max, 6);
            cusp_func->SetParameter(1,  1.61422);   // a
            cusp_func->SetParameter(2,  0.0601015); // b
            cusp_func->SetParameter(3, -2.58405);   // r_re
            cusp_func->SetParameter(4, -0.0814372); // r_im

            if (save_name.Contains("Kp")) {
                cusp_func->SetParameter(0, 18.6321);  // amp
                cusp_func->SetParameter(5,  6.04046); // theta
            } else if (save_name.Contains("K0n")) {
                cusp_func->SetParameter(0, 7.65647); // amp
                cusp_func->SetParameter(5, 2.98705); // theta
            } else if (save_name.Contains("piPlusSigmaMinus")) {
                cusp_func->SetParameter(0, 13.8160);   // amp
                cusp_func->SetParameter(5,  0.156233); // theta
            }
            // 等間隔に分割する点
            Double_t num_points = 100000.0;
            Double_t step = (x_max - x_min) / (num_points - 1.0);

            for (Int_t i = 0; i < num_points; i++) {
                Double_t mom = x_min + i*step;
                Double_t sqrt_s = FF::cal_sqrt_s(mom);
                g_coeff[order]->SetPoint(i, mom, cusp_func->Eval(sqrt_s));
            }

        } else {
            Int_t i_point = 0;
            for (const auto &mom : mom_kaons) {
                g_coeff[order]->SetPoint(i_point, mom, legendre_coeff.at(mom)[order]);
                i_point++;
            }
        }
        g_coeff[order]->SetMarkerSize(1.);
        g_coeff[order]->SetMarkerStyle(20);

        // -- prepare spline ------
        spline[order] = new TSpline3(Form("A%d_spline", order), g_coeff[order]);
        spline[order]->SetName(Form("A%d_spline", order));
    }
    

    // +------------------------------+
    // | For python prepare tree data |
    // +------------------------------+
    TTree* tree = new TTree("tree", "");

    std::vector<std::vector<Double_t>> x_values(n_coeff);
    std::vector<std::vector<Double_t>> y_values(n_coeff);
    std::vector<std::vector<Double_t>> mom_values(n_coeff);
    std::vector<std::vector<Double_t>> coeff_values(n_coeff);
    
    for (Int_t order = 0; order < n_coeff; order++) {
        Double_t x_min = spline[order]->GetXmin();
        Double_t x_max = spline[order]->GetXmax();
        
        // 等間隔に分割する点
        Double_t num_points = 10000.0;
        Double_t step = (x_max - x_min) / (num_points - 1.0);

        // x, y の値を計算して保存
        std::vector<Double_t> x_vec, y_vec;
        for (Double_t i = 0.0; i < num_points; ++i) {
            Double_t x = x_min + i * step;
            Double_t y = spline[order]->Eval(x);
            x_vec.push_back(x);
            y_vec.push_back(y);
        }

        // ブランチ用の vector に追加
        x_values[order] = x_vec;
        y_values[order] = y_vec;

        tree->Branch(Form("A%d_spline_x", order), &x_values[order]);
        tree->Branch(Form("A%d_spline_y", order), &y_values[order]);

        // -- for coeff data itself -----
        for (const auto &mom : mom_kaons) {
            mom_values[order].push_back(mom);
            coeff_values[order].push_back(legendre_coeff.at(mom)[order]);
        }
        tree->Branch(Form("A%d_coeff_x", order), &mom_values[order]);
        tree->Branch(Form("A%d_coeff_y", order), &coeff_values[order]);
    }
    tree->Fill();

    // +------------+
    // | Write data |
    // +------------+
    fout.cd(); 
    tree->Write();
    for (Int_t order = 0; order < n_coeff; order++) {
        spline[order]->Write();
        g_coeff[order]->Write();
    }
    fout.Close(); 

    for (Int_t order = 0; order < n_coeff; order++) {
        delete spline[order];
        delete g_coeff[order];
    }
}

Int_t main(int argc, char** argv) {

    // +----+
    // | Kp |
    // +----+
    TString save_name = "Kp";
    mom_kaons = DiffCrossSection::mom_kaons_Kp_for_cusp;
    legendre_coeff = DiffCrossSection::legendre_coeff_Kp_for_cusp;
    analyze(save_name);

    // +-----+
    // | K0n |
    // +-----+
    save_name = "K0n";
    mom_kaons = DiffCrossSection::mom_kaons_K0n_for_cusp;
    legendre_coeff = DiffCrossSection::legendre_coeff_K0n_for_cusp;
    analyze(save_name);
    
    // +------------------+
    // | piPlusSigmaMinus |
    // +------------------+
    save_name = "piPlusSigmaMinus";
    mom_kaons = DiffCrossSection::mom_kaons_piPlusSigmaMinus_for_cusp;
    legendre_coeff = DiffCrossSection::legendre_coeff_piPlusSigmaMinus_for_cusp;
    analyze(save_name);

    return 0;
}