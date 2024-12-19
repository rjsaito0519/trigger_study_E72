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
#include <TSpline.h>     
#include <TGraphErrors.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "DiffCrossSection.hh"

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
    TString output_path = Form("%s/data/spline/%s_spline.root", WORK_DIR.Data(), save_name.Data());
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
        Int_t nth_order = 0;
        for (const auto &mom : mom_kaons) {
            g_coeff[order]->SetPoint(nth_order, mom, legendre_coeff.at(mom)[order]);
            nth_order++;
        }
        g_coeff[order]->SetMarkerSize(1.);
        g_coeff[order]->SetMarkerStyle(20);

        // -- prepare spline ------
        spline[order] = new TSpline3(Form("A%d_spline", order), g_coeff[order]);
        spline[order]->SetName(Form("A%d_spline", order));
    }
    
    // +------------+
    // | Write data |
    // +------------+
    fout.cd(); 
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
    TString save_name = "Kp_bubble1970";
    mom_kaons = DiffCrossSection::mom_kaons_Kp_bubble1970;
    legendre_coeff = DiffCrossSection::legendre_coeff_Kp_bubble1970;
    analyze(save_name);

    // +-----+
    // | K0n |
    // +-----+
    save_name = "K0n_bubble1970";
    mom_kaons = DiffCrossSection::mom_kaons_K0n_bubble1970;
    legendre_coeff = DiffCrossSection::legendre_coeff_K0n_bubble1970;
    analyze(save_name);
    
    // +-----------+
    // | pi0Sigma0 |
    // +-----------+
    save_name = "pi0Sigma0_bubble1970";
    mom_kaons = DiffCrossSection::mom_kaons_pi0Sigma0_bubble1970;
    legendre_coeff = DiffCrossSection::legendre_coeff_pi0Sigma0_bubble1970;
    analyze(save_name);

    // +-----------+
    // | pi0Lambda |
    // +-----------+
    save_name = "pi0Lambda_bubble1970";
    mom_kaons = DiffCrossSection::mom_kaons_pi0Lambda_bubble1970;
    legendre_coeff = DiffCrossSection::legendre_coeff_pi0Lambda_bubble1970;
    analyze(save_name);

    // +------------------+
    // | piMinusSigmaPlus |
    // +------------------+
    save_name = "piMinusSigmaPlus_bubble1970";
    mom_kaons = DiffCrossSection::mom_kaons_piMinusSigmaPlus_bubble1970;
    legendre_coeff = DiffCrossSection::legendre_coeff_piMinusSigmaPlus_bubble1970;
    analyze(save_name);

    // +------------------+
    // | piPlusSigmaMinus |
    // +------------------+
    save_name = "piPlusSigmaMinus_bubble1970";
    mom_kaons = DiffCrossSection::mom_kaons_piPlusSigmaMinus_bubble1970;
    legendre_coeff = DiffCrossSection::legendre_coeff_piPlusSigmaMinus_bubble1970;
    analyze(save_name);

    // +-----------+
    // | etaLambda |
    // +-----------+
    save_name = "etaLambda_CB";
    mom_kaons = DiffCrossSection::mom_kaons_etaLambda_CB;
    legendre_coeff = DiffCrossSection::legendre_coeff_etaLambda_CB;
    analyze(save_name);

    return 0;
}