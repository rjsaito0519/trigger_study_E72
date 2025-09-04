// Standard C++
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <sys/stat.h>

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

// +------------------------+
// | Main Analysis Function |
// +------------------------+
void analyze(TString path) {
    Config& conf = Config::getInstance();

    // +---------+
    // | Styling |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.1216, 0.4667, 0.7059);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980, 0.0549);
    gROOT->GetColor(kGreen)->SetRGB(44.0 / 256, 160.0 / 256, 44.0 / 256);

    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x");
    gStyle->SetTitleSize(0.06, "y");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // +-----------+
    // | Load File |
    // +-----------+
    TFile* f = new TFile(path.Data());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << path << std::endl;
        return;
    }

    TTreeReader reader("g4hyptpc", f);
    TTreeReaderValue<Int_t> generator(reader, "generator");
    TTreeReaderValue<std::vector<TParticle>> HTOF(reader, "HTOF");
    Int_t total_entry = reader.GetEntries();

    // +--------------------------+
    // | Prepare Output Root File |
    // +--------------------------+
    TString save_name;
    Int_t dot_index = path.Last('.');
    Int_t slash_index = path.Last('/');
    for (Int_t i = slash_index + 1; i < dot_index; i++) save_name += path[i];

    TString output_path = Form("%s/root/htof_%s_hit_position.root", OUTPUT_DIR.Data(), save_name.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "RECREATE");

    // +------------------------------+
    // | Prepare Segment-wise Storage |
    // +------------------------------+
    const Int_t num_of_seg = 34;
    std::vector<std::vector<Double_t>> vx(num_of_seg), vy(num_of_seg);
    std::vector<std::vector<Double_t>> vx_high(num_of_seg), vy_high(num_of_seg);
    std::vector<std::vector<Double_t>> vx_others(num_of_seg), vy_others(num_of_seg);

    // +----------------------+
    // | Event Loop & Filling |
    // +----------------------+
    Int_t evnum = 0;
    reader.Restart();
    while (reader.Next()) { displayProgressBar(++evnum, total_entry);
        if (*generator != conf.beam_generator) {
            for (const auto& item : *HTOF) {
                Int_t seg = item.GetMother(1);
                if (seg < 0 || seg >= num_of_seg) continue;

                Double_t x = item.Vx();
                Double_t y = item.Vy();
                Double_t z = item.Vz();

                // Segment-dependent rotation
                Double_t theta = (seg < 6) ? 0.0 : TMath::Pi() / 4.0 * ((seg - 2) / 4);
                Double_t x_prime = x * TMath::Cos(theta) - z * TMath::Sin(theta);
                Double_t z_prime = x * TMath::Sin(theta) + z * TMath::Cos(theta);
                if (seg >= 6) {
                    x_prime -= 280.0 * ((seg - 2) / 4);
                }

                if (item.GetPdgCode() == 2212) {
                    // Store in segment vector
                    vx[seg].push_back(x_prime);
                    vy[seg].push_back(y);

                    if (item.GetWeight() > 3.0) {
                        vx_high[seg].push_back(x_prime);
                        vy_high[seg].push_back(y);
                    }
                } else {
                    if (item.GetWeight() > 3.0) {
                        vx_others[seg].push_back(x_prime);
                        vy_others[seg].push_back(y);
                    }
                }
            }
        }
    }

    // +----------------+
    // | Write to TTree |
    // +----------------+
    TTree output_tree("tree", "");
    std::vector<Double_t> x, y, x_high, y_high, x_others, y_others;
    Int_t seg;

    output_tree.Branch("seg", &seg, "seg/I");
    output_tree.Branch("x", &x);
    output_tree.Branch("y", &y);
    output_tree.Branch("x_high", &x_high);
    output_tree.Branch("y_high", &y_high);
    output_tree.Branch("x_others", &x_others);
    output_tree.Branch("y_others", &y_others);

    for (Int_t i = 0; i < num_of_seg; i++) {
        seg = i;
        x = vx[i];
        y = vy[i];
        x_high = vx_high[i];
        y_high = vy_high[i];
        x_others = vx_others[i];
        y_others = vy_others[i];
        output_tree.Fill();
    }

    output_tree.Write();
    fout.Close();
    delete f;
}

// +---------------------+
// | Program Entry Point |
// +---------------------+
Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // -- Check CLI arguments --
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input ROOT file>" << std::endl;
        return 1;
    }

    TString path = argv[1];

    // -- Beam configuration --
    if (path.Contains("beam")) {
        conf.beam_initialize();
        std::cout << " -- beam config loaded -- " << std::endl;
    }

    analyze(path);
    return 0;
}
