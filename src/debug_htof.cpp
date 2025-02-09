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

#include <TCanvas.h>
#include <TGraph.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TMath.h>

void draw_htof() {
    // パラメータ設定
    const Double_t l = 687.264 / 2;  // HTOFの中心座標
    const Double_t htof_width = 70;  // HTOFの幅
    const Double_t htof_thick = 10;  // HTOFの厚み
    const Double_t htof_gap = 1;     // HTOFの間隔
    const Int_t num_rows = 8;
    const Int_t num_cols = 4;
    const Double_t pi = TMath::Pi();
    const Double_t angle_step = 45.0;  // 各回転角度

    // キャンバス作成
    TCanvas *c1 = new TCanvas("c1", "HTOF Drawing", 800, 800);

    // まずTGraphで座標軸を設定
    Double_t dummy_x[4] = {-400, -400, 400, 400};
    Double_t dummy_y[4] = {-400, 400, 400, -400};
    TGraph *graph = new TGraph(4, dummy_x, dummy_y);
    graph->SetMarkerStyle(1);
    graph->Draw("AP");  // "A" で座標系を確保

    for (Int_t i = 0; i < num_rows; ++i) {
        Double_t theta = angle_step * i * pi / 180.0;  // ラジアン変換
        Double_t cos_theta = TMath::Cos(theta);
        Double_t sin_theta = TMath::Sin(theta);

        for (Int_t j = 0; j < num_cols; ++j) {
            Double_t x0 = htof_width * (j - 2) + htof_gap * (j - 2);
            Double_t y0 = l; // 初期のy座標

            if (j == 0) {
                // 左端の多角形
                Double_t x[5] = {
                    x0,
                    x0 + htof_width,
                    x0 + htof_width,
                    x0 - htof_thick / TMath::Tan(pi * 3 / 8),
                    x0
                };
                Double_t y[5] = {
                    y0,
                    y0,
                    y0 + htof_thick,
                    y0 + htof_thick,
                    y0
                };

                // 回転を適用
                for (Int_t k = 0; k < 5; ++k) {
                    Double_t x_rot = x[k] * cos_theta - y[k] * sin_theta;
                    Double_t y_rot = x[k] * sin_theta + y[k] * cos_theta;
                    x[k] = x_rot;
                    y[k] = y_rot;
                }

                TPolyLine *polygon = new TPolyLine(5, x, y);
                polygon->SetLineColor(kBlack);
                polygon->SetLineWidth(2);
                polygon->SetFillStyle(0);
                polygon->Draw("LSAME");  // SAME を追加
            }
            else if (j == 3) {
                // 右端の多角形
                Double_t x[5] = {
                    x0,
                    x0 + htof_width,
                    x0 + htof_width + htof_thick / TMath::Tan(pi * 3 / 8),
                    x0,
                    x0
                };
                Double_t y[5] = {
                    y0,
                    y0,
                    y0 + htof_thick,
                    y0 + htof_thick,
                    y0
                };

                // 回転を適用
                for (Int_t k = 0; k < 5; ++k) {
                    Double_t x_rot = x[k] * cos_theta - y[k] * sin_theta;
                    Double_t y_rot = x[k] * sin_theta + y[k] * cos_theta;
                    x[k] = x_rot;
                    y[k] = y_rot;
                }

                TPolyLine *polygon = new TPolyLine(5, x, y);
                polygon->SetLineColor(kBlack);
                polygon->SetLineWidth(2);
                polygon->SetFillStyle(0);
                polygon->Draw("LSAME");
            }
            else {
                // 通常の長方形 (TPolyLine に変更)
                Double_t x1 = x0, y1 = y0;
                Double_t x2 = x0 + htof_width, y2 = y0;
                Double_t x3 = x0 + htof_width, y3 = y0 + htof_thick;
                Double_t x4 = x0, y4 = y0 + htof_thick;

                // 回転を適用
                Double_t x_rot[5], y_rot[5];
                x_rot[0] = x1 * cos_theta - y1 * sin_theta;
                y_rot[0] = x1 * sin_theta + y1 * cos_theta;

                x_rot[1] = x2 * cos_theta - y2 * sin_theta;
                y_rot[1] = x2 * sin_theta + y2 * cos_theta;

                x_rot[2] = x3 * cos_theta - y3 * sin_theta;
                y_rot[2] = x3 * sin_theta + y3 * cos_theta;

                x_rot[3] = x4 * cos_theta - y4 * sin_theta;
                y_rot[3] = x4 * sin_theta + y4 * cos_theta;

                // 最初の点に戻る
                x_rot[4] = x_rot[0];
                y_rot[4] = y_rot[0];

                // TPolyLine を使って描画
                TPolyLine *polygon = new TPolyLine(5, x_rot, y_rot);
                polygon->SetLineColor(kBlack);
                polygon->SetLineWidth(2);
                polygon->SetFillStyle(0);
                polygon->Draw("LSAME");

            }
        }
    }

    // 更新処理
    gPad->Modified();
    gPad->Update();
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

    auto *h = new TH2D("test", "", 100, 0., 150.0, 100, 0., 10.0);

    // +----------------------+
    // | check and fill event |
    // +----------------------+
    Int_t evnum = 0;
    std::vector<Int_t> forward_seg = conf.forward_seg_wide;
    reader.Restart();
    while (reader.Next()){ displayProgressBar(++evnum, total_entry);

        // -- HTOF -----
        std::cout << "----" << std::endl;
        for(const auto& item : (*HTOF)) {
            if (std::binary_search(forward_seg.begin(), forward_seg.end(), item.GetMother(1))) {
                std::cout << item.GetMother(1) << ", " << item.GetWeight() << ", " << item.GetPdgCode() << std::endl;
            }
        }
    }

}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file path> [focus_pdg_code]" << std::endl;
        return 1;
    }
    TString path = argv[1];

    TApplication *theApp = new TApplication("App", &argc, argv);    


    // analyze(path);
    draw_htof();

    theApp->Run();
    
    return 0;

}