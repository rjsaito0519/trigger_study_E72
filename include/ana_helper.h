#ifndef ANA_HELPER_
#define ANA_HELPER_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>

// ROOTヘッダー
#include <Rtypes.h>
#include <TCanvas.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
#include <TLine.h>
#include <TBox.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TH2Poly.h>

#include "config.h"
#include "pad_helper.h"

// -- fitting result container -----
struct FitResult {
    std::vector<Double_t> par;
    std::vector<Double_t> err;
    Double_t chi_square;
    Int_t ndf;
    Int_t migrad_stats;
    std::vector<Double_t> additional; // 何か追加で値を返したいときのコンテナ
};

// raw と trig をペアで管理する構造体
struct HistPair {
    TH1D* raw;  // raw histogram
    TH1D* trig; // trigged histogram

    HistPair(const TString& object_name, const TString& title, Int_t bins, Double_t range_min, Double_t range_max) {
        raw  = new TH1D(object_name+"_raw",  title, bins, range_min, range_max);
        trig = new TH1D(object_name+"_trig", title, bins, range_min, range_max);
    }

    // // デストラクタでメモリを解放 (これをするとグラフが消えてしまう。ちゃんと動作はするので、グラフ化しないんだったらOK)
    // ~HistPair() {
    //     delete raw;
    //     delete trig;
    // }
};

namespace ana_helper {
    TCanvas* add_tab(TGTab *tab, const char* tabName);
    Int_t get_index(Double_t mom_value);
    std::vector<std::vector<Double_t>> load_data(TString path);
}

#endif  // ANA_HELPER_