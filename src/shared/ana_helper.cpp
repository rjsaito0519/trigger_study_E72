#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    TCanvas* add_tab(TGTab *tab, const char* tabName) {
        // タブを作成し、キャンバスを埋め込む
        TGCompositeFrame *tf = tab->AddTab(tabName);
        TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas(tabName, tf, 1000, 800);
       tf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        return embeddedCanvas->GetCanvas();
    }

    // ____________________________________________________________________________________________
    Int_t get_index(Double_t mom_value) {
        Config& conf = Config::getInstance();

        // 指定範囲外のチェック
        if (mom_value < conf.mom_start || mom_value >= conf.mom_start + conf.n_mom_points * conf.mom_step_size) {
            return -1; // 範囲外の場合 -1 を返す
        }
        
        // インデックス計算
        return static_cast<Int_t>((mom_value - conf.mom_start) / conf.mom_step_size);
    }

}