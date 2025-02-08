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

    // ____________________________________________________________________________________________
    std::vector<std::vector<Double_t>> load_data(TString path) {
        std::ifstream file(path.Data());
        std::string line;

        std::vector<std::vector<Double_t>> data;

        // Skip the header line
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string value;

            std::vector<Double_t> buf;
            while (std::getline(ss, value, ',')) {
                buf.push_back(std::stod(value));
            }
            data.push_back(buf);
        }
        return data;
    }

    // ____________________________________________________________________________________________
    void TPC_pad_template(TH2Poly *h)
    {
        Double_t X[5];
        Double_t Y[5];
        for (Int_t l=0; l<pad_helper::NumOfLayersTPC; ++l) {
            Double_t pLength = pad_helper::padParameter[l][5];
            Double_t st      = (180.-(360./pad_helper::padParameter[l][3]) *
                            pad_helper::padParameter[l][1]/2.);
            Double_t sTheta  = (-1+st/180.)*TMath::Pi();
            Double_t dTheta  = (360./pad_helper::padParameter[l][3])/180.*TMath::Pi();
            Double_t cRad    = pad_helper::padParameter[l][2];
            Int_t    nPad    = pad_helper::padParameter[l][1];
            for (Int_t j=0; j<nPad; ++j) {
                X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
                X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
                X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
                X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
                X[0] = X[4];
                Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
                Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
                Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
                Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
                Y[0] = Y[4];
                for (Int_t k=0; k<5; ++k) X[k] += -143.;
                h->AddBin(5, X, Y);
            }
        }
        h->SetMaximum(0x1000);
    }


}