#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TList.h>
#include <TString.h>
#include <TClass.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

// C++の実行可能ファイルとしてコンパイルするため、main関数を定義
int main(int argc, char* argv[]) {
    // --- 0. コマンドライン引数のチェック ---
    if (argc < 2) {
        std::cerr << "エラー: 入力ROOTファイルを引数に指定してください。" << std::endl;
        std::cerr << "使用法: ./combine_hists <input_file.root>" << std::endl;
        return 1; // 異常終了
    }
    // 最初の引数を入力ファイル名として受け取る
    std::string inFileName = argv[1];

    // --- 1. ファイルを開く ---
    TFile *inFile = TFile::Open(inFileName.c_str());
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "エラー: 入力ファイル " << inFileName << " を開けません。" << std::endl;
        return 1;
    }
    std::cout << "入力ファイル: " << inFileName << " を開きました。" << std::endl;

    // --- 出力ファイル名を自動生成 ---
    std::string outFileName = inFileName;
    size_t pos = outFileName.rfind(".root");
    if (pos != std::string::npos) {
        outFileName.replace(pos, 5, "_combined.root");
    } else {
        outFileName += "_combined.root";
    }

    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "エラー: 出力ファイル " << outFileName << " を作成できません。" << std::endl;
        inFile->Close();
        return 1;
    }
    std::cout << "出力ファイル: " << outFileName << " を作成します。" << std::endl;

    // --- 2. 運動量の範囲とステップを定義 ---
    const double p_start = 724.0;
    const double p_end = 770.0;
    const double p_step = 2.0;

    // --- 3. 合算用のヒストグラムを準備 ---
    std::map<int, TH1D*> sum_raw_hists;
    std::map<int, TH1D*> sum_trig_hists;

    // --- 4. ファイル内の全ヒストグラムを一度だけループ ---
    TIter next(inFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D")) {
            continue;
        }

        TH1D *h = (TH1D*)key->ReadObj();
        TString h_name = h->GetName();
        TString h_title = h->GetTitle();

        float current_p, dummy;
        if (sscanf(h_title.Data(), "%f , %f", &current_p, &dummy) != 2) {
            delete h;
            continue;
        }

        if (current_p < p_start || current_p >= p_end) {
            delete h;
            continue;
        }

        int bin_index = static_cast<int>((current_p - p_start) / p_step);
        int bin_p_start = p_start + bin_index * p_step;

        if (h_name.Contains("cos_theta_raw")) {
            if (sum_raw_hists.find(bin_index) == sum_raw_hists.end()) {
                TString new_name = TString::Format("cos_theta_raw_%d_%d", bin_p_start, bin_p_start + (int)p_step);
                sum_raw_hists[bin_index] = (TH1D*)h->Clone(new_name);
                sum_raw_hists[bin_index]->SetDirectory(nullptr); // メモリ管理を outFile に移すため
                sum_raw_hists[bin_index]->Reset();
            }
            sum_raw_hists[bin_index]->Add(h);
        } else if (h_name.Contains("cos_theta_trig")) {
            if (sum_trig_hists.find(bin_index) == sum_trig_hists.end()) {
                TString new_name = TString::Format("cos_theta_trig_%d_%d", bin_p_start, bin_p_start + (int)p_step);
                sum_trig_hists[bin_index] = (TH1D*)h->Clone(new_name);
                sum_trig_hists[bin_index]->SetDirectory(nullptr); // メモリ管理を outFile に移すため
                sum_trig_hists[bin_index]->Reset();
            }
            sum_trig_hists[bin_index]->Add(h);
        }
        
        delete h;
    }

    // --- 5. 結果をファイルに書き込む ---
    outFile->cd();
    std::cout << "\n--- 結果をファイルに書き込みます ---" << std::endl;
    for (auto const& [key, val] : sum_raw_hists) {
        val->Write();
        std::cout << "書き込み完了: " << val->GetName() << std::endl;
    }
    for (auto const& [key, val] : sum_trig_hists) {
        val->Write();
        std::cout << "書き込み完了: " << val->GetName() << std::endl;
    }

    // --- 6. クリーンアップ ---
    inFile->Close();
    outFile->Close();
    delete inFile;
    delete outFile;

    std::cout << "\n処理が完了しました。結果は " << outFileName << " に保存されました。" << std::endl;
    return 0; // 正常終了
}