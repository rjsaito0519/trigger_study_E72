// config.h
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    Int_t max_htof_ch = 34;
    Int_t max_bh2_ch  =  8;
    Int_t max_kvc_ch  =  8;

    Int_t beam_generator = 7201;
    Double_t refractive_index_kvc = 1.46;

    Double_t edep_threshold = 0.2; // MeV

    // -- for histogram range -----
    Double_t htof_seg_min = -0.5;
    Double_t htof_seg_max = 33.5;

    Int_t edep_bin_num = 120;
    Double_t edep_min = -0.5;
    Double_t edep_max = 11.5;




    std::unordered_map<Int_t, std::vector<Int_t>> focus_decay_particle = {
        { 7202, {2212, -211} }, // eta Lambda
        { 7203, {2212, -211} }, // pi0 Lambda
        { 7204, {+211, -211} }, // pi+ Sigma-
        { 7205, {2212, -211} }, // pi0 Sigma0
        { 7206, {-211, +211, 2212} }, // pi- Sigma+
        { 7207, {-321, 2212} }, // K p
        { 7208, {+211, -211} }  // k0 n
    };


    Int_t n_bin = 100;
    Int_t n_mom_points = 600;
    Double_t mom_step_size = 0.5; // MeV/c
    Double_t mom_start = 600.0; // MeV

    Double_t density_LH2 = 0.07085; // g/cm^3
    Double_t Na = 6.02214076 * TMath::Power(10.0, 23); // /mol
    Double_t d = 6.898; // cm
    Double_t daq_eff = 0.9;


    // void bac_initialize(Int_t tmp_npe_bin_num = 525) {
    //     adjust_adc_bin_num = 1024;
    //     sumadc_bin_num = 1024;
    //     npe_bin_num = tmp_npe_bin_num;
    //     threshold_fit_range_max = 100.0;
    // }
    

private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
