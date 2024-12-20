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

    Int_t hitpos_bin_num = 800;
    Double_t hitpos_min = -400.0;
    Double_t hitpos_max =  400.0;
    

    // -- for Forward Proton trigger -----
    std::vector<Int_t> forward_seg_narrow{17, 18, 19, 20, 21};
    std::vector<Int_t> forward_seg_wide{10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    std::vector<Int_t> forward_seg_all{6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};


    // std::unordered_map<Int_t, std::vector<Int_t>> focus_decay_particle = {
    //     { 7202, {2212, -211} }, // eta Lambda
    //     { 7203, {2212, -211} }, // pi0 Lambda
    //     { 7204, {+211, -211} }, // pi+ Sigma-
    //     { 7205, {2212, -211} }, // pi0 Sigma0
    //     { 7206, {-211, +211, 2212} }, // pi- Sigma+
    //     { 7207, {-321, 2212} }, // K p
    //     { 7208, {+211, -211} }  // k0 n
    // };


    Int_t cos_theta_n_bin = 100;

    // -- for momentum binning -----
    Int_t n_mom_points = 600;
    Double_t mom_step_size = 0.5; // MeV/c
    Double_t mom_start = 600.0; // MeV
    Int_t cos_theta_bin_num = 50;

    // -- for yield estimation -----
    Double_t spill_length = 4.24; // second
    Double_t density_LH2 = 0.07085; // g/cm^3
    Double_t Na = 6.02214076 * TMath::Power(10.0, 23.0); // /mol
    Double_t d = 6.898; // cm
    Double_t daq_eff = 0.9;

    // -- for dcs calc -----
    Int_t n_smoothing = 5;

    void beam_initialize() {
        beam_generator = -1;
    }
    

private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
