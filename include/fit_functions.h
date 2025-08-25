#ifndef FIT_FUNCTIONS_H_
#define FIT_FUNCTIONS_H_

#include <Math/SpecFuncMathMore.h>
#include <vector>

namespace fit_functions
{
    static const Double_t hbarc = 197.3269804;
    static const Double_t mass_Lambda  = 1115.683;
    static const Double_t mass_eta     =  547.862;
    static const Double_t mass_kaon_pm =  493.677;
    static const Double_t mass_kaon_0  =  497.611;
    static const Double_t mass_proton  =  938.27208816;
    static const Double_t mass_neutron =  939.56542052;
    static const Double_t mass_pion_pm =  139.57039;
    static const Double_t mass_pion_0  =  134.9768;
    static const Double_t mass_Sigma_p = 1189.37;
    static const Double_t mass_Sigma_0 = 1192.642;
    static const Double_t mass_Sigma_m = 1197.449;

    // static const TString f_amplitude_str = "[0]*ROOT::Math::legendre(0, x) + (2.0*[1] + [2])*ROOT::Math::legendre(1, x) + (3.0*[3] + 2.0*[4])*ROOT::Math::legendre(2, x) + (4.0*[5] + 3.0*[6])*ROOT::Math::legendre(3, x)";
    // static const TString g_amplitude_str = "([0] - [1])*ROOT::Math::assoc_legendre(1, 1, x) + ([2] - 2.0*[3])*ROOT::Math::assoc_legendre(2, 1, x) + ([4] - 3.0*[5])*ROOT::Math::assoc_legendre(3, 1, x)";

    inline Double_t cal_sqrt_s(Double_t mom_kaon)
    {
        Double_t energy_kaon = TMath::Sqrt(mom_kaon * mom_kaon + mass_kaon_pm * mass_kaon_pm);
        return TMath::Sqrt(TMath::Power(energy_kaon + mass_proton, 2.0) - mom_kaon * mom_kaon);
    }

    inline Double_t cal_sqrt_s_err(Double_t mom_kaon, Double_t mom_kaon_err) {
        Double_t energy_kaon = TMath::Sqrt(mom_kaon * mom_kaon + mass_kaon_pm * mass_kaon_pm);
        Double_t numerator = (2.0 * mom_kaon * (energy_kaon + mass_proton)) / energy_kaon - 2.0 * mom_kaon;
        Double_t denominator = 2.0 * TMath::Sqrt(TMath::Power(energy_kaon + mass_proton, 2.0) - mom_kaon * mom_kaon);
        Double_t derivative = numerator / denominator;
        return TMath::Sqrt(derivative * derivative * mom_kaon_err * mom_kaon_err);
    }

    inline Double_t f_cusp(Double_t *x, Double_t *par)
    {
        Double_t E     = x[0];
        Double_t a     = par[1];
        Double_t b     = par[2];
        Double_t r     = par[3];
        Double_t theta = par[4];
        Double_t m1    = mass_Lambda;
        Double_t m2    = mass_eta;

        TComplex A_c = TComplex(a/hbarc, b/hbarc);
        TComplex i_c = TComplex(0, 1);
        TComplex const_term = r * ( TMath::Cos(theta) + i_c*TMath::Sin(theta));
        
        Double_t k_insqrt = ( TMath::Power(E, 2.0) - TMath::Power(m1+m2, 2.0) ) * ( TMath::Power(E, 2.0) - TMath::Power(m1-m2, 2.0) );
        
        TComplex k;
        if (k_insqrt >= 0.) {
            k = TMath::Sqrt(k_insqrt) / (2.0*E);
        } else {
            k = i_c * TMath::Sqrt(-k_insqrt) / (2.0*E);
        }

        TComplex f = TMath::Sqrt(b/hbarc) / (1.0 - i_c * k * A_c);
        TComplex f_tot = f + const_term;

        return par[0]*f_tot.Rho2();
    }

    inline TF1Convolution f_cusp_conv(Double_t sigma, Double_t range_min, Double_t range_max, Int_t NofPointsFFT = 10000)
    {
        auto *gauss = new TF1("gauss", Form("TMath::Gaus(x, 0, %f, true)", sigma));
        auto *cusp  = new TF1("cusp", f_cusp, range_min, range_max, 5);
        auto *cusp_conv = new TF1Convolution(cusp, gauss, range_min-100.0, range_max+100.0, true);
        cusp_conv->SetNofPointsFFT(NofPointsFFT);
        return *cusp_conv;
    }
    
    
    inline Double_t f_cusp_with_range(Double_t *x, Double_t *par)
    {
        Double_t E     = x[0];
        Double_t a     = par[1];
        Double_t b     = par[2];
        Double_t r_re  = par[3];
        Double_t r_im  = par[4];
        Double_t theta = par[5];
        Double_t m1    = mass_Lambda;
        Double_t m2    = mass_eta;

        TComplex A_c = TComplex(a/hbarc, b/hbarc);
        TComplex r_c = TComplex(r_re/hbarc, r_im/hbarc);
        TComplex i_c = TComplex(0, 1.0);
        TComplex phase_pos = TComplex(TMath::Cos(theta), TMath::Sin(theta));
        TComplex phase_neg = TComplex(TMath::Cos(theta), -1.0*TMath::Sin(theta));
        
        Double_t k_insqrt = ( TMath::Power(E, 2.0) - TMath::Power(m1+m2, 2.0) ) * ( TMath::Power(E, 2.0) - TMath::Power(m1-m2, 2.0) );
        
        TComplex k;
        if (k_insqrt >= 0) {
            k = TMath::Sqrt(k_insqrt) / (2.0*E);
        } else {
            k = i_c * TMath::Sqrt(-k_insqrt) / (2.0*E);
        }

        // TComplex numerator = TMath::Sin(theta) + i_c*(b*TMath::Cos(theta) - a*TMath::Sin(theta))*k/hbarc - 1/2*( (a*TMath::Cos(theta) + b*TMath::Sin(theta))*r_im + (b*TMath::Cos(theta) - a*TMath::Sin(theta))*r_re)*k*k/(hbarc*hbarc);
        TComplex numerator = TMath::Sin(theta) + i_c*(phase_neg*A_c).Im()*k - (phase_neg*A_c*r_c).Im()*k*k/2.0;
        TComplex denominator = 1.0 - i_c*k*A_c + A_c*r_c*k*k/2.0;
        TComplex f = phase_pos * numerator / denominator;

        return par[0]*f.Rho2();
    }

    inline TF1Convolution f_cusp_with_range_conv(Double_t sigma, Double_t range_min, Double_t range_max, Int_t NofPointsFFT = 10000)
    {
        auto *gauss = new TF1("gauss", Form("TMath::Gaus(x, 0, %f, true)", sigma));
        auto *cusp  = new TF1("cusp", f_cusp_with_range, range_min, range_max, 6);
        auto *cusp_conv = new TF1Convolution(cusp, gauss, range_min-100.0, range_max+100.0, true);
        cusp_conv->SetNofPointsFFT(NofPointsFFT);
        return *cusp_conv;
    }

    inline TComplex amplitude_f(Double_t *x, Double_t *par)
    {
        Double_t cos_theta = x[0];
        Double_t W         = par[0];
        Double_t m1        = par[1];
        Double_t m2        = par[2];
        // Double_t S1_re    = par[3];
        // Double_t S1_im    = par[4];
        // Double_t P1_re    = par[5];
        // Double_t P1_im    = par[6];
        // ...

        TComplex f = TComplex(par[3], par[4]) * ROOT::Math::legendre(0, cos_theta);
        for (Int_t order = 1; order < 4; order++) {
            f += ( 
                   static_cast<Double_t>(order)        * TComplex(par[4*order+1], par[4*order+2]) 
                + (static_cast<Double_t>(order) + 1.0) * TComplex(par[4*order+3], par[4*order+4]) 
            ) * ROOT::Math::legendre(order, cos_theta);
        }

        Double_t k = TMath::Power( 
            (  TMath::Power(W, 2.0) - TMath::Power(m1+m2, 2.0) ) 
            *( TMath::Power(W, 2.0) - TMath::Power(m1-m2, 2.0) ),
            0.5
        ) / (2.0 * W);

        return f/ k;
    }


    inline TComplex amplitude_g(Double_t *x, Double_t *par)
    {
        Double_t cos_theta = x[0];
        Double_t W         = par[0];
        Double_t m1        = par[1];
        Double_t m2        = par[2];
        // Double_t S1_re    = par[3]; not use
        // Double_t S1_im    = par[4]; not use
        // Double_t P1_re    = par[5];
        // Double_t P1_im    = par[6];
        // ...

        TComplex g = TComplex(0., 0.);
        for (Int_t order = 1; order < 4; order++) {
            g += ( 
                - TComplex(par[4*order+1], par[4*order+2]) 
                + TComplex(par[4*order+3], par[4*order+4]) 
            ) * ROOT::Math::assoc_legendre(order, 1, cos_theta);
        }
        
        Double_t k = TMath::Power( 
            (  TMath::Power(W, 2.0) - TMath::Power(m1+m2, 2.0) ) 
            *( TMath::Power(W, 2.0) - TMath::Power(m1-m2, 2.0) ),
            0.5
        ) / (2.0 * W);

        return g / k;
    }

    
    inline Double_t f_Kp(Double_t *x, Double_t *par)
    {
        Double_t W         = par[0];
        Double_t scale     = par[1];
        // __ isospin 0 ______
        // Double_t S01_re    = par[3];
        // Double_t S01_im    = par[4];
        // Double_t P01_re    = par[5];
        // Double_t P01_im    = par[6];
        // ...
        // Double_t F07_re    = par[15];
        // Double_t F07_im    = par[16];
        // __ isospin 1 ______
        // Double_t S11_re    = par[17];
        // Double_t S11_im    = par[18];
        // ...
        // Double_t F07_re    = par[29];
        // Double_t F07_im    = par[30];

        
        std::vector<Double_t> new_par0;
        new_par0.push_back(W);
        new_par0.push_back(mass_kaon_pm);
        new_par0.push_back(mass_proton);
        for (Int_t i = 3; i < 17; i++) new_par0.push_back(par[i]);

        std::vector<Double_t> new_par1;
        new_par1.push_back(W);
        new_par1.push_back(mass_kaon_pm);
        new_par1.push_back(mass_proton);
        for (Int_t i = 17; i < 31; i++) new_par1.push_back(par[i]);
    
        TComplex f0 = amplitude_f(x, new_par0.data()); // vector から pointer に変換
        TComplex f1 = amplitude_f(x, new_par1.data());
        TComplex g0 = amplitude_g(x, new_par0.data());
        TComplex g1 = amplitude_g(x, new_par1.data());

        TComplex f = ( 1.0/TMath::Power(2.0, 0.5)*f0 + 1.0/TMath::Power(2.0, 0.5)*f1 );
        TComplex g = ( 1.0/TMath::Power(2.0, 0.5)*g0 + 1.0/TMath::Power(2.0, 0.5)*g1 );

        // std::vector<Double_t> new_par;
        // new_par.push_back(W);
        // new_par.push_back(mass_kaon_pm);
        // new_par.push_back(mass_proton);
        // for (Int_t i = 3; i < 17; i++) {
        //     new_par.push_back( 1.0/TMath::Power(2.0, 0.5)*par[i] + 1.0/TMath::Power(2.0, 0.5)*par[i+14] );
        // }
        // TComplex f = amplitude_f(x, new_par.data());
        // TComplex g = amplitude_g(x, new_par.data());

        return scale * TMath::Power(hbarc, 2.0) * 10.0 * (f.Rho2() + g.Rho2()) / 2.0;  // |f|^2 + |g|^2
    }

    inline Double_t f_K0n(Double_t *x, Double_t *par)
    {
        Double_t W         = par[0];
        Double_t scale     = par[2];
        // __ isospin 0 ______
        // Double_t S01_re    = par[3];
        // Double_t S01_im    = par[4];
        // Double_t P01_re    = par[5];
        // Double_t P01_im    = par[6];
        // ...
        // Double_t F07_re    = par[15];
        // Double_t F07_im    = par[16];
        // __ isospin 1 ______
        // Double_t S11_re    = par[17];
        // Double_t S11_im    = par[18];
        // ...
        // Double_t F07_re    = par[29];
        // Double_t F07_im    = par[30];
        
        std::vector<Double_t> new_par0;
        new_par0.push_back(W);
        new_par0.push_back(mass_kaon_0);
        new_par0.push_back(mass_neutron);
        for (Int_t i = 3; i < 17; i++) new_par0.push_back(par[i]);

        std::vector<Double_t> new_par1;
        new_par1.push_back(W);
        new_par1.push_back(mass_kaon_0);
        new_par1.push_back(mass_neutron);
        for (Int_t i = 17; i < 31; i++) new_par1.push_back(par[i]);
    
        TComplex f0 = amplitude_f(x, new_par0.data()); // vector から pointer に変換
        TComplex f1 = amplitude_f(x, new_par1.data());
        TComplex g0 = amplitude_g(x, new_par0.data());
        TComplex g1 = amplitude_g(x, new_par1.data());

        TComplex f = ( - 1.0/TMath::Power(2.0, 0.5)*f0 + 1.0/TMath::Power(2.0, 0.5)*f1 );
        TComplex g = ( - 1.0/TMath::Power(2.0, 0.5)*g0 + 1.0/TMath::Power(2.0, 0.5)*g1 );
        
        return scale * TMath::Power(hbarc, 2.0) * 10.0 * (f.Rho2() + g.Rho2()) / 2.0;  // |f|^2 + |g|^2
    }


    inline Double_t f_pipSigmam(Double_t *x, Double_t *par)
    {
        Double_t W         = par[0];
        Double_t scale     = par[1];
        // __ isospin 0 ______
        // Double_t S01_re    = par[3];
        // Double_t S01_im    = par[4];
        // Double_t P01_re    = par[5];
        // Double_t P01_im    = par[6];
        // ...
        // Double_t F07_re    = par[15];
        // Double_t F07_im    = par[16];
        // __ isospin 1 ______
        // Double_t S11_re    = par[17];
        // Double_t S11_im    = par[18];
        // ...
        // Double_t F07_re    = par[29];
        // Double_t F07_im    = par[30];
        
        std::vector<Double_t> new_par0;
        new_par0.push_back(W);
        new_par0.push_back(mass_pion_pm);
        new_par0.push_back(mass_Sigma_m);
        for (Int_t i = 3; i < 17; i++) new_par0.push_back(par[i]);

        std::vector<Double_t> new_par1;
        new_par1.push_back(W);
        new_par1.push_back(mass_pion_pm);
        new_par1.push_back(mass_Sigma_m);
        for (Int_t i = 17; i < 31; i++) new_par1.push_back(par[i]);
    
        TComplex f0 = amplitude_f(x, new_par0.data()); // vector から pointer に変換
        TComplex f1 = amplitude_f(x, new_par1.data());
        TComplex g0 = amplitude_g(x, new_par0.data());
        TComplex g1 = amplitude_g(x, new_par1.data());

        TComplex f = ( 1.0/TMath::Power(2.0, 0.5)*f0 + 1.0/TMath::Power(3.0, 0.5)*f1 );
        TComplex g = ( 1.0/TMath::Power(2.0, 0.5)*g0 + 1.0/TMath::Power(3.0, 0.5)*g1 );
        
        return scale * TMath::Power(hbarc, 2.0) * 10.0 * (f.Rho2() + g.Rho2());  // |f|^2 + |g|^2
    }

    inline Double_t f_pimSigmap(Double_t *x, Double_t *par)
    {
        Double_t W         = par[0];
        Double_t scale     = par[2];
        // __ isospin 0 ______
        // Double_t S01_re    = par[3];
        // Double_t S01_im    = par[4];
        // Double_t P01_re    = par[5];
        // Double_t P01_im    = par[6];
        // ...
        // Double_t F07_re    = par[15];
        // Double_t F07_im    = par[16];
        // __ isospin 1 ______
        // Double_t S11_re    = par[17];
        // Double_t S11_im    = par[18];
        // ...
        // Double_t F07_re    = par[29];
        // Double_t F07_im    = par[30];
        
        std::vector<Double_t> new_par0;
        new_par0.push_back(W);
        new_par0.push_back(mass_pion_pm);
        new_par0.push_back(mass_Sigma_p);
        for (Int_t i = 3; i < 17; i++) new_par0.push_back(par[i]);

        std::vector<Double_t> new_par1;
        new_par1.push_back(W);
        new_par1.push_back(mass_pion_pm);
        new_par1.push_back(mass_Sigma_p);
        for (Int_t i = 17; i < 31; i++) new_par1.push_back(par[i]);
    
        TComplex f0 = amplitude_f(x, new_par0.data()); // vector から pointer に変換
        TComplex f1 = amplitude_f(x, new_par1.data());
        TComplex g0 = amplitude_g(x, new_par0.data());
        TComplex g1 = amplitude_g(x, new_par1.data());

        TComplex f = ( 1.0/TMath::Power(2.0, 0.5)*f0 - 1.0/TMath::Power(3.0, 0.5)*f1 ) ;
        TComplex g = ( 1.0/TMath::Power(2.0, 0.5)*g0 - 1.0/TMath::Power(3.0, 0.5)*g1 ) ;
        
        return scale * TMath::Power(hbarc, 2.0) * 10.0 * (f.Rho2() + g.Rho2());  // |f|^2 + |g|^2
    }


}


#endif  // FIT_FUNCTIONS_H_