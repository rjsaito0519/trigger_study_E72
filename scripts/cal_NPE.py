import numpy as np
from scipy.integrate import quad
from scipy.stats import poisson
import scipy.stats as stats
def normal_pdf(x, mu=0, sigma=1):
    return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

alpha = 7.2973525693*10**-3
m_e = 0.51099895000
m_k = 493.677
m_pi = 139.57039

def frank_tamm(n, beta, L): # L: cm
    return 2*np.pi*alpha * (1 - 1/(n**2 * beta**2)) * L * (1/320 - 1/900)*10**7

def beta(mom, mass):
    return mom / np.sqrt(mass**2 + mom**2)

def cal_poisson(lamb, thre):
    prob = 1.0
    for i in range(thre):
        prob -= poisson.pmf(i, lamb)
    return prob

def cal_gauss(mu, sigma, thre):
    below_thre = stats.norm.cdf(thre, loc=mu, scale=sigma)
    prob = 1.0 - below_thre
    return prob

beta_kekar = beta(2010, m_e)
beta_e72_kaon = beta(735, m_k)
beta_e72_kaon_low = beta(660, m_k)
beta_e72_kaon_low = beta(600, m_k)
beta_e72_kaon_high = beta(790, m_k)
beta_e72_pion = beta(735, m_pi)
beta_e72_pion_low = beta(660, m_pi)
beta_e72_pion_high = beta(790, m_pi)

BAC_npe_th_kekar = 20
BAC_npe_th_jparc = 15
KVC_npe_th_kekar = 20
KVC_npe_th_jparc = 10

n_BAC = 1.115
n_KVC = 1.46

# +------------------------+
# | E72 BAC 3-layer (0, 0) |
# +------------------------+
factor = frank_tamm(n_KVC, beta_e72_kaon_low, 2.0)/frank_tamm(n_KVC, beta_e72_kaon, 2.0)
print(factor)
mu    = 84.987344 * factor
sigma = 15.667632 * factor
eff_poisson = cal_poisson(mu, BAC_npe_th_jparc)
eff_gauss   = cal_gauss(mu, sigma, BAC_npe_th_jparc)
print("\n-- E72 BAC 3-layer (0, 0) -----")
print(mu)
print(f"poisson : {eff_poisson*100.0}")
print(f"gauss   : {eff_gauss*100.0}")
