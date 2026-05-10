import math

def gut_normalized_alphas_at_MZ(alpha_em, sin2thetaW, alpha3):
    """
    Convert (alpha_em, sin^2 theta_W, alpha3) at MZ into GUT-normalized (alpha1, alpha2, alpha3).
    alpha1 here uses SU(5) normalization: alpha1 = (5/3) * alpha_Y.
    """
    alpha2 = alpha_em / sin2thetaW              # alpha2 = g2^2 / (4π)
    alphaY  = alpha_em / (1.0 - sin2thetaW)     # hypercharge
    alpha1 = (5.0/3.0) * alphaY                 # GUT-normalized U(1)_Y
    return alpha1, alpha2, alpha3

def run_one_loop(alpha_inv_mu0, b, mu0, mu):
    """1-loop running of 1/alpha: 1/alpha(mu) = 1/alpha(mu0) - (b/2π) ln(mu/mu0)."""
    two_pi = 2.0*math.pi
    return alpha_inv_mu0 - (b / two_pi) * math.log(mu/mu0)

def match_SM_to_LR_with_ratio(alpha1_inv_SM, alpha2_inv_SM, alpha3_inv_SM, r_ratio=1.0):
    """
    At MI, identify alpha2_SM = alphaL_LR and alpha3_SM = alpha3_LR.
    For U(1): alpha1^{-1} = (3/5) alphaR^{-1} + (2/5) alphaBL^{-1}.
    Use a ratio r = alphaR^{-1} / alphaBL^{-1} at MI to split alpha1.
    Given A = alpha1^{-1}, solve:
        A = (3/5) x + (2/5) y,   with x/y = r  ->  x = r y
        => y = A / ((3/5) r + (2/5)),   x = r*y
    """
    A = alpha1_inv_SM
    denom = (3.0/5.0)*r_ratio + (2.0/5.0)
    y = A / denom
    x = r_ratio * y
    alphaL_inv = alpha2_inv_SM
    alpha3_inv = alpha3_inv_SM
    alphaR_inv = x
    alphaBL_inv = y
    return alpha3_inv, alphaL_inv, alphaR_inv, alphaBL_inv

def find_MGUT_by_alpha2_eq_alpha3(alpha2_inv_at_MI, alpha3_inv_at_MI, b2_above, b3_above, MI, log10_bracket=(14.5,17.5), max_iters=60):
    """
    Analytical solution: ln(μ/MI) = [alpha2_inv_at_MI - alpha3_inv_at_MI] / [(b2-b3)/(2π)].
    """
    two_pi = 2.0*math.pi
    delta_b = (b2_above - b3_above) / two_pi
    if abs(delta_b) < 1e-12:
        return None, None  # no intersection if slopes equal
    f_MI = alpha2_inv_at_MI - alpha3_inv_at_MI
    ln_mu_over_MI = f_MI / delta_b
    mu = MI * math.exp(ln_mu_over_MI)
    log10_mu = math.log10(mu)
    lo, hi = log10_bracket
    if log10_mu < lo or log10_mu > hi:
        mu = 10**max(min(log10_mu, hi), lo)
        log10_mu = math.log10(mu)
    return mu, log10_mu

def mismatch_metrics(alpha1_inv, alpha2_inv, alpha3_inv):
    """
    Return (max_frac_mismatch, spread, mean), where max_frac_mismatch compares each alpha^{-1} to the mean.
    """
    vals = [alpha1_inv, alpha2_inv, alpha3_inv]
    mean = sum(vals)/3.0
    spread = max(vals) - min(vals)
    max_frac = max(abs(v-mean)/mean for v in vals)
    return max_frac, spread, mean
