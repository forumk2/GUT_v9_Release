import math
import numpy as np

def step_rge_gauge_two_loop(alphas, b_vec, bij_mat):
    """
    Compute d alpha_i / d t (t = ln mu) at gauge two-loop:
      d alpha_i / dt = (b_i / (2π)) * alpha_i^2 + (1 / (8π^2)) * sum_j b_ij * alpha_i^2 * alpha_j
    Inputs:
      alphas: np.array shape (n,)
      b_vec:  np.array shape (n,)
      bij_mat: np.array shape (n,n)
    Returns:
      dalphas_dt: np.array shape (n,)
    """
    two_pi = 2.0*math.pi
    term1 = (b_vec / two_pi) * (alphas**2)
    term2 = np.zeros_like(alphas)
    if bij_mat is not None:
        # alpha_i^2 * sum_j b_ij * alpha_j
        term2 = (1.0/(8.0*math.pi**2)) * (alphas**2) * (bij_mat @ alphas)
    return term1 + term2

def run_segment_two_loop(alphas_init, b_vec, bij_mat, mu0, mu1, nsteps=2000):
    """
    Integrate two-loop RGEs from mu0 to mu1 using RK4 on t = ln mu.
    alphas_init at mu0; returns alphas_out at mu1.
    """
    import math
    alphas = np.array(alphas_init, dtype=float)
    t0 = math.log(mu0)
    t1 = math.log(mu1)
    dt = (t1 - t0) / nsteps
    t = t0
    for _ in range(nsteps):
        k1 = step_rge_gauge_two_loop(alphas, b_vec, bij_mat)
        k2 = step_rge_gauge_two_loop(alphas + 0.5*dt*k1, b_vec, bij_mat)
        k3 = step_rge_gauge_two_loop(alphas + 0.5*dt*k2, b_vec, bij_mat)
        k4 = step_rge_gauge_two_loop(alphas + dt*k3, b_vec, bij_mat)
        alphas = alphas + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    return alphas

def run_segment_one_loop(alpha_inv_init, b_vec, mu0, mu1):
    """
    1-loop analytic running for inverse couplings: alpha^{-1}(mu) = alpha^{-1}(mu0) - (b/2π) ln(mu/mu0)
    """
    import math
    two_pi = 2.0*math.pi
    ln = math.log(mu1/mu0)
    return alpha_inv_init - (b_vec/ two_pi) * ln
