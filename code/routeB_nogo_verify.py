#!/usr/bin/env python3
"""
Route B: locally verify GPT's no-go claims for leg 1 (Nowosad) and the
S0 beta-system setup for leg 2.

GPT claims (must verify):
 (A) Fourier curvature identity: for DFT modes m=1,2,3, theta_m=2*pi*m/7,
     gamma_m=cos theta_m, T_i=x_i/d_i, alpha_i=p*x_{i+1}/d_i, beta_i=q*x_{i+2}/d_i,
     M0=Sum T_i=P, M1=Sum T_i*alpha_i, M2=Sum T_i*alpha_i^2:
       Q_m := Q_x(c^(m))+Q_x(s^(m)) = (7/4)(1-gamma_m)[2M2-(3+2 gamma_m)M1+(2+2 gamma_m)M0]
     and  Sum_m Q_m = 4M2-4M1+2M0 = 2 Sum_i T_i (alpha_i^2+beta_i^2) > 0.   (NOT c*(P-7))
 (B) Stronger obstruction: at p=1/4, S2 minimizer xbar=(0,1,b,c,d,0,t) with
     b=1.29267618722160523516, c=1.20199648972223459363,
     d=0.37292344231504565761, t=1.67101172500978658986,
     P(xbar)=6.9560827...<7.  Lift x(eps)=(eps,1,b,c,d,eps,t), eps=1e-3:
     P(x(eps))=6.9574199036486...<7, and log-Hessian (mod scale) is PSD
     (eigenvalues ~ 0, 3.28e-4, 1.25e-3, 6.56e-2, 2.32, 3.39, 6.18).
     -> "P<7 => exists negative curvature direction" is FALSE.
 (C) S0 beta-system: at a non-uniform S0 stationary point (interior, x>0),
     r_i=x_i/x_{i+1}, prod r_i=1, beta_i=q*r_{i+1}/(p+q*r_{i+1});
     A_i=1/[r_i(p+q r_{i+1})]; stationarity A_i+beta_{i-1} A_{i-1}=C (const);
     a_i=A_i/C satisfies cyclic recurrence a_{i+1}=1-beta_i a_i (mod 7);
     a_0 formula (13); closure K0=prod(1-beta_i)/beta_i = (p/q)^7  (i.e. rho=p/q);
     stationary value P_S0^stat = rho(1+rho)*(a0/h0)*sum_i a_i.   (17)
     Verify all against a numerical S0 stationary point at p=1/4.
"""
import numpy as np
from scipy.optimize import root

n = 7

def P_val(x, p):
    q = 1.0 - p
    return sum(x[i] / (p * x[(i+1) % n] + q * x[(i+2) % n]) for i in range(n))

def grad_P(x, p):
    q = 1.0 - p
    g = np.zeros(n)
    den = np.array([p * x[(i+1) % n] + q * x[(i+2) % n] for i in range(n)])
    for j in range(n):
        g[j] += 1.0 / den[j]
    for i in range(n):
        g[(i+1) % n] += -x[i] * p / den[i]**2
        g[(i+2) % n] += -x[i] * q / den[i]**2
    return g

def hessian_P(x, p, h=1e-5):
    H = np.zeros((n, n))
    for i in range(n):
        xp = x.copy(); xp[i] += h
        xm = x.copy(); xm[i] -= h
        H[i, i] = (P_val(xp, p) - 2 * P_val(x, p) + P_val(xm, p)) / (h * h)
    for i in range(n):
        for j in range(i+1, n):
            xpp = x.copy(); xpp[i] += h; xpp[j] += h
            xpm = x.copy(); xpm[i] += h; xpm[j] -= h
            xmp = x.copy(); xmp[i] -= h; xmp[j] += h
            xmm = x.copy(); xmm[i] -= h; xmm[j] -= h
            val = (P_val(xpp, p) - P_val(xpm, p) - P_val(xmp, p) + P_val(xmm, p)) / (4 * h * h)
            H[i, j] = val; H[j, i] = val
    return H

def log_hessian(x, p):
    """H_log = X H_x X + Diag(x_i * grad_i P). Quadratic form = log 2nd variation."""
    H = hessian_P(x, p)
    X = np.diag(x)
    g = grad_P(x, p)
    return X @ H @ X + np.diag(x * g)

def Qx_logvar(v, x, p):
    """GPT (1): Q_x(v) = Sum_i T_i[(v_i-a_i v_{i+1}-b_i v_{i+2})^2 - a_i b_i (v_{i+1}-v_{i+2})^2]."""
    q = 1.0 - p
    d = np.array([p * x[(i+1) % n] + q * x[(i+2) % n] for i in range(n)])
    T = x / d
    al = p * np.roll(x, -1) / d          # alpha_i = p x_{i+1}/d_i
    be = q * np.roll(x, -2) / d          # beta_i  = q x_{i+2}/d_i
    s = 0.0
    for i in range(n):
        vi = v[i]; vi1 = v[(i+1) % n]; vi2 = v[(i+2) % n]
        s += T[i] * ((vi - al[i]*vi1 - be[i]*vi2)**2 - al[i]*be[i]*(vi1-vi2)**2)
    return s

# -------------------- (A) Fourier identity --------------------
def check_fourier_identity():
    print("="*70)
    print("(A) Fourier curvature identity  Sum_m Q_m = 2 Sum_i T_i(alpha_i^2+beta_i^2)")
    print("="*70)
    rng = np.random.default_rng(42)
    for p in [0.25, 0.27, 0.30]:
        for trial in range(3):
            x = np.exp(rng.normal(0, 0.6, n))
            q = 1.0 - p
            d = np.array([p*x[(i+1)%n]+q*x[(i+2)%n] for i in range(n)])
            T = x/d
            al = p*np.roll(x,-1)/d
            be = q*np.roll(x,-2)/d
            M0 = T.sum(); M1 = (T*al).sum(); M2 = (T*al**2).sum()
            total = 0.0
            for m in range(1,4):
                th = 2*np.pi*m/n
                gm = np.cos(th)
                cc = np.sqrt(2.0/n)*np.cos(th*np.arange(n))
                ss = np.sqrt(2.0/n)*np.sin(th*np.arange(n))
                # CORRECTED prefactor 4/n (GPT wrote n/4 -- inverted; sum identity (5) unaffected)
                Qm_formula = (4.0/n)*(1-gm)*(2*M2-(3+2*gm)*M1+(2+2*gm)*M0)
                Qm_direct = Qx_logvar(cc,x,p)+Qx_logvar(ss,x,p)
                total += Qm_direct
                assert abs(Qm_formula-Qm_direct)<1e-6*max(1,abs(Qm_direct)), (m,Qm_formula,Qm_direct)
            rhs = 2*(T*(al**2+be**2)).sum()
            alt = 4*M2-4*M1+2*M0
            print(f" p={p} t{trial}: Sum Q_m={total:.6f}  2Sum T(a^2+b^2)={rhs:.6f}  4M2-4M1+2M0={alt:.6f}  P-7={P_val(x,p)-7:+.4f}")
            assert abs(total-rhs)<1e-6*max(1,abs(rhs))
            assert abs(total-alt)<1e-6*max(1,abs(alt))
    print(" -> (A) VERIFIED: per-mode Q_m = (4/n)(1-g_m)*B_m [GPT's n/4 prefactor was inverted; corrected to 4/n],")
    print("    and sum identity (5) holds: Sum Q_m = 2Sum T_i(a_i^2+b_i^2) = trace(H_log) > 0, NOT c*(P-7).")
    print()

# -------------------- (B) epsilon-lift obstruction --------------------
def check_eps_lift():
    print("="*70)
    print("(B) epsilon-lift obstruction: P<7 with log-Hessian PSD (mod scale)")
    print("="*70)
    p = 0.25
    b=1.29267618722160523516; c=1.20199648972223459363
    d=0.37292344231504565761; t=1.67101172500978658986
    xbar = np.array([0.0, 1.0, b, c, d, 0.0, t])
    # P at the boundary S2 point (two zeros -> terms with zero numerator drop, but
    # denominators must stay nonzero; x0=0,x5=0). Compute directly.
    Pbar = P_val(xbar, p)
    print(f" P(xbar) = {Pbar:.10f}  (GPT: 6.9560827...)  <7? {Pbar<7}")
    for eps in [1e-2, 1e-3, 1e-4]:
        xe = xbar.copy(); xe[0]=eps; xe[5]=eps
        Pe = P_val(xe, p)
        Hlog = log_hessian(xe, p)
        w = np.linalg.eigvalsh(Hlog)
        wsort = np.sort(w)
        # scale-kernel eigenvalue ~0; check the rest >= -tol
        min_nonzero = wsort[1]  # smallest after the ~0 mode
        print(f" eps={eps:.0e}: P={Pe:.10f} <7? {Pe<7}  "
              f"eig(H_log)~{[f'{e:.4e}' for e in wsort]}  PSD(mod scale)? {min_nonzero>-1e-6}")
    print(" -> (B) check: does a small lift keep P<7 AND H_log PSD (mod scale)?")
    print()

# -------------------- (C) S0 beta-system --------------------
def find_s0_stationary(p, seed=0):
    """Find a non-uniform interior stationary point of P at given p (full 7-dim KKT).
    Use log-coords y, x=exp(y), normalize x0=1 (scale). Solve grad_log P = 0 on y1..y6."""
    rng = np.random.default_rng(seed)
    def F(y):
        x = np.empty(n); x[0]=1.0; x[1:]=np.exp(y)
        g = grad_P(x, p)
        # log-gradient = x * g  (since dP/dy_i = x_i dP/dx_i); fix scale by dropping i=0
        return np.array([x[i]*g[i] for i in range(1,n)])
    best=None
    for _ in range(400):
        y0 = rng.normal(0,0.7,6)
        sol = root(F, y0, method='hybr', tol=1e-13)
        if sol.success and np.max(np.abs(sol.fun))<1e-9:
            x = np.empty(n); x[0]=1.0; x[1:]=np.exp(sol.x)
            # reject uniform (all ones) and near-periodic duplicates
            if np.max(np.abs(x-np.mean(x)))<1e-4:
                continue
            P = P_val(x,p)
            # dedup
            key = tuple(sorted(np.round(x,4)))
            if best is None or abs(P-best[1])>1e-6:
                if best is None or abs(P-best[1])>1e-4:
                    best=(x,P,key)
            break
    return best

def check_s0_beta_system():
    print("="*70)
    print("(C) S0 beta-system: recurrence / closure / stationary-value formula (17)")
    print("="*70)
    p = 0.25; q = 1.0-p
    found = None
    for s in range(40):
        r = find_s0_stationary(p, seed=s)
        if r is not None:
            found = r; break
    if found is None:
        print(" (no non-uniform S0 stationary point found by nsolve -- skip)")
        return
    x, Pnum, _ = found
    print(f" S0 stationary @ p={p}: P={Pnum:.8f}  x={np.round(x,5)}")
    print(f" grad (log) max = {max(abs(x[i]*grad_P(x,p)[i]) for i in range(n)):.2e}")
    # ratios r_i = x_i/x_{i+1}, prod = 1
    r = np.array([x[i]/x[(i+1)%n] for i in range(n)])
    print(f" prod r_i = {r.prod():.2e}  (should be 1)")
    # beta_i = q*r_{i+1}/(p+q*r_{i+1})
    beta = np.array([q*r[(i+1)%n]/(p+q*r[(i+1)%n]) for i in range(n)])
    print(f" beta = {np.round(beta,5)}  all in (0,1)? {np.all((beta>0)&(beta<1))}")
    # A_i = 1/[r_i (p+q r_{i+1})]
    A = np.array([1.0/(r[i]*(p+q*r[(i+1)%n])) for i in range(n)])
    # stationarity: A_i + beta_{i-1} A_{i-1} = C (const)
    Cvec = np.array([A[i]+beta[(i-1)%n]*A[(i-1)%n] for i in range(n)])
    C = Cvec.mean()
    print(f" A_i+beta_{{i-1}}A_{{i-1}} = {np.round(Cvec,6)}  spread={Cvec.max()-Cvec.min():.2e}  (const C={C:.6f})")
    a = A/C
    # recurrence a_{i+1} = 1 - beta_i a_i (mod 7)
    rec = np.array([1-beta[i]*a[i] for i in range(n)])
    a_next = np.array([a[(i+1)%n] for i in range(n)])
    print(f" recurrence 1-beta_i*a_i  vs  a_{{i+1}}:  max diff = {np.max(np.abs(rec-a_next)):.2e}")
    # closure K0 = prod (1-beta_i)/beta_i  vs (p/q)^7
    K0 = np.prod((1-beta)/beta)
    print(f" K0=prod(1-beta)/beta = {K0:.6f}   (p/q)^7 = {(p/q)**7:.6f}   match? {abs(K0-(p/q)**7)<1e-6}")
    rho = K0**(1/7)
    print(f" rho=K0^(1/7)={rho:.6f}  p/q={p/q:.6f}   p_recover=1/(1+rho)={1/(1+rho):.6f}")
    # a_0 formula (13): a0 = (1+prod beta)/(1 - b6 + b6 b5 - b6 b5 b4 + b6 b5 b4 b3 - ... + b6..b1)
    pb = 1.0
    for i in range(7): pb *= beta[i]
    num = 1 + pb
    # denominator: 1 - b6 + b6 b5 - b6 b5 b4 + b6 b5 b4 b3 - b6 b5 b4 b3 b2 + b6 b5 b4 b3 b2 b1
    den = 1.0; term=1.0; sgn=-1.0
    for k in range(6):  # multiply beta_6, beta_5, ..., beta_1
        term *= beta[6-k]
        den += sgn*term
        sgn *= -1
    a0_formula = num/den
    print(f" a0 formula (13) = {a0_formula:.6f}  vs a0 actual = {a[0]:.6f}  match? {abs(a0_formula-a[0])<1e-5}")
    # h_i = beta_{i-1}/[(1-beta_{i-1})(1-beta_i)]  -- try BOTH readings, report which matches
    h_div = np.array([beta[(i-1)%n]/((1-beta[(i-1)%n])*(1-beta[i])) for i in range(n)])
    h_prod = np.array([beta[(i-1)%n]*(1-beta[(i-1)%n])*(1-beta[i]) for i in range(n)])
    # s = a0/h0 ; check a_i * h_i const (= a0 h0 / ... ) actually GPT: a_i h_i independent of i, and s=a0/h0
    # From A_i = (p^2/q) s a_i  and also A_i = (p^2/q) h_i  => a_i h_i = s? Let's test: a_i * h_i const?
    for name,h in [("h=beta/((1-beta)(1-beta'))",h_div),("h=beta*(1-beta)*(1-beta')",h_prod)]:
        ah = a*h
        print(f"  [{name}]  a_i*h_i = {np.round(ah,5)} spread={ah.max()-ah.min():.2e}")
    # stationary value (17): P = rho(1+rho) * (a0/h0) * sum a_i
    for name,h in [("h_div",h_div),("h_prod",h_prod)]:
        s = a[0]/h[0]
        Pformula = rho*(1+rho)*s*np.sum(a)
        print(f"  [{name}] P_formula(17)={Pformula:.8f}  vs P_actual={Pnum:.8f}  match? {abs(Pformula-Pnum)<1e-4}")
    print()

if __name__ == "__main__":
    check_fourier_identity()
    check_eps_lift()
    check_s0_beta_system()
