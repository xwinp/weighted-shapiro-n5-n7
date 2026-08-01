# Supplementary material

**Companion to:** *Exact holding region of the Tuan–Thuong weighted Shapiro cyclic inequality for $n=5$ and $n=7$.*
**Author:** Weipeng Xue, Sun Yat-sen University — xuewp5@mail2.sysu.edu.cn

This document records the computer-assisted verification details supporting the main paper: the software environment, the machine-checkable certificates and the scripts producing them, the machine records, and the reproducibility commands. The mathematical statements and proofs live in the main paper; this supplement makes every computational claim auditable.

## 1. Software environment

All certificates run under **Python 3.13.9** with **SymPy 1.14** (exact rational algebra), **mpmath 1.3** (interval arithmetic via `mpmath.iv`, with native outward-rounded `sqrt`/$n$th-root), **python-flint 0.9** (Arb, the validated-numerics backend for the $H_C$ cover), and **NumPy/SciPy** (display/sanity only, never load-bearing). Exact package versions are pinned in `code/requirements_lock.txt`.

Two independent arithmetic backends are used, both outward-directed (rigorous):
- **exact rational** — `Poly.resultant`, `Poly.count_roots`/Sturm sequences, exact `exquo` division with zero remainder, finite-field irreducibility;
- **interval** — `mpmath.iv` (monomial-wise bounds, native outward-rounded `sqrt`/$n$th-root) and Arb/`python-flint` (a hand-written lo/hi interval layer with outward-rounded $\sqrt{\cdot}$ and $\rho^{1/7}$, factored dependency-reduced evaluators). The $H_C$ cover is verified by both backends independently.

## 2. The $S_1$ crossing polynomial $\Phi_{35}(\rho)$

Leading coefficient $2^{18}$, trailing coefficient $191102976$, degree $35$:
$$\scriptsize\begin{aligned}(c_0,\ldots,c_{35})=\bigl(&191102976,\,2193551360,\,11883053056,\,41692980224,\,107322284272,\,221073893824,\\&410624553696,\,823732168256,\,1952186204080,\,4815489002744,\,10772897763040,\,20701541228760,\\&33774085926224,\,46846168825232,\,55494399759599,\,56982453166940,\,51848980402968,\,43526037225048,\\&34334035751596,\,25279161341952,\,15843550970240,\,7848268705920,\,2306869918304,\,88546441088,\\&-543801267968,\,-327985403648,\,-166554596992,\,-34693234688,\,-9171843072,\,2185738240,\\&547314432,\,444107776,\,57874432,\,19453952,\,1211392,\,262144\bigr).\end{aligned}$$

Sturm sign-variation counts: $V(0)-V(\infty)=2$; isolation $V(0)-V(2)=0$, $V(2)-V(5/2)=1$, $V(5/2)-V(3)=0$, $V(3)-V(7/2)=1$, $V(7/2)-V(\infty)=0$ (exactly two positive roots $\rho_1\in(2,5/2)$, $\rho_2\in(3,7/2)$).

## 3. Certificate scripts

The $S_1$/$S_2$ crossing and determinant certificates are univariate (Sylvester determinants of size $\le 12$). The $H_C$ critical-event enumeration forms bivariate resultants (eliminate $v$ via $E_2^u$, then $w$ via $H_C$), but each factor has degree $\le 3$ in the eliminated variable, so every intermediate has bounded degree and the scripts run in bounded memory. Scripts are grouped by theorem piece; the authoritative three-tier role map (load-bearing core vs. formal certificates vs. diagnostics) is in `code/README.md`.

**$n=5$.**
- `n5_s1_certificate.py` — stationary curve $q^3a^6-p^3a^2-p^2q=0$, $\operatorname{Res}_a=r^{17}Q_{11}(r)$ (degree-11, all-positive), $M(1)>0$.
- `n5_s1_m1_interval.py` — $a(1)\in(1.150,1.151)$, $M(1)\in[0.5662,0.5821]>\tfrac12$ by rational interval arithmetic.

**$n=7$, $S_2$.**
- `n7_resultant.py` / `n7_roots.py` — exact quotient $\operatorname{Res}_t(R,B)=p^{15}(p-1)^6F$ (exact polynomial division, remainder $0$); $F$ (degree 15) squarefree; $F\bmod 23$ irreducible (finite-field certificate $\Rightarrow$ irreducible over $\mathbb Q$ by Gauss's lemma); exactly two $(0,1)$-roots $=a_7,b_7$.
- `verify_endpoints.py` — exact Sturm isolation confirming the four cited decimals to their stated precision ($a_7,b_7$ to 11 digits, $p_1,p_2$ to 9 digits, each within $0.5$ ulp) and $\{p_1,p_2\}\subset(a_7,b_7)$.
- `n7_s2_rational_signs.py` — three-sign rational-interval samples of $B$ at $p=1/5,1/4,1/3$ (width-$4{\times}10^{-5}$ $t$-brackets), $D_0,D_5>0$ at $p=1/4$.
- `n7_inactive_sturm.py` — $G$ (degree 21), Sturm (0 roots in $(1/5,1/3)$).

**$n=7$, $S_1$.**
- `n7_s1_crossing_resultant.py` — builds $A,B$, computes $\operatorname{Res}_z=896\rho^{13}(\rho{+}1)^7(8\rho^2{+}8\rho{+}7)^6\Phi_{35}$ (exact `exquo`, remainder $0$; all 36 coefficients match), Sturm counts.
- `_verify_crossings_in_band.py` — exact Sturm isolation of $a_7,b_7,p_1,p_2$ certifying $\{p_1,p_2\}\subset(a_7,b_7)$ and $p_1,p_2<1/3<p_0$.
- `_verify_phi35_squarefree.py` — $\gcd(\Phi_{35},\Phi_{35}')=1$ (squarefree, so both crossings are simple) + Sturm zero-counts; corroboration of the three-segment sign portrait certified by `_verify_s1_three_samples.py`.
- `n7_s1_rigorous_certs.py` — `mpmath.iv`: $P_{S_1}^{\rm stat}(2/5)>7$; $\det$ sign at $z=17/20,7/8,9/10,19/20$; $w_{\rm det}(z_7)<0$; $p_0\in(3/8,2/5)$.
- `_verify_s1_three_samples.py` — three rigorous interval samples ($p=2/5,1/4,1/5$) fixing the $H_B$ three-segment sign portrait (Krawczyk-unique $z$, native outward-rounded `sqrt`).
- `n7_s1_rho_monotone_Lcert.py` — $L$-certificate: $\operatorname{Res}_w(H_B,L)=z(z{-}1)Q_7(z)$, rational-interval $w_L(z_7)<0<w_+(z_7)$, $L(w_+(1/2),1/2)<0$, hence $\rho(z)$ strictly increasing.
- `n7_s1_rho_monotone_check.py` — $F'(z)\ne0$ on $(0,1)$ (Sturm on $G=A^2-B^2\Delta$, deg 18; sole $(0,1)$-root $z_7$ extraneous).
- `n7_s1_correct_morse_trace.py` — finite-difference Hessian (diag/off formulas, three step sizes), Morse index.
- `n7_s1_all_branches_correct_hess.py` — all-branch enumeration.
- `n7_s1_hc_definitive_trace.py` — $H_C$ trace (8 admissible points, all $P>7$).
- `n7_s1_hc_kkt_check.py` — all-five-KKT-residual check at $H_C$ lifts; proves the $\beta$-variety is a superset (spurious branch fails $g_4,g_5$).
- `n7_s1_hc_Cstruct.py` — probes $P=C+2\sqrt{AB}$ ($C<7$, $2\sqrt{AB}\in(2.95,3.47)$).
- `n7_s1_hc_critical_events.py` — enumerates the five critical $z$-events (resultants $E_1$–$E_9$, sizes $\le 8$) partitioning $(0,1)$.
- `n7_s1_hc_rigorous_cert.py` — the full rigorous $H_C$ certificate: $(c,s)$ desingularisation, Krawczyk root isolation on $G$ and $E_2^{\rm red}$, mean-value-form $P$ bound, global $P_C\ge L_C=141/20$ (`mpmath.iv`, native outward-rounded `sqrt`, ivprec $130$).
- `n7_s1_hc_exhaustiveness_cert.py` — the load-bearing **branch-level exhaustiveness certificate**: exact `CRootOf` critical events; a $2\times2$ Krawczyk graph joining boxes only when their unique roots are certified to coincide at an exact rational overlap parameter; two connected root components with gap-free exact $s$-coverage; exact `CRootOf` enumeration at one sample of each of six event-free cells; and certified critical-fiber enumeration using rational isolating intervals, full-parameter interval-Newton $c$-tubes, outward quadratic $v$-enclosures, and $3\times3$ boundary Krawczyk certificates. Every admissible root is matched to one active graph component (`all_cells_branch_bijection=True`, `all_critical_admissible_roots_boxed=True`). It uses exact `CRootOf` isolation and rational interval arithmetic throughout (no numerical root-finding, number-field factorisation, or residual thresholds).
- `n7_s1_hc_completeness.py` — a numerical cross-check superseded by `n7_s1_hc_exhaustiveness_cert.py`.
- `n7_s1_hc_s1collar.py` — the $s\to1$ **second-blow-up collar** certificate: $\delta=1{-}s,\ c=\delta\bar c,\ v=\delta\bar v$ with $\bar G,\bar E_2$ regular at $\delta=0$ ($\bar G(\bar c,0)=1{-}\bar c$, $\bar E_2(\bar v,0)=\bar v{-}1$), rescaled admissibility ($\bar c,\bar v,\hat u,\widetilde K,\tilde\rho_7,a_5>0$, $O(1)$ as $\delta\to0$), a two-regime $P$ bound (direct interval $P$ for $\delta\ge10^{-10}$, monotone $\delta^{-1/7}$ crude bound below), full-parameter interval-Newton uniqueness on $\bar c,\bar v$ for every rational $\delta$-interval (including $[0,10^{-10}]$), exact-`Fraction` $\delta$-abutment; produces the $60$-piece collar `code/_hc_s1collar.json` tiling $[s_{\max},1]$ with exact seam/endpoints and re-verified $\min P_{\rm lo}^{\rm re}=7.3440544147\ldots>L_C$.
- `n7_s1_hc_cover_dump.py` — dumps the machine-checkable cover to `code/_hc_cover.json` (rational $s/c/v$ boxes, admissibility lower bounds, $L_C=141/20$).
- `n7_s1_hc_arb_refine.py` — refines the few tight real-arc $s$-pieces (whose mean-value lower bound, before refinement, lay in $(7.0487,7.05)$ — below $L_C$ due to interval dependency, not because the true value is below $L_C$; the minimum of the certified piecewise lower bounds after refinement is $\min_j P_{\rm MV}^{(j)}=7.050005023\ldots$, only $\approx 5{\times}10^{-6}$ above $L_C$) by $s$-bisection until the rigorous Arb mean-value bound clears $L_C$; produces the canonical $2604$-piece cover.
- `n7_s1_hc_cover_checker.py` — short deterministic checker (`mpmath.iv`): re-derives per piece the Krawczyk inclusion, strict admissibility, and $P_{\rm MV}>L_C$; verifies $2604/2604$ pieces.
- `n7_s1_hc_arb_checker.py` — independent Arb / `python-flint` checker (prec $150$): a hand-written lo/hi interval layer (monotone-endpoint $\sqrt{\cdot}$, $\rho^{1/7}$; 4-product signed mul) independently re-derives Krawczyk uniqueness, admissibility, and the mean-value bound; verifies $2604/2604$ pieces, with minimum certified lower bound $\min_j P_{\rm MV}^{(j)}=7.050005023\ldots>141/20$.
- `n7_s1_hc_mv_verify.py` — validity check $P_{\rm MV}\le P_{\rm true}$ at sampled lifts.

**Exact symbolic verifiers (Lemma 2.2 bridge + $S_1$ structure).**
- `verify_ny_spectral.py` — exact SymPy proof of the N–Y spectral bridge (Lemma 2.2): $\mu_k=2|\lambda_k|^2-2\operatorname{Re}\lambda_k\equiv 2(1-\cos\theta_k)[2q\cos\theta_k+2p^2-3p+2]$; $\operatorname{sign}\mu_k=\operatorname{sign}(|\lambda_k-\tfrac12|^2-\tfrac14)$, i.e. $\mu_k>0\Leftrightarrow\lambda_k$ outside the closed disk $|\lambda-\tfrac12|\le\tfrac12$; $\Delta_n=4c_n^2-4c_n-7<0$ for $n=5,7$ with the dangerous-mode bracket discriminant $=\Delta_n$ and leading coefficient $2>0$, hence every $\mu_k>0$.
- `verify_hc_closedform.py` — exact SymPy proof of the $S_1$ stationary closed form $P=A/\rho_1+B\rho_1+C$ (rational identity, numerator $\equiv0$), $g_1=\rho_1(B-A/\rho_1^2)$, and the squared check $(A/\rho_1+B\rho_1)^2|_{\rho_1^2=A/B}=4AB$, giving the exact identity $P=C+2\sqrt{AB}$.
- `verify_s1_elimination.py` — the exact $S_1$ elimination-completeness certificate: stepwise resultant elimination yielding $\operatorname{Res}_v(F_1,F_2)=w^2z^2(w-1)^2H_BH_C$ as an exact identity (`expand(R-target)==0`), proving the $S_1$ KKT locus is contained in $\{H_B=0\}\cup\{H_C=0\}$ with no third projected component; saturation prefactors $w,z,w-1$ localised away at $\beta\in(0,1)^4$.
- `verify_hb_exact_backsub.py` — a forward-construction record on $H_B$ (the load-bearing forward containment is `verify_s1_elimination.py`); reconfirms $H_B\mid\operatorname{Res}_v(F_1,F_2)$ by exact univariate polynomial division (`sp.div` in $w$, coefficients in $\mathbb Q[z]$, remainder $0$), with exact quotient $R/H_B=w^2(w-1)^2H_C$.

## 4. Machine records

- `code/_hc_cover.json` — the canonical $2604$-piece $H_C$ cover (rational $s/c/v$ boxes, admissibility lower bounds, $L_C=141/20$).
- `code/_hc_s1collar.json` — the $60$-piece second-blow-up collar tiling $[s_{\max},1]$.
- `code/_hc_exhaustiveness.json` — branch graph, six cell bijections, and five critical-fiber records.
- `code/_archival_run.txt` — captured verification output of every Tier 1 + Tier 2 certificate script (21 certificates, each `[exit 0]`).

## 5. Reproducibility

Run the full Tier 1 + Tier 2 certificate chain from the repository root:

```bash
python code/run_all_certificates.py      # runs every Tier 1 + Tier 2 certificate, asserts exit 0
```

or the Tier 1 core individually:

```bash
python code/verify_ny_spectral.py
python code/verify_hc_closedform.py
python code/verify_s1_elimination.py
python code/n7_s1_hc_s1collar.py
python code/n7_s1_hc_exhaustiveness_cert.py
python code/n7_s1_hc_cover_checker.py
python code/n7_s1_hc_arb_checker.py
```

Every load-bearing artifact is committed to a self-contained git repository with a standard SHA-256 manifest (`sha256_manifest.txt`); on a fresh extraction `sha256sum -c sha256_manifest.txt` verifies all files. The repository snapshot is identified by `git rev-parse HEAD`; a Zenodo DOI will be added upon deposition. The authoritative script role map (load-bearing core vs. formal certificates vs. diagnostics) is in `code/README.md`.
