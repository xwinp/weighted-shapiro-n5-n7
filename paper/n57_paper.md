# Exact holding region of the Tuan–Thuong weighted Shapiro cyclic inequality for $n=5$ and $n=7$

**薛炜鹏 (Weipeng Xue)** — 中山大学 (Sun Yat-sen University) — xuewp5@mail2.sysu.edu.cn

**Keywords:** Shapiro cyclic inequality; weighted cyclic inequality; holding region; resultant; Sturm sequence.

**2020 Mathematics Subject Classification:** 26D15, 26D20, 11C08.

## Abstract

Tuan and Thuong introduced the weighted Shapiro cyclic inequality $P_{n,p,q}(x)=\sum x_i/(p\,x_{i+1}+q\,x_{i+2})\ge n/(p+q)$ and asked, as Open Question (b), for the exact $(p,q)$-region where it holds. Both sides scale alike under $(p,q)\mapsto c(p,q)$, so we normalize $p+q=1$ (bound $n$) and determine the holding region $H_n=\{p\in(0,1):P_{n,p,1-p}(x)\ge n\ \forall\,x\}$ completely for $n=5$ and $n=7$. We prove $H_5=(0,1)$: the inequality holds for every $p$ at $n=5$, so the Tuan–Thuong sufficient interval is strict, not exact. For $n=7$, $H_7=(0,a_7]\cup[b_7,1)$ where $a_7,b_7$ are the two $(0,1)$-roots of an explicit irreducible degree-$15$ integer polynomial ($a_7\approx0.21427$, $b_7\approx0.32863$); failure occurs exactly on $(a_7,b_7)$. This first failure band, absent at $n=5$, proves the $n=5$ sufficient condition does not extend to $n=7$. Every face is settled by exact algebraic certificates — resultant elimination, rational Sturm sequences, a positive-coefficient crossing resultant, a negative-determinant saddle certificate, and the Nowosad–Yamagami interior-uniqueness theorem — each independently verified by exact division (remainder zero).

## 1. Introduction

The classical Shapiro cyclic inequality [8] $\sum x_i/(x_{i+1}+x_{i+2})\ge n/2$ is false in general; the first counterexample is $n=20$ (Lighthill [5], 1956). The modern picture is that the inequality holds for even $4\le n\le12$ and odd $3\le n\le23$, and fails for even $n\ge14$ and odd $n\ge25$ [6]; Ando [1] gives a new proof for the even cases $4\le n\le12$ and treats $n=23$. Drinfel'd [3] determined the asymptotic constant. Tuan and Thuong [9], under the title *On an Extension of Shapiro's Cyclic Inequality*, introduced the two-parameter weighted form above, completely classified $n=4$, gave the sufficient interval $[(5-\sqrt5)/10,(5+\sqrt5)/10]$ for $n=5$, and proved failure for every even $n\ge4$ when $p<q$. The exact $(p,q)$-region for odd $n\ge5$ was left open.

We normalize $p+q=1$ (bound $n$) and write $P_{n,p}:=P_{n,p,1-p}$. The value function $m_n(p)=\inf_x P_{n,p}(x)$ is continuous on $(0,1)$ (Lemma 2.1 below); the holding region is $H_n=\{p:m_n(p)\ge n\}$. The two smallest odd cases are settled here.

**Main results.**
- **Theorem 1 ($n=5$, rigorous).** $H_5=(0,1)$; equivalently $m_5(p)\equiv5$ on $(0,1)$. The Tuan–Thuong sufficient interval is strict: $[(5-\sqrt5)/10,(5+\sqrt5)/10]\subsetneq(0,1)=H_5$.
- **Theorem 2 ($n=7$, rigorous).** $H_7=(0,a_7]\cup[b_7,1)$ where $a_7<b_7$ are the two $(0,1)$-roots of the irreducible degree-$15$ integer polynomial $F$ in §4, pinned down by exact Sturm isolation as $a_7\in(0.214273512,\,0.214273525)$ and $b_7\in(0.328627665,\,0.328627686)$ (so $a_7,b_7$ are defined algebraically, not as experimental decimals). Here $H_n\subset(0,1)$ throughout, so the endpoints $0$ and $1$ are excluded by the notation. Equivalently $m_7(p)=\min\{7,P_{S_2}^{\rm curve}(p)\}$, and $m_7(p)<7\iff p\in(a_7,b_7)$.

Theorems 1–2 give the exact answer to Open Question (b) of Tuan–Thuong for $n=5,7$. The appearance of the first failure band at $n=7$ (absent at $n=5$) is a boundary-induced transition: since $a_7<\frac{5-\sqrt5}{10}\approx0.2764<b_7<\frac{5+\sqrt5}{10}$, the band $(a_7,b_7)$ satisfies $\bigl(\frac{5-\sqrt5}{10},b_7\bigr)\subset(a_7,b_7)$, so the $n=5$ sufficient condition provably does **not** extend to $n=7$.

## 2. Preliminaries

**Lemma 2.1 (continuity of $m_n$).** For each fixed odd $n$, $m_n(p)=\inf_{x\in\Delta_n}P_{n,p}(x)$ is continuous on $(0,1)$.

*Proof.* Extend $P_{n,p}$ to the closed simplex $\Delta_n=\{x\ge0:\sum x_i=1\}$ with values in $\mathbb R\cup\{+\infty\}$, setting a summand $x_i/(p\,x_{i+1}+q\,x_{i+2})$ to $+\infty$ whenever $x_i>0$ and the denominator is $0$ (and $0/0:=0$). Each summand is lower-semicontinuous (lsc) in $(x,p)$ on $\Delta_n\times(0,1)$ — it is continuous wherever the denominator is positive and jumps to $+\infty$ only when a denominator vanishes under a positive numerator — so $P$ is jointly lsc. Let $\Omega:=\{(x,p):P_{n,p}(x)<+\infty\}$ be the **admissible interior** (all relevant denominators $>0$); on $\Omega$, $P$ is $C^\infty$ (rational). Since $P_{n,p}(\mathbf1)=n$ for every $p$, we have $m_n(p)\le n$.

*Attainment (coercivity).* For fixed $p$, the sublevel $\{x\in\Delta_n:P_{n,p}(x)\le n\}$ is closed (lsc) in the compact $\Delta_n$, hence compact, and nonempty (it contains $\mathbf1$); so the infimum is attained at some $x^\star$ with $P_{n,p}(x^\star)=m_n(p)\le n<+\infty$, i.e. $x^\star\in\Omega_p$. This is the boundary blow-up / coercivity: $P_{n,p}(x)\to+\infty$ at every *infeasible* boundary point (one with $x_i>0$ but $p\,x_{i+1}+q\,x_{i+2}\to0$), uniformly in $p$ on compact subintervals (since $p,q\ge\alpha>0$), so no minimizer can sit on the infeasible boundary; the minimizer lies in the admissible interior.

*Continuity.* Let $p_k\to p$ and let $x_k$ attain $m_n(p_k)$ (so $P_{n,p_k}(x_k)\le n$). By compactness of $\Delta_n$, after passing to a subsequence $x_k\to x^\star$; by joint lsc, $m_n(p)\le P_{n,p}(x^\star)\le\liminf_k m_n(p_k)$ (lower semicontinuity). Conversely, let $x^\star$ attain $m_n(p)$; then $x^\star\in\Omega_p$, and since $\Omega$ is open and $P$ is jointly continuous on $\Omega$, $P_{n,p_k}(x^\star)\to P_{n,p}(x^\star)=m_n(p)$, whence $m_n(p_k)\le P_{n,p_k}(x^\star)$ gives $\limsup_k m_n(p_k)\le m_n(p)$ (upper semicontinuity). Thus $m_n$ is continuous. $\square$

(Equivalence of finiteness and admissibility: $P_{n,p}(x)<+\infty$ iff every denominator $p\,x_{i+1}+q\,x_{i+2}$ is strictly positive. Necessity is the lsc blow-up above. For sufficiency, if no denominator vanishes then each summand is finite and continuous, so $P$ is finite; the only way a denominator can vanish while $x\in\Delta_n$ is a maximal consecutive zero-block producing a positive numerator over a zero denominator — exactly the $+\infty$ case. Hence $\Omega=\{P<+\infty\}$ is precisely the admissible domain, and the boundary-inheritance recursion of Lemma 2.3 ranges only over viable independent zero-sets.)

Consequently $H_n=\{p:m_n(p)\ge n\}$ is closed in $(0,1)$, so its complement (the failure set) is open — a union of open bands.

**Faces and supports.** $P$ is degree-$0$ homogeneous in $x$, so we minimize on the simplex $\sum x_i=1$. A support $S=\{i:x_i>0\}$ is **viable** only if its zero-set contains no two cyclically adjacent indices (otherwise a denominator vanishes). Hence viable zero-sets are **independent sets of the cycle $C_n$**. By degree-$0$ homogeneity and Euler's theorem, at an interior minimizer on a support, $\nabla_S P=0$ (the Lagrange multiplier $\lambda$ in $\nabla P=\lambda\mathbf 1$ vanishes since $\sum x_i\partial_iP=0$). Each orbit's infimum is either an interior KKT point (a **stable** support) or is attained on its boundary (degenerating to a smaller orbit).

**Gap structure.** A zero-set decomposes the cycle into gaps; up to dihedral symmetry a zero-set is a **gap word** $(g_1,\ldots,g_m)$ with $\sum g_i=n$ and every $g_i\ge2$. The boundary of an orbit is the union of orbits obtained by enlarging the zero-set (adding one more zero while keeping the set independent).

> **Lemma 2.3 (boundary inheritance).** For a precise-support stratum $\Sigma_Z:=\{x\in\Delta_n:\{i:x_i=0\}=Z\}$ with $Z$ an independent set of $C_n$,
> $$\inf_{\Sigma_Z} P \;=\; \min\Bigl\{\,\inf\bigl\{P(x):x\in\Sigma_Z\text{ is an interior local-minimum KKT point}\bigr\},\ \ \min_{Z'\supsetneq Z,\ Z'\text{ independent}}\inf_{\Sigma_{Z'}}P\,\Bigr\},$$
> with the convention $\inf\varnothing=+\infty$, and $\inf_{\Sigma_{Z'}}P=+\infty$ when $Z'$ is non-viable (contains adjacent zeros).
>
> *Proof.* Extend $P$ to the closed face $\overline{\Sigma_Z}=\{x\in\Delta_n:\mathrm{supp}(x)^c\supseteq Z\}$ as in Lemma 2.1 (lsc, $+\infty$ on infeasible points); $P<+\infty$ iff every denominator is strictly positive, so the admissible domain equals the union of precise-support strata $\bigcup_{Z'\supseteq Z,\,Z'\text{ ind.}}\Sigma_{Z'}$. Take a minimizing sequence $x_k\in\Sigma_Z$ with $P(x_k)\to\inf_{\Sigma_Z}P$; by compactness of $\Delta_n$, after passing to a subsequence $x_k\to x^\ast\in\overline{\Sigma_Z}$. If $\mathrm{supp}(x^\ast)^c=Z$ (the support did not shrink), then $x^\ast\in\Sigma_Z$ and, being an interior minimizer on $\Sigma_Z$, satisfies $\nabla_{\Sigma_Z}P=0$ (Euler, $\lambda=0$) — an interior KKT value. Otherwise the zero-set enlarged to an independent $Z'\supsetneq Z$, so $x^\ast\in\Sigma_{Z'}$ and $P(x^\ast)\ge\inf_{\Sigma_{Z'}}P$ (lsc). If the enlargement ever creates two cyclically adjacent zeros, take a maximal consecutive zero-block: the positive component just before it has a numerator bounded away from $0$ while its denominator (a sum involving a zero neighbour) tends to $0$, so $P\to+\infty$ — hence non-viable $Z'$ contributes $+\infty$ and never attains the infimum. Conversely, each viable child $Z'$ is approachable from $\Sigma_Z$ by lifting its new zero to $\varepsilon\downarrow0$ (independence keeps all denominators $>0$ for small $\varepsilon$), so $\inf_{\Sigma_{Z'}}P$ is indeed inherited from below. Taking the infimum over the interior local-minimum KKT values and the minimum over all viable children gives the displayed formula. $\square$

This is the stratum recursion used throughout §3–§4: $\inf_{S_0}P=\min\{7,\inf_{S_1}P,\ldots\}$, etc. (Here $S_0$ denotes the *open* full-support stratum $\{x_i>0:\sum x_i=1\}$, distinct from the closed simplex $\Delta_n$; the minimizer on $\Delta_n$ is either the interior point $\mathbf1$ or lies on a viable boundary stratum, which is exactly what the lemma formalises.) Yamagami [11, §5, pp.524–527] likewise notes that interior Hessian uniqueness does not control the boundary, which must be analysed separately — the stratum recursion is that analysis.

**Nowosad–Yamagami theorem** (Nowosad [7], CPAM 21 (1968), Thm 1.8; Yamagami [11], Proc. AMS 118 (1993), 521–527, Thm 1, Lemma 3, Cor. 5). We state the finite-dimensional specialisation we need, in the commutative C$^\ast$-algebra $A=\mathbb C^n$ (entrywise product, state $\varphi=\sum_i$, self-adjoint part $A_s=\mathbb R^n$). For a linear $T:A\to A$ write $X_T(y)=\varphi(y^{-1}T(y))=\sum_i T(y)_i/y_i$ on $A_{++}$.

(NY-1) *Constancy on the closed subgroup* (Yamagami Thm 1, restating Nowosad Thm 1.8). For $a\in A_{++}$, let $[a]\subset A_s^\times$ denote Yamagami's closed subgroup consisting of those self-adjoint invertibles $x$ for which both $x$ and $x^{-1}$ are norm-limits of Laurent polynomials in $a$; it contains $a$ (so it is *not* the geometric segment $\{(1{-}t)\mathbf1+ta:0\le t\le1\}$, and in general is neither compact nor contained in the positive cone). If $T$ is $\ast$-preserving ($T(x^\ast)=T(x)^\ast$) and $X_T$ has local maxima at $\mathbf1$ and at $a\in A_s\setminus\mathbb R\mathbf1$, then $X_T$ is **constant on $[a]$**. In finite dimension $a^t=\exp(t\log a)\in[a]$ for every $t\in\mathbb R$ (by polynomial interpolation on the finite spectrum of $a$); constancy gives $\tfrac{d^2}{dt^2}X_T(a^t)\big|_{t=0}=0$, so the tangent direction $\log a$ (mod the scale direction $\mathbf1$, nonzero when $a$ is non-scalar) lies in $\ker\nabla^2X_T(\mathbf1)$. Consequently, if $\ker\nabla^2X_T(\mathbf1)=\mathrm{span}\{\mathbf1\}$, no such non-scalar $a$ exists: $X_T$ has at most one isolated interior local maximum, up to the scalar ray.

(NY-2) *Spectral non-degeneracy* (Yamagami Lemma 2 + Lemma 3; Cor. 5). Yamagami's Hessian formula (Lemma 2) gives, for $P=f_S$ with $S\mathbf1=\mathbf1$ ($s=1$), the Hessian on $\mathbf1^\perp$ diagonalised by the Fourier modes; the same Fourier calculation as in Cor. 5 yields the eigenvalue $\mu_k=2|\lambda_k|^2-2\operatorname{Re}\lambda_k$ (identically $2(1-\cos\theta_k)[2q\cos\theta_k+2p^2-3p+2]$, verified below), where $\lambda_k=p\omega^k+q\omega^{2k}$ ($k\ne0$). Cor. 5 itself states the **interior-disk** condition $|\lambda_k-\tfrac12|<\tfrac12$ for the $f_S$-*maximum* case ($\mu_k<0$, the negative Hessian $s(S+S^\top)-2S^\top S$ being positive semidefinite); our $P=f_S$ *minimum* case is the **symbol-dual**: $\mu_k>0$ iff $|\lambda_k-\tfrac12|>\tfrac12$ (exterior disk), equivalently $\operatorname{sign}\mu_k=\operatorname{sign}(|\lambda_k-\tfrac12|^2-\tfrac14)$ (verified in `verify_ny_spectral.py`). Thus $\mu_k>0$ for all $k\ne0$ makes $\mathbf1$ a **strict local minimum** with $\ker\nabla^2P(\mathbf1)=\mathrm{span}\{\mathbf1\}$, and by (NY-1) the **unique** interior local-minimum ray.

> **Lemma 2.2 (interior uniqueness, finite-dimensional N–Y).** Let $S=pC+qC^2$ ($p,q>0$, $p+q=1$, $n$ odd) and $P(x)=\sum_{i=0}^{n-1}x_i/(Sx)_i$ on the positive cone. If $\nabla^2P(\mathbf1)\succeq0$ with $\ker\nabla^2P(\mathbf1)=\mathrm{span}\{\mathbf1\}$ (equivalently every $\mu_k>0$; for $n=5,7$ this follows uniformly in $p$ from $\Delta_n<0$, proved in §3.1 and §4.1), then $\mathbf1$ is the **unique interior local-minimum ray** of $P$ on the full-support face. The inequality $P\ge n$ on the whole face is *not* a consequence of this lemma; it requires a separate boundary analysis (Yamagami [11], §5, pp.524–527), performed here by stratum recursion (boundary inheritance) in §3.4 ($n=5$) and §4.5 ($n=7$).

*Proof.* (1) **$S$ is invertible.** $S$ is circulant with Fourier eigenvalues $\lambda_k=\omega^k(p+q\omega^k)$; $\lambda_k=0$ would force $\omega^k=-p/q$ of modulus $1$, hence $p=q$ and $\omega^k=-1$, impossible for odd $n$. (2) **$S\mathbf1=\mathbf1$ and $S^\top\mathbf1=\mathbf1$** (each column sums to $p+q=1$), hence $S^{-1}\mathbf1=\mathbf1$. (3) **$T=-S^{-1}$ is $\ast$-preserving, not self-adjoint.** $S$ is real, hence so is $S^{-1}$, so $T(x^\ast)=T(x)^\ast$ (complex conjugation commutes with the real matrix $T$). This is the structure N-Y requires. ($S$ itself is *not* self-adjoint — $S^\top=pC^{-1}+qC^{-2}\ne S$ for $n\ge3$ — and self-adjointness of $T$ is neither needed nor true in general; only $\ast$-preservation is.) (4) **Sign bridge.** The linear map $S:A_{++}\to S(A_{++})$ is a diffeomorphism onto its image (by (1), $S$ is invertible); $S(A_{++})$ is open in $\mathbb R^n$ and contains a neighbourhood of $S\mathbf1=\mathbf1$. Since $P(x)=\sum_i(S^{-1}y)_i/y_i=X_{S^{-1}}(y)$ and $X_{-S^{-1}}=-X_{S^{-1}}$ by linearity of $\varphi$, a local minimum of $P$ at $\mathbf1$ corresponds to a local maximum of $X_{-S^{-1}}$ at $\mathbf1$ *within $S(A_{++})$*; because $S(A_{++})$ is a neighbourhood of $\mathbf1$, this is also a local maximum in $A_{++}$, with $X_{-S^{-1}}(\mathbf1)=-n$. (We do not claim $S^{-1}$ preserves the positive cone — only the local correspondence at $\mathbf1$ is needed.) (5) **Hessian transfer (congruence).** Differentiating $X_T(y)=-P(S^{-1}y)$ twice at $y=\mathbf1$ (with $T=-S^{-1}$) gives, in the same coordinate basis, the **congruence**
$$\nabla_y^2 X_T(\mathbf1)\;=\;-S^{-\top}\,\nabla_x^2 P(\mathbf1)\,S^{-1}.$$
Since $S\mathbf1=\mathbf1$, this congruence flips the inertia sign ($H\succeq0\iff \nabla_y^2X_T(\mathbf1)\preceq0$) and preserves the kernel: $\ker\nabla_y^2X_T(\mathbf1)=S\,\ker\nabla_x^2P(\mathbf1)=\mathrm{span}\{\mathbf1\}$. Thus $\mathbf1$ is a non-degenerate strict local maximum of $X_T$ modulo scale — the Hessian hypothesis of (NY-1). (6) **Uniqueness via (NY-1).** Suppose $P$ had a second non-scalar interior local-minimum ray $x^\ast$. Then $y^\ast=Sx^\ast$ is a second local maximum of $X_T$; by (NY-1) $X_T$ is constant on the closed subgroup $[y^\ast]$, so the curve $(y^\ast)^t=\exp(t\log y^\ast)\in[y^\ast]$ has $\tfrac{d^2}{dt^2}X_T((y^\ast)^t)\big|_{t=0}=0$, placing the non-scale tangent $\log y^\ast$ in $\ker\nabla_y^2X_T(\mathbf1)$. This contradicts $\ker\nabla_y^2X_T(\mathbf1)=\mathrm{span}\{\mathbf1\}$ from (5). Hence no second interior local minimum exists: $\mathbf1$ is the unique one. (7) This is an *interior* statement on the full-support face; the global infimum may lie on a boundary stratum and is closed separately by stratum recursion. $\square$

**Uniform Hessian spectrum.** At $x=\mathbf 1$, the Hessian of $P_{n,p}$ on $\mathbf 1^\perp$ is diagonalized by the Fourier modes $\theta_k=2\pi k/n$ ($k=1,\dots,n-1$):
$$\mu_k=2(1-\cos\theta_k)\bigl[\,2q\cos\theta_k+2p^2-3p+2\,\bigr].$$
The most dangerous mode is $k=(n-1)/2$; writing $c_n=-\cos(\pi/n)$, the bracket there is a quadratic in $p$ with discriminant $\Delta_n=4c_n^2-4c_n-7$. For $n=5,7$, $\Delta_n<0$ (proved rigorously in §3.1 and §4.4), so the bracket is strictly positive for all $p\in(0,1)$, hence every $\mu_k>0$: the uniform point is a **strict** local minimum with $\ker\nabla^2P(\mathbf 1)=\mathrm{span}\{\mathbf 1\}$, and Lemma 2.2 applies.

**Certificate toolbox.** All algebraic assertions are obtained by: (i) eliminate the auxiliary variable(s) via a **resultant**; (ii) count real roots in an interval by an exact rational **Sturm sequence** (recording the sign-variation counts at the endpoints); (iii) certify sign-constancy by a **positive-coefficient** resultant factor (a polynomial with all coefficients of one sign is sign-definite on $(0,\infty)$, hence has no positive root), or by resultant+Sturm on a degeneracy locus; (iv) isolate each real root in a rational interval and certify signs by rational interval arithmetic. Every polynomial below is locally reconstructed and verified by exact division (quotient and remainder recorded, remainder zero). The certificate scripts are listed in the supplement.

## 3. The $n=5$ classification (Theorem 1, rigorous)

The independence number of $C_5$ is $\alpha(C_5)=2$, so viable zero-sets have at most two elements. Up to the dihedral action of $C_5$ there are exactly **three** orbits:

| orbit | zero-set | $\lvert$support$\rvert$ | gap word | type |
|---|---|---|---|---|
| $S_0$ | $\emptyset$ | 5 | — | full (uniform) |
| $S_1$ | $\{0\}$ | 4 | $(5)$ | 1-zero, one-gap $L=5$ |
| $S_2$ | $\{0,2\}$ | 3 | $(2,3)$ | 2-zero, one-gap $L=3$ (**AM–GM**) |

We show $P>5$ on $S_1,S_2$ and $P=5$ only at the uniform point of $S_0$; boundary degeneration then gives $m_5(p)=5$.

### 3.1 $S_0$ — Nowosad–Yamagami
At $x=\mathbf 1$, $P=5$. The uniform Hessian spectrum (§2) has $\Delta_5=4c_5^2-4c_5-7$ with $c_5=-\cos(\pi/5)=-\frac{1+\sqrt5}{4}$. Then $4c_5^2=\frac{3+\sqrt5}{2}$ and $4c_5=-(1+\sqrt5)$, so
$$\Delta_5=\tfrac{3+\sqrt5}{2}+(1+\sqrt5)-7=\tfrac{3(\sqrt5-3)}{2}\approx-1.146<0.$$
Hence every $\mu_k>0$ on $\mathbf 1^\perp$, so $\nabla^2P(\mathbf 1)\succeq0$ with kernel $\mathrm{span}\{\mathbf 1\}$. By Lemma 2.2, $\mathbf 1$ is the **unique interior local-minimum ray** of $P$ on the full-support stratum $S_0$ for every $p\in(0,1)$, with value $P(\mathbf1)=5$. This is a *local* statement; the face infimum $\inf_{S_0}P$ is closed by **boundary inheritance** (Lemma 2.3): the minimizer on the closed face $\overline{S_0}=\Delta_5$ is either the interior critical point $\mathbf1$ (value $5$) or lies on a viable boundary stratum, so $\inf_{S_0}P=\min\{5,\inf_{S_1}P,\inf_{S_2}P\}$.

### 3.2 $S_2$ — AM–GM closed form
With $x_0=x_2=0$ and support $\{1,3,4\}$,
$$P=\frac{x_1}{q\,x_3}+\frac{x_3}{p\,x_4}+\frac{x_4}{q\,x_1},$$
whose three summands have product $1/(p\,q^2)$. By AM–GM,
$$P\ge 3\,(p\,q^2)^{-1/3}=:M_{5,3}(p),\qquad\text{with equality iff }x_1/(qx_3)=x_3/(px_4)=x_4/(qx_1).$$
The minimum of $M_{5,3}$ over $p\in(0,1)$ is at $p=1/3$ (where $p\,q^2$ is maximized):
$$\min_p M_{5,3}=\frac{3}{\bigl(\tfrac13(\tfrac23)^2\bigr)^{1/3}}=\frac{9}{\sqrt[3]{4}}\approx5.6696>5,$$
rigorously since $9^3=729>500=5^3\cdot4$. Hence **$S_2$ never fails**: $\inf_{S_2}P\ge9/\sqrt[3]{4}>5$.

### 3.3 $S_1$ — stationary branch, positive-coefficient certificate
This is the only non-closed-form face. Set $x_0=0$, $x_4=1$ (homogeneity), $a=x_1,b=x_2,c=x_3$, $q=1-p$. Then
$$P=\frac{a}{p\,b+q\,c}+\frac{b}{p\,c+q}+\frac{c}{p}+\frac{1}{q\,a}.$$
The KKT equations $\partial_aP=\partial_bP=\partial_cP=0$ (Euler, $\lambda=0$) reduce exactly to
$$p\,b+q\,c=q\,a^2,\qquad p\,c+q=\frac{q^2 a^3}{p},\qquad\boxed{\,q^3 a^6-p^3 a^2-p^2 q=0\,}\quad(\text{stationary curve}).$$
(The first two give $c=(q^2a^3-pq)/p^2$, $b=q(a^2-c)/p$; substituting into $\partial_cP=0$ and cancelling the $pq^2a^3$ pair yields the curve.)

**Dehomogenization.** Write $\widetilde P_r:=P_{5,r,1}$ (parameters $p:q=r:1$); since both sides of the inequality scale identically under $(p,q)\mapsto c(p,q)$, the normalized value at $p=r/(1+r),q=1/(1+r)$ is $(1+r)\widetilde P_r$. The curve becomes $a^6-r^3a^2-r^2=0$, which for each $r>0$ has a **unique positive root** $a(r)$ (with $u=a^2$: $u^3-r^3u-r^2=0$ has $f(0)<0$, $f\to+\infty$, and a negative minimum, hence one positive crossing). The branch is a smooth connected curve over $r>0$.

**Positivity of the support point.** With $p:q=r:1$ the recovery formulas are $c=(a^3-r)/r^2$ and $b=(a^2-c)/r$. Put $y:=a^3/r>0$; the curve gives $r a^2=y^2-1$. Since $a,r>0$ we have $ra^2>0$, hence $y^2>1$ and, as $y>0$, $y>1$. Therefore $c=(y-1)/r>0$, and $a^2-c=(y^2-1)/r-(y-1)/r=y(y-1)/r>0$, so $b=(a^2-c)/r>0$. The eliminated branch thus lies entirely in the positive support.

**No-crossing certificate.** The degree-$0$ target is $M:=(1+r)\widetilde P_r-5$. Let $N(a,r)=\operatorname{num}((1+r)\widetilde P_r-5)$ on the branch. Eliminating $a$:
$$\operatorname{Res}_a\!\bigl(a^6-r^3a^2-r^2,\;N(a,r)\bigr)=r^{17}\,Q_{11}(r),$$
where (coefficients high-to-low)
$$\begin{aligned}
Q_{11}(r)={}&1024r^{11}+3644r^{10}+6360r^{9}+11580r^{8}+15895r^{7}+28772r^{6}\\
&+23892r^{5}+21120r^{4}+9360r^{3}+7680r^{2}+2368r+1728.
\end{aligned}$$
**Every coefficient of $Q_{11}$ is positive**, so $Q_{11}(r)>0$ for all $r>0$. Since $r>0$ on the branch and $Q_{11}$ is (up to the prefactor $r^{17}$) the crossing resultant, $M$ **never vanishes** on the stationary curve: $P\neq5$ there. The sign is constant on the connected branch, and it suffices to evaluate at one point.

**Sign certificate at $r=1$.** At $r=1$ ($p=q=1/2$) the curve is $a^6-a^2-1=0$; its unique positive root is isolated in $a\in(1.150,\,1.151)$ ($f(1.150)=-604111/64000000<0<f(1.151)=354174281094401/10^{15}$, and $f'(a)=6a^5-2a>0$ on this interval, so the root is unique there). On this interval rational interval arithmetic gives
$$M(1)=\tfrac{6}{a}+\tfrac{2}{a^3}+2a^3-9\;\in\;\bigl[\,0.5662,\;0.5821\,\bigr]\;>\;\tfrac12>0$$
(the simplification uses $b+c=a^2$, $c+1=a^3$ on the branch). Hence **$P>5$ on the entire $S_1$ stationary branch**. The boundary $\partial S_1=S_2$ has $P>5$ (§3.2), so $\inf_{S_1}P>5$.

### 3.4 Conclusion of Theorem 1
The three orbits exhaust all viable supports ($\alpha(C_5)=2$). **Boundary inheritance:** for each face, $P$ is lsc on the compact face-closure and blows up ($+\infty$) on its *infeasible* boundary, so its infimum is attained either at an interior critical point or on a viable boundary stratum. Thus
$$\inf_{S_0}P=\min\{\,P(\mathbf1),\,\inf_{S_1}P,\,\inf_{S_2}P\,\}=\min\{5,\,>5,\,>5\}=5,$$
using Lemma 2.2 (unique interior local min on $S_0$, value $5$), $\inf_{S_1}P>5$ (§3.3), and $\inf_{S_2}P\ge9/\sqrt[3]4>5$ (§3.2). The recursion $S_0\to S_1\to S_2\to\varnothing$ (each $\partial S_k$ the next viable stratum, $S_2$ terminal) closes every face. Therefore
$$m_5(p)=5\quad\forall p\in(0,1),\qquad H_5=(0,1).$$
The Tuan–Thuong interval $[(5-\sqrt5)/10,(5+\sqrt5)/10]\approx[0.2764,0.7236]$ is a strict subset of $H_5$: sufficient but not exact. $\blacksquare$

## 4. The $n=7$ classification (Theorem 2, rigorous)

The independent sets of $C_7$ fall into five dihedral orbits:

| orbit | zero-set | $\lvert$support$\rvert$ | gap word | type |
|---|---|---|---|---|
| $S_0$ | $\emptyset$ | 7 | — | full (uniform) |
| $S_1$ | $\{0\}$ | 6 | $(7)$ | 1-zero, $L=7$ |
| $S_2$ | $\{0,2\}$ | 5 | $(2,5)$ | 2-zero, $L=5$ (**main / failing**) |
| $S_3$ | $\{0,3\}$ | 5 | $(3,4)$ | 2-zero |
| $S_4$ | $\{0,2,4\}$ | 4 | $(2,2,3)$ | 3-zero, $L=3$ |

The boundary of an orbit is obtained by enlarging the zero-set: $\partial S_0=S_1$, $\partial S_1=S_2\cup S_3$ (a second zero is either $\{0,2\}\sim S_2$ or $\{0,3\}\sim S_3$), $\partial S_2\subset S_4$, $\partial S_3\subset S_4$. Since $\alpha(C_7)=3$, no viable zero-set has more than three elements; the three-zero orbit is $S_4$ (unique, gap word $(2,2,3)$, the only partition of $7$ into three parts $\ge2$), and these **five orbits exhaust all viable supports**. $S_4$ is *terminal*: adding a fourth zero to $\{0,2,4\}$ always creates two cyclically adjacent zeros (positions $1,3,5,6$ each border an existing zero), so $\partial S_4$ is entirely non-viable ($P\to+\infty$). Thus the boundary recursion $S_0\to S_1\to\{S_2,S_3\}\to S_4\to\varnothing$ terminates, and $m_7(p)=\min\{\inf_{S_k}P:0\le k\le4\}$ is controlled by the five face certificates below.

### 4.1 Closed forms: $S_0,S_4,S_3$
- **$S_0$ (uniform).** $x_i\equiv c$ gives $P=7$. Write $u_7=\cos(\pi/7)\in(9/10,91/100)$ (it is the root of $8u^3-4u^2-4u+1=0$ in $(0,1)$; $f(9/10)=-1/125<0<f(91/100)$), and $c_7=-u_7$ as in §2. The discriminant is $\Delta_7=4u_7^2+4u_7-7=4c_7^2-4c_7-7$; the function $g(u)=4u^2+4u-7$ is increasing ($g'(u)=8u+4>0$) with $g(9/10)=-4/25<0$ and $g(91/100)<0$, so $\Delta_7=g(u_7)<0$. Hence every uniform-Hessian eigenvalue $\mu_k>0$, Lemma 2.2 applies, and $\mathbf1$ is the **unique interior local-minimum ray** on $S_0$ (value $7$). This is a *local* statement; the face infimum is closed by **boundary inheritance** (Lemma 2.3): the minimizer on the closed face $\overline{S_0}=\Delta_7$ is either $\mathbf1$ (value $7$) or lies on a viable boundary stratum, so $\inf_{S_0}P=\min\{7,\inf_{S_1}P,\inf_{S_2}P,\inf_{S_3}P,\inf_{S_4}P\}$ — *not* pointwise $P\ge7$ on the open stratum $S_0$.
- **$S_4$ — AM–GM.** With $x=(0,1,0,b,0,c,d)$, $q=1-p$, the four summands are $\frac1{qb},\frac b{qc},\frac c{pd},\frac dq$, whose product is $\frac1{p\,q^3}$. AM–GM gives
$$P_{S_4}\ge 4\,(p\,q^3)^{-1/4},\qquad\text{equality iff }\tfrac1{qb}=\tfrac b{qc}=\tfrac c{pd}=\tfrac dq,$$
i.e. $b=(p/q)^{1/4},c=(p/q)^{1/2},d=(q/p)^{1/4}$ (the KKT solution). The minimum over $p$ is at $p=1/4$ (maximizing $p\,q^3$):
$$\min_p P_{S_4}=\frac{16}{3^{3/4}}\approx7.0190614>7\qquad(16^4=65536>64827=7^4\cdot3^3).$$
**$S_4$ never fails.** All finite boundaries of $S_2,S_3$ lie in $S_4$, hence are $>7$.
- **$S_3$ — exact sign certificate.** With $r=(q/p)^{1/5}$, the KKT solution is $b=r^2,c=r^{-1},e=r^3,d=r(1-r^7)$; $d>0\iff r<1\iff p>1/2$. For $p\le1/2$ (the whole failure band $(a_7,b_7)\subset(0,1/3)$) there is **no** positive interior stationary point, so $\inf_{S_3}P$ is attained on $\partial S_3\subset S_4$ and is $>7$. For $p>1/2$ (so $0<r<1$) the stationary value is
$$P_{S_3}^{\rm stat}(r)=\frac{(1+r^5)(5-r^7)}{r^2}.$$
The **exact identity**
$$\begin{aligned}
\bigl(1+r^5\bigr)\bigl(5-r^7\bigr)-8r^2&=(1-r)\,H(r),\\
H(r)&=r^{11}+r^{10}+r^9+r^8+r^7+2r^6+2r^5-3r^4-3r^3-3r^2+5r+5.
\end{aligned}$$
holds (verified by expansion). For $0<r<1$, $H(r)\ge 5+5r-9r^2\ge1$: indeed $H(r)-(5+5r-9r^2)=r^2\!\bigl(r^9+r^8+r^7+r^6+r^5+2r^4+2r^3-3r^2-3r+6\bigr)$ and the bracket is $\ge 6-3r^2-3r=3(2-r^2-r)>0$ on $(0,1)$; and $5+5r-9r^2$ is concave with endpoint values $5,1$, hence $\ge1$. Therefore $P_{S_3}^{\rm stat}-8=\frac{(1-r)H(r)}{r^2}>0$, i.e.
$$P_{S_3}^{\rm stat}(r)>8>7\qquad(0<r<1).$$
**$S_3$ never fails.**

### 4.2 $S_2$ (main) — the only failing orbit
$x=(0,1,b,c,d,0,e)$, $e=t$. Euler ($\lambda=0$) gives $d=t^2$, $c=q(q-pt^4)/(p^2t^2)$, $b=q/(pt)-q^2(q-pt^4)/(p^3t^2)$, with the **stationary curve**
$$R(p,t)=q^3-p^3t^5-p^2qt^8=0.$$
**Positive lift.** Put $\rho=q/p$ (so $q=\rho p$); $R=0$ reads $\rho^3=t^5+\rho t^8$. As a function of $t>0$ the left-minus-right side is strictly decreasing ($-5t^4-8\rho t^7<0$), with value $\rho^3>0$ at $t=0$ and value $-\rho^{5/4}<0$ at $t=\rho^{1/4}$; hence the unique positive root satisfies $t<\rho^{1/4}$, i.e. $t^4<\rho$. From $R=0$, $\rho(\rho-t^4)(\rho+t^4)=\rho^3-\rho t^8=t^5$. Therefore
$$c=\frac{\rho(\rho-t^4)}{t^2}>0,\qquad d=t^2>0,\qquad b=\frac{\rho}{t^2}\Bigl[t-\rho(\rho-t^4)\Bigr]=\frac{\rho}{t^2}\Bigl[t-\frac{t^5}{\rho+t^4}\Bigr]=\frac{\rho^2}{t(\rho+t^4)}>0,$$
so the stationary point is strictly positive for every $p\in(0,1)$.

On $R=0$, $P_{\rm curve}-7=\frac{B(p,t)}{q\,p^2}+\frac{R\cdot(\cdots)}{q\,p^3t^4}$, so $P_{\rm curve}=7\iff B=5p^2t+2pqt^4-2q^2-7qp^2=0$. Eliminating $t$:
$$\operatorname{Res}_t(R,B)=p^{15}(p-1)^6\,F(p),$$
where $F$ is the irreducible degree-$15$ integer polynomial ($5764801=7^8$ leading, $65536=2^{16}$ trailing):
$$\begin{aligned}
F(p)={}&5764801 p^{15}-47765494 p^{14}+190003135 p^{13}-486209703 p^{12}\\
&+901678743 p^{11}-1287828143 p^{10}+1464952167 p^{9}-1351039522 p^{8}\\
&+1017028633 p^{7}-624621984 p^{6}+310300032 p^{5}-122238368 p^{4}\\
&+36836352 p^{3}-7952896 p^{2}+1073408 p-65536.
\end{aligned}$$
**Irreducibility:** $F\bmod23\in\mathbb F_{23}[p]$ is irreducible (leading coeff $7^8\not\equiv0$), so $F$ is irreducible over $\mathbb Q$ (Gauss's lemma [4]). **Sturm:** $F$ has exactly two roots in $(0,1)$ — these are
$$a_7\approx0.21427352090984097\in(1/5,1/4),\qquad b_7\approx0.32862767791659197\in(1/4,1/3)$$
(a third real root $\approx1.3266$ lies outside).

**Three-sign certificate (rational samples).** Isolating the unique positive root $t(p)$ of $R(p,\cdot)$ in a rational bracket (width $4\!\times\!10^{-5}$, $R$ monotone) and evaluating $B$ by exact rational interval arithmetic gives the sign of $P_{\rm curve}-7$ at three rational parameters spanning the three regions:
$$\begin{array}{c|c|c}
p & \text{region} & B(p,t(p))\\\hline
1/5 & p<a_7 & B\in\!\left[\tfrac{30\,116\,709\,200\,321\,381\,334\,555\,841}{31\,250\,000\,000\,000\,000\,000\,000\,000},\,\tfrac{34\,746\,273\,526\,146\,252\,534\,263\,041}{31\,250\,000\,000\,000\,000\,000\,000\,000}\right]>0\\[2pt]
1/4 & a_7<p<b_7 & B\in\!\left[\tfrac{-170\,398\,618\,635\,119\,584\,219\,204\,477}{8\!\cdot\!10^{28}},\,\tfrac{-159\,030\,250\,817\,069\,792\,431\,615\,677}{8\!\cdot\!10^{28}}\right]<0\\[2pt]
1/3 & p>b_7 & B\in\!\left[\tfrac{130\,644\,355\,439\,570\,080\,856\,803}{263\,671\,875\,000\,000\,000\,000\,000\,000},\,\tfrac{164\,831\,571\,916\,353\,743\,995\,603}{263\,671\,875\,000\,000\,000\,000\,000\,000}\right]>0
\end{array}$$
Together with $\operatorname{sign}(P_{\rm curve}-7)=\operatorname{sign}B$ on $R=0$ and the Sturm count (exactly two $P=7$ crossings, at $a_7,b_7$), this fixes $P_{\rm curve}>7$ on $(0,a_7)\cup(b_7,1)$ and $P_{\rm curve}<7$ on $(a_7,b_7)$ (minimum $\approx6.951$ at $p\approx0.27$).

**Inactive KKT (zeros stay zero).** For $S_2$ to give the global minimum the inactive derivatives $D_0=\partial_{x_0}P|_{x_0=x_5=0}$, $D_5=\partial_{x_5}P|_{x_0=x_5=0}$ must be $\ge0$. The curve $R=0$ is a single branch ($\partial R/\partial t<0$ for $t>0$; one positive root $t(p)$ per $p$). Eliminating $t$:
$$\operatorname{Res}_t(R,\mathrm{num}\,D_0)=-p^{23}(p-1)^{11}G(p),\quad \operatorname{Res}_t(R,\mathrm{num}\,D_5)=p^{19}(p-1)^{18}G(p),$$
sharing the degree-$21$ factor $G$ (listed in §4.5). **Sturm:** $G$ has $0$ roots in $(1/5,1/3)\supset(a_7,b_7)$; its sole $(0,1)$-root is $\approx0.39387>1/3>b_7$. Hence $D_0,D_5$ do not vanish on the curve in the band. Their sign is fixed by the rational sample $p=1/4$ (same $t$-isolation as above): $D_0\in[0.29319,0.29464]>0$, $D_5\in[1.04553,1.04895]>0$. Thus $D_0,D_5$ are **strictly positive on all of $(a_7,b_7)$**. Since $R(p,\cdot)$ has exactly one positive root, $S_2$ has a **unique** interior stationary point, and every finite boundary of $S_2$ lies in $S_4$ where $P>7$ (§4.1). By boundary inheritance (Lemma 2.3), whenever $P_{S_2}^{\rm curve}<7$ — i.e. on $(a_7,b_7)$ — this unique stationary point attains $\inf_{S_2}P$: the boundary values exceed $7$ while the interior stationary value is $<7$, so the infimum cannot lie on $\partial S_2$, and the only interior candidate is this point. The strict inequalities $D_0,D_5>0$ then certify **stability against activating the two inactive coordinates** $x_0,x_5$ (moving into the full positive cone only raises $P$), confirming the attained minimum is genuinely on the $S_2$ face.

### 4.3 $S_1$ — full $p\in(0,1)$ certificate

The $S_1$ KKT system (5 equations in the support ratios) reduces, via the cyclic-ratio / $\beta$-coordinate elimination, to a $p$-free real algebraic curve. The elimination ideal has two components: the **symmetric main branch** $H_B$ (on which $x_2=x_5$, $x_3=x_4$) and an **asymmetric component** $H_C$ (treated at the end of this subsection).

**Completeness of the elimination — no third component.** We record the exact resultant certificate (`code/verify_s1_elimination.py`) that the two components above are *exhaustive*. Use the path-ratio / $\beta$ coordinates $\rho_1,\ldots,\rho_5$ with $\rho_{i+1}=\beta_i/(\rho(1-\beta_i))$, $(u,v,w,z)=\beta_1,\ldots,\beta_4$, $\rho=q/p$. As proved in the closed-form certificate (`code/verify_hc_closedform.py`, §4.3(i)),
$$P=A/\rho_1+B\rho_1+C,\qquad A=(1-u)(1+\rho),\quad B=\tfrac{(1+\rho)uvwz}{\rho^5(1-u)(1-v)(1-w)(1-z)},\quad C=\rho(1+\rho)S,$$
with $A,B,C$ independent of $\rho_1$ ($S=(1-v)(1-u)/u+\cdots+(1-z)/z$). The Euler/KKT system ($\lambda=0$, $\beta\in(0,1)^4$) is $g_1:=\rho_1\partial_{\rho_1}P=\rho_1(B-A/\rho_1^2)=0$ (hence $\rho_1^2=A/B$) and $\partial_{\beta_j}P=0$ for $j=1..4$. Each $L_j:=B\rho_1\,\partial_{\beta_j}C+B\,\partial_{\beta_j}A+A\,\partial_{\beta_j}B$ is *linear* in $\rho_1$, write $L_j=a_j\rho_1+b_j$. Elimination chain (all exact, resultants / linear cross-multiplication, no floating-point arithmetic):
- *(u-elimination)* cross-multiply $L_u,L_v$ $\Rightarrow$ $E_{uv}=\mathrm{(prefactor)}\cdot E_3$ with $E_3=uv^2+uw-u-v^2+v$; solve $E_3=0$ for $u=v(1-v)/((1-w)-v^2)$.
- *(substitute u)* in $E_{uw},E_{uz}$, strip the common prefactor $G=\rho v^3 w z(\rho+1)^4(v-1)^2(v+w-1)$ (admissibility / denominator-zero boundaries) $\Rightarrow$ two $\rho$-free $v$-relations $F_1,F_2$.
- *(v-elimination)* $\operatorname{Res}_v(F_1,F_2)=w^2\,z^2\,(w-1)^2\,H_B\,H_C$ — an **exact identity** (SymPy `expand(R-target)==0`; $H_B$ as below, $H_C$ as in the asymmetric subsection), verified with remainder-zero exact division.

On the admissible set $w,z\in(0,1)$ the prefactors $w^2,z^2,(w-1)^2$ are nonzero ($\beta_3,\beta_4\in(0,1)$), so
$$\pi(V_{\rm KKT}^{\rm adm})\;\subseteq\;V(H_B)\cup V(H_C),$$
i.e. **every admissible $S_1$ KKT point projects to one of the two resultant factors $H_B,H_C$; no third projected component is possible.** (We write containment, not equality: a resultant captures *necessary* conditions — common zeros, leading-coefficient degenerations, and lifts failing the remaining KKT equations may enlarge $V(H_B)\cup V(H_C)$; see the superset remark below. This is exactly what completeness requires.)

**Localization (every stripped factor is nonzero on the admissible set).** The elimination divides out, in succession: $\rho,u,v,w,z$ and $\rho{+}1$ (strictly positive: $\rho>0$, $\beta\in(0,1)^4$, $\rho{+}1>0$); $1{-}v,1{-}w,1{-}z$ (nonzero, $\beta\in(0,1)$); and two non-obvious factors:
- *$v+w-1$.* If $v+w-1=0$ then $w=1-v$ and $v^2+w-1=v^2-v=-v(1-v)$, so $E_3=u(v^2+w-1)+v(1-v)=(-v(1-v))u+v(1-v)=v(1-v)(1-u)=0$ forces $u=1$ (as $v(1-v)>0$), contradicting admissibility $u\in(0,1)$. Hence $v+w-1\ne0$ on admissible KKT points.
- *$(1-w)-v^2$ (the $u$-recovery denominator).* If $1-w-v^2=0$ then $v^2+w-1=0$, so $E_3=u\cdot0+v(1-v)=v(1-v)>0$, impossible for $E_3=0$. Hence the denominator is nonzero.
Thus every factor removed in the chain is rigorously nonzero on $V_{\rm KKT}^{\rm adm}$, justifying the localisation.

*Forward direction is what the theorem uses.* The theorem requires only the **forward** containment proved above — every genuine $S_1$ KKT point projects into $V(H_B)\cup V(H_C)$ — together with bounds on the *whole* positive $H_B$-candidate graph (monotonicity, determinant, crossing portrait below) and on the whole $H_C$-superset (Proposition 4.3). The reverse implication ("$H_B=0\Rightarrow$ a genuine KKT lift exists") is **not** needed and is not claimed: a vanishing specialised resultant can also reflect leading-coefficient degeneration (a projective root at infinity), a complex common root, or a real root outside $(0,1)$, so $H_B\mid\operatorname{Res}_v(F_1,F_2)$ does not by itself exhibit an admissible $v$-root. This is harmless: even if the $H_B$ resultant locus contained extraneous points, analysing the entire positive $H_B$-candidate graph is *stronger* than analysing only the true KKT subset, and the value/determinant/crossing certificates below apply to every true KKT point that does project to $H_B$. For $H_C$, the $\beta$-variety $\{H_C,E_2^u,\overline K,g_1\}$ is a **superset** of the true $H_C$ stationary set (its spurious lifts violate $g_4,g_5$), and this inclusion is what matters, since Proposition 4.3 bounds $P$ on the whole superset — stronger than on the true set. (The reconfirmation of $H_B\mid\operatorname{Res}_v(F_1,F_2)$ by exact univariate polynomial division is recorded in `code/verify_hb_exact_backsub.py`; see the supplement for the divisibility detail.)

On $H_B$ introduce
$$w=\beta_3,\quad z=\beta_4,\qquad H_B(w,z)=zw^2+(1-z^2)w+z^2-z=0,\quad D=1-z^2+wz^2,$$
with the positive $w$-root $w_+(z)=\bigl(-(1-z^2)+\sqrt{\Delta(z)}\bigr)/(2z)$, $\Delta(z)=z^4-4z^3+2z^2+1>0$ on $(0,1)$. Setting $\rho=q/p$ ($p=1/(1+\rho)$), the closure and stationary value are
$$\rho^7(1-z)D^3-wz^5=0,\qquad P_{S_1}^{\rm stat}=\frac{2\rho(1+\rho)}{z}\bigl[3-2z-z^2+wz(1+z)\bigr].\tag{4.1}$$

**Monotonicity of the branch.** The closure determines $\rho$ as a function of $z$ alone: $\rho(z)^7=K(z)=w_+(z)z^5/((1-z)D(z)^3)$. Since $H_w=2zw+1-z^2>0$ on the positive-$w$ branch (the two $H_B$-roots in $w$ have opposite sign, so the positive root has positive derivative), implicit differentiation gives the exact logarithmic derivative
$$\frac{K'(z)}{K(z)}=-\frac{2\,L(w,z)}{w\,z\,D\,H_w},\qquad L(w,z)=\bigl(2z^5{-}3z^4{-}z^3{-}2z^2{+}2\bigr)w-2z^5{+}3z^4{-}z^3{+}5z^2{-}5z,$$
whose denominator is strictly positive on the branch. Thus $K'=0\iff L=0$. Because $L$ is *linear* in $w$, its joint zeros with $H_B$ are governed by a single resultant,
$$\operatorname{Res}_w(H_B,L)=z(z-1)\,Q_7(z),$$
where $Q_7$ is exactly the determinant polynomial below (unique $(0,1)$-root $z_7\approx0.87618$). The unique candidate $L=0$ lift at $z_7$ is $w_L(z_7)=-b(z_7)/a(z_7)$; rational interval arithmetic gives $w_L(z_7)\in[-0.50859,-0.50858]<0$, whereas the branch root is $w_+(z_7)\approx0.24345>0$, so $L\neq0$ on the positive-$w$ branch. A single sign sample fixes the sign: $L(w_+(1/2),1/2)=(-35+5\sqrt{17})/16<0$ (as $\sqrt{17}<7$). Hence $L<0$ throughout, so $K'(z)>0$ and $\rho(z)$ is **strictly increasing**; equivalently $z(p)$ is strictly decreasing and the positive-$w$ branch is a graph over $p\in(0,1)$.

**Determinant sign.** The reduced Hessian determinant on $H_B$ factors as
$$\begin{aligned}
\operatorname{Res}_w(H_B,\mathrm{num\_red})&=4\,z\,(z-1)^9\,Q_5(z)\,Q_7(z),\\
Q_5&=2z^5+2z^3-2z^2-1,\\
Q_7&=8z^7-24z^6+20z^5-9z^4+30z^3-15z^2-6.
\end{aligned}$$
**Sturm:** $Q_5$ has a unique root $z_0\approx0.89756$ in $(0,1)$; $Q_7$ has a unique root $z_7\approx0.87618$ in $(0,1)$. Since $\mathrm{num\_red}=P_v(z)\,w+Q_v(z)$ is *linear* in $w$, the determinant vanishes at a single $w$-value $w_{\rm det}=-Q_v/P_v$. At $z_7$, rational interval arithmetic gives $w_{\rm det}(z_7)\in[-0.866,-0.248]<0$, whereas the branch root is $w_+(z_7)\approx0.243>0$; thus $z_7$ is **not** a determinant zero on the positive-$w$ branch. At $z_0$, $w_{\rm det}=w_+(z_0)>0$ (the genuine transition). Hence on the positive-$w$ branch $\det H_{S_1}=0\iff z=z_0$. Sign samples (rational interval, $D>0$):
$$z=\tfrac{17}{20},\tfrac78\;(<z_0):\ \det>0;\qquad z=\tfrac{9}{10},\tfrac{19}{20}\;(>z_0):\ \det<0.$$
So $\det H_{S_1}<0$ for $z>z_0$ and $\det>0$ for $z<z_0$. **Load-bearing use:** only the $\det<0$ side is rigorous for the theorem — a negative determinant of a $5\times5$ symmetric Hessian forces an *odd* (hence $\ge1$) number of negative eigenvalues, so the stationary point is a **saddle** for $z>z_0$; this is what defers the in-band infimum to the boundary ($\partial S_1=S_2\cup S_3$) below. The $\det>0$ side ($z<z_0$) is *not* used as a Morse certificate: $\det>0$ permits $0,2,$ or $4$ negative eigenvalues and does not by itself prove a local minimum; the out-of-band conclusion for $p>p_0$ instead rests directly on the stationary *value* $P_{S_1}^{\rm stat}>7$ (Proposition 4.2) together with Proposition 4.3 and the boundary, independent of Hessian index (see "Conclusion for $H_B$" below). By monotonicity $z>z_0\iff p<p_0$ and $z<z_0\iff p>p_0$, where $p_0=1/(1+\rho(z_0))=0.388528131361\ldots$ and rational interval arithmetic gives
$$p_0\in(\tfrac38,\tfrac25)\subset(\tfrac13,\tfrac25)\quad\Longrightarrow\quad \rho_0:=\rho(z_0)\in(\tfrac32,2).\tag{4.2}$$
*Numerical remark.* A finite-difference inertia trace (see the supplement) over the branch finds Morse index $1$ for $p<p_0$ (one negative eigenvalue, consistent with the saddle above) and index $0$ for $p>p_0$; this is an independent numerical corroboration, not a certificate the theorem relies on.

**Crossing certificate (region $p>p_0$).** From (4.1), $P_{S_1}^{\rm stat}=7$ is linear in $w$; solving for $w$, substituting into $H_B$ and the closure, and taking $\operatorname{Res}_z$ of the two resulting bivariate polynomials $A(z,\rho)$ ($\deg_z4,\deg_\rho4$), $B(z,\rho)$ ($\deg_z8,\deg_\rho11$) yields
$$\operatorname{Res}_z(A,B)=896\,\rho^{13}(\rho+1)^7(8\rho^2+8\rho+7)^6\,\Phi_{35}(\rho),$$
where the first three factors are strictly positive for $\rho>0$ and $\Phi_{35}$ is the degree-$35$ polynomial (coefficient vector in the supplement). **Sturm:** $\Phi_{35}$ has exactly two positive roots,
$$\rho_1=2.3312465633\ldots\in(2,\tfrac52),\qquad \rho_2=3.2394376235\ldots\in(3,\tfrac72),$$
both $>2>\rho_0$. Hence every $P_{S_1}^{\rm stat}=7$ crossing has $\rho>\rho_0$, i.e. $z>z_0$ (saddle region), and corresponds to
$$p_1=\tfrac1{1+\rho_1}=0.300187927\ldots,\qquad p_2=\tfrac1{1+\rho_2}=0.235880343\ldots,\qquad \{p_1,p_2\}\subset(a_7,b_7).$$
By monotonicity the branch over $p$ is split by $\{p_2,p_1\}$ into three intervals. Since the resultant has no zero except at $\rho_1,\rho_2$, on each open interval the sign of $P_{S_1}^{\rm stat}-7$ is constant; three rigorous interval samples (Krawczyk-unique $z$, native outward-rounded $\sqrt{\cdot}$) fix it:
$$P_{S_1}^{\rm stat}(p)>7\ \text{on }(0,p_2)\cup(p_1,1),\qquad P_{S_1}^{\rm stat}(p)<7\ \text{on }(p_2,p_1)\subset(a_7,b_7),$$
with $P_{S_1}^{\rm stat}(2/5)\in[7.157602,7.157603]>7$ ($\rho=3/2\in(0,\rho_1)$), $P_{S_1}^{\rm stat}(1/4)\in[6.989976,6.989977]<7$ ($\rho=3\in(\rho_1,\rho_2)$), $P_{S_1}^{\rm stat}(1/5)\in[7.059316,7.059317]>7$ ($\rho=4\in(\rho_2,\infty)$). This is formalised as Proposition 4.2 below.

**Proposition 4.2 (out-of-band $S_1$; in particular $p>p_0$).** *On the symmetric branch $H_B$, the only $p\in(0,1)$ at which $P_{S_1}^{\rm stat}=7$ are $p_1=0.300187927\ldots$ and $p_2=0.235880343\ldots$, both lying in $(a_7,b_7)$. Consequently*
$$P_{S_1}^{\rm stat}(p)>7\ \ \text{on}\ \ (0,p_2)\cup(p_1,1)\supset(0,a_7]\cup[b_7,1),\qquad P_{S_1}^{\rm stat}(p)<7\ \ \text{only on}\ \ (p_2,p_1)\subset(a_7,b_7).$$
*In particular, for every $p>p_0$ (where $p_0\in(\tfrac38,\tfrac25)\subset(\tfrac13,1)$ is the unique Hessian-determinant transition, by (4.2)), $P_{S_1}^{\rm stat}(p)>7$.*

*Proof.* The three certificates above are mutually independent and each is rigorous. (a) *Monotonicity:* $\operatorname{Res}_w(H_B,L)=z(z-1)Q_7(z)$ and the sign analysis $w_L(z_7)<0<w_+(z_7)$, $L(w_+(1/2),1/2)<0$ give $K'(z)>0$, so $\rho(z)$ is strictly increasing and the positive-$w$ branch is a single graph over $p\in(0,1)$. (b) *Crossing set:* $\operatorname{Res}_z(A,B)=896\rho^{13}(\rho+1)^7(8\rho^2+8\rho+7)^6\Phi_{35}(\rho)$ with the first three factors strictly positive for $\rho>0$; Sturm gives $\Phi_{35}$ exactly two positive roots $\rho_1\in(2,5/2)$, $\rho_2\in(3,7/2)$, both $>\rho_0$ (since $\rho_0\in(3/2,2)$ by (4.2)), i.e. $p_1,p_2<1/3<b_7$. Thus no crossing has $\rho\le\rho_0$ (equivalently $p\ge p_0$): $P_{S_1}^{\rm stat}\neq7$ on $[p_0,1)$. (c) *Sign of the three components.* Since the resultant has no zero except at $\rho_1,\rho_2$, on each open interval $(0,\rho_1)$, $(\rho_1,\rho_2)$, $(\rho_2,\infty)$ the sign of $P_{S_1}^{\rm stat}-7$ is constant; three rigorous interval samples (Krawczyk-unique $z$ on the closure, native outward-rounded $\sqrt{\cdot}$, `mpmath.iv`; `_verify_s1_three_samples.py`) fix it:
$$\begin{array}{c|c|c}
\text{component} & \text{sample }(\rho,p) & P_{S_1}^{\rm stat}\\\hline
(0,\rho_1)\leftrightarrow(p_1,1) & (3/2,\,2/5) & [7.157602,7.157603]>7\\
(\rho_1,\rho_2)\leftrightarrow(p_2,p_1) & (3,\,1/4) & [6.989976,6.989977]<7\\
(\rho_2,\infty)\leftrightarrow(0,p_2) & (4,\,1/5) & [7.059316,7.059317]>7
\end{array}$$
In particular $[p_0,1)\subset(p_1,1)$ has $P_{S_1}^{\rm stat}>7$. (Corroboration: $\Phi_{35}$ is squarefree, $\gcd(\Phi_{35},\Phi_{35}')=1$, so both crossings are simple — `_verify_phi35_squarefree.py`.) The inclusions $\{p_1,p_2\}\subset(a_7,b_7)$ follow from exact Sturm isolation of all four roots: $a_7\in(0.214,0.215)$, $b_7\in(0.328,0.329)$ (the two $(0,1)$-roots of $F$), $p_2\in(0.235,0.236)$, $p_1\in(0.300,0.301)$ (from $\rho_1\in(2,5/2)$, $\rho_2\in(3,7/2)$); hence $a_7<0.215<0.235<p_2<0.236<0.328<b_7$ and $a_7<0.215<0.300<p_1<0.301<0.328<b_7$ (`n7_s1_crossing_resultant.py`, `_verify_crossings_in_band.py`). $\square$

**Conclusion for $H_B$.** Combining the certificates (and recalling that $\inf_{S_1}P=\min\{$interior local-min values, $\inf_{\partial S_1}P\}$, with interior local minima a subset of the stationary points $H_B\cup H_C$):
- *In-band* $p\in(a_7,b_7)\subset(0,p_0)$: $\det H_{S_1}<0\Rightarrow$ saddle, so the $H_B$ stationary point is *not* a local minimum and does not supply the infimum; $\inf_{S_1}P$ is attained on $\partial S_1=S_2\cup S_3$, where $S_2<7$ (the failure, §4.2) and $S_3>7$ (§4.1). Hence $\inf_{S_1}P=\inf_{S_2}P<7$ in-band.
- *Out-of-band* $p\notin(a_7,b_7)$: every interior stationary value is $>7$ — the $H_B$ value by Proposition 4.2 ($P_{S_1}^{\rm stat}>7$ on $(0,a_7]\cup[b_7,1)$), the $H_C$ value by Proposition 4.3 ($\ge L_C>7$) — and $\inf_{\partial S_1}P=\min\{\inf_{S_2}P,\inf_{S_3}P\}\ge7$ out-of-band (§4.2, §4.1). Therefore $\inf_{S_1}P\ge7$ out-of-band, *regardless* of whether the $H_B$ stationary point is a local minimum or a saddle (its value is $>7$ either way).

**Asymmetric component $H_C$.** $H_C$ is the cubic-in-$w$ factor of the elimination ideal,
$$H_C(w,z)=zw^3+z^3w^2-zw^2+z^4w-3z^3w+2z^2w+zw-w-z^4+3z^3-3z^2+z,$$
appearing as the second factor of the exact completeness identity $\operatorname{Res}_v(F_1,F_2)=w^2z^2(w-1)^2H_BH_C$ recorded above. Its admissible positive lifts (solving $H_C=0$ and the two remaining $p$-free $\beta$-recurrences, then the closure) exist only for $z\in(\zeta,1)$, $\zeta$ the unique $(0,1)$-root of $z^3+2z^2-z-1$, and correspond to $p\in(0.3885,0.3930)\subset(b_7,1)$ (out-of-band). The rigorous bound on $H_C$ is **not** a pointwise Morse classification but a validated cover of the *entire admissible $H_C$-superset* (Proposition 4.3 below): every admissible point has $P\ge L_C=141/20>7$, certified by interval-Newton uniqueness, strict admissibility, and a mean-value $P$-bound, independently re-verified on two interval backends (mpmath.iv and Arb). The would-be second lift branch near $p\approx0.52$ gives KKT residual $>2$ and is not a stationary point (script `n7_s1_hc_kkt_check.py`), so the $\beta$-variety is a genuine *superset* of the true $H_C$ stationary set; bounding $P$ on the superset is stronger than on the true set.

*Numerical remark (representative samples, not a rigorous count).* A finite-difference trace over $z\in\{0.80,\ldots,0.99\}$ (reconstructing $x$ from the 1D free-height equation $g_1=0$ and checking all five KKT residuals $\sim10^{-7}$) finds eight representative real $H_C$ lifts, all with $P>7$ (smallest sampled $P\approx7.1334$ near $z\approx0.90$); these are illustrative samples of the validated arc, not a rigorous enumeration or Hessian-inertia certificate. The theorem relies solely on the superset bound $P\ge L_C>7$.

**Proposition 4.3 (validated exclusion of $H_C$).** *Every admissible point of the reduced $H_C$-superset satisfies $P\ge L_C>7$, where $L_C=\tfrac{141}{20}=7.05$ is an explicit rational. Consequently every true $H_C$ stationary point has value $>7$.*

*Proof.* As noted above, the reduced system $H_C=0\,\wedge\,E_2^u=0\,\wedge\,\overline K=0\,\wedge\,g_1=0$ (closure $K=u v w z^3 a_5^2/((1{-}v)(1{-}w)(1{-}z)^3)$, $u$ eliminated by $E_2^u$) is a *superset* of the true $H_C$ stationary set (its spurious branch violates $g_4,g_5$), so bounding $P$ on the whole admissible superset is *stronger* than on the true set.

**(i) Closed form (exact).** Write the $S_1$ path ratios $\rho_1,\dots,\rho_5$ ($x=(0,1,\rho_1,\rho_1\rho_2,\dots,\rho_1\cdots\rho_5)$) with $\rho_{i+1}=\beta_i/(\rho(1-\beta_i))$ ($\beta_i\in(0,1)$, $\rho=q/p$, $p=1/(1+\rho)$), and set $(u,v,w,z)=(\beta_1,\beta_2,\beta_3,\beta_4)$. Expanding $P=\sum_{i=1}^4\frac1{\rho_i(p+q\rho_{i+1})}+\frac1{p\rho_5}+\frac{\prod\rho_i}{q}$ in these coordinates gives the **rational identity** (verified by exact SymPy expansion, numerator identically zero — `verify_hc_closedform.py`)
$$P \;=\; \frac{A}{\rho_1} + B\,\rho_1 + C,$$
with
$$A=(1{-}u)(1{+}\rho),\quad B=\tfrac{(1+\rho)uvwz}{\rho^5(1{-}u)(1{-}v)(1{-}w)(1{-}z)},\quad C=\rho(1{+}\rho)S,\quad S=\tfrac{(1-v)(1-u)}{u}+\tfrac{(1-w)(1-v)}{v}+\tfrac{(1-z)(1-w)}{w}+\tfrac{1-z}{z},$$
where $A,B,C$ are **independent of $\rho_1$**. The first KKT equation is $g_1:=\rho_1\,\partial P/\partial\rho_1=\rho_1(B-A/\rho_1^2)$, so on the admissible set ($\rho_1>0$)
$$g_1=0 \iff \rho_1^2=A/B \iff \rho_1=+\sqrt{A/B}\quad(A,B>0\text{ on the admissible set}).$$
Substituting $\rho_1=\sqrt{A/B}$: $\;A/\rho_1+B\rho_1=A\sqrt{B/A}+B\sqrt{A/B}=\sqrt{AB}+\sqrt{AB}$, whence
$$\boxed{\,P \;=\; C + 2\sqrt{AB}\,}$$
as an **exact algebraic identity** on the $H_C$ stationary set — no numerical approximation is involved (the squared check $(A/\rho_1+B\rho_1)^2\big|_{\rho_1^2=A/B}=4AB$ is also verified identically). Numerically $C\in(3.7,5.2)<7$ throughout and $2\sqrt{AB}\in(2.95,3.47)$. On the admissible set all denominators and radicands are strictly positive, so $P$ is $C^1$ there.

**(ii) Desingularisation (both endpoints).** The admissible arc terminates at $(w,z)=(0,1)$, where $H_C(w,1)=w^3$ and $\partial_w H_C,\partial_z H_C\to0$ (triple root). With $s=1{-}z,\ c=w/s$ the factor divides out: $H_C(cs,1{-}s)=s^3\,G(c,s)$, $G=(1{-}s)c^3+(1{-}s)(s{-}2)c^2+(s^2{-}s{-}1)c+(1{-}s)$, **regular** at $s=0$ ($G(c,0)=c^3-2c^2-c+1$ has the arc root $c_0\approx2.24698$, $\partial_cG(c_0,0)\approx5.15\ne0$); the closure and $E_2^u$ lose their $s$-power degeneracy and give regular reduced equations $G(c,s)=0$, $E_2^{\rm red}(v,c,s)=0$. The arc also approaches the opposite endpoint $s=1$ ($z=0$: $c=w/s\to0$, $v\to0$, $u\to0$, $K\to0$ — an inadmissible point, but the admissible arc tends to it). There the raw $c,v,w,u,K,\rho_7$ all vanish like powers of $\delta:=1{-}s=z$, which would defeat interval arithmetic (underflow to $0$); a **second blow-up** $\delta=1{-}s,\ c=\delta\bar c,\ v=\delta\bar v$ desingularises it. Setting $\bar G(\bar c,\delta)=G(\delta\bar c,1{-}\delta)/\delta$ and $\bar E_2(\bar v,\bar c,\delta)=E_2^{\rm red}(\delta\bar v,\delta\bar c,1{-}\delta)/\delta$ (both polynomials, denominator $1$) gives quantities **regular at $\delta=0$**: $\bar G(\bar c,0)=1{-}\bar c$ (simple root $\bar c=1$, $\partial_{\bar c}\bar G|_{0}=-1\ne0$) and $\bar E_2(\bar v,0)=\bar v{-}1$ (simple root $\bar v=1$). Admissibility is certified in **rescaled** variables that stay $O(1)$ and bounded away from $0$ as $\delta\to0^+$ — $\bar c>0,\ \bar v>0,\ w=\bar c(1{-}\delta)>0,\ \hat u=u/\delta>0,\ \widetilde K=K/\delta^6>0,\ \tilde\rho_7=\rho_7/\delta^6>0,\ a_5>0$ — each raw quantity being $\delta^k$ times a positive regular factor, so rescaled admissibility for $\delta>0$ is *equivalent* to raw admissibility (no underflow, no $1/\delta$ singularity). The $s\to1$ tail is thus a regular desingularised region, not a gap.

**(iii) Validated cover (exhaustiveness).** The admissible-lift set can change only at the root-multiplicity and admissibility-boundary events $E_1$--$E_9$ generated in `n7_s1_hc_critical_events.py`. Their roots in $(0,1)$ give five exact algebraic critical values and hence six event-free $s$-cells. The load-bearing exhaustiveness checker is `n7_s1_hc_exhaustiveness_cert.py`; it does not use `nroots`, `polyroots`, `nsimplify`, residual tolerances, or a count-only midpoint argument.

The checker combines the $2604$ original boxes with the $60$ second-blow-up collar boxes and constructs a **root-branch graph**. Two boxes are joined only after a full $2\times2$ Krawczyk operator, evaluated at an exact rational parameter in their overlap, proves that their Krawczyk-unique roots are the same root. Raw $(c,v)$ coordinates are used for original-cover edges and at the cover--collar seam; $(\bar c,\bar v)$ coordinates are used inside the collar. All $2823$ candidate overlaps are certified. The graph has exactly two connected components: the principal branch, represented by $1435$ original boxes and all $60$ collar boxes, covers $s\in[0,1]$; the secondary branch, represented by $1169$ original boxes, covers $s\in[0,s_\dagger]$, where $s_\dagger=55749571201593/281474976710656\approx0.1980622642$. For each component the union of its exact rational $s$-intervals is checked to be gap-free. This is a branch-level statement, not merely coverage of the total $s$-projection.

At a certified rational sample in each of the six event-free cells, all roots are enumerated algebraically as follows. The cubic $G(c,s_0)$ is solved by exact `CRootOf`; for each isolated $c$-root, the factor $Q(v,c,s)=E_2^{\rm red}(v,c,s)/(v-1)$ is quadratic in $v$, so all of its real roots are enclosed by an outward-rounded interval quadratic formula. Strict admissibility is then decided by outward interval signs. The admissible roots are in bijection with the branch-graph components active throughout that cell: two roots/components in the first cell and one in each of the other five.

The five critical fibers are treated without floating root recovery or number-field factorisation. Each critical value $s_\ast$ is retained as an exact `CRootOf` and enclosed in a certified rational isolating interval containing no other event. After removing the endpoint factor $(s-1)^2$, the exact resultant $\operatorname{Res}_c(G,\partial_cG)$ has no root in $(0,1)$; hence full-parameter interval Newton produces disjoint tubes containing every real root of $G(c,s_\ast)=0$. For each tube, the factor $Q(v,c,s)=E_2^{\rm red}/(v-1)$ is quadratic in $v$, so all real $v$-roots are enclosed by the outward-rounded quadratic formula. Strictly admissible roots are interval-classified directly. Boundary tubes are excluded by certified equations: generic boundary events use a $3\times3$ Krawczyk operator for $(G,Q,B)=0$, while the terminating $u=1$ branch is identified exactly by the last nonzero linear subresultant of $H_C(w,z)$ and $E_2^u(1-w,w,z)$; its coefficient is interval-certified nonzero, giving the exact algebraic root $w=-B(z)/A(z)$ and therefore $v=1-w$. Thus every critical-fiber admissible root is enclosed by exactly one branch component, while the secondary branch terminates at the first critical fiber on the boundary $u=1$ and is not admissible there. The machine record `code/_hc_exhaustiveness.json` reports `n_components=2`, `all_cells_branch_bijection=True`, and `all_critical_admissible_roots_boxed=True`. Since no multiplicity or admissibility-boundary event occurs inside a cell, the certified graph and the cell/fiber enumeration exhaust every admissible lift on $(0,1)$.

**(iv) Validated cover (bounding).** Each admissible lift is continued over its $s$-cell by interval-Newton (Krawczyk) contraction on $G$ (in $c$) and $E_2^{\rm red}$ (in $v$) with factored derivative evaluators: $N(W)\subset\operatorname{int}W$ certifies a *unique* root in the box. The stationary value is bounded over the 4-box $(v,c,s,\rho)$ by the **mean-value form** $P(\mathrm{box})\subseteq P(\mathrm{mid})+\sum_i(\partial P/\partial x_i)(\mathrm{box})\cdot(\mathrm{box}_i{-}\mathrm{mid}_i)$, the four partials being closed-form factored (no symbolic simplification, avoiding interval blow-up from fractional powers with sign-changing bases); inputs are outward-rounded. The $s$-partition is adaptively bisected whenever the lower bound falls below $L_C=141/20$.

**(v) Result.** The branch graph has two admissible components on the first event-free cell and one thereafter. At the terminating fiber $s_\dagger\approx0.1980622642$, the secondary component reaches the exact boundary $u=1$ and is excluded; the principal component remains strictly admissible. At each of the five critical fibers, every strictly admissible root is enclosed by exactly one certified component box. The endpoint $s=0$ ($z=1$, degenerate limit $w\to0$, $K\to0$, not admissible) is **covered, not merely approached**: the desingularisation $G(c,s)$ is regular at $s=0$ ($G(c,0)=c^3-2c^2-c+1$, simple arc root $c_0\approx2.24698$), so the cover includes two pieces with closed lower endpoint $s_{\rm lo}=0$ (namely $s\in[0,\,9.67\!\times\!10^{-5}]$, real-arc and a second branch), certified with $P_{\rm MV}\ge7.1139$ and $\ge8.2036$ respectively — both $>7$, so $(0,10^{-7})\subset(0,9.67\!\times\!10^{-5})$ is rigorously covered (real-arc limit $P\to7.146684>7$). The opposite endpoint $s=1$ ($z=0$: $c=v=u=0$, $K=0$, an inadmissible point) is likewise covered, not gapped. The original $2604$ boxes reach only $s_{\max}=8840966100154093/9007199254740992\approx0.9815444124$. The remaining interval $[s_{\max},1]$ is covered by `n7_s1_hc_s1collar.py`, which performs the second blow-up of §4.3(ii). Its $60$ exact-rational pieces are certified by **parametric interval Newton on the full parameter interval** $D=[\delta_a,\delta_b]$ of every piece, not at the midpoint alone:
$$N_{\bar c}(C;D)\subset\operatorname{int}C,\qquad N_{\bar v}(V;C,D)\subset\operatorname{int}V.$$
The same full-$D$ inclusions, rescaled admissibility inequalities, and $P$ bounds are independently re-evaluated from the recorded rational boxes. The final closed tail $D=[0,10^{-10}]$ is also handled by parametric Newton at the regular endpoint $(\bar c,\bar v,\delta)=(1,1,0)$; admissibility is restricted to $\delta>0$, while the excluded point $\delta=0$ has $K=0$. The pieces abut exactly, with first endpoint $s_{\max}$ and final endpoint $1$, and the rechecked lower bound is
$$\min_j P_{{\rm collar},j}^{\rm lo}=7.3440544147\ldots>141/20.$$
In particular, the previously uncovered admissible point at $s=491/500$ is inside a parametric-Newton-certified collar piece. The canonical cover ($2604$ rational lift-records, produced by `n7_s1_hc_arb_refine.py`) is certified by the **formal validated-numerics backend** `n7_s1_hc_arb_checker.py` (Arb / `python-flint`, prec $150$): it re-reads the dumped boxes and re-derives, per piece, Krawczyk uniqueness, strict admissibility, and the mean-value bound in a hand-written lo/hi interval layer (monotone-endpoint enclosures for $\sqrt{\cdot}$ and $\rho^{1/7}$, every operation rigorously outward-rounded, dependency reduced by factored evaluators), verifying $P_{\rm MV}>L_C$ for **all $2604/2604$** pieces; the minimum of these certified piecewise lower bounds is $\min_j P_{\rm MV}^{(j)}=7.050005023\ldots>141/20$. (This is a *rigorous lower bound on $P$* over the superset, obtained as the minimum of the per-piece mean-value lower bounds; we make no claim that it equals the exact infimum of $P$ on the $H_C$-superset — only that it exceeds $141/20>7$, which is all the proof requires.) As an independent defence-in-depth, the cover is **re-verified by a second interval backend** — `n7_s1_hc_cover_checker.py` (`mpmath.iv`, native outward-rounded `sqrt`) — which independently re-derives the same Krawczyk inclusion, admissibility, and mean-value bound and confirms $P_{\rm MV}>L_C$ for **all $2604/2604$** pieces. (The few tight real-arc pieces, whose mean-value lower bound lay within $0.002$ of $L_C$ and were widened below it by interval dependency, were resolved by $s$-bisection — `n7_s1_hc_arb_refine.py` — not by any relaxation of rigor; both backends are outward-directed.) Hence
$$P\ge L_C=\tfrac{141}{20}=7.05>7\qquad\text{on the entire admissible }H_C\text{-superset},$$
and in particular on every true $H_C$ stationary point. (Full per-piece data and checkers are in the supplement.) $\square$

Thus every positive $S_1$ stationary point either has $P\ge7$, or lies in $(a_7,b_7)$ with a non-PSD Hessian: on $H_B$ by Proposition 4.2 (out-of-band $P_{S_1}^{\rm stat}>7$ for $p>p_0$; saddle for $p<p_0$, deferring to the boundary), on $H_C$ by Proposition 4.3 ($P\ge L_C>7$). Consequently $\inf_{S_1}P\ge7$ for $p\notin(a_7,b_7)$ and $=\inf_{S_2}P<7$ for $p\in(a_7,b_7)$.

### 4.4 $S_0$ — Nowosad–Yamagami
The uniform point is handled by Lemma 2.2: the Hessian eigenvalues $\mu_k=2(1-\cos\theta_k)[2q\cos\theta_k+2p^2-3p+2]$ ($\theta_k=2\pi k/7$) are strictly positive because $\Delta_7<0$ (§4.1), and $S$ is invertible (no Fourier eigenvalue $\lambda_k=p\omega^k+q\omega^{2k}$ vanishes for odd $n$). Hence $\mathbf1$ is the **unique interior local-minimum ray** on $S_0$ for every $p\in(0,1)$, with value $7$. The face infimum is *not* $7$ pointwise: by boundary inheritance (§4.5) $\inf_{S_0}P=\min\{7,\inf_{S_1}P,\inf_{S_2}P,\inf_{S_3}P,\inf_{S_4}P\}$, so in the band $(a_7,b_7)$ where $\inf_{S_2}P<7$, one has $\inf_{S_0}P<7$ — the failure is inherited from the $S_2$ boundary stratum, not contradicted by anything at the uniform point.

### 4.5 Conclusion of Theorem 2
Assembling the five face certificates for every $p\in(0,1)$ by **boundary inheritance** ($\inf_{S_0}P=\min\{7,\inf_{S_1}P,\inf_{S_2}P,\inf_{S_3}P,\inf_{S_4}P\}$, and each $\inf_{S_k}$ is itself either an interior stationary value or inherited from $\partial S_k$): $S_0$ contributes the attained uniform value $7$ (unique interior local min, Lemma 2.2, §4.4); $S_4$ gives $P\ge\tfrac{16}{3^{3/4}}>7$ (§4.1); $S_3$ gives $P>8$ ($p>1/2$) or defers to $S_4$ ($p\le1/2$) (§4.1); $S_1$ gives $\inf_{S_1}P\ge7$ out-of-band and $=\inf_{S_2}P$ in-band (Propositions 4.2–4.3, §4.3); and $S_2$ gives $P_{S_2}^{\rm curve}>7$ out-of-band, $<7$ in-band (§4.2). Since $S_4$ is terminal and every boundary recursion lands in one of these faces, the only face that can drop below $7$ is $S_2$, hence
$$m_7(p)=\min\{7,P_{S_2}^{\rm curve}(p)\},\qquad H_7=(0,a_7]\cup[b_7,1).$$
The degree-$21$ polynomial
$$\begin{aligned}G(p)=&\,13p^{21}-238p^{20}+1841p^{19}-8295p^{18}+25690p^{17}-60200p^{16}\\&+114261p^{15}-184911p^{14}+263550p^{13}-333970p^{12}+372736p^{11}\\&-359996p^{10}+295750p^9-203770p^8+116300p^7-54264p^6\\&+20349p^5-5985p^4+1330p^3-210p^2+21p-1\end{aligned}$$
certifies the inactive KKT (Sturm: $0$ roots in $(1/5,1/3)\supset(a_7,b_7)$). $\blacksquare$

*Remark (sharpness vs. $n=5$).* The Tuan–Thuong sufficient interval $I_{\rm TT}=[(5-\sqrt5)/10,(5+\sqrt5)/10]\approx[0.2764,0.7236]$ (a strict subset of $H_5=(0,1)$, Theorem 1) does **not** extend to $n=7$: since $a_7\approx0.2143<(5-\sqrt5)/10\approx0.2764<b_7\approx0.3286$, the whole interval $\bigl((5-\sqrt5)/10,\,b_7\bigr)$ lies in $H_5=(0,1)$ but inside the $n=7$ failure band $(a_7,b_7)$. The transition from "$n=5$ safe everywhere" to "$n=7$ first failure band" is boundary-induced, exactly as the gap-concentration picture predicts (the $L=5$ one-gap face $S_2$ first drops below $7$).

## 5. Computer-assisted proof and code availability

The proofs are computer-assisted in the following sense. Resultant eliminations, exact rational Sturm sign-variation counts [2], positive-coefficient crossing factors, Hessian-determinant factors, and rational interval sign samples are all produced by exact symbolic computation, and every asserted identity is checked by exact polynomial division (quotient and remainder recorded, remainder zero). Inequality and root-isolation certificates use rigorous interval arithmetic with outward rounding [10]. Two independent arithmetic backends are used, both outward-directed: exact rational arithmetic (SymPy) and interval arithmetic (`mpmath.iv` and Arb via `python-flint`); the $H_C$ cover is verified by both backends independently. The $S_1$/$S_2$ crossing and determinant certificates are univariate (Sylvester determinants of size $\le 12$); the $H_C$ critical-event enumeration forms bivariate resultants, but each factor has degree $\le 3$ in the eliminated variable, so every intermediate has bounded degree and the scripts run in bounded memory.

The complete software environment, the $\Phi_{35}$ coefficient vector and its Sturm isolation counts, the per-file certificate script listings, the machine records (`_hc_cover.json`, `_hc_s1collar.json`, `_hc_exhaustiveness.json`, `_archival_run.txt`), and the reproducibility commands are recorded in the supplementary material (`SUPPLEMENT.md`); the authoritative script role map is in `code/README.md`, and pinned dependencies in `code/requirements_lock.txt`. The proof artifact is committed to a self-contained git repository with a SHA-256 manifest (`sha256_manifest.txt`); a Zenodo DOI will be added upon deposition and is not a proof ingredient.

*AI assistance disclosure.* Some of the computer-assisted verification scripts were developed with the assistance of a large language model (ChatGPT) and were adopted only after the author's independent audit and exact recomputation; all mathematical content, theorems, and proofs are the author's own.

## References

[1] T. Ando, *A new proof of Shapiro inequality*, Math. Inequal. Appl. **16** (2013), 611–632.

[2] S. Basu, R. Pollack, and M.-F. Roy, *Algorithms in Real Algebraic Geometry*, Algorithms Comput. Math. **10**, Springer-Verlag, Berlin, 2006.

[3] V. G. Drinfel'd, *A cyclic inequality*, Math. Notes Acad. Sci. USSR **9** (1971), 68–71.

[4] D. Dummit and R. Foote, *Abstract Algebra*, 3rd ed., John Wiley & Sons, 2004.

[5] M. J. Lighthill, *An invalid inequality*, Amer. Math. Monthly **63** (1956), 191–192.

[6] D. S. Mitrinović, J. E. Pečarić, and A. M. Fink, *Classical and New Inequalities in Analysis*, Kluwer Academic Publishers, Dordrecht, 1993.

[7] P. Nowosad, *Isoperimetric eigenvalue problems in algebras*, Comm. Pure Appl. Math. **21** (1968), 401–465.

[8] H. S. Shapiro, Problem 4603, Amer. Math. Monthly **61** (1954), 571.

[9] N. M. Tuan and L. Q. Thuong, *On an extension of Shapiro's cyclic inequality*, J. Inequal. Appl. **2009**, Article ID 491576.

[10] W. Tucker, *Validated Numerics: A Short Introduction to Rigorous Computations*, Princeton University Press, Princeton, 2011.

[11] S. Yamagami, *Cyclic inequalities*, Proc. Amer. Math. Soc. **118** (1993), 521–527.
