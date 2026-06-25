# Group-algebra (Kani–Rosen / Paulhus–Rojas) attack on genus 56

Magma scripts applying the **non-modular group-action method** to the search for a
genus-56 curve with completely decomposable Jacobian. This complements the
modular / Atkin–Lehner quotient work in the repo root (`harvest*.m`, `SUMMARY.md`).

Full write-up: [`../group-algebra/Summary-third-run-claude.md`](../group-algebra/Summary-third-run-claude.md).

**Result: no genus-56 example found by any of these searches, but every engine
reproduces published Ekedahl–Serre / Paulhus–Rojas examples** (A₄→g4, PGL(2,7)→g36,
`(1344,11289)`→g57, ES genus-47 = `X₀(600)/w₂₄`, PR's genus-54/58/67/69 quotients).
The systematic constructions *straddle* 56 (hit 54,55,57,58,59,61,… never 56).

Run scripts with `magma -b <file>` from this directory.

## Files

| file | purpose |
|---|---|
| `engine.m` | Reusable engine: group-algebra factor-dimension formula `dim Bᵢ = ½⟨ψᵢ,χ_V⟩` (PR eq. 4) + Frobenius structure-constant test for existence of a generating vector. `load`ed by the drivers. |
| `galg.m` | Standalone validation of the dimension formula (A₄ on genus 4 → `E×E′³`; PGL(2,7) on genus 36 → 5 elliptic factors). |
| `gsearch.m` | Engine with built-in positive controls (PGL(2,7) g36, `SmallGroup(1344,11289)` g57). |
| `gdrive.m` | **§2.1 direct sweep** — all candidate group orders ≤ 1980, `g₀∈{0,1}`, `r≤4`, incl. the uncovered `\|G\|≤220` regime. Result: 0 hits. (loads `engine.m`) |
| `gquot.m` | **§2.2 intermediate covers** from completely-decomposable sources of genus 57–91 (too small — only trivial quotient ≥56). (loads `engine.m`) |
| `gquot2.m` | **§2.2** from large genus-109 / 145 `SmallGroup` sources. Reproduces PR's 54, 67, 69. No 56. (loads `engine.m`) |
| `region1b.m` | **§2.2** reconstructs genus-121 (order 2880) and genus-163 (order 2592) sources via `LowIndexNormalSubgroups`, scans order-2/3 quotients. No 56. |
| `region1c.m` | **§2.2** genus-193 (order 5760) and genus-244 (order 11664) sources via Paulhus–Rojas's explicit presentations, scans small quotients. Reproduces PR's 55, 58. No 56. |
| `quot.m` | **Modular**: Atkin–Lehner quotient genera of every completely-decomposable `J₀(N)` (86 good levels). Reproduces `X₀(600)/w₂₄=47`. No 56. |
| `search56.m` | **Modular**: partial-decomposability AL-quotient search over all `N` with `g(X₀(N))≥56` (invariant part all-elliptic test). No 56 up to N≈1100. |
| `orders56.txt` | Candidate group orders for a genus-56, `g₀=0` action (Riemann–Hurwitz). |
| `goodlevels.txt` | The 86 levels with `J₀(N)` completely decomposable (max N=1200), reconstructed from LMFDB. |
