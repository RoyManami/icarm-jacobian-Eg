# Hunting a genus-56 curve with completely decomposable Jacobian

*Research note — state of the problem, work done, and future approaches.*

## 1. The problem

A principally polarized abelian variety over **C** is **completely decomposable** if it is
isogenous to a product of elliptic curves. Ekedahl and Serre (1993) asked:

> **Q1.** For every genus *g*, is there a curve of genus *g* whose Jacobian is completely decomposable?
> **Q2.** Are the genera of such curves bounded?

Both are still open in characteristic 0. Ekedahl–Serre exhibited examples up to genus 1297 but
with many gaps. Over the decades the gaps have shrunk (Yamauchi 2007; Paulhus–Rojas 2017;
Rodríguez–Rojas 2024; Paulhus–Sutherland 2025).

**The current frontier: genus 56 is the smallest genus for which *no* completely decomposable
Jacobian is known** (over any number field, or even over **C**). Paulhus–Sutherland (2025) resolved
genus 38, which pushed the smallest open case up to 56.

Known genera in char. 0 (Paulhus–Sutherland 2025):
`0–55, 57–58, 61–65, 67–69, 71–73, 75–82, 85, 89, 91, 93, 95, 97, 101, 103, 105–107, 109, 113,
118, 121, 125, 129, 135, 137, 142, 145, 154, 157, 159, 161, 163, 169, 193, 199, 205, 211, 213,
217, 244, 257, 325, 409, 433, 649, 1297.`
→ Smallest gap = **56**; then 59, 60, 66, 70, 74, ...

### Key references
- Ekedahl–Serre, *Exemples de courbes … à jacobienne complètement décomposable*, C.R. Acad. Sci. 1993.
- Yamauchi, *On Q-simple factors of Jacobian varieties of modular curves*, 2007 (classifies N with J₀(N) completely decomposable).
- Paulhus–Rojas, *Completely decomposable Jacobian varieties in new genera*, Exp. Math. 2017 (arXiv:1603.00331) — group-action method.
- Rodríguez–Rojas 2024 — genus 101 via simple-factor decomposition algorithm.
- **Paulhus–Sutherland, *Completely decomposable modular Jacobians*, 2025 (arXiv:2502.16007)** — modular X_H method; states 56 is the smallest open genus.

### The two systematic methods
- **Modular.** J₀(N), and more generally Jac(X_H) for X_H a modular curve, decompose into modular
  abelian varieties A_f. Completely decomposable ⇔ every A_f is an elliptic curve (weight-2 rational
  newform). PS2025 searched all X_H of level **N < 240** exhaustively.
- **Group actions (Kani–Rosen / Paulhus–Rojas).** For a curve X with finite group G acting, the
  group algebra **Q**[G] induces an isogeny `JX ~ B₁^{n₁} × … × B_r^{n_r}`, one factor per rational
  irreducible character. Completely decomposable (via this method) ⇔ every `dim Bᵢ ∈ {0,1}`, computed
  by `dim Bᵢ = ½⟨ψᵢ, χ_V⟩` (Paulhus–Rojas eq. 4). Two sub-methods:
  - **§2.1 direct:** find a single curve whose own decomposition is all-elliptic.
  - **§2.2 intermediate covers:** if X is completely decomposable, *every* quotient X/H is too — so a
    subgroup quotient of genus 56 would settle it.

## 2. Work done (this investigation)

Tools: **Magma** (`/Applications/Magma/magma`) — the same tool the papers use — plus the **LMFDB**
database. Every search engine was **validated against published examples before use**.

### A. Modular Atkin–Lehner quotients
- Independently reconstructed from LMFDB the **86 completely-decomposable levels** N (J₀(N) a product
  of elliptic curves), max **N = 1200** — matches Yamauchi.
- Reproduced the Ekedahl–Serre genus-47 example `X₀(600)/w₂₄` ✓.
- Computed all Atkin–Lehner quotient genera of every decomposable J₀(N): reachable set =
  `{2..55 (most), 57, 61, 73, 77, 81, 97, 103, 121, 161, 205}`. **56 is not reachable.**
- Broader *partial-decomposability* quotient search over **all N ≤ 1100** (any Atkin–Lehner subgroup,
  testing whether the invariant part is all-elliptic): **no genus-56 completely-decomposable quotient.**

### B. Group-algebra method — implemented & validated
Reproduced exactly: A₄ on genus 4 → `E×E′³`; PGL(2,7) on genus 36 → 5 elliptic factors `[0;2,6,8]`;
`SmallGroup(1344,11289)` on genus 57 → 6 elliptic factors `[0;2,4,6]`.

- **§2.1 direct sweep** — all candidate group orders ≤ 1980, both `g₀=0` and `g₀=1`, branch points
  `r ≤ 4`, **including the uncovered `|G| ≤ 220` regime** (no complete classification exists for genus
  56, since Breuer stops at genus 48): **0 hits.**

- **§2.2 intermediate covers** — reconstructed the large completely-decomposable source curves and
  scanned every subgroup quotient that could reach genus 56:

  | source (order, signature) | genus | quotient genera near 56 | 56? |
  |---|---|---|---|
  | 2880, `[0;2,4,6]` | 121 | `\|H\|=2`: 49, **55,57,58**, 59, 61 | ✗ |
  | 2592, `[0;2,4,8]` | 163 | `\|H\|=2`: 73, 79–82 · `\|H\|=3`: **55** | ✗ |
  | 5760, `[0;2,3,10]` | 193 | `\|H\|=2`: 89, 93, 97 · `\|H\|=3`: 63, 65 | ✗ |
  | 11664, `[0;2,3,8]` | 244 | `\|H\|=3`: 79, 82 · `\|H\|=4`: **55, 58** | ✗ |
  | 1296 (×3) | 109 | 52, 53, 54, 55 | ✗ |
  | 1728 (×3) | 145 | 55, 65, 67, 69, 71 | ✗ |

  Correctness anchored by reproducing PR's *published* quotients: genus **54** (genus-109 source),
  **58** (genus-244 source), 67/69 (genus-145 sources).

### Central finding: **the "straddling" phenomenon**
Every systematic construction lands on 54, 55, 57, 58 (and 59, 61, 63, 65…) but **jumps over 56**.
The genus-121 source is the sharpest: its order-2 quotients give 55, 57, 58, 59, 61 — bracketing 56
on both sides — yet never 56. This is *why* 56 is the smallest open genus: the elliptic-factor
multiplicities that must sum to exactly 56 under Riemann–Hurwitz simply don't occur for the available
`(group, signature)` data. **All searches negative — but all tooling reproduces known results.**

## 3. Future approaches (prioritized)

**Group-theoretic (extends what's done):**
1. **Genus-433 source** (order 5184, `[0;2,6,6]`) — the richest untested source; needs subgroups of
   order ≤ 7. Blocked only by `LowIndexNormalSubgroups` being slow at index 5184 and PR publishing no
   presentation. Reconstruct via targeted coset enumeration / RWS, or pull Conder's group data file.
   *Best single remaining shot.*
2. Other large sources: genus **257** (order 12288), **325** (15552), **129** (10752).
3. **Genuinely uncovered region:** group actions on genus 56 with small `|G| ≤ 220` and/or `g₀ ≥ 1`
   and **`r ≥ 5` branch points** — not covered by the `r ≤ 4` sweep and with no existing
   classification. Combinatorially large but searchable in slices.

**Modular (extends Paulhus–Sutherland):**
4. **The explicit Remark-2 program PS2025 did *not* run:** take their 100,000+ completely-decomposable
   X_H of genus > 100 and quotient by *extra* automorphisms (beyond Atkin–Lehner) to look for genus 56.
5. Modular curves X_H at **level 240 ≤ N ≤ 707** — PS2025 ran the full algorithm only for N < 240.
6. **Nearly-decomposable** modular Jacobians (PS2025 Remark 2): J_H with one small non-elliptic factor
   killed by a quotient.

**Other constructions (untouched here):**
7. Generalized Fermat curves and their isogenous decompositions (Hidalgo et al.).
8. Shimura curves / fake elliptic curves; fiber products of elliptic surfaces.
9. Prym / Recillas–Rodríguez intermediate-cover algorithm (route by which Rodríguez–Rojas reached
   genus 101) applied with genus-56 targets.

## 4. Reusable assets produced
Validated Magma engines (group-algebra dimension formula + structure-constant generating-vector test;
Atkin–Lehner quotient genus + all-elliptic rationality checker; source reconstruction via
`LowIndexNormalSubgroups` and via explicit presentations). LMFDB query reconstructing the 86
completely-decomposable levels. All cross-checked against published Ekedahl–Serre / Paulhus–Rojas /
Paulhus–Sutherland results.
