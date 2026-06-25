# Investigation summary: hunting new Ekedahl–Serre genera by quotienting completely decomposable modular Jacobians

Goal (from Remark 2 of Paulhus–Sutherland, *Completely decomposable modular
Jacobians*, arXiv:2502.16007): take quotients of the >100,000 completely
decomposable modular curves they found (and of related modular curves) and look
for a quotient whose genus is one of the **missing** ("gap") genera — a genus
for which no completely decomposable Jacobian is currently known. The smallest
gap is **56**; gaps run throughout `[56, 1297]`.

All work is over the data/code of
`github.com/AndrewVSutherland/CompletelyDecomposableModularJacobians`
(`grpdata.txt` = 561,077 completely decomposable `X_H` of level `< 240`) using
Magma 2.29-8. "Realized genus" = the set encoded in `KnownRealizedGenera()` in
`harvest.m`; the negative results below are relative to that set.

---

## 0. The bug that blocked everything (`harvest.m` would not run)

`harvest.m` Part A crashed with `Undefined reference 'atoi'` inside the GL2
package. Root cause: the `Magma` **submodule** of the CDMJ repo was pinned at an
old commit (`eed5c96`) whose `genus1.m` had `end intrinsic` **missing a
semicolon**, so `utils.m`/`genus1.m` failed to attach and `atoi` (hence
`GL2FromGenerators`, all of Part A) was undefined.

**Fix:** update the submodule to the upstream commit that fixes it:
```bash
git submodule update --init --depth 1
git -C Magma fetch --depth 50 origin main && git -C Magma checkout 71a8fdd  # "fixing semicolons"
```
After this `AttachSpec("spec")` is clean and `CheckAnchor()` prints
`g(X_0(600)/w_24) = 47`.

---

## 1. Modular (`GL2(Ẑ)`-internal) quotients of census curves — `harvest.m` Part A

`ModularQuotientGenera(H)` computes `GL2Genus(B)` for every intermediate group
`H ≤ B ≤ N(H)` (normalizer **inside** `GL2(ℤ/Nℤ)`): diamond operators and
`GL2`-internal automorphisms.

**Rigorous result — provably no new genus.** If `H ≤ B` then `B` is defined mod
the same `N < 240`, has surjective determinant, contains `−I`, and
`Jac(X_B)` is an isogeny factor of `Jac(X_H)` (Poincaré) hence completely
decomposable. Those are *exactly* the census membership conditions, so **every
such `B` is already in `grpdata.txt`**. Verified directly: the set of distinct
census genus values minus the realized set is **empty**.

**Performance note.** The shipped enumeration `Subgroups(Nz : OrderMultipleOf)`
times out at level 48/60 (ambient ~10⁶). Replacing it with "quotient
`Q = N(H)/H`, enumerate `Subgroups(Q)`, pull back via `@@ pi`" makes it
millisecond-fast (used in `harvestD.m`).

---

## 2. Atkin–Lehner quotients of `X_0(N)`, good levels — `harvest.m` Part B & `harvestB.m`

The crucial limitation of §1: the `GL2(ℤ/N)`-normalizer **misses Atkin–Lehner /
Fricke involutions** (their matrices have non-unit determinant mod `N`, so they
are *not* in `GL2(Ẑ)`). E.g. for `X_0(30)` the normalizer gives `|N(H)/H| = 1` —
Part A sees **zero** quotients, though `X_0(30)` has AL group `(ℤ/2)³`. AL
quotients escape the census and *can* give new genera (genus 47 = `X_0(600)/w_24`
is one).

- **`harvest.m` `HarvestX0AtkinLehner`** — single involutions `w_d`, 86 good
  levels `N ≤ 1200`: **0 new genera**.
- **`harvestB.m`** — *every* subgroup of the full AL group `(ℤ/2)^k`, 86 good
  levels: **0 new gaps**. Proper-quotient genera `≥56` = `{57,61,77,81,97,103}`,
  all realized. Anchors reproduced: `X_0(600)→…,47,49`; `X_0(1200)→25,49,52,97,103`.

---

## 3. Atkin–Lehner CD-quotients of `X_0(N)`, **non-CD sources allowed** — `harvestC.m`

This is the "nearly but not completely decomposable" idea: take `X_0(N)` which
need *not* be CD, quotient by an AL subgroup `U`, and **keep it iff the quotient
Jacobian is CD** (its `U`-invariant part contains only rational/elliptic newform
orbits — sound test, may miss CM cases). Genus = `dim` of the `U`-invariant
cusp-form subspace.

Reproduces `g(X_0(600)/w_24)=47` (CD) and `g(X_0(43)/w_43)=1` (CD).

Scans (all composite levels with `g(X_0(N)) ≥ 56`):
```
N ∈ [400,1500]  : CD AL-quotient genera ≥56 = {57,61,65,73,77,81,97,103,129}
N ∈ [1500,2400] : adds {76,79,157}   (source genus up to ~370)
NEW gap genera  : NONE.  Every genus found is already realized.
```

---

## 4. Quotients of the census curves `X_H` themselves, genus `> 112` — `harvestD.m`

The headline of Remark 2: quotients of the >100,000 genus-`>100` census curves.
Computed `FastQuotientGenera(H)` (the modular `GL2(Ẑ)`-internal quotients, via
the fast `N(H)/H` method) for **all 61,129** census curves of genus `> 112`
(levels 48, 60, 96, 120, 144, 180, 192 — all `< 240`).

```
processed 61,129 curves of genus > 112
Union of ALL quotient genera:
  1–5,7–19,21–25,27–29,32–38,40,41,43,45,49,51,52,53,55,57,61,65,68,69,73,
  75,76,77,78,79,81,82,85,97,103,105,109,113,121,129,135,137,157,159,161,
  163,169,193,205,217,325,409
Quotient genera that are GAPS (new): []   <-- NONE
```
Every quotient genus is already realized — exactly as §1 proves (these are
`GL2(Ẑ)`-internal quotients, so each `X_B` is itself a census curve).

Examples: `120.1440.113 → {25,53,57,113}`, `60.1728.121 → {13,25,29,57,61,121}`,
`48.2304.161 → {17,33,37,41,73,81,161}`, `60.5760.409 → {49,97,103,193,205,409}`.

---

## 5. Structural findings (why the modular-quotient route is stuck)

- **Census curves are not diamond curves.** Of 6000 sampled census curves, **0**
  contain the full unipotent — none is `X_0(N)`/`X_1(N)`/`X_Δ(N)`. They are
  general (Cartan-type / partial-Borel) subgroups, so there is no character-space
  shortcut for their Atkin–Lehner quotients, and explicit models are infeasible
  at genus `>100`.
- **AL is the only escape, and it is the hard "second pass."** Any quotient that
  becomes CD *inside* `GL2(Ẑ)` is already a census row (already realized). The
  only way to a new genus is a quotient by an automorphism **outside** `GL2(Ẑ)`
  (Atkin–Lehner / Fricke). For the general high-genus census `H` this needs a
  custom realization of `S_2(Γ_H)` + the AL operator, or an explicit model +
  `AutomorphismGroup` — the pass the paper's authors state they "have not
  attempted."
- **On diamond curves, AL only lives over `ℚ(ζ_N)`** (verified:
  `AtkinLehnerOperator(JOne(N), N)` errors over `ℚ`), while "completely
  decomposable" is a statement over `ℚ`; and a "degree-1 newform = elliptic"
  test is wrong for character spaces (those forms have non-rational
  `q`-coefficients — verified). A correct CD-check must detect elliptic-over-`ℚ`
  factors including CM/`ℚ`-curve cases.
- **LMFDB cross-check.** `gps_gl2zhat_fine` (modular curves well beyond the
  level-`<240` census) contains **no completely decomposable curve** (all
  `dims = 1`) at any gap genus `56,59,60,62,66,70,74,80,83,84`. The gaps are real
  across all known modular curves. The genus-56 curves it does have are exactly
  "nearly CD" — a few elliptic factors plus dimension-2,3,4,8 pieces.

---

## 6. Bottom line

Within everything tractable — every `GL2(Ẑ)`-internal quotient of the census
(incl. all genus-`>112` curves), and every Atkin–Lehner CD-quotient of `X_0(N)`
for `N ≤ 2400` (incl. nearly-CD sources) — **no new Ekedahl–Serre genus is
found.** The `GL2(Ẑ)`-internal route is *provably* empty; the `X_0(N)` AL route
is empirically empty and LMFDB shows no CD modular curve at any gap genus.

The remaining, genuinely uncomputed avenue is **Atkin–Lehner / full-automorphism
quotients of the general census `X_H`** (and intermediate diamond curves) with a
CM-aware over-`ℚ` CD-check — a real research build, not a code tweak. The
historically successful alternative for *new* genera is the **non-modular
automorphism / group-algebra method** (Paulhus, arXiv:1603.00331).

---

## Files

| file | what it does | result |
|---|---|---|
| `harvest.m` | original: Part A (`GL2(Ẑ)` quotients of census), Part B (`X_0(N)` single-involution AL, good levels) | A provably empty; B = 0 |
| `harvestB.m` | all AL subgroups of `X_0(N)`, 86 good levels | 0 new gaps |
| `harvestC.m` | AL CD-quotients of `X_0(N)`, non-CD sources allowed, scan `N` | 0 new gaps (`N ≤ 2400`) |
| `harvestD.m` | modular quotient genera of all census `X_H` with genus `> 112` | 0 new gaps |
| `RESULTS-claude-code.md` | detailed running write-up | — |
| `SUMMARY.md` | this file | — |

### Reproduce
```bash
git clone --recurse-submodules \
  https://github.com/AndrewVSutherland/CompletelyDecomposableModularJacobians
cd CompletelyDecomposableModularJacobians
git -C Magma fetch --depth 50 origin main && git -C Magma checkout 71a8fdd
cp ../harvest.m ../harvestB.m ../harvestC.m ../harvestD.m .
magma -b
> AttachSpec("spec"); load "harvest.m"; CheckAnchor();             // 47
> HarvestX0AtkinLehner("x0.txt");                                  // 0
> load "harvestB.m";                                               // 0
> load "harvestC.m"; ScanX0(400,1500,"scanC.txt");                 // 0 new gaps
> load "harvestD.m";
> RunD("grp_g112.txt","harvestD_gaps.txt");  // grp_g112.txt = awk -F: '$3>112' grpdata.txt ; 0 new gaps
```
