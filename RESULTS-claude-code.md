# Results: automorphism-quotient harvest (rerun + fix)

Rerun of the harvest described in `README-claude-code.md`, against
`AndrewVSutherland/CompletelyDecomposableModularJacobians`, with Magma 2.29-8.

## 1. The bug that broke `harvest.m`

`harvest.m` Part A failed not in `harvest.m` itself but in the **GL2 package it
depends on**. The symptom chain:

```
GL2FromGenerators(...)  ->  Runtime error: Undefined reference 'atoi'
```

`atoi` lives in `Magma/utils.m`. That file (and `Magma/genus1.m`) failed to
attach with `syntax error` / `bad syntax` — `genus1.m` line ~224 had
`end intrinsic` **missing its semicolon**, which Magma 2.29 rejects. So the
whole `utils.m`/`genus1.m` compilation aborted, `atoi` never got defined, and
every `GL2…` call that parses strings (i.e. all of Part A) blew up.

**Root cause:** the `Magma` submodule was pinned in CDMJ at an old commit
(`eed5c96`) that predates the upstream fix.

**Fix:** update the submodule to the commit that fixes it:

```bash
git submodule update --init --depth 1
git -C Magma fetch --depth 50 origin main
git -C Magma checkout 71a8fdd      # "Merge … main"; includes a46c13b "fixing semicolons"
```

After this, `AttachSpec("spec")` is clean, `GL2FromGenerators`, `GL2Genus`,
`GL2DeterminantIndex`, the normalizer, etc. all work, and `CheckAnchor()`
prints `g(X_0(600)/w_24) = 47`.

(If you cloned without `--recurse-submodules`, the `Magma/` dir is empty — that
alone also breaks everything; the first command above fixes that.)

## 2. Part A (census quotients) is provably empty — no need to run it

A single genus-161/genus-409 source already runs >2 min: `Normalizer` +
`Subgroups` inside `GL2(Z/48Z)` / `GL2(Z/60Z)` (order ~1–2 million) does not
scale to 88,034 high-genus rows. But it doesn't matter, because of this:

> **Every modular quotient of a census curve is itself in the census.**
> If `H` is a census group and `H ≤ B` (any quotient curve `X_B` that `X_H`
> covers, in particular any automorphism quotient), then `B` is defined mod the
> same `N`, has surjective determinant, contains `-I`, has level `< 240`, and
> `Jac(X_B)` is an isogeny factor of `Jac(X_H)` hence completely decomposable
> (Poincaré). Those are **exactly** the census membership conditions, so
> `B ∈ grpdata.txt`.

Therefore every Part-A quotient genus is a genus already present in
`grpdata.txt`. Computed directly:

```
distinct census genera: 1..29,31..41,43,45,49..53,55,57,61,65,68,69,72,73,
  75..79,81,82,85,89,97,103,105,109,113,121,129,135,137,157,159,161,163,169,
  193,205,217,325,409
census genera NOT already realized: []     <-- empty
```

So **Part A cannot produce a new Ekedahl–Serre genus.** The 88k normalizer
computations are unnecessary. (This is the rigorous form of the README's note
that "within the N<240 census, quotient genera are already-realized small
genera.")

## 3. Part B (X_0(N) Atkin–Lehner), strengthened to all AL subgroups — also 0

The shipped Part B only quotients by single involutions `w_d`. `harvestB.m`
strengthens this to **every subgroup `U` of the full Atkin–Lehner group
`W(N) ≅ (Z/2)^k`** for all 86 good levels `N ≤ 1200`, with
`g(X_0(N)/U) = dim` of the simultaneous `(+1)`-eigenspace of `U` on
`S_2(Γ_0(N))`.

```
AL-subgroup harvest: 0 NEW gap genera: []
all proper-quotient genera in [56,1297]: [57, 61, 77, 81, 97, 103]
of those, how many are gaps: 0
```

Full spreads for the four largest sources (47 reproduces Ekedahl–Serre):

```
X_0(600)  g=97   -> 10,20,21,22,23,24,25,43,45,47,49
X_0(720)  g=121  -> 13,27,29,31,57,61
X_0(1152) g=161  -> 37,77,81
X_0(1200) g=205  -> 25,49,52,97,103
```

Every quotient genus ≥ 56 ({57,61,77,81,97,103}) is already realized.

## 4. IMPORTANT CORRECTION — §2/§3 only cover quotients *inside* `GL2(Ẑ)`

The "Part A is provably empty" statement above is correct **only for the
quotients `harvest.m` actually computes**: those coming from the normalizer of
`H` *inside* `GL2(ℤ/Nℤ)` (diamond operators and `GL2`-internal automorphisms).
It does **NOT** cover **Atkin–Lehner / Fricke** involutions, and that is exactly
the gap Remark 2 of the paper points at.

Why the census-closure argument fails for Atkin–Lehner: AL/Fricke matrices
(e.g. `(0,-1;N,0)`) have determinant `N`, not a unit — **they do not live in
`GL2(Ẑ)`**. So `X_H/⟨w⟩` is the modular curve for `⟨H, w⟩`, which is *not a
subgroup of `GL2(Ẑ)`*, hence *not* in the census, hence its genus *can* be new.
The Ekedahl–Serre genus 47 = `X_0(600)/w_24` is precisely this kind.

Concrete proof the normalizer misses AL: for `X_0(30)` (= Borel mod 30) the
`GL2(ℤ/30)`-normalizer gives `|N(H)/H| = 1` — `harvest.m` Part A finds **zero**
quotients of `X_0(30)`, even though its Atkin–Lehner group is `(ℤ/2)^3`.

## 5. Remark 2, done correctly: Atkin–Lehner CD-quotients (`harvestC.m`)

Remark 2 has two ideas. We implemented the tractable one over `ℚ`:

> take an `X_0(N)` (which need NOT be completely decomposable — "nearly CD"),
> quotient by an Atkin–Lehner subgroup `U`, and keep it if the quotient
> Jacobian *is* completely decomposable (its `U`-invariant part contains only
> rational/elliptic newform orbits). This is the genus-47 mechanism generalized
> to non-CD sources.

`harvestC.m` computes, for `X_0(N)` and **every** AL subgroup `U`:
`genus = dim(U-invariant cusp forms)` and `CD? = (no degree>1 newform orbit is
U-invariant)`. It reproduces both anchors: `g(X_0(600)/w_24)=47` (CD) and
`g(X_0(43)/w_43)=1` (CD).

Scans (all composite levels with `g(X_0(N)) ≥ 56`):

```
N ∈ [400,1500] : CD AL-quotient genera ≥56 = {57,61,65,73,77,81,97,103,129}
N ∈ [1500,2400]: adds {76,79,157}   (source genus up to ~370)
NEW gap genera in either range: NONE.   Every genus found is already realized.
```

So the over-`ℚ`, trivial-character branch of Remark 2 — *including* non-CD
sources whose quotients become CD — yields **no new Ekedahl–Serre genus**
through `N ≈ 2400`.

## 6. Why the diamond curves / general `X_H` AL-quotients are the hard "second pass"

Both natural extensions hit the same wall:

- **Quotients of the census `X_H` themselves (general `H` from `grpdata.txt`).**
  The `GL2(Ẑ)`-internal quotients are provably already realized (§2). The
  promising ones are again **Atkin–Lehner**, but for a *general* `H` there is no
  off-the-shelf `ModularSymbols(Γ_H)` + `AtkinLehner` in Magma — it needs an
  explicit model + `AutomorphismGroup`, or a hand-built realization of
  `S_2(Γ_H)` inside `S_2(Γ(N))` with the Fricke operator. `X_0(N)` was used as
  the tractable proxy that exhibits the identical mechanism.
- **Intermediate diamond curves `X_Δ(N)`** (`Γ_1 ⊆ Γ_H ⊆ Γ_0`). These DO carry
  extra content vs `X_0(N)` (non-trivial-character newforms). But: (i) on the
  character/diamond Jacobian, `AtkinLehnerOperator` is only defined over
  `ℚ(ζ_N)` (verified: it errors over `ℚ`), while "completely decomposable" is a
  statement over `ℚ`; and (ii) a naive "degree-1 newform = elliptic" test is
  WRONG for character spaces — those degree-1-over-`ℚ(χ)` forms have
  non-rational `q`-coefficients (verified) and are dimension `≥2` over `ℚ`. A
  correct CD-check must detect elliptic-over-`ℚ` factors including the
  **CM / `ℚ`-curve** character newforms (the RSZB cases the census's `cmfdata`
  files use) — i.e. Frobenius-trace matching, not a degree count.

`Decomposition` of `JOne(N)` does cleanly expose elliptic vs higher-dim factors
(e.g. `J_1(33)` → dims `[1,1,1,8,4,4,2]`), so a `ModularAbelianVariety`-based
pipeline (diamond Jacobian → AL quotient → `Decomposition` → flag CD iff every
factor is elliptic/`ℚ`) is the right tool — but doing AL over the correct field
and the CM-aware over-`ℚ` test is a genuine build, exactly the pass the paper's
authors say they "have not attempted."

## 7. Bottom line (updated)

- `GL2(Ẑ)`-internal quotients of census curves: **provably already realized**
  (§2) — Part A as written is empty *for the right reason*, but only for these.
- Atkin–Lehner CD-quotients of `X_0(N)`, `N ≤ 2400`, including nearly-CD
  sources (the genus-47 mechanism): computed, **0 new genera**.
- Atkin–Lehner CD-quotients of the general census `X_H` and of the diamond
  curves `X_Δ(N)`: **not yet computed** — the real "second pass" (explicit
  models / character-aware AL over `ℚ` + CM-aware CD-check). This is the open
  next step.

All conclusions assume the realized-genus set `KnownRealizedGenera()`
(smallest gap = 56) is correct.

## Reproduce

```bash
git clone --recurse-submodules \
  https://github.com/AndrewVSutherland/CompletelyDecomposableModularJacobians
cd CompletelyDecomposableModularJacobians
git -C Magma fetch --depth 50 origin main && git -C Magma checkout 71a8fdd
cp ../harvest.m ../harvestB.m ../harvestC.m .
magma -b
> AttachSpec("spec");
> load "harvest.m";  CheckAnchor();                          // 47
> HarvestX0AtkinLehner("x0_al_newgenera.txt");               // 0 (single involutions)
> load "harvestB.m";                                         // 0 (all AL subgroups, good levels)
> load "harvestC.m";  ScanX0(400,1500,"scanC.txt");          // 0 new gaps (CD AL quotients of X_0(N))
```
