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

## 4. Bottom line

Within the scope of this pipeline — census curves of level `< 240`, plus
`X_0(N)` for the 86 good levels `N ≤ 1200`, quotiented by their **modular**
(diamond / Atkin–Lehner) automorphisms — **no new Ekedahl–Serre genus is
found, and Part A is mathematically guaranteed to find none.**

### Where new genera could still come from (honest forward directions)

The census-closure argument in §2 only kills *modular* quotients. Genuinely
open avenues, none of which this pipeline covers:

1. **Non-modular automorphisms of census curves.** A quotient `X_H/B` by extra
   automorphisms need not be a modular curve `X_{B'}`, so it escapes the census
   and its genus can be new. Requires `AutomorphismGroup` on an explicit model
   of each `X_H` (the README's "second, larger pass") — expensive.
2. **`X_0(N)` with `N > 1200`** whose `J_0(N)` is completely decomposable
   (extend Yamauchi's good-level list), then their AL quotients.
3. **Intermediate curves strictly between `Γ_1(N)` and `Γ_0(N)`** for good `N`
   (level can be ≥ 240, so not in the census, and not `X_0(N)`): genuinely new
   sources with AL + diamond quotient genera worth scanning.

All conclusions assume the realized-genus set encoded in
`KnownRealizedGenera()` (smallest gap = 56) is correct; the negative result is
relative to that set.

## Reproduce

```bash
git clone --recurse-submodules \
  https://github.com/AndrewVSutherland/CompletelyDecomposableModularJacobians
cd CompletelyDecomposableModularJacobians
git -C Magma fetch --depth 50 origin main && git -C Magma checkout 71a8fdd
cp ../harvest.m ../harvestB.m .
magma -b
> AttachSpec("spec");
> load "harvest.m";  CheckAnchor();                         // 47
> HarvestX0AtkinLehner("x0_al_newgenera.txt");              // 0
> load "harvestB.m";                                        // 0 (all AL subgroups)
```
