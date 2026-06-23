# Automorphism-quotient harvest for completely decomposable modular Jacobians

Turning the "obvious next step" of Paulhus–Sutherland (arXiv:2502.16007) —
taking quotients of the >100,000 completely decomposable modular curves they
found — into concrete new Ekedahl–Serre genera.

## The mechanism (and why every quotient works for free)

If `J_H` is isogenous to a product of elliptic curves, then by Poincaré
reducibility **every** isogeny factor is too. For any subgroup `B ≤ Aut(X_H)`,
`Jac(X_H/B)` is an isogeny factor of `J_H`, so `X_H/B` automatically has a
completely decomposable Jacobian. The only question is *which genus* `g(X_H/B)`
each quotient has. So the harvest is: enumerate quotients, read off genera,
keep the ones that land on an unrealized genus.

**Validated on real LMFDB data.** The completely decomposable curve
`60.480.37.bg.1` (genus 37) covers exactly five modular-quotient curves —
genera 3, 9, 15, 17, 19 — and all five are themselves completely decomposable
(`parent_is_cd = true` in `gps_gl2zhat_fine`). That is the mechanism in action.

**Where new genera come from.** Within the `N<240` census, quotient genera are
already-realized small genera. The *new* genera come from quotients of sources
*outside* the census — high-level `X_0(N)` with `J_0(N)` completely
decomposable, via Atkin–Lehner involutions. This is exactly how Ekedahl–Serre
got genus 47 = `X_0(600)/w_24`. The modular automorphisms of `X_H` are the
normalizer quotient `N(H)/H` in `GL2(Ẑ)` (Atkin–Lehner, diamond, Fricke all
live there), and the quotient of `X_H` by a subgroup of `N(H)/H` is the modular
curve `X_B` for `H ≤ B ≤ N(H)`, whose genus is just `GL2Genus(B)`.

## What is in this folder

| file | contents |
|---|---|
| `targets.json` | the gap list (genera not yet realized) + the 86 completely decomposable `X_0(N)` source levels + their genera, all derived from LMFDB |
| `harvest.m` | production Magma harvest: Part A over the census `grpdata.txt`, Part B over the `X_0(N)` Atkin–Lehner family; uses the repo's `GL2Genus` so quotient genera are computed exactly |
| `orchestrate.py` | prints the targets; ingests Magma output and reports which gaps closed |
| `good_levels.json` | the 86 good `N ≤ 1200` and `g(X_0(N))` for each |

## Source data (computed here, validated against the paper)

The 86 levels `N ≤ 1200` with `J_0(N)` completely decomposable were computed
from LMFDB classical-newform data: `N` is good iff every divisor `M | N` carries
only rational weight-2 newforms (Yamauchi's criterion). The resulting `X_0(N)`
genera reproduce Yamauchi's published list **exactly**:

```
0,1,2,3,4,5,6,7,8,9,10,11,13,17,19,21,25,29,33,37,43,49,53,55,57,61,73,97,121,161,205
```

and `g(X_0(1200)) = 205` matches Remark 1 of the paper. The high-level sources
(outside the `N<240` census) whose Atkin–Lehner quotients drive the harvest:

```
X_0(600):97   X_0(720):121   X_0(1152):161   X_0(1200):205
X_0(240):37   X_0(288):33    X_0(300):43     X_0(336):53
X_0(360):57   X_0(384):49    X_0(396):61     X_0(400):43   X_0(432):55   X_0(576):73
```

The four largest (97, 121, 161, 205) each have Atkin–Lehner group `(Z/2)^3`
(8 quotient curves; more once the `X_H` refinements between `Γ_1` and `Γ_0` are
included), and their quotient genera interpolate down through the gaps in
`[56, 205]`.

## Targets

- Smallest unrealized genus: **56**.
- 131 gaps in `[56, 240]`; 1181 gaps in `[56, 1297]`.
- The small ones to aim at first: 56, 59, 60, 66, 70, 74, 83, 84, 86–88, 90, …

## Running it (needs Magma + the repo)

```bash
git clone --recurse-submodules \
  https://github.com/AndrewVSutherland/CompletelyDecomposableModularJacobians
cd CompletelyDecomposableModularJacobians
# drop harvest.m here, then:
magma -b
> AttachSpec("spec");
> load "harvest.m";
> CheckAnchor();                                  // must print genus 47
> HarvestX0AtkinLehner("x0_al_newgenera.txt");    // Part B (fast)
> HarvestCensus("grpdata.txt","census_newgenera.txt": MinGenus:=100);  // Part A
```

then

```bash
python3 orchestrate.py targets.json x0_al_newgenera.txt census_newgenera.txt
```

Part B is cheap (a few hundred small groups). Part A is embarrassingly parallel
over the rows of `grpdata.txt`; split the file and run shards. Only `GL2Genus`
of intermediate groups is needed — no traceform/decomposition recomputation,
since complete decomposability is inherited automatically.

## Notes / honest boundaries

- The genus-`GL2Genus` step and the normalizer enumeration are delegated to the
  repo's GL2 package because it computes Atkin–Lehner fixed points / quotient
  genera exactly; a hand-rolled CM class-number formula was deliberately avoided
  to prevent convention bugs.
- A quotient genus equal to a gap is a *guaranteed* new realization (the curve
  exists and its Jacobian is completely decomposable by Poincaré). No further
  verification is needed beyond confirming `J_H` (the source) is completely
  decomposable, which for the `X_0(N)` family is Yamauchi's criterion.
- Beyond Atkin–Lehner, the census sources may have extra (non-modular)
  automorphisms; those quotients need `AutomorphismGroup(X_H)` on an explicit
  model and are a second, larger pass.
