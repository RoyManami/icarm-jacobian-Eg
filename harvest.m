// ============================================================================
//  harvest.m  --  automorphism-quotient harvest for completely decomposable
//                 modular Jacobians (Ekedahl-Serre gap filling)
//
//  Depends on the GL2 package shipped with
//    github.com/AndrewVSutherland/CompletelyDecomposableModularJacobians
//  (which bundles AndrewVSutherland/Magma as a submodule and gl2split.m).
//
//  Run:   magma -b
//         > AttachSpec("spec");
//         > load "harvest.m";
//
//  KEY FACT (Poincare reducibility): if J_H is isogenous to a product of
//  elliptic curves, then so is every isogeny factor, hence J(X_H / B) is
//  completely decomposable for EVERY B <= Aut(X_H).  So every quotient curve
//  we produce automatically has a completely decomposable Jacobian; the only
//  question is which genera arise.  Validated on LMFDB data: the five modular
//  quotients of 60.480.37.bg.1 (genera 3,9,15,17,19) are all completely
//  decomposable.
//
//  The modular automorphisms of X_H are exactly N(H)/H, where N(H) is the
//  normalizer of H in GL2(Zhat) (Atkin-Lehner + diamond + Fricke live here).
//  A quotient of X_H by a subgroup of N(H)/H is the modular curve X_B for the
//  intermediate group H <= B <= N(H); its genus is GL2Genus(B).  So the
//  quotient-genus harvest reduces to: enumerate B between H and N(H) (keeping
//  det(B) = Zhat^x and -I in B so X_B is nice over Q) and record GL2Genus(B).
// ============================================================================

// ---- targets: genera not yet realized (from the LMFDB-derived gap list) ----
GAPS := { g : g in [56..1297] } diff KnownRealizedGenera();  // see helper below

// KnownRealizedGenera(): the consolidated Paulhus-Sutherland set.
function KnownRealizedGenera()
    R := [<0,37>,<38,38>,<39,55>,<57,58>,<61,65>,<67,69>,<71,73>,<75,78>,
          <79,82>,<85,85>,<89,89>,<91,91>,<93,93>,<95,95>,<97,97>,<101,101>,
          <103,103>,<105,107>,<109,109>,<113,113>,<118,118>,<121,121>,
          <125,125>,<129,129>,<135,135>,<137,137>,<142,142>,<145,145>,
          <154,154>,<157,157>,<159,159>,<161,161>,<163,163>,<169,169>,
          <193,193>,<199,199>,<205,205>,<211,211>,<213,213>,<217,217>,
          <244,244>,<257,257>,<325,325>,<409,409>,<433,433>,<649,649>,
          <1297,1297>];
    return {Integers()| g : g in [a..b], t in R | true where a:=t[1] where b:=t[2]};
end function;

// ---------------------------------------------------------------------------
// Part A.  Quotient harvest over the completely decomposable census.
//   Input  : grpdata.txt rows (label N.i.g.c.n) with J_H completely decomposable.
//   Output : every quotient genus GL2Genus(B) for H <= B <= N(H), flagged if
//            it lands in GAPS.  Focus on the >100,000 sources of genus > 100.
// ---------------------------------------------------------------------------

function ModularQuotientGenera(H)
    // H : subgroup of GL(2,Z/NZ) of level N, det index 1, -I in H.
    N := GL2Level(H);
    G := GL2Ambient(N);
    Nz := Normalizer(G, H);                 // N(H) in GL2(Z/NZ)
    // intermediate groups H <= B <= Nz, up to conjugacy in Nz
    L := [ B`subgroup : B in Subgroups(Nz : OrderMultipleOf := #H) |
                        #B`subgroup mod #H eq 0 and H subset B`subgroup ];
    gen := {Integers()|};
    for B in L do
        if GL2DeterminantIndex(B) eq 1 and (-G!1) in B then
            Include(~gen, GL2Genus(B));     // genus of the quotient curve X_B
        end if;
    end for;
    return gen;
end function;

procedure HarvestCensus(infile, outfile : MinGenus := 100)
    fp := Open(outfile, "w");
    new := {Integers()|};
    for line in Split(Read(infile), "\n") do
        if #line eq 0 or line[1] eq '#' then continue; end if;
        H := GroupFromGrpdataRow(line);     // build H <= GL(2,Z/NZ) from the row
        if GL2Genus(H) lt MinGenus then continue; end if;
        for g in ModularQuotientGenera(H) do
            if g in GAPS and not (g in new) then
                Include(~new, g);
                fprintf fp, "%o\t%o\t(source %o)\n", g, "NEW", line;
            end if;
        end for;
    end for;
    delete fp;
    printf "Census harvest: %o new genera filled: %o\n", #new, Sort(SetToSequence(new));
end procedure;

// ---------------------------------------------------------------------------
// Part B.  Atkin-Lehner harvest over the X_0(N) family (N <= 1200 with
//   J_0(N) completely decomposable).  This is where ES got genus 47 from
//   X_0(600)/w_24, and it reaches sources OUTSIDE the N<240 census.
//
//   GOOD_N (86 levels, computed from LMFDB classical-newform data; the X_0(N)
//   genera reproduce Yamauchi's published list exactly):
// ---------------------------------------------------------------------------

GOOD_N := [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,
  28,30,32,33,34,36,37,38,40,42,44,45,48,49,50,52,54,56,57,60,64,66,72,75,76,80,
  84,90,96,99,100,108,112,114,120,121,128,132,144,150,168,180,192,198,200,216,
  240,288,300,336,360,384,396,400,432,576,600,720,1152,1200];

procedure HarvestX0AtkinLehner(outfile)
    fp := Open(outfile, "w");
    new := {Integers()|};
    for N in GOOD_N do
        H := GL2Borel(N);                   // Gamma_0(N) image in GL(2,Z/NZ)
        Nz := Normalizer(GL2Ambient(N), H); // contains all Atkin-Lehner elements
        for S in Subgroups(Nz : OrderMultipleOf := #H) do
            B := S`subgroup;
            if #B mod #H ne 0 or not (H subset B) then continue; end if;
            if GL2DeterminantIndex(B) eq 1 and (-GL2Ambient(N)!1) in B then
                g := GL2Genus(B);
                if g in GAPS and not (g in new) then
                    Include(~new, g);
                    fprintf fp, "%o\tNEW\tX_0(%o) Atkin-Lehner quotient\n", g, N;
                end if;
            end if;
        end for;
    end for;
    delete fp;
    printf "X_0(N) AL harvest: %o new genera: %o\n", #new, Sort(SetToSequence(new));
end procedure;

// ---------------------------------------------------------------------------
// Sanity anchor: X_0(600)/w_24 must come out as genus 47 (Ekedahl-Serre).
// ---------------------------------------------------------------------------
procedure CheckAnchor()
    assert 47 in ModularQuotientGenera(GL2Borel(600));
    printf "Anchor OK: genus 47 realized as an X_0(600) Atkin-Lehner quotient.\n";
end procedure;

// Entry point:
//   CheckAnchor();
//   HarvestX0AtkinLehner("x0_al_newgenera.txt");
//   HarvestCensus("grpdata.txt", "census_newgenera.txt" : MinGenus := 100);
