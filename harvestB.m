// harvestB.m -- exhaustive Atkin-Lehner SUBGROUP quotient harvest over the
// X_0(N) family (the 86 levels N<=1200 with J_0(N) completely decomposable).
//
// For each good N: the Atkin-Lehner group W(N) ~ (Z/2)^k (k = #primes | N).
// For EVERY subgroup U <= W(N), the quotient X_0(N)/U is a modular curve whose
// Jacobian is an isogeny factor of J_0(N) (Poincare), hence completely
// decomposable.  Its genus = dim of the simultaneous (+1)-eigenspace of the AL
// involutions generating U on S_2(Gamma_0(N)) = n - Rank(stack of (M_w - I)).
//
// This generalizes the shipped Part B (single involutions only) to all 2^k
// subgroups, which is where the rank-2/rank-3 (smaller-genus) quotients live.

load "harvest.m";   // for KnownRealizedGenera / GAPS

GOOD_N := [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,
  28,30,32,33,34,36,37,38,40,42,44,45,48,49,50,52,54,56,57,60,64,66,72,75,76,80,
  84,90,96,99,100,108,112,114,120,121,128,132,144,150,168,180,192,198,200,216,
  240,288,300,336,360,384,396,400,432,576,600,720,1152,1200];

KNOWN := KnownRealizedGenera();
GAPSET := {Integers()| g : g in [56..1297]} diff KNOWN;

procedure HarvestX0ALsubgroups(outfile)
    fp := Open(outfile, "w");
    newgap := {Integers()|};
    allquot := {Integers()|};   // all proper-quotient genera seen (any size)
    for N in GOOD_N do
        ppeds := [p^Valuation(N, p) : p in PrimeDivisors(N)];  // AL generators
        k := #ppeds;
        if k eq 0 then continue; end if;
        MS := ModularSymbols(N, 2, +1);
        S := CuspidalSubspace(MS);
        n := Dimension(S);
        if n eq 0 then continue; end if;
        // AL involution matrices for the k prime-power generators.
        Mgen := [Matrix(Rationals(), AtkinLehner(S, pp)) : pp in ppeds];
        I_n := IdentityMatrix(Rationals(), n);
        // Enumerate all subgroups of the AL group (Z/2)^k.
        A := AbelianGroup([2 : i in [1..k]]);
        for srec in Subgroups(A) do
            U := srec`subgroup;
            if #U eq 1 then continue; end if;   // trivial -> the curve itself
            // genus of quotient = dim of common (+1)-eigenspace of generating involutions.
            rows := [];
            for v in Generators(U) do
                c := Eltseq(A!v);
                M := I_n;
                for i := 1 to k do if c[i] ne 0 then M := M * Mgen[i]; end if; end for;
                Append(~rows, M - I_n);
            end for;
            Big := VerticalJoin(rows);
            g := n - Rank(Big);
            Include(~allquot, g);
            if g in GAPSET and not (g in newgap) then
                Include(~newgap, g);
                fprintf fp, "%o\tNEW-GAP\tX_0(%o)/U  (AL subgroup, dim %o)\n", g, N, Dimension(U);
            end if;
        end for;
    end for;
    delete fp;
    printf "AL-subgroup harvest: %o NEW gap genera: %o\n", #newgap, Sort(SetToSequence(newgap));
    printf "all proper-quotient genera in [56,1297]: %o\n",
        Sort([g : g in allquot | g ge 56 and g le 1297]);
    printf "of those, how many are gaps: %o\n", #({g : g in allquot | g ge 56} meet GAPSET);
end procedure;

HarvestX0ALsubgroups("x0_al_subgroup_newgenera.txt");
