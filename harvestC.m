// harvestC.m -- Atkin-Lehner quotient harvest realizing Remark 2 of Paulhus-Sutherland.
//
// For X_0(N) and EVERY subgroup U of its Atkin-Lehner group (Z/2)^k, compute:
//   * genus(X_0(N)/U)  = dim of the U-invariant cusp-form subspace;
//   * whether J(X_0(N)/U) is COMPLETELY DECOMPOSABLE, i.e. its U-invariant part
//     contains only rational (degree-1) newform orbits  ==> isogenous /Q to a
//     product of elliptic curves.  (Sufficient, sound; may miss CM cases.)
//
// This is the genus-47 = X_0(600)/w_24 mechanism, generalized: the SOURCE need
// NOT be completely decomposable, only the quotient.  We flag any CD quotient
// whose genus is an unrealized (gap) genus -- that would be a NEW Ekedahl-Serre
// genus.

load "harvest.m";                       // KnownRealizedGenera
KNOWN := KnownRealizedGenera();
GAPSET := {Integers()| g : g in [56..1297]} diff KNOWN;

// genus + CD? of all AL-subgroup quotients of X_0(N).  Returns list of <genus, isCD, Udesc>.
function ALQuotients(N)
    S := CuspidalSubspace(ModularSymbols(N, 2, 1));
    g0 := Dimension(S);
    if g0 eq 0 then return [], 0; end if;
    D := NewformDecomposition(S);
    degs := [Dimension(AssociatedNewSpace(P)) : P in D];
    ppeds := [p^Valuation(N, p) : p in PrimeDivisors(N)];   // AL generators
    k := #ppeds;
    // Precompute, per piece, the AL operator for each prime-power generator.
    pieceOps := [* [ AtkinLehnerOperator(D[j], pp) : pp in ppeds ] : j in [1..#D] *];
    pieceId  := [* IdentityMatrix(Rationals(), Dimension(D[j])) : j in [1..#D] *];
    res := [];
    A := AbelianGroup([2 : i in [1..k]]);
    for srec in Subgroups(A) do
        U := srec`subgroup;
        if #U eq 1 then continue; end if;          // trivial -> the curve itself
        basis := [Eltseq(A!v) : v in Generators(U)];
        gtot := 0; gbad := 0;
        for j := 1 to #D do
            // common (+1)-eigenspace of U's generators on piece D[j]
            rows := [];
            for c in basis do
                M := pieceId[j];
                for i := 1 to k do if c[i] ne 0 then M := M * pieceOps[j][i]; end if; end for;
                Append(~rows, M - pieceId[j]);
            end for;
            inv := Dimension(D[j]) - Rank(VerticalJoin(rows));
            gtot +:= inv;
            if degs[j] gt 1 then gbad +:= inv; end if;
        end for;
        Append(~res, <gtot, gbad eq 0>);
    end for;
    return res, g0;
end function;

procedure ScanX0(Nlo, Nhi, outfile)
    fp := Open(outfile, "w");
    newgaps := {Integers()|};
    cdgenera := {Integers()|};            // all CD quotient genera >= 56 seen
    for N := Nlo to Nhi do
        if #PrimeDivisors(N) lt 2 then continue; end if;   // need a nontrivial AL group
        quots, g0 := ALQuotients(N);
        if g0 lt 56 then continue; end if;                 // can't reach a gap
        for q in quots do
            g := q[1]; cd := q[2];
            if cd and g ge 56 then
                Include(~cdgenera, g);
                if g in GAPSET and not (g in newgaps) then
                    Include(~newgaps, g);
                    fprintf fp, "%o\tNEW-GAP\tcompletely-decomposable AL quotient of X_0(%o) (g0=%o)\n", g, N, g0;
                    Flush(fp);
                    printf "*** NEW GAP GENUS %o from X_0(%o) (source genus %o) ***\n", g, N, g0;
                end if;
            end if;
        end for;
        printf "N=%o g0=%o  CD-quotient genera>=56 so far: %o\n", N, g0, Sort(SetToSequence(cdgenera));
    end for;
    delete fp;
    printf "DONE [%o,%o]. NEW gap genera: %o\n", Nlo, Nhi, Sort(SetToSequence(newgaps));
    printf "All CD AL-quotient genera >=56 found: %o\n", Sort(SetToSequence(cdgenera));
end procedure;
