// harvestD.m -- quotient genera of all census curves X_H with genus > 112.
// Uses the fast normalizer/quotient-group method.  Accumulates the union of
// all quotient genera and flags any that lands on a missing (gap) genus.

load "harvest.m";                       // KnownRealizedGenera
KNOWN  := KnownRealizedGenera();
GAPSET := {Integers()| g : g in [56..1297]} diff KNOWN;

function FastQuotientGenera(H)
    N := GL2Level(H); G := GL2Ambient(N);
    Nz := Normalizer(G, H);
    Q, pi := quo< Nz | H >;
    gen := {Integers()|};
    for srec in Subgroups(Q) do
        B := srec`subgroup @@ pi;        // H <= B <= N(H)
        if GL2DeterminantIndex(B) eq 1 and (-G!1) in B then
            Include(~gen, GL2Genus(B));
        end if;
    end for;
    return gen;
end function;

procedure RunD(infile, gapfile)
    fp := Open(gapfile, "w");
    allq := {Integers()|};
    gaphits := {Integers()|};
    nrows := 0;
    for line in Split(Read(infile), "\n") do
        if #line eq 0 then continue; end if;
        r := Split(line, ":");
        if #r lt 4 or StringToInteger(r[3]) le 112 then continue; end if;
        nrows +:= 1;
        H := GL2FromGenerators(r[1], r[2], r[4]);
        for g in FastQuotientGenera(H) do
            if not (g in allq) then
                Include(~allq, g);
                if g in GAPSET then
                    Include(~gaphits, g);
                    fprintf fp, "%o\tGAP-HIT\tquotient of %o\n", g, r[5]; Flush(fp);
                    printf "*** GAP GENUS %o from quotient of %o ***\n", g, r[5];
                end if;
            end if;
        end for;
        if nrows mod 1000 eq 0 then
            printf "[%o rows] union of quotient genera so far: %o\n", nrows, Sort(SetToSequence(allq));
            printf "          gap hits: %o\n", Sort(SetToSequence(gaphits));
        end if;
    end for;
    delete fp;
    printf "\nDONE: processed %o curves of genus>112.\n", nrows;
    printf "Union of ALL quotient genera: %o\n", Sort(SetToSequence(allq));
    printf "Quotient genera that are GAPS (new): %o\n", Sort(SetToSequence(gaphits));
end procedure;
