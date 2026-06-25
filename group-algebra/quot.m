// Genera of Atkin-Lehner quotients X_0(N)/W for completely-decomposable levels N.
// For such N, J_0(N) ~ product of elliptic curves, so EVERY quotient X_0(N)/W
// has completely decomposable Jacobian. We search for quotient genus = 56.

goodlevels := [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28,30,32,33,34,36,37,38,40,42,44,45,48,49,50,52,54,56,57,60,64,66,72,75,76,80,84,90,96,99,100,108,112,114,120,121,128,132,144,150,168,180,192,198,200,216,240,288,300,336,360,384,396,400,432,576,600,720,1152,1200];

// Hall divisors q || N (gcd(q, N/q)=1), these index the Atkin-Lehner group.
function HallDivisors(N)
  ds := [];
  for d in Divisors(N) do
    if Gcd(d, N div d) eq 1 then Append(~ds, d); end if;
  end for;
  return ds;
end function;

allreached := {};
hit56 := [];

for N in goodlevels do
  if N le 10 then continue; end if;  // genus 0, nothing
  M := ModularSymbols(N, 2, +1);
  S := CuspidalSubspace(M);
  g := Dimension(S);
  if g lt 28 then
     // a single involution at best halves genus; need g >= 2*56-? to reach 56.
     // quotient genus <= g, and to reach 56 need g >= 56. Skip small g.
     if g lt 56 then continue; end if;
  end if;
  halls := HallDivisors(N);           // includes 1
  // trace of each Atkin-Lehner operator on S (weight-2 cusp forms, sign +1 space)
  tr := AssociativeArray();
  tr[1] := g;
  for q in halls do
    if q eq 1 then continue; end if;
    w := AtkinLehnerOperator(S, q);
    tr[q] := Trace(Matrix(w));
  end for;
  // The Atkin-Lehner group: elements are Hall divisors, op: q1*q2 -> q1*q2/gcd(q1,q2)^2
  grp := halls;
  op := function(a,b) c := a*b div Gcd(a,b)^2; return c; end function;
  // enumerate all subgroups W (genss of grp closed under op, containing 1)
  // grp is small (<= 8 elements)
  n := #grp;
  Q := {};  // set of quotient genera for this N
  for mask in [0..2^n-1] do
    gens := [];
    for i in [0..n-1] do
      if (mask div 2^i) mod 2 eq 1 then Append(~gens, grp[i+1]); end if;
    end for;
    if #gens eq 0 then continue; end if;
    // close under op to form subgroup W
    W := {1} join Set(gens);
    repeat
      changed := false;
      for a in W do for b in W do
        c := op(a,b);
        if c notin W then Include(~W, c); changed := true; end if;
      end for; end for;
    until not changed;
    // dim of W-fixed cusp forms = (1/|W|) sum_{w in W} tr(w)
    s := 0;
    for q in W do s +:= tr[q]; end for;
    gq := Integers()!(s / #W);          // genus of X_0(N)/W
    Include(~Q, gq);
  end for;
  allreached := allreached join Q;
  if 56 in Q then
    Append(~hit56, N);
    printf "*** N=%o (g=%o): quotient of genus 56 FOUND ***\n", N, g;
  end if;
  printf "N=%o g=%o quotient-genera=%o\n", N, g, Sort(SetToSequence(Q));
end for;

printf "\n==== SUMMARY ====\n";
printf "Levels giving a genus-56 quotient: %o\n", hit56;
printf "All quotient genera reached (completely decomposable): %o\n", Sort(SetToSequence(allreached));
printf "Is 56 reachable this way? %o\n", 56 in allreached;
quit;
