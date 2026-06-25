// Search for genus-56 completely decomposable quotients X_0(N)/W,
// W an Atkin-Lehner subgroup, over all levels with g(X_0(N)) >= 56.
// A quotient counts if dim(W-fixed cusp forms) = 56 AND that invariant
// subspace is all-elliptic (every Hecke eigenform in it is rational).

SetColumns(0);
Nlo := 11; Nhi := 1500;
TARGET := 56;

function HallDivisors(N)
  return [ d : d in Divisors(N) | Gcd(d, N div d) eq 1 ];
end function;

// rationality test: B basis (rows) of Hecke-stable subspace in S's g-dim basis
function IsAllElliptic(S, N, B)
  d := Nrows(B);
  if d eq 0 then return true; end if;
  for p in PrimesUpTo(100) do
    if N mod p eq 0 then continue; end if;
    Tp := Matrix(HeckeOperator(S, p));
    C := Solution(B, B*Tp);
    for fa in Factorization(CharacteristicPolynomial(C)) do
      if Degree(fa[1]) gt 1 then return false; end if;
    end for;
  end for;
  return true;
end function;

op := func< a,b | a*b div Gcd(a,b)^2 >;

nhits := 0;
for N in [Nlo..Nhi] do
  if Genus(Gamma0(N)) lt TARGET then continue; end if;
  M := ModularSymbols(N, 2, +1);
  S := CuspidalSubspace(M);
  g := Dimension(S);
  if g lt TARGET then continue; end if;
  I := IdentityMatrix(Rationals(), g);
  halls := HallDivisors(N);
  // AL operator matrices for every Hall divisor
  ALmat := AssociativeArray();
  ALmat[1] := I;
  for q in halls do
    if q eq 1 then continue; end if;
    ALmat[q] := Matrix(AtkinLehnerOperator(S, q));
  end for;
  // enumerate all subgroups W of the elementary-abelian AL group
  n := #halls;  // includes 1
  seen := {};
  for mask in [1..2^n-1] do
    gen := [];
    for i in [0..n-1] do
      if (mask div 2^i) mod 2 eq 1 then Append(~gen, halls[i+1]); end if;
    end for;
    W := {1} join Set(gen);
    repeat
      ch := false;
      for a in W do for b in W do
        c := op(a,b); if c notin W then Include(~W,c); ch := true; end if;
      end for; end for;
    until not ch;
    if W in seen then continue; end if;  Include(~seen, W);
    // Fix_W = intersection of (+1)-eigenspaces of all q in W
    Fix := VectorSpace(Rationals(), g);
    for q in W do
      if q eq 1 then continue; end if;
      Fix := Fix meet Kernel(ALmat[q] - I);
    end for;
    if Dimension(Fix) eq TARGET then
      B := BasisMatrix(Fix);
      if IsAllElliptic(S, N, B) then
        nhits +:= 1;
        printf "*** HIT: N=%o W=%o  genus(X_0(N)/W)=56 and ALL-ELLIPTIC ***\n", N, Sort(SetToSequence(W));
      end if;
    end if;
  end for;
  if N mod 25 eq 0 then printf "  ..scanned up to N=%o (g=%o), hits so far=%o\n", N, g, nhits; end if;
end for;
printf "DONE. total genus-56 completely-decomposable quotient hits: %o\n", nhits;
quit;
