// Region (1), lean: reconstruct PR completely-decomposable sources of genus
// 121 and 163, scan order-2 and order-3 subgroup quotients for genus 56.
SetColumns(0);

function GADims(G, RT, clsfull, g0, reps)
  ordG := Order(G); dims := [];
  for psi in RT do
    deg := Integers()!(psi(Id(G)));
    mtriv := (&+[ clsfull[i][2]*psi(clsfull[i][3]) : i in [1..#clsfull] ]) / ordG;
    val := 2*mtriv + 2*(g0-1)*deg;
    for x in reps do
      m := Order(x); val +:= deg - (&+[ psi(x^j) : j in [0..m-1] ])/m;
    end for;
    Append(~dims, Rationals()!(val/2));
  end for;
  return dims;
end function;

function chiV(G, gv)
  reg := PermutationCharacter(G, sub<G|>);
  cv := 2*PrincipalCharacter(G) - 2*reg;
  for x in gv do cv +:= reg - PermutationCharacter(G, sub<G|x>); end for;
  return cv;
end function;

sources := [ <[2,4,8],2592,163>, <[2,4,6],2880,121> ];

for s in sources do
  ms := s[1]; ord := s[2]; gx := s[3];
  F := FreeGroup(3);
  K := quo< F | F.1^ms[1], F.2^ms[2], F.3^ms[3], F.1*F.2*F.3 >;
  printf "SOURCE [0;%o] order %o genus %o: LINS...\n", ms, ord, gx;
  t0 := Cputime();
  L := LowIndexNormalSubgroups(K, ord);
  cands := [ N : N in L | N`Index eq ord ];
  printf "  LINS %o s, %o subgroups at index %o\n", Cputime(t0), #cands, ord;
  for N in cands do
    f, Qp := CosetAction(K, N`Group);
    gv := [ f(K.1), f(K.2), f(K.3) ];
    if [Order(x):x in gv] ne ms then continue; end if;
    cv := chiV(Qp, gv);
    if Integers()!(cv(Id(Qp))/2) ne gx then continue; end if;
    cl := Classes(Qp); clsfull := [<c[1],c[2],c[3]>:c in cl];
    d := GADims(Qp, RationalCharacterTable(Qp), clsfull, 0, gv);
    if not forall{z:z in d|z in {0,1}} then continue; end if;
    printf "  >> COMPLETELY DECOMPOSABLE genus-%o source reconstructed (%o elliptic factors). Scanning order-2,3 subgroups...\n",
           gx, #[z:z in d|z eq 1];
    qgen := {};
    for kk in [2,3] do
      for H in Subgroups(Qp : OrderEqual:=kk) do
        gq := Integers()!((1/2)*InnerProduct(cv, PermutationCharacter(Qp, H`subgroup)));
        Include(~qgen, <kk,gq>);
        if gq eq 56 then
          printf "     *** GENUS 56 QUOTIENT: H order %o -> COMPLETELY DECOMPOSABLE genus 56 ***\n", kk;
        end if;
      end for;
    end for;
    printf "     (order,quotient-genus) seen: %o\n", Sort(Setseq(qgen));
  end for;
end for;
printf "DONE region1b\n";
quit;
