// Region (1), large sources via explicit Paulhus-Rojas presentations:
// genus 193 (order 5760, [0;2,3,10]) and genus 244 (order 11664, [0;2,3,8]).
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

procedure DoSource(G, gx, Horders)
  // faithful permutation representation (regular rep)
  fa, Gp := CosetAction(G, sub<G|>);
  gv := [ fa(G.1), fa(G.2), fa(G.3) ];
  printf "  order=%o  gen orders=%o\n", Order(Gp), [Order(x):x in gv];
  cv := chiV(Gp, gv);
  g := Integers()!(cv(Id(Gp))/2);
  printf "  computed genus = %o (expected %o)\n", g, gx;
  cl := Classes(Gp); clsfull := [<c[1],c[2],c[3]>:c in cl];
  d := GADims(Gp, RationalCharacterTable(Gp), clsfull, 0, gv);
  printf "  completely decomposable? %o (%o elliptic factors, max factor dim %o)\n",
         forall{z:z in d|z in {0,1}}, #[z:z in d|z eq 1], Maximum(d);
  if not forall{z:z in d|z in {0,1}} then return; end if;
  qgen := {};
  for kk in Horders do
    for H in Subgroups(Gp : OrderEqual:=kk) do
      gq := Integers()!((1/2)*InnerProduct(cv, PermutationCharacter(Gp, H`subgroup)));
      Include(~qgen, <kk,gq>);
      if gq eq 56 then
        printf "     *** GENUS 56 QUOTIENT: H order %o -> COMPLETELY DECOMPOSABLE genus 56 ***\n", kk;
      end if;
    end for;
  end for;
  printf "     (Horder, quotient genus) seen: %o\n", Sort(Setseq(qgen));
end procedure;

// ---- genus 193, order 5760, [0;2,3,10] ----
printf "SOURCE genus 193 (order 5760, [0;2,3,10]):\n";
F1 := FreeGroup(3); x:=F1.1; y:=F1.2; z:=F1.3;
G193 := quo< F1 |
  x^2, y^3, z^10, z^-1*y^-1*x,
  x*z^2*y*z^-1*x*z*y^-1*z^-2*x*z*y^-1*z^-2,
  y*z^-1*x*z^4*y*z^-1*x*y^-1*z^-1*x*y^-1*z^-2*x*z^4*y*z^-1*x*y^-1*z^-1*x*z,
  z^2*y^-1*z^-4*x*z*y*z^-1*x*z^2*y*x*z*y*x*z^-1*x*z*y^-1*z^-3*x >;
DoSource(G193, 193, [2,3]);

// ---- genus 244, order 11664, [0;2,3,8] ----
printf "SOURCE genus 244 (order 11664, [0;2,3,8]):\n";
F2 := FreeGroup(3); x:=F2.1; y:=F2.2; z:=F2.3;
G244 := quo< F2 |
  x^2, y^3, z^8, z^-1*y^-1*x,
  z*y*x*z*y*x*z*y*x*y^-1*z^-1*x*y^-1*z^-1*x*y^-1*z^-1*x,
  z^2*y*x*z^2*y*x*z^2*y*x*z^2*y*x*y^-1*x*y^-1*z^-1*x*z*y^-1*z^-1*x >;
DoSource(G244, 244, [2,3,4]);

printf "DONE region1c\n";
quit;
