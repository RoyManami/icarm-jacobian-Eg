// Group-algebra (Kani-Rosen / Paulhus-Rojas) decomposition of Jacobians.
// For a curve X of genus g with G acting, signature [g0; m_1..m_r] and
// generating vector gv = (g_1..g_r), the Jacobian decomposes as
//   JX ~ prod_i B_i^{n_i}
// with dim B_i = (1/2) <psi_i, chi_V>,  psi_i = rational irreducible char,
//   chi_V = 2 chi_triv + 2(g0-1) rho_reg + sum_k (rho_reg - rho_<g_k>).
// Completely decomposable (via this method) iff every dim B_i in {0,1}.

SetColumns(0);

// dims of the factors B_i over the rational irreducible characters
function GADims(G, g0, gv)
  RT := RationalCharacterTable(G);
  cls := Classes(G);
  ordG := Order(G);
  dims := [];
  for psi in RT do
    deg := Integers()!(psi(Id(G)));
    // multiplicity of trivial char in psi = (1/|G|) sum_g psi(g)
    mtriv := (&+[ cls[i][2]*psi(cls[i][3]) : i in [1..#cls] ]) / ordG;
    val := 2*mtriv + 2*(g0-1)*deg;
    for x in gv do
      m := Order(x);
      s := &+[ psi(x^j) : j in [0..m-1] ];   // = m * dim of <x>-fixed space
      val +:= deg - s/m;
    end for;
    Append(~dims, Rationals()!(val/2));
  end for;
  return dims;
end function;

// Riemann-Hurwitz genus check
function RHgenus(ordG, g0, ms)
  mu := 2*g0 - 2 + &+[ 1 - 1/m : m in ms ];
  return 1 + ordG*mu/2;
end function;

// random search for a generating vector with prescribed element orders, product 1
function FindGV(G, orders : ntries := 200000)
  r := #orders;
  for t in [1..ntries] do
    gv := [];
    ok := true;
    for i in [1..r-1] do
      cnt := 0; x := Id(G);
      repeat x := Random(G); cnt +:= 1; until (Order(x) eq orders[i]) or (cnt gt 500);
      if Order(x) ne orders[i] then ok := false; break; end if;
      Append(~gv, x);
    end for;
    if not ok then continue; end if;
    last := (&*gv)^-1;
    if Order(last) ne orders[r] then continue; end if;
    Append(~gv, last);
    if sub<G | gv> eq G then return true, gv; end if;
  end for;
  return false, [];
end function;

// ---- validation 1: A_4 on genus 4, signature [0;2,3,3,3] -> E x E'^3 ----
G := AlternatingGroup(4);
ok, gv := FindGV(G, [2,3,3,3]);
printf "A4 [0;2,3,3,3]: found gv? %o\n", ok;
if ok then
  d := GADims(G, 0, gv);
  printf "  genus(RH)=%o   factor dims B_i = %o\n", RHgenus(12,0,[2,3,3,3]), Sort(d);
  printf "  all dims in {0,1}? %o\n", forall{x : x in d | x in {0,1}};
end if;

// ---- validation 2: PGL(2,7) on genus 36, signature [0;2,6,8] -> all elliptic ----
G := PGL(2,7);
ok, gv := FindGV(G, [2,6,8]);
printf "PGL(2,7) [0;2,6,8]: found gv? %o\n", ok;
if ok then
  d := GADims(G, 0, gv);
  printf "  genus(RH)=%o   nonzero factor dims = %o\n", RHgenus(336,0,[2,6,8]), Sort([x : x in d | x ne 0]);
  printf "  all dims in {0,1}? %o   #elliptic factors=%o\n",
         forall{x : x in d | x in {0,1}}, #[x : x in d | x eq 1];
end if;
quit;
