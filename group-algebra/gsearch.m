// Group-algebra search engine for completely decomposable Jacobians of a target genus.
// Key fact: dim B_i depends on the generating vector only via the conjugacy CLASSES
// of the branch generators, so we enumerate class-multisets directly. A genuine
// generating vector (=> a curve) exists when the Frobenius structure constant > 0
// and the chosen classes generate G.
SetColumns(0);

function GADims(G, RT, clsdata, g0, reps)
  // reps: list of branch-point class representatives
  ordG := Order(G);
  dims := [];
  for psi in RT do
    deg := Integers()!(psi(Id(G)));
    mtriv := (&+[ clsdata[i][2]*psi(clsdata[i][3]) : i in [1..#clsdata] ]) / ordG;
    val := 2*mtriv + 2*(g0-1)*deg;
    for x in reps do
      m := Order(x);
      val +:= deg - (&+[ psi(x^j) : j in [0..m-1] ])/m;
    end for;
    Append(~dims, Rationals()!(val/2));
  end for;
  return dims;
end function;

function StructConst(CT, clsdata, ordG, idxs)
  // number of tuples (g_1..g_r), g_i in class idxs[i], with product = 1
  r := #idxs;
  reps := [ clsdata[k][3] : k in idxs ];
  sz := &*[ clsdata[k][2] : k in idxs ];
  s := 0;
  for chi in CT do
    num := &*[ chi(reps[i]) : i in [1..r] ];
    s +:= num / chi(Id(Parent(reps[1])))^(r-2);
  end for;
  return Rationals()!((sz/ordG)*s);
end function;

// enumerate non-decreasing class-index tuples whose orders satisfy the RH sum
procedure rec(clsdata, orders, target, ~combos, cur, start, rem, rmax)
  if rem eq 0 then
    if #cur ge 3 then Append(~combos, cur); end if;
    return;
  end if;
  if #cur ge rmax then return; end if;
  for k in [start..#clsdata] do
    o := clsdata[k][1];
    t := 1 - 1/o;
    if t le rem + 1/1000000 then
      // allow if can still reach (prune: remaining slots * max term >= rem)
      rec(clsdata, orders, target, ~combos, Append(cur,k), k, rem - t, rmax);
    end if;
  end for;
end procedure;

function SearchGroup(G, targetg : g0:=0, rmax:=6, verbose:=false)
  ordG := Order(G);
  munum := 2*(targetg-1);            // = ordG * (Sum(1-1/m_i) + 2g0 - 2)
  if munum mod 1 ne 0 then return []; end if;
  target := munum/ordG - 2*g0 + 2;   // = Sum(1-1/m_i)
  if target le 0 then return []; end if;
  cl := Classes(G);                  // [<order,size,rep>]
  clsdata := [ <c[1], c[2], c[3]> : c in cl ];
  RT := RationalCharacterTable(G);
  CT := CharacterTable(G);
  combos := [];
  rec(clsdata, [], target, ~combos, [], 2, target, rmax);  // skip class 1 (identity, order1)
  found := [];
  for idxs in combos do
    reps := [ clsdata[k][3] : k in idxs ];
    if sub<G | reps> ne G then continue; end if;       // must generate
    sc := StructConst(CT, clsdata, ordG, idxs);
    if sc le 0 then continue; end if;                  // need a product-1 tuple
    dims := GADims(G, RT, clsdata, g0, reps);
    if forall{x : x in dims | x in {0,1}} then
      orders := [ clsdata[k][1] : k in idxs ];
      Append(~found, <orders, [x : x in dims | x ne 0]>);
      printf "   COMPLETELY DECOMPOSABLE: signature [%o; %o]  ->  %o elliptic factors (dims %o)\n",
             g0, orders, #[x:x in dims|x eq 1], dims;
    end if;
  end for;
  return found;
end function;

// ---------- positive control: genus 36 PGL(2,7) should be found ----------
printf "== control PGL(2,7), target genus 36 ==\n";
_ := SearchGroup(PGL(2,7), 36 : rmax:=4);

// ---------- positive control: genus 57 via SmallGroup(1344,11289) ----------
printf "== control SmallGroup(1344,11289), target genus 57 ==\n";
_ := SearchGroup(SmallGroup(1344,11289), 57 : rmax:=4);
quit;
