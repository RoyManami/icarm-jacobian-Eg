// Group-algebra search engine (optimized): cheap class-combo enumeration first,
// character tables only when genus is reachable.

function GADims(G, RT, clsdata, g0, reps)
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
  r := #idxs;
  reps := [ clsdata[k][3] : k in idxs ];
  sz := &*[ clsdata[k][2] : k in idxs ];
  s := 0;
  for chi in CT do
    s +:= (&*[ chi(reps[i]) : i in [1..r] ]) / chi(Id(Parent(reps[1])))^(r-2);
  end for;
  return Rationals()!((sz/ordG)*s);
end function;

procedure rec(clsdata, ~combos, cur, start, rem, rmax)
  if rem eq 0 then
    if #cur ge 3 then Append(~combos, cur); end if;
    return;
  end if;
  if #cur ge rmax then return; end if;
  // prune: even using the largest available (1-1/o) for all remaining slots, can we reach rem?
  maxterm := 1 - 1/clsdata[#clsdata][1];
  if (rmax - #cur)*maxterm lt rem - 1/1000000 then return; end if;
  for k in [start..#clsdata] do
    t := 1 - 1/clsdata[k][1];
    if t le rem + 1/1000000 then
      rec(clsdata, ~combos, Append(cur,k), k, rem - t, rmax);
    end if;
  end for;
end procedure;

function SearchGroup(G, targetg : g0:=0, rmax:=6, tag:="")
  ordG := Order(G);
  num := 2*(targetg-1);
  target := num/ordG - 2*g0 + 2;          // = Sum(1-1/m_i)
  if target le 0 then return []; end if;
  cl := Classes(G);
  clsfull := [ <c[1], c[2], c[3]> : c in cl ];               // includes identity (for GADims mtriv)
  clsdata := [ <c[1], c[2], c[3]> : c in cl | c[1] gt 1 ];   // drop identity (for combos)
  // sort by order ascending so pruning maxterm = last is correct
  Sort(~clsdata, func<a,b | a[1]-b[1]>);
  combos := [];
  rec(clsdata, ~combos, [], 1, target, rmax);
  if #combos eq 0 then return []; end if;
  RT := RationalCharacterTable(G);
  CT := CharacterTable(G);
  found := [];
  for idxs in combos do
    reps := [ clsdata[k][3] : k in idxs ];
    if sub<G | reps> ne G then continue; end if;
    if StructConst(CT, clsdata, ordG, idxs) le 0 then continue; end if;
    dims := GADims(G, RT, clsfull, g0, reps);
    if forall{x : x in dims | x in {0,1}} then
      orders := [ clsdata[k][1] : k in idxs ];
      Append(~found, <orders, g0>);
      printf "*** HIT genus %o: %o  signature [%o; %o] -> %o elliptic factors\n",
             targetg, tag, g0, orders, #[x:x in dims|x eq 1];
    end if;
  end for;
  return found;
end function;
