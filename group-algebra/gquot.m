// Intermediate-cover (Paulhus-Rojas section 2.2) search for genus 56.
// Take known completely-decomposable source curves (group + signature, genus>56);
// every subgroup quotient X/H is then completely decomposable. Look for g(X/H)=56.
// g(X/H) = (1/2) <chi_V, Ind_H^G 1>.
SetColumns(0);
load "engine.m";

// build the H^1 character chi_V for action (G,g0,gv)
function chiV(G, g0, gv)
  triv := PrincipalCharacter(G);
  reg  := PermutationCharacter(G, sub<G|>);          // regular representation character
  cv := 2*triv + 2*(g0-1)*reg;
  for x in gv do
    cv +:= reg - PermutationCharacter(G, sub<G|x>);
  end for;
  return cv;
end function;

// find a generating vector for signature (orders) that is COMPLETELY DECOMPOSABLE
function FindCDGV(G, g0, orders : ntries:=400000)
  cl := Classes(G); clsfull := [<c[1],c[2],c[3]>: c in cl];
  RT := RationalCharacterTable(G);
  r := #orders;
  for t in [1..ntries] do
    gv := [];
    ok := true;
    for i in [1..r-1] do
      cnt:=0; x:=Id(G);
      repeat x:=Random(G); cnt+:=1; until (Order(x) eq orders[i]) or cnt gt 800;
      if Order(x) ne orders[i] then ok:=false; break; end if;
      Append(~gv,x);
    end for;
    if not ok then continue; end if;
    last := (&*gv)^-1;
    if Order(last) ne orders[r] then continue; end if;
    Append(~gv,last);
    if sub<G|gv> ne G then continue; end if;
    d := GADims(G, RT, clsfull, g0, gv);
    if forall{z: z in d | z in {0,1}} then return true, gv; end if;
  end for;
  return false, [];
end function;

sources := [
  <1344,11289, 0, [2,4,6],  57>,
  <1440, 4605, 0, [2,4,6],  61>,
  <1728,46270, 0, [2,4,6],  73>,
  <1152, 5806, 0, [2,4,8],  73>,
  < 432,  682, 0, [2,2,2,6],73>,
  <1152,157853,0, [2,4,9],  81>,
  < 432,  686, 0, [2,2,2,12],91>
];

for s in sources do
  id1:=s[1]; id2:=s[2]; g0:=s[3]; orders:=s[4]; gexp:=s[5];
  G := SmallGroup(id1,id2);
  okgv, gv := FindCDGV(G, g0, orders);
  if not okgv then
    printf "SOURCE (%o,%o) g=%o: no completely-decomposable gv found in budget\n", id1,id2,gexp;
    continue;
  end if;
  cv := chiV(G, g0, gv);
  gcheck := Integers()!(cv(Id(G))/2);
  printf "SOURCE (%o,%o) signature[%o;%o] genus=%o : enumerating subgroups...\n", id1,id2,g0,orders,gcheck;
  subs := Subgroups(G);
  hit := false;
  for H in subs do
    gq := (1/2)*InnerProduct(cv, PermutationCharacter(G, H`subgroup));
    if gq eq 56 then
      printf "   *** genus-56 quotient: H of order %o (index %o) -> COMPLETELY DECOMPOSABLE genus 56 ***\n",
             Order(H`subgroup), Order(G) div Order(H`subgroup);
      hit := true;
    end if;
  end for;
  if not hit then
    qs := Sort(Setseq({ Integers()!((1/2)*InnerProduct(cv, PermutationCharacter(G, H`subgroup))) : H in subs }));
    printf "   no genus-56 quotient. quotient genera available: %o\n", qs;
  end if;
end for;
printf "DONE gquot\n";
quit;
