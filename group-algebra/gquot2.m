// Intermediate-cover search for genus 56 from LARGE completely-decomposable
// sources (genus 109-145), where a proper subgroup quotient can land on 56.
SetColumns(0);
load "engine.m";

function chiV(G, g0, gv)
  triv := PrincipalCharacter(G);
  reg  := PermutationCharacter(G, sub<G|>);
  cv := 2*triv + 2*(g0-1)*reg;
  for x in gv do cv +:= reg - PermutationCharacter(G, sub<G|x>); end for;
  return cv;
end function;

function FindCDGV(G, g0, orders : ntries:=600000)
  cl := Classes(G); clsfull := [<c[1],c[2],c[3]>: c in cl];
  RT := RationalCharacterTable(G);
  r := #orders;
  for t in [1..ntries] do
    gv := []; ok := true;
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
  <1728,13293, 0, [2,6,6],   145>,
  <1728,32233, 0, [2,6,6],   145>,
  <1728,46119, 0, [2,2,2,3], 145>,
  <1296, 2945, 0, [2,6,6],   109>,
  <1296, 3498, 0, [2,4,12],  109>,
  <1296, 2940, 0, [2,2,2,3], 109>
];

for s in sources do
  id1:=s[1]; id2:=s[2]; g0:=s[3]; orders:=s[4]; gexp:=s[5];
  G := SmallGroup(id1,id2);
  okgv, gv := FindCDGV(G, g0, orders);
  if not okgv then
    printf "SOURCE (%o,%o) g=%o: no completely-decomposable gv found\n", id1,id2,gexp; continue;
  end if;
  cv := chiV(G, g0, gv);
  gg := Integers()!(cv(Id(G))/2);
  printf "SOURCE (%o,%o) [%o;%o] genus=%o : scanning subgroups (this may take a bit)...\n", id1,id2,g0,orders,gg;
  subs := Subgroups(G);
  qs := {};
  for H in subs do
    gq := Integers()!((1/2)*InnerProduct(cv, PermutationCharacter(G, H`subgroup)));
    Include(~qs, gq);
    if gq eq 56 then
      printf "   *** genus-56 quotient FOUND: H order %o (index %o) -> COMPLETELY DECOMPOSABLE genus 56 ***\n",
             Order(H`subgroup), Order(G) div Order(H`subgroup);
    end if;
  end for;
  printf "   quotient genera: %o\n", Sort(Setseq(qs));
end for;
printf "DONE gquot2\n";
quit;
