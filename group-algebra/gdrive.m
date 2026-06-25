// Driver: group-algebra search for a completely decomposable genus-56 Jacobian,
// over all candidate group orders (g0 in {0,1}), Conder regime AND the
// uncovered small-group regime. Skips orders with explosively many groups.
SetColumns(0);
load "engine.m";

orders := [29,30,32,33,34,35,36,37,38,39,40,42,44,45,48,50,51,52,54,55,56,58,60,63,64,65,66,70,72,75,76,77,78,80,84,87,88,90,91,96,99,100,102,104,105,108,110,112,115,120,121,125,126,130,132,135,140,143,144,150,152,153,154,160,165,168,174,175,176,180,182,189,192,195,198,200,210,220,225,231,240,242,252,260,264,270,273,275,280,286,288,300,306,308,320,330,336,352,360,363,375,385,396,400,420,429,440,450,495,520,528,550,600,616,660,700,770,780,792,825,840,858,880,924,990,1020,1056,1100,1155,1320,1452,1650,1980];

CAP := 500;
nhit := 0;
for ord in orders do
  ng := NumberOfSmallGroups(ord);
  if ng gt CAP then
    printf "ORDER %o: SKIP (%o groups > cap)\n", ord, ng;
    continue;
  end if;
  for k in [1..ng] do
    G := SmallGroup(ord, k);
    if IsAbelian(G) then continue; end if;
    for g0 in [0,1] do
      f := SearchGroup(G, 56 : g0:=g0, rmax:=4, tag:=Sprintf("SmallGroup(%o,%o) g0=%o", ord,k,g0));
      nhit +:= #f;
    end for;
  end for;
  printf "ORDER %o done (%o groups), running hit count=%o\n", ord, ng, nhit;
end for;
printf "ALL DONE. total genus-56 completely-decomposable group-algebra hits: %o\n", nhit;
quit;
