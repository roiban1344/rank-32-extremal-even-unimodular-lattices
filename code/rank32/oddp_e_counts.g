# GAP functions of rank32/oddp_e_counts.py: e(h), and the direct count from the conjugacy classes of Aut L.
# e(h): g in C(h), n = ord(g). When the order of h is an odd prime p, the smallest odd prime dividing n is p and g^(n/p) = h.
#
# rank32/oddp_e_counts.py の GAP の関数：e(h) と、Aut L の共役類からの直接の数え上げ。
# e(h)：g in C(h)、n = ord(g)。h の位数が奇素数 p のとき、n を割る最小の奇素数が p で g^(n/p) = h。

SmallestOddPrime := function(n)
  local f;
  f := Filtered(Set(FactorsInt(n)), q -> q <> 2);
  if Length(f) = 0 then return 0; fi;
  return f[1];
end;

# count by conjugacy classes / 共役類で数える
EofH := function(C, h, p)
  local tot, cl, g, n;
  tot := 0;
  for cl in ConjugacyClasses(C) do
    g := Representative(cl);
    n := Order(g);
    if n mod p = 0 and SmallestOddPrime(n) = p and g^(n/p) = h then
      tot := tot + Size(cl);
    fi;
  od;
  return tot;
end;

# the number of elements of G whose smallest odd prime divisor of the order is p (directly), and sum |G| e(h) / |C(h)| (h over the classes of order p)
# G の中で、位数の最小の奇素因数が p である元の個数（直接）と、sum |G| e(h) / |C(h)|（h は位数 p の共役類）
CheckGroup := function(G, primes)
  local classes, res, p, direct, cl, n, reps, x, C, e, viaE, row;
  classes := ConjugacyClasses(G);
  res := [];
  for p in primes do
    direct := 0;
    for cl in classes do
      n := Order(Representative(cl));
      if n > 1 and SmallestOddPrime(n) = p then
        direct := direct + Size(cl);
      fi;
    od;
    viaE := 0;
    row := [];
    for cl in classes do
      x := Representative(cl);
      if Order(x) = p then
        C := Centralizer(G, x);
        e := EofH(C, x, p);
        viaE := viaE + Size(G) * e / Size(C);
        Add(row, [Size(C), e]);
      fi;
    od;
    Add(res, [p, direct, viaE, direct = viaE, row]);
  od;
  return res;
end;

# reduce the degree: add orbits, smallest first, until the action is faithful (on all 146880 points one lattice needed more than 3 GB)
# 次数を落とす：小さい軌道から順に足し、作用が忠実になったところで止める（146880 点のままだと、ある格子で 3 GB を超えた）
ReduceDegree := function(G)
  local orbs, pts, o, A;
  orbs := ShallowCopy(Orbits(G, MovedPoints(G)));
  SortBy(orbs, Length);
  pts := [];
  for o in orbs do
    Append(pts, o);
    A := Action(G, pts);
    if Size(A) = Size(G) then
      return A;
    fi;
  od;
  return G;
end;

# [|C(h)|, e(h), number of fixed points] per conjugacy class of elements of order p, without all conjugacy classes;
# representatives of the classes of a Sylow p-subgroup are fused in G. Fixed points are counted on all shortest vectors
# (to tell the type of the fixed lattice; 0 means f = 0).
# 全体の共役類を作らずに、位数 p の元の共役類ごとの [|C(h)|, e(h), 不動点の数] を求める。Sylow p の共役類の代表を G で融合する。
# 不動点は最小ベクトルの全部の上で数える（固定格子の型の判別用。0 なら f = 0）。
EKnown := function(G, p)
  local n, orbs, pts, o, phi, A, S, reps, c, x, res, C;
  n := NrMovedPoints(G);
  orbs := ShallowCopy(Orbits(G, MovedPoints(G)));
  SortBy(orbs, Length);
  pts := [];
  for o in orbs do
    Append(pts, o);
    phi := ActionHomomorphism(G, pts, "surjective");
    A := Image(phi);
    if Size(A) = Size(G) then
      break;
    fi;
  od;
  S := SylowSubgroup(A, p);
  reps := [];
  for c in ConjugacyClasses(S) do
    x := Representative(c);
    if Order(x) = p and ForAll(reps, y -> not IsConjugate(A, x, y)) then
      Add(reps, x);
    fi;
  od;
  res := [];
  for x in reps do
    C := Centralizer(A, x);
    Add(res, [Size(C), EofH(C, x, p), n - NrMovedPoints(PreImagesRepresentative(phi, x))]);
  od;
  return res;
end;
