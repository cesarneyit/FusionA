# Readme

FusionA is written in [GAP](https://www.gap-system.org/) and computes the modular data and the fusion algebra of the modular categories $SU(N)_k$ or equivalently the semisimplification of $U_q(sl_N)$ for $q$ a primitive root of order $2(N+1+k)$.

## Using FusionA
The simple objects of $SU(N)_k$ are displayed as N-tuples of non-negative integers. The following function displays information on labels and twists.
```GAP
gap> Display_labels(4,3);

L.   Grnd.   Twist
(1)   [ 0, 0, 0, 0 ]   0
(2)   [ 1, 0, 0, 0 ]   3/10
(3)   [ 2, 0, 0, 0 ]   7/10
(4)   [ 3, 0, 0, 0 ]   6/5
(5)   [ 0, 1, 0, 0 ]   9/20
(6)   [ 1, 1, 0, 0 ]   33/40
(7)   [ 2, 1, 0, 0 ]   13/10
(8)   [ 0, 2, 0, 0 ]   21/20
(9)   [ 1, 2, 0, 0 ]   3/2
(10)   [ 0, 3, 0, 0 ]   9/5
(11)   [ 0, 0, 1, 0 ]   9/20
(12)   [ 1, 0, 1, 0 ]   4/5
(13)   [ 2, 0, 1, 0 ]   5/4
(14)   [ 0, 1, 1, 0 ]   1
(15)   [ 1, 1, 1, 0 ]   57/40
(16)   [ 0, 2, 1, 0 ]   17/10
(17)   [ 0, 0, 2, 0 ]   21/20
(18)   [ 1, 0, 2, 0 ]   29/20
(19)   [ 0, 1, 2, 0 ]   17/10
(20)   [ 0, 0, 3, 0 ]   9/5
(21)   [ 0, 0, 0, 1 ]   3/10
(22)   [ 1, 0, 0, 1 ]   5/8
(23)   [ 2, 0, 0, 1 ]   21/20
(24)   [ 0, 1, 0, 1 ]   4/5
(25)   [ 1, 1, 0, 1 ]   6/5
(26)   [ 0, 2, 0, 1 ]   29/20
(27)   [ 0, 0, 1, 1 ]   33/40
(28)   [ 1, 0, 1, 1 ]   6/5
(29)   [ 0, 1, 1, 1 ]   57/40
(30)   [ 0, 0, 2, 1 ]   3/2
(31)   [ 0, 0, 0, 2 ]   7/10
(32)   [ 1, 0, 0, 2 ]   21/20
(33)   [ 0, 1, 0, 2 ]   5/4
(34)   [ 0, 0, 1, 2 ]   13/10
(35)   [ 0, 0, 0, 3 ]   6/5
gap> time;
5
```

The modular data is computed as follows:
```GAP
gap> A12:=ModularData_TipoA(1,2);;
gap> S:=A12.Smatrix;;
gap> T:=A12.Tmatrix;;
gap> Display(S);
[ [             1,   E(8)-E(8)^3,             1 ],
  [   E(8)-E(8)^3,             0,  -E(8)+E(8)^3 ],
  [             1,  -E(8)+E(8)^3,             1 ] ]
gap> Display(T);
[ [        1,        0,        0 ],
  [        0,  E(16)^3,        0 ],
  [        0,        0,       -1 ] ]
gap> A45:=ModularData_TipoA(4,5);; # 126 simples
gap> time;
8521
```
