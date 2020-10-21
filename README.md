# Readme

FusionA is written in [GAP](https://www.gap-system.org/) and computes the modular data and the fusion algebra of the modular categories $SU(N)_k$ or equivalently the semisimplification of tilting representation of  $U_q(sl_N)$ for $q$ a primitive root of order $2(N+1+k)$.

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

The modular data is computed as follows: (Recall that $E(n)$ means the complex number $Exp(2 \pi i/n)$ )

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

## Fusion algebra

The following example shows how to compute the fusion algebra of the Ising using FusionA: 

```GAP
gap> Ising:=FusionAlgebra(1,2);; # Ising corresponds to SU(2)_2
gap> v:=Ising.Simples;; # v will be a register with the simples in the same orden obstained with the function Display_labels
gap> Display_labels(1,2);

L.   Grnd.   Twist
(1)   [ 0 ]   0
(2)   [ 1 ]   3/16
(3)   [ 2 ]   1/2

gap> v.2*v.2;
v.1+v.3
gap> v.3*v.3;
v.1
gap> v.2*v.3;
v.2

```
As a bigger example let us compute some fusion products of $SU(5,5)$ that has 126 simples. 

```GAP
Size(Labels_A(4,5)); # Number of simple of SU(5,5)
126
gap> A45:=FusionAlgebra(4,5);;time; # it takes some time
814672
gap> v:=A45.Simples;; #record with the simples
```
In general, v.1 is the unit object, v.2 is the standard the representation and v.126 (last label) is the generator of the invertibles.


```GAP
# the first power of v.2 where v.1 appears should be 5
gap> v.2^2;
v.3+v.7
gap> v.2^3;
v.4+(2)*v.8+v.22
gap> v.2^4;
v.5+(3)*v.9+(2)*v.12+(3)*v.23+v.57
gap> v.2^5;
v.1+v.6+(4)*v.10+(5)*v.13+(6)*v.24+(5)*v.27+(4)*v.58
# v.126 should generates a cyclic group of order 5
gap> v.126^2;
v.56
gap> v.126^3;
v.21
gap> v.126^4;
v.6
gap> v.126^5;
v.1
```
## Zesting SU(N)_k

FusionA allows to zest the fusion rules and modular data of SU(N)_k

The next fusion provides some information of the pairs that you can use to zest SU(N)_k

```GAP
gap> Display_Info_Zesting(2,3);

The generator of the transparent objects of PGL(3,3)
is the power 1 of the simple object v.10
and its twist is 1
The pairs for computing braided zestings are
[ 0, 0 ]
[ 1, 2 ]
[ 2, 1 ]
gap> 
```
If we use the pair [1,2] we can construct the modular data  as follows:
```GAP
gap> ZestA23:=Zesting_Modular_Data(2,3,[1,2]);;
gap> ZestA23.Smatrix;
[ [ 1, 2, 2, 1, 2, 3, 2, 2, 2, 1 ], 
  [ 2, -2*E(9)^7, -2*E(9)^2, 2*E(3), 2*E(9)^2+2*E(9)^5, 0, 2*E(9)^4+2*E(9)^7, 
      -2*E(9)^4, -2*E(9)^5, 2*E(3)^2 ], 
  [ 2, -2*E(9)^2, -2*E(9)^7, 2*E(3)^2, 2*E(9)^4+2*E(9)^7, 0, 
      2*E(9)^2+2*E(9)^5, -2*E(9)^5, -2*E(9)^4, 2*E(3) ], 
  [ 1, 2*E(3), 2*E(3)^2, 1, 2*E(3)^2, 3, 2*E(3), 2*E(3), 2*E(3)^2, 1 ], 
  [ 2, 2*E(9)^2+2*E(9)^5, 2*E(9)^4+2*E(9)^7, 2*E(3)^2, -2*E(9)^4, 0, 
      -2*E(9)^5, -2*E(9)^2, -2*E(9)^7, 2*E(3) ], 
  [ 3, 0, 0, 3, 0, -3, 0, 0, 0, 3 ], 
  [ 2, 2*E(9)^4+2*E(9)^7, 2*E(9)^2+2*E(9)^5, 2*E(3), -2*E(9)^5, 0, -2*E(9)^4, 
      -2*E(9)^7, -2*E(9)^2, 2*E(3)^2 ], 
  [ 2, -2*E(9)^4, -2*E(9)^5, 2*E(3), -2*E(9)^2, 0, -2*E(9)^7, 
      2*E(9)^4+2*E(9)^7, 2*E(9)^2+2*E(9)^5, 2*E(3)^2 ], 
  [ 2, -2*E(9)^5, -2*E(9)^4, 2*E(3)^2, -2*E(9)^7, 0, -2*E(9)^2, 
      2*E(9)^2+2*E(9)^5, 2*E(9)^4+2*E(9)^7, 2*E(3) ], 
  [ 1, 2*E(3)^2, 2*E(3), 1, 2*E(3), 3, 2*E(3)^2, 2*E(3)^2, 2*E(3), 1 ] ]
gap> ZestA23.Tmatrix;
[ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, -E(9)^4-E(9)^7, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, -E(9)^4-E(9)^7, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, E(9)^7, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, -1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, E(9)^7, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, E(9)^4, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, E(9)^4, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ]
```
We can construct the fusion algebra using:

```GAP
gap> Zesting_FusionAlgebra(2,3,[1,2]);; # record of fusion algebra and canonical basis
gap> v:=last.Simples;; # stracting the basis
```
Let us check that the fusion algebra of SU(3)_3 and its zesting are non-isomorphic. For this we will see that the exponent of v.2 (the standard representation) in the zesting has order 9 and not 3

```GAP
gap> v.2^3;
v.4+(2)*v.6+v.10
gap> v.2^4;
(2)*v.2+(3)*v.7+(3)*v.8
gap> v.2^5;
(5)*v.3+(5)*v.5+(6)*v.9
gap> v.2^6;
(6)*v.1+(5)*v.4+(16)*v.6+(5)*v.10
gap> v.2^7;
(22)*v.2+(21)*v.7+(21)*v.8
gap> v.2^8;
(43)*v.3+(43)*v.5+(42)*v.9
gap> v.2^9;
(42)*v.1+(43)*v.4+(128)*v.6+(43)*v.10
```

