# Readme

FusionA is written in [GAP](https://www.gap-system.org/) and computes the modular data and the fusion algebra of the modular categories $SU(N)_k$ or equivalently the semisimplification of the category of tilting representations of  $U_q(sl_N)$ for $q=Exp(pi i/ (N+1+k))$.

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
If we use the pair [1,2] we can construct the modular data of the three zestings as follows:
```GAP
# zesting 1 
gap> ZestA23:=Zesting_Modular_Data(2,3,[1,2],0);;
gap> ZestA23.Smatrix;
[ [ 1, 2, 2, 1, 2, 3, 2, 2, 2, 1 ], 
  [ 2, -2, -2*E(3)^2, 2*E(3), -2*E(3), 0, -2*E(3), -2*E(3)^2, -2, 2*E(3)^2 ], 
  [ 2, -2*E(3)^2, -2*E(3)^2, 2*E(3)^2, -2, 0, -2*E(3), -2, -2*E(3), 2*E(3) ], 
  [ 1, 2*E(3), 2*E(3)^2, 1, 2*E(3)^2, 3, 2*E(3), 2*E(3), 2*E(3)^2, 1 ], 
  [ 2, -2*E(3), -2, 2*E(3)^2, -2*E(3), 0, -2, -2*E(3)^2, -2*E(3)^2, 2*E(3) ], 
  [ 3, 0, 0, 3, 0, -3, 0, 0, 0, 3 ], 
  [ 2, -2*E(3), -2*E(3), 2*E(3), -2, 0, -2*E(3)^2, -2, -2*E(3)^2, 2*E(3)^2 ], 
  [ 2, -2*E(3)^2, -2, 2*E(3), -2*E(3)^2, 0, -2, -2*E(3), -2*E(3), 2*E(3)^2 ], 
  [ 2, -2, -2*E(3), 2*E(3)^2, -2*E(3)^2, 0, -2*E(3)^2, -2*E(3), -2, 2*E(3) ], 
  [ 1, 2*E(3)^2, 2*E(3), 1, 2*E(3), 3, 2*E(3)^2, 2*E(3)^2, 2*E(3), 1 ] ]
gap> ZestA23.Tmatrix;
[ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, E(3)^2, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, E(3), 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, -1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, E(3)^2, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, E(3), 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ]
# Zesting 2
gap> ZestA23:=Zesting_Modular_Data(2,3,[1,2],1);;
gap>  ZestA23.Smatrix;
[ [ 1, 2, 2, 1, 2, 3, 2, 2, 2, 1 ], 
  [ 2, -2*E(3)^2, -2, 2*E(3), -2*E(3)^2, 0, -2, -2*E(3), -2*E(3), 2*E(3)^2 ], 
  [ 2, -2, -2*E(3), 2*E(3)^2, -2*E(3)^2, 0, -2*E(3)^2, -2*E(3), -2, 2*E(3) ], 
  [ 1, 2*E(3), 2*E(3)^2, 1, 2*E(3)^2, 3, 2*E(3), 2*E(3), 2*E(3)^2, 1 ], 
  [ 2, -2*E(3)^2, -2*E(3)^2, 2*E(3)^2, -2, 0, -2*E(3), -2, -2*E(3), 2*E(3) ], 
  [ 3, 0, 0, 3, 0, -3, 0, 0, 0, 3 ], 
  [ 2, -2, -2*E(3)^2, 2*E(3), -2*E(3), 0, -2*E(3), -2*E(3)^2, -2, 2*E(3)^2 ], 
  [ 2, -2*E(3), -2*E(3), 2*E(3), -2, 0, -2*E(3)^2, -2, -2*E(3)^2, 2*E(3)^2 ], 
  [ 2, -2*E(3), -2, 2*E(3)^2, -2*E(3), 0, -2, -2*E(3)^2, -2*E(3)^2, 2*E(3) ], 
  [ 1, 2*E(3)^2, 2*E(3), 1, 2*E(3), 3, 2*E(3)^2, 2*E(3)^2, 2*E(3), 1 ] ]
gap>  ZestA23.Tmatrix;
[ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, E(3)^2, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, E(3), 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, -1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, E(3), 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, E(3)^2, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ]
# Zesting 2
gap> ZestA23:=Zesting_Modular_Data(2,3,[1,2],2);;
gap>  ZestA23.Smatrix;
[ [ 1, 2, 2, 1, 2, 3, 2, 2, 2, 1 ], 
  [ 2, -2*E(3), -2*E(3), 2*E(3), -2, 0, -2*E(3)^2, -2, -2*E(3)^2, 2*E(3)^2 ], 
  [ 2, -2*E(3), -2, 2*E(3)^2, -2*E(3), 0, -2, -2*E(3)^2, -2*E(3)^2, 2*E(3) ], 
  [ 1, 2*E(3), 2*E(3)^2, 1, 2*E(3)^2, 3, 2*E(3), 2*E(3), 2*E(3)^2, 1 ], 
  [ 2, -2, -2*E(3), 2*E(3)^2, -2*E(3)^2, 0, -2*E(3)^2, -2*E(3), -2, 2*E(3) ], 
  [ 3, 0, 0, 3, 0, -3, 0, 0, 0, 3 ], 
  [ 2, -2*E(3)^2, -2, 2*E(3), -2*E(3)^2, 0, -2, -2*E(3), -2*E(3), 2*E(3)^2 ], 
  [ 2, -2, -2*E(3)^2, 2*E(3), -2*E(3), 0, -2*E(3), -2*E(3)^2, -2, 2*E(3)^2 ], 
  [ 2, -2*E(3)^2, -2*E(3)^2, 2*E(3)^2, -2, 0, -2*E(3), -2, -2*E(3), 2*E(3) ], 
  [ 1, 2*E(3)^2, 2*E(3), 1, 2*E(3), 3, 2*E(3)^2, 2*E(3)^2, 2*E(3), 1 ] ]
gap>  ZestA23.Tmatrix;
[ [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 0, E(3), 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, E(3)^2, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, -1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0, 0, E(3)^2, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, E(3), 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ]
```
The next function test that in fact the pair of matrices is a modar data. (It test: (1) that normalized S and T are unitary, (2) S_{0a}>0, (3) (ST)^3=Gauus_Sum(S,T)*S^2 )

```GAP
IsModularData(ZestA23.Smatrix,ZestA23.Tmatrix);
[ true, true, true ]
```
We can construct the fusion algebra using:

```GAP
gap> Zesting_FusionAlgebra(2,3,[1,2,0]);; # record of fusion algebra and canonical basis
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
As a bigger example let us compute the braided zestings of SU(5)_5

```GAP
gap> Display_Info_Zesting(4,5);

The generator of the transparent objects of PGL(5,5)
is the power 1 of the simple object v.126
and its twist is 1
The pairs for computing braided zestings are
[ 0, 0 ]
[ 1, 3 ]
[ 2, 1 ]
[ 3, 4 ]
[ 4, 2 ]

gap> zA45:=Zesting_Modular_Data(4,5,[1,3]);;time;
8039
gap> IsModularData(zA45.Smatrix,zA45.Tmatrix);
[ true, true, true ]
```
The fusion algebra 

```GAP
gap> zA45:=Zesting_FusionAlgebra(4,5,[1,3],0);;time;
779499
gap> v:=zA45.Simples;;
```
We can see again that the exponent of v.2 changed from 5 to 
15

```GAP
gap> v.2^5;
v.6+(4)*v.10+v.21+(4)*v.35+(5)*v.43+(6)*v.68+(5)*v.75
gap> v.2^7;
(15)*v.15+(20)*v.17+(14)*v.26+(35)*v.29+(14)*v.46+(14)*v.51+(14)*v.60+(
15)*v.71+(35)*v.80+(21)*v.84+(21)*v.98
gap> v.2^8;
(14)*v.4+(35)*v.18+(35)*v.19+(64)*v.30+(90)*v.32+(56)*v.39+(28)*v.52+(
14)*v.54+(28)*v.61+(70)*v.64+(64)*v.81+(70)*v.86+(56)*v.100+(42)*v.104
gap> v.2^9;
(42)*v.5+(84)*v.9+(70)*v.20+(189)*v.33+(189)*v.34+(120)*v.40+(216)*v.42+(
42)*v.55+(162)*v.65+(216)*v.67+(168)*v.74+(162)*v.87+(84)*v.89+(120)*v.101+(
168)*v.106+(42)*v.114
gap> v.2^10;
(210)*v.13+(42)*v.21+(210)*v.27+(288)*v.35+(300)*v.43+(567)*v.44+(525)*v.48+(
70)*v.56+(252)*v.68+(450)*v.69+(768)*v.77+(450)*v.82+(448)*v.90+(252)*v.96+(
525)*v.107+(567)*v.109+(300)*v.116+(288)*v.118+(42)*v.126
gap> v.2^11;
(462)*v.14+(660)*v.16+(1188)*v.28+(330)*v.36+(660)*v.37+(1155)*v.45+(
825)*v.49+(1540)*v.50+(462)*v.62+(990)*v.70+(1320)*v.78+(2310)*v.79+(
2310)*v.83+(518)*v.91+(1320)*v.97+(990)*v.102+(1540)*v.110+(825)*v.117+(
1155)*v.119+(330)*v.122
gap> v.2^12;
(462)*v.7+(462)*v.15+(2112)*v.17+(2970)*v.29+(4158)*v.31+(4158)*v.38+(
1485)*v.46+(3520)*v.51+(2058)*v.53+(2970)*v.63+(1320)*v.71+(2112)*v.72+(
5775)*v.80+(4455)*v.84+(7700)*v.85+(2640)*v.98+(4455)*v.99+(5775)*v.103+(
2058)*v.111+(1320)*v.112+(3520)*v.120+(1485)*v.123
gap> v.2^13;
(3432)*v.8+(2574)*v.18+(3432)*v.19+(2574)*v.22+(3432)*v.30+(15015)*v.32+(
11583)*v.39+(16016)*v.41+(5005)*v.52+(5578)*v.54+(8580)*v.64+(11583)*v.66+(
15015)*v.73+(8580)*v.81+(21450)*v.86+(11816)*v.88+(3432)*v.92+(12870)*v.100+(
12870)*v.104+(21450)*v.105+(8580)*v.113+(5578)*v.121+(5005)*v.124
gap> v.2^14;
(12012)*v.9+(15015)*v.12+(6006)*v.20+(21021)*v.23+(21021)*v.33+(27027)*v.34+(
15015)*v.40+(64064)*v.42+(27832)*v.47+(10583)*v.55+(6006)*v.57+(12012)*v.65+(
48048)*v.67+(48048)*v.74+(64064)*v.76+(35035)*v.87+(38844)*v.89+(27027)*v.93+(
21450)*v.101+(68640)*v.106+(38844)*v.108+(21450)*v.114+(35035)*v.115+(
10583)*v.125
gap> v.2^15;
(10583)*v.1+(6006)*v.6+(54054)*v.10+(96525)*v.13+(126126)*v.24+(125125)*v.27+(
24024)*v.35+(75075)*v.43+(81081)*v.44+(96525)*v.48+(6006)*v.56+(84462)*v.58+(
81081)*v.68+(75075)*v.69+(100100)*v.75+(292864)*v.77+(125125)*v.82+(
54054)*v.90+(130740)*v.94+(146328)*v.96+(100100)*v.107+(126126)*v.109+(
130740)*v.116+(84462)*v.118+(10583)*v.126
```


# Frobenius-Schur Indicator

The following example shows how to compute the Frobenius-Schur Indicator of an object

```GAP
gap> Display_labels(2,3);

L.   Grnd.   Twist
(1)   [ 0, 0 ]   0
(2)   [ 1, 0 ]   2/9
(3)   [ 2, 0 ]   5/9
(4)   [ 3, 0 ]   1
(5)   [ 0, 1 ]   2/9
(6)   [ 1, 1 ]   1/2
(7)   [ 2, 1 ]   8/9
(8)   [ 0, 2 ]   5/9
(9)   [ 1, 2 ]   8/9
(10)   [ 0, 3 ]   1
```

Now we can compute the list of all n-th FS indicators ( in the following example for SU(3,3) )

```GAP
gap> FS_indicator(2,3,2);
[ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0 ]
gap> FS_indicator(2,3,3);
[ 1, 1, 1, 1, 1, 2, 1, 1, 1, 1 ]
gap> FS_indicator(2,3,4);
[ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0 ]
```

# TODO:
- Extend the program to fractional leveles.
- Extend to all simple lie algebras.
