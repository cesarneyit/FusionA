




#####
# Calcula la matriz de la forma bilineal en el producto, 
# de tal forma que (r,r)=2 para las raices mas bajas
########

BilinearFormMat_Normalized:= function(R)
# R es root system
local B, D, d;
B:=BilinearFormMat(R);
D:=ShallowCopy((DiagonalOfMat(B))); # matriz diagonal / usando ShallowCopy para no cambiar B
SortBy(D, x -> x); # ordenandola para encontrar el valor de la raiz mas pequeña
d:=D[1]; 
return 2*d^(-1)*B;
end;




################
# Calcula el producto interno normalizado 
# R es el sistema de raices y wt1 y wt2 son dos pesos
#########
Inner_Product := function (R,wt1,wt2) 
local C, B, B_P;
C:=SimpleSystem(R)^(-1); 
B:=BilinearFormMat_Normalized(R);
B_P:=TransposedMat(C)*B*C;
return wt1*B_P*wt2;
end;

######
# Matriz del producto interno normalizado
########
MBP:=function(R)
local B,D,d,C,B_P;
C:=SimpleSystem(R)^(-1);
B:=BilinearFormMat_Normalized(R);
B_P:=TransposedMat(C)*B*C;
return B_P;
end;


######
# producto interno mejorado
#####
Inner_Product_m := function (B_P,wt1,wt2) 
return wt1*B_P*wt2;
end;


#######
# Función recursiva que calcula los simples en el caso SU(N)_k
# Es decir suma de N enteros no negativos cuya suma es k
####
Labels_SU:=function(length, total)
local lista, i,t,Length;
lista:=[];
if length=1 then 
        return [[total]];
else
    for i in [0..total] 
        do 
            for t in Labels_SU(length-1,total-i)
                do
                    Append(t,[i]);
                    Append(lista,[t]);
                od;
        od;
fi;
return lista;
end;


#######
# función recursiva que calcula los simples en el caso A_n
# Todo??: crearlo como un iterador 
####
Labels_A:=function(length, total)
local lista, i;
lista:=Labels_SU(length+1,total);
for i in lista
    do
        Remove(i,1);
    od; 
return lista;
end;








######
# action ciclica
# 
########
act_cyclic:=function(wt,g)
return Permuted( wt, CycleFromList([1..Size(wt)]) );
end;
######
##
## Calcula numero de orbitas
###
Numero_orbitas:=function(n,k)
local labels,G,orb;
labels:=Labels_A_ext(n,k);
G:=Group(CycleFromList([1..n+1]));
orb:=Orbits(G,labels,act_cyclic);
return [Size(Labels_A_ext(n,k)),Size(orb),Size(orb)/Size(Labels_A_ext(n,k))];
end;


########
# contruye el peso rho para A_n
# en este caso es una fila de tamaño n de 1's
#
#########
rho_max:=function(n)
local rho,i;
rho:=[];
for i in [1..n]
    do
        Append(rho,[1]);
    od;
return rho;
end;



######
# 
# Calcula el numerador en la entrada S_ij 
# Una parte importante y confusa es el 
# significado de q^a con a un racional
# este es si: a=p/q y  q=E(2*x), entonces q^a= E(2*x*q)^p
########

S_ij:=function(R,k,wt1,wt2)
local x,K, n ,w, W,rho,a;
K:=0;
n:=Size(SimpleSystem(R)); # rango de R
x:=k+n+1; # nivel 
W:= WeylGroup( R );
rho:= rho_max(n);
for w in W do
a:=Inner_Product(R,wt1+rho ,(rho+wt2)*w );
K:=K+DeterminantMat(w)*E(2*x*DenominatorRat(a))^(2*NumeratorRat(a));
od;
return K;
end;




S_ij_m:=function(R,k,wt1,wt2)
local W, K,x,w,rho,B_P;
K:=0;
n:=Size(SimpleSystem(R)); # rango de R
x:=k+n+1; # nivel
rho:=rho_max(n);
W:= WeylGroup( R ); # grupo de Weyl que es  S_n
B_P:=MBP(R);;
for w in W 
    do
        a:=Inner_Product_m(B_P,wt1+rho ,(rho+wt2)*w );
        K:=K+DeterminantMat(w)*E(2*x*DenominatorRat(a))^(2*NumeratorRat(a));
    od;
return K;
end;

SS:=function(n,k)
local R,wt;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
wt:=0*rho_max(n);;
return S_ij_m(R,k,wt,wt);
end;


#############
# S-matriz (CONJUGADA RESPECTO A LA FORMULA DE BAKI)
# no rmalizada en el sentido que S_00=1
#TODO: CAMBIAR EL NOMBRE A 
###################

SMatrix:=function(R,k,labels)
local x,y,S,r;
r:=Size(labels);
S:=NullMat(r,r); 
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij_m(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return ComplexConjugate( (1/S[1][1])*S);
end;
#####
# S-matrix mejorada
#####
SMatrix_m:=function(R,k,labels)
local x,y,S,r;
r:=Size(labels);
S:=NullMat(r,r); 
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij_m(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return  (1/S[1][1])*S;
end;

SMatrix2:=function(R,k,labels)
local x,y,S,r;
r:=Size(labels);
S:=NullMat(r,r); 
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij2(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return ComplexConjugate( (1/S[1][1])*S);
end;

#####
# Antigua function 
###

SMatrix3:=function(n,k)
local x,y,S,r,R,labels;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
labels:=Labels_A(n,k);;
r:=Size(labels);
S:=NullMat(r,r); # identidad pero como matriz immutable
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return ComplexConjugate( (1/S[1][1])*S);
end;

#####
# Nueva function 
###

SMatrix4:=function(n,k)
local x,y,S,r,R,labels;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
labels:=Labels_A(n,k);;
r:=Size(labels);
S:=NullMat(r,r); # identidad pero como matriz immutable
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij2(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return ComplexConjugate( (1/S[1][1])*S);
end;


SMatrix5:=function(n,k)
local x,y,S,r,R,labels;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
labels:=Labels_A(n,k);;
r:=Size(labels);
S:=NullMat(r,r); # identidad pero como matriz immutable
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij3(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return ComplexConjugate( (1/S[1][1])*S);
end;


SMatrix6:=function(n,k)
local x,y,S,r,R,labels;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
labels:=Labels_A(n,k);;
r:=Size(labels);
S:=NullMat(r,r); # identidad pero como matriz immutable
for x in [1..r] do
    for y in [x..r] do
        S[x][y]:=S_ij4(R,k,labels[x],labels[y]);
        S[y][x]:=S[x][y];
    od;
od;
return ComplexConjugate( (1/S[1][1])*S);
end;


#########
# S.matrix normalizada
# Normalizada en el sentido que S^2=C
# con C_ij=\delta_{i,j^*}
##########

SMatrix_Normalized:=function(R,k,labels)
local D,S;
S:=SMatrix(R,k,labels);
D:=Sqrt(S[1]*S[1]);
return D^(-1)*S;
end;


#############
#
# Entrada ii como racional
##########

Twist_Racional:=function(R,k,wt)
local x,rho,a,theta,n;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
rho:= DiagonalOfMat(IdentityMat(n)); 
a:=Inner_Product(R,wt ,wt+2*rho )*(2*x)^(-1);
return a;
end;



#############
#
# Entrada ii de la T-matrix
##########

Twist:=function(R,k,wt)
local x,rho,a,theta,n,B_P;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
rho:= rho_max(n);
B_P:=MBP(R);
a:=Inner_Product_m(B_P,wt ,wt+2*rho );
theta:=E(2*x*DenominatorRat(a))^NumeratorRat(a);
return theta;
end;

########
# T-matrix
# 
########


Tmatrix:=function(R,k,labels)
local x,rho,a,theta,n,T,r;
r:=Size(labels);
T:=IdentityMat(r);
for x in [1..r] do
    T[x][x]:=Twist(R,k,labels[x]);
od;
return T;
end;

########
# T-matrix como enteros
########


Tmatriz_Racional:=function(n,k,labels)
local x,rho,a,theta,T,r,R, L;
r:=Size(labels);
L:= SimpleLieAlgebra( "A", n, Rationals );; R:= RootSystem( L );;
T:=IdentityMat(r);
for x in [1..r] do
    T[x][x]:=Twist_Racional(R,k,labels[x]);
od;
return T;
end;

#########
# Calcula la
# Multiplicative and additive Central charges
# el output es un par de c la aditiva y la cmultiplicativa
#########

CentralCharge:=function(R,k)
local c,n,x;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
c:= k*(n+1)*x^(-1); # additive central charge
return [c,E(12*DenominatorRat(c))^NumeratorRat(c)];
end;


####
# Transpuesta conjugado
#########
DagaMatrix:=function(U)
return TransposedMat( ComplexConjugate(U) );
end;

########
# comprueba si U es unitaria
#########

IsUnitaryMatrix:=function(U)
local v,n;
n:=Size(U);
v:=false;
    if U*DagaMatrix(U)=IdentityMat(n) then 
        v:=true;
    fi;
return v;
end;



#########
# Suma de Gauss de un dato modular
#  S y T es el dato modular
# signo debe ser 1 o -1
#############

Gauss_Sum:=function(S,T,signo)
local k,s,x;
s:=Size(S);
k:=0;
for x in [1..s] do
    k:=k+(S[1][x])^2*T[x][x]^(signo);
od;
return k;
end;



#######
# Calcula dato modular de A_n en nivel principal k
#  En otras palabras el dato modular de SU(N+1)_k
#########
ModularData_TipoA:=function(n,k)
local R,labels,S,T,Labels, i;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
labels:=Labels_A(n,k);
S:=SMatrix(R,k,labels);
T:=Tmatrix(R,k,labels);
Labels:=rec();
for i in [1..Size(labels)] do
        Labels.(i):=labels[i];
    od;
return rec(Smatrix:=S,Tmatrix:=T,labels:=Labels);
end;


######
# Nueva probando
#####
ModularData_TipoA2:=function(n,k)
local R,labels,S,T,Labels, i;
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
labels:=Labels_A(n,k);
S:=SMatrix2(R,k,labels);
T:=Tmatrix(R,k,labels);
Labels:=rec();
for i in [1..Size(labels)] do
        Labels.(i):=labels[i];
    od;
return rec(Smatrix:=S,Tmatrix:=T,labels:=Labels);
end;



#######
# Velinde N^z_{xy}
# Calcula la formula Dim(C)*N^z_{xy}=\Sum_w \frac{S_{xw}S_{yw}\overline{S_{zw}}}{S_{0w}}
#######
Velinde_xyz:=function(S,x,y,z)
local k,w,n,D;
D:=S[1]*S[1];
n:=Size(S);
k:=0;
for w in [1..n] do
    k:=k+S[x][w]*S[y][w]*ComplexConjugate(S[z][w])*S[1][w]^(-1);
od;
return D^(-1)*k;
end;
######
# Formula de verlinde 
# Dada S constuye una lista
# N cuya entrada i-esima contien
# la matriz (N_i)_{jk}=N_{i,j}^k
#######
Verlinde_formula_mala:=function(S)
local N,tmp,x,y,z,n;
n:=Size(S);
N:=[];
tmp:=NullMat(n,n);
for x in [1..n] do
    for y in [1..n] do
        for z in [1..n] do
            tmp[y][z]:=Velinde_xyz(S,x,y,z);
        od;
    od;
    Append(N,[tmp]);
od;
return N;
end;


######
# Formula de verlinde MEJORADA
# Dada S constuye una lista
# N cuya entrada i-esima contien
# la matriz (N_i)_{jk}=N_{i,j}^k
#######
Verlinde_formula:=function(S)
local N,tmp,x,y,z,n,D;
n:=Size(S);
N:=[];
D:=S[1]*S[1];
tmp:=[];
for x in [1..n] do
    tmp:=DiagonalMat(S[1])^(-1)*DiagonalMat(S[x]); # matrix diagonal  Gamma_a, dada por  s_{aj}/s_{0j}
    Append(N,[D^(-1)*S*tmp*ComplexConjugate(S)]);
od;
return N;
end;



######
# Formula de verlinde MEJORADA
# Dada S y un label i construye
# N_i cuya entrada i-esima contiene
# la matriz (N_i)_{jk}=N_{i,j}^k
#######
Verlinde_N:=function(S,x)
local N,tmp,y,z,n,D;
n:=Size(S);
N:=[];
D:=S[1]*S[1];
tmp:=DiagonalMat(S[1])^(-1)*DiagonalMat(S[x]); # matrix diagonal  Gamma_a, dada por  s_{aj}/s_{0j}
return D^(-1)*S*tmp*ComplexConjugate(S);
end;

#######ç
# revisa si el par es un dato modular
#####
IsModularData:=function(S,T)
local v1,v2,v3,D;
# revisa unitariedad salvo por la dimesion
if S*ComplexConjugate(S)=(S[1]*S[1])*IdentityMat(Size(S)) then 
    v1:=true;
    else v1:=false;
fi;
# revisa que S_{0a}>0 para todo a
D:=ShallowCopy(S[1]);
Sort(D);
if D[1]>0 then
    v2:=true;
    else v2:=false;
fi;
# revisa si (ST)^3=Gauus_Sum(S,T)*S^2
if (S*T)^3=Gauss_Sum(S,T,1)*S^2 then 
    v3:=true;
    else v3:=false;
fi;
return [v1,v2,v3];
end;

######
#
# Test de modular
########
Test_modular:=function(n,k)
local MD,test;
MD:=ModularData_TipoA(n,k);
test:=IsModularData(MD.Smatrix,MD.Tmatrix);
return test;
end;

######
#
# Test de verlinde
########
Test_verlinde:=function(n,k)
local MD,test;
MD:=ModularData_TipoA(n,k);
return Verlinde_formula(MD.Smatrix);
end;


########
#
# Inveritble
# Calcula los labels de los invertibles y N_a
######

invertibles:=function(MD)
local n,x,L,S,labels,tmp,N;
S:=MD[1];
labels:=MD[3];
n:=Size(S);
L:=[];
N:=[];
tmp:=[];
for x in [1..n] do 
    if S[1][x]=1 then 
        Append(N,[Verlinde_N(S,x)]);
        Append(L,[labels[x]]);
    fi;
od;
return [N,L];
end;

######
# Divide labels en n+1 sublistas de los gradings Z/n+1
# la primera corresponde a E(n+1), la segunda a E(n+1)^2 ...
#########

Grading_A:=function(MD)
local S,labels,k,zeta,labels_graded,ftw,n,x,w,tmp,KK,raices,r;
S:=MD[1];
labels:=MD[3];
KK:=[]; #borrar
k:=MD[4]; # nivel
labels_graded:=[];
n:=Size(labels[1]); # dimension del sistema de raices
raices:=[];
r:=Size(labels); #  numero de labels
#ftw:=0*labels[1]; # peso cero
#ftw[1]:=k; # primer peso fundamental multiplcado por k (corresponde a invertible generador)
#w:=Position( labels, ftw ); # posicion de peso fundamental
 
# lista con n+1 sublistas
for x in [1..n+1] do
    Append(labels_graded,[[]]);
od;
###
# Creo lista con las raices de la unidad de orden n+1
# las cuales usaré para ordenar
####
for x in [1..n+1] do
    Append(raices,[E(n+1)^x]);
od;
###
# Organiza los pesos en labels_graded dependiendo S_{w,x}/S_{1,x}
###
for x in [1..r] do 
    tmp:=Position(raices, S[r][x]*(S[1][x])^(-1));
    Append(labels_graded[tmp],[labels[x]]);
od;
#return labes_graded;
return labels_graded;
end;


#####
# Funcion de grado
# 
#####

grado_label:=function(wt)
local x, r,k;
r:=Size(wt);
k:=0;
#if wt<>0*wt  then 
    for x in [1..r] 
        do
            k:=k+x*wt[x];
        od;
#fi;
return RemInt(k, r+1);
end;


#####
# PArticion associada un peso
# 
###
Particion_wt:=function(wt)
local r,x,P,tmp,PP;
r:=Size(wt);
P:=[wt[1]];
# crea la particion de la forma a1,a1+a2,a1+a2+a3
for x in [1..r-1] do
    Append(P,[P[x]+wt[x+1]]);
od;
## cambio el orden para que quede como una particion normal
PP:=[];
for x in [1..r] do
    Append(PP,[P[r+1-x]]);
od;
return PP;
end;

##
# Cambio de base a particiones
# no muy util
###
Cambio_base_particiones:=function(wt)
local n,L,R,M,x,C;
n:=Size(wt);
L:= SimpleLieAlgebra( "A", n, Rationals );; R:= RootSystem( L );;
M:=NullMat(n,n+1);
for x in [1..n] do
    M[x][x]:=1;
    M[x][x+1]:=-1;
od;
C:=SimpleSystem(R);
return wt*C^(-1)*M;
end;

##
# a particiones
# 
###
Cambio_base_particiones:=function(wt)
local n,L,R,M,x,C;
n:=Size(wt);
#L:= SimpleLieAlgebra( "A", n, Rationals );; R:= RootSystem( L );;
M:=NullMat(n,n);
M[1][1]:=1;
for x in [1..n-1] do
    M[x+1][x+1]:=1;
    M[x+1]:=M[x]+M[x+1];
od;
return wt*M;
end;

###
#  Grado
# Calcula la altura de la particion asociada al peso wt
####
ht:=function(wt)
return Sum(Cambio_base_particiones(wt)) mod Size(wt);
end;


#######
#######
# Canonical zesting, se debe verificar que (N,k)\neq 1
#
#######
######
Indice_modularidad:=function(n,k)
local gcd,v,h;
h:=Gcd(n,k);
v:=0;
if h mod 2 =0 and n/h mod 2 =1 and k/h mod 2 =1 then 
    v:=1;
fi;
return v;
end;
######
# Dada o si no hay zesting
#####
Obstruccion:=function(n,k)
local m, epsilon,v,x;
m:=n/Gcd(n,k);
epsilon:=Indice_modularidad(n,k);
v:=0;
for x in [1..n] do 
    if 2*x+epsilon=n/m mod n then
        v:=1;
    fi;
od;
return v;
end;
########
zesting:=function(n,k)
local gcd, MD, Tz, x, labels, T;
gcd:=Gcd(n,k);
MD:=ModularData_TipoA(n,k);
labels:=MD[3];
T:=MD[2];
Tz:=ShallowCopy(T);
for x in [1..n] do
    Tz[x][x]:=E(2)*T[x][x];
od;
end;

Display_labels:=function(n,k)
local labels, x,T,R;
labels:=Labels_A(n,k);
R:= RootSystem( SimpleLieAlgebra( "A", n, Rationals ) );;
Print("\n");
Print("L.", "   " , "Grnd.", "   ",  "Twist", "\n");
for x in [1..Size(labels)] do
    Print("(",x,")", "   " , labels[x],"   " ,Twist_Racional(R,k,labels[x]),  "\n");
od;
end;

######################
# Funcion usada para entrar 
# los coeficientes de estructura en la funcion StructureConstantsTable
# Nij es una fila con Nij[k]=N_{ij}^k
############
FormatoTabla:=function(Nij)
local d,x,posicion,coeficientes;
d:=Size(Nij);
posicion:=[];
coeficientes:=[];
for x in [1..d] 
    do
        if  Nij[x] <>  0 then
            Append(posicion,[x]);
            Append(coeficientes,[Nij[x]]);
        fi; 
    od;
return [posicion,coeficientes];
end;

FusionAlgebra:=function(n,k)
local T,x,y,d,Alg,base,Base,N;
N:=Test_verlinde(n,k);; # lista de matrices (N_i)_{jk}=N_{ij}^k
d:=Size(N);; # dimension del algebra de fusion
T:=NullMat(d,d);;
for x in [1..d]
    do
        for y in [1..d]
            do 
                T[x][y]:=FormatoTabla(N[x][y]);
            od;
    od;
Append(T,[1,0]); # datos adicionales pedidos por GAP
Alg:=AlgebraByStructureConstants( Rationals, T );;
base:=BasisVectors( CanonicalBasis( Alg ) );;
Base:=rec();
for x in [1..Size(base)] 
    do
        Base.(x):=base[x];
    od;
return rec(FusionAlgebra:=Alg,Simples:=Base);
end;





epsilon:=function(n,k)
local e,l;
e:=0;
l:=Gcd(n+1,k);
if RemInt(n+1,2)=0 and RemInt((n+1)/l,2)=1 and RemInt(k/l,2) =1 then
    e:=1;
fi;
return e;
end;


######
# función usada en el proxima funcion
######
temporal:=function(a)
local i;
i:=1;
if a=0 then i:=0;
fi;
return i;
end;

Display_Info_Zesting:=function(n,k)
local e,l,d,A,B,Zestings,a,b;
l:=Gcd(n+1,k);
d:=Size(Labels_A(n,k));
e:=epsilon(n,k);
Print("\n");
Print("The generator of the transparent objects of PGL(", n+1,",",k,")");
Print("\n");
Print("is the power ", (n+1)/l, " of the simple object ", "v.", d);
Print("\n");
Print("and its twist is ", (-1)^e);
Print("\n");
Print("The pairs for computing braided zestings are");
Print("\n");
A:=[0..(l-1)];
B:=[0..n];
Zestings:=[];
for a in A 
    do
        for b in B 
        do 
            if RemInt(a*((n+1)/l)-temporal(a)*e^a-2*b,n+1)=0 then
                 Append(Zestings,[[a,b]]);
                 Print([a,b]);
                 Print("\n");
            fi;    
        od;
    od;
end;



Zesting_Modular_Data:=function(n,k,pair,c)
local MD,S,L,d,x,y,Tz,MDz,Szz,s,a,b,i,ij;
MD:=ModularData_TipoA(n,k);
S:=MD.Smatrix;
L:=MD.labels;
d:=Size(S);
a:=2*pair[2]+epsilon(n,k);
s:=E(2*(n+1)^2 )^a*E(n+1)^c;
Szz:=NullMat(d,d);
Tz:=ShallowCopy(MD.Tmatrix);
MDz:=ShallowCopy(MD);
for x in [1..d]
    do 
        for y in [x..d]
            do 
                ij:=grado_label(L.(x))*grado_label(L.(y));
                Szz[x][y]:=s^(2*ij)*S[x][y];
                 Szz[y][x]:=Szz[x][y];
            od;
    od;
for x in [1..d]
    do 
        i:=grado_label(L.(x));
        Tz[x][x]:=s^(-i^2)*Tz[x][x];
    od;
MDz.Smatrix:=Szz;
MDz.Tmatrix:=Tz;
return MDz;
end;


######
#
# verlinde zesting
########
zesting_verlinde:=function(n,k,par)
local MD,test;
MD:=Zesting_Modular_Data(n,k,par,0);
return Verlinde_formula(MD.Smatrix);
end;


Zesting_FusionAlgebra:=function(n,k,par)
local T,x,y,d,Alg,base,Base,N;
N:=zesting_verlinde(n,k,par);; # lista de matrices (N_i)_{jk}=N_{ij}^k
d:=Size(N);; # dimension del algebra de fusion
T:=NullMat(d,d);;
for x in [1..d]
    do
        for y in [1..d]
            do 
                T[x][y]:=FormatoTabla(N[x][y]);
            od;
    od;
Append(T,[1,0]); # datos adicionales pedidos por GAP
Alg:=AlgebraByStructureConstants( Rationals, T );;
base:=BasisVectors( CanonicalBasis( Alg ) );;
Base:=rec();
for x in [1..Size(base)] 
    do
        Base.(x):=base[x];
    od;
return rec(FusionAlgebra:=Alg,Simples:=Base);
end;


######
# Computes the list of n-th FS-indicators of SU(N+1)_k
####
FS_indicator:=function(N,k,n)
local MD,S,T,NN,FS,R,dimension,list_FS;
MD:=ModularData_TipoA(N,k);
S:=MD.Smatrix;
T:=MD.Tmatrix;
NN:=Test_verlinde(N,k);; # lista de matrices (N_i)_{jk}=N_{ij}^k
#FS(x):= 1/dim(C)sum_{i,j} N_{ij}^x S_{i,0} S_{j,0} (T_i/T_j)^2 
R:=Size(Labels_A(N,k));
FS:=0;
dimension:=0;
list_FS:=0*[1..R];
for i in [1..R] do
        dimension:=S[i][1]^2+dimension;
    od;

for x in [1..R] 
    do
        for i in [1..R] do 
             for j in [1..R] do 
        list_FS[x]:=NN[i][j][x]*S[i][1]*S[j][1]*(T[i][i]/T[j][j] )^n+list_FS[x];
    od;
od;
    od;

return (1/dimension)*list_FS;
end;
