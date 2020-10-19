

#######################################
# Calcula la
# Multiplicative and additive Central charges
# el output es un par de c la aditiva y la cmultiplicativa
#######################################

CentralCharge:=function(R,k)
local c,n,x;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
c:= k*(n+1)*x^(-1); # additive central charge
return [c,E(12*DenominatorRat(c))^NumeratorRat(c)];
end;


#######################################
# Transpuesta conjugado
#######################################
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


######
#
# Test de verlinde
########
Test_verlinde:=function(n,k)
local MD,test;
MD:=ModularData_TipoA(n,k);
return Verlinde_formula(MD.Smatrix);
end;




