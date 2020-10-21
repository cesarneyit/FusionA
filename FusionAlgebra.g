

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

#######################################
# comprueba si U es unitaria
#######################################

IsUnitaryMatrix:=function(U)
local v,n;
n:=Size(U);
v:=false;
    if U*DagaMatrix(U)=IdentityMat(n) then 
        v:=true;
    fi;
return v;
end;



#######################################
# Suma de Gauss de un dato modular
#  S y T es el dato modular
# signo debe ser 1 o -1
#######################################

Gauss_Sum:=function(S,T,signo)
local k,s,x;
s:=Size(S);
k:=0;
for x in [1..s] do
    k:=k+(S[1][x])^2*T[x][x]^(signo);
od;
return k;
end;


#######################################
# Formula de verlinde
# Dada S constuye una lista
# N cuya entrada i-esima contien
# la matriz (N_i)_{jk}=N_{i,j}^k
#######################################
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



#######################################
# Formula de verlinde MEJORADA
# Dada S y un label i construye
# N_i cuya entrada i-esima contiene
# la matriz (N_i)_{jk}=N_{i,j}^k
#######################################
Verlinde_N:=function(S,x)
local N,tmp,y,z,n,D;
n:=Size(S);
N:=[];
D:=S[1]*S[1];
tmp:=DiagonalMat(S[1])^(-1)*DiagonalMat(S[x]); # matrix diagonal  Gamma_a, dada por  s_{aj}/s_{0j}
return D^(-1)*S*tmp*ComplexConjugate(S);
end;


#######################################
# Test de verlinde
#######################################
Test_verlinde:=function(n,k)
local MD,test;
MD:=ModularData_TipoA(n,k);
return Verlinde_formula(MD.Smatrix);
end;


#######################################
# Funcion usada para entrar 
# los coeficientes de estructura en la funcion StructureConstantsTable
# Nij es una fila con Nij[k]=N_{ij}^k
#######################################
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


#######################################
# Fusion principal
# Estacion produce un registro que 
# contiene el algebra de fusion y por comodidad
# la base canonica como una lista ordenada siguiendo.
#######################################
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

