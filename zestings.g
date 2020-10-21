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

epsilon:=function(n,k)
local e,l;
e:=0;
l:=Gcd(n+1,k);
if RemInt(n+1,2)=0 and RemInt((n+1)/l,2)=1 and RemInt(k/l,2) =1 then
    e:=1;
fi;
return e;
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
            if RemInt(a*((n+1)/l)-e-2*b,n+1)=0 then
                 Append(Zestings,[[a,b]]);
                 Print([a,b]);
                 Print("\n");
            fi;    
        od;
    od;
end;



Zesting_Modular_Data:=function(n,k,pair)
local MD,S,L,d,x,y,Tz,MDz,Szz,s,a,b,i,ij;
MD:=ModularData_TipoA(n,k);
S:=MD.Smatrix;
L:=MD.labels;
d:=Size(S);
a:=2*pair[2]+epsilon(n,k);
s:=E(2*(n+1)^2 )^a;
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
MD:=Zesting_Modular_Data(n,k,par);
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
