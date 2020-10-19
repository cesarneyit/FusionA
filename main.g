########################################################
########################################################
## Primera Parte:
## Labels
## Calculamos los Labels
########################################################
########################################################

############################
# Función recursiva que calcula los simples en el caso SU(N)_k
# Es decir N-tuplas de enteros no negativos cuya suma es k
############################
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







############################
## Segunda Parte:
## Producto interno normalizado
############################

############################
# Calcula la matriz de la forma bilineal en el producto, 
# de tal forma que (r,r)=2 para las raices mas bajas
############################
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
############################
Inner_Product := function (R,wt1,wt2) 
local C, B, B_P;
C:=SimpleSystem(R)^(-1); 
B:=BilinearFormMat_Normalized(R);
B_P:=TransposedMat(C)*B*C;
return wt1*B_P*wt2;
end;

############################
# Matriz del producto interno normalizado
############################
MBP:=function(R)
local B,D,d,C,B_P;
C:=SimpleSystem(R)^(-1);
B:=BilinearFormMat_Normalized(R);
B_P:=TransposedMat(C)*B*C;
return B_P;
end;

############################
# producto interno mejorado
############################
Inner_Product_m := function (B_P,wt1,wt2) 
return wt1*B_P*wt2;
end;




########################################################
########################################################
## Tercera Parte:
## Construyendo la T-matriz
########################################################
########################################################

############################
# Entrada ii de la T-matrix
############################
Twist:=function(R,k,wt)
local x,rho,a,theta,n,B_P;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
rho:= rho_max(n);
B_P:=MBP(R)
a:=Inner_Product_m(B_P,wt ,wt+2*rho );
theta:=E(2*x*DenominatorRat(a))^NumeratorRat(a);
return theta;
end;

############################
# T-matrix
############################
Tmatrix:=function(R,k,labels)
local x,rho,a,theta,n,T,r;
r:=Size(labels);
T:=IdentityMat(r);
for x in [1..r] do
    T[x][x]:=Twist(R,k,labels[x]);
od;
return T;
end;

############################
# Twist como racional
############################
Twist_Racional:=function(R,k,wt)
local x,rho,a,theta,n;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
rho:= DiagonalOfMat(IdentityMat(n)); 
a:=Inner_Product(R,wt ,wt+2*rho )*(2*x)^(-1);
return a;
end;

############################
# Funcion que imprime los labels y los twist como racional
############################
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


########################################################
########################################################
## Cuarta Parte:
## S-matriz
########################################################
########################################################

############################
# contruye el peso rho para A_n
# en este caso es una fila de tamaño n de 1's
############################
rho_max:=function(n)
local rho,i;
rho:=[];
for i in [1..n]
    do
        Append(rho,[1]);
    od;
return rho;
end;



############################
# Calcula la entrada S_ij 
# Una parte importante y confusa es el 
# significado de q^a con a un racional en la formula de BaKi
# esto signnifica que si: a=p/q y  q=E(2*x), entonces q^a= E(2*x*q)^p
############################
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

############################
# S-matriz (CONJUGADA RESPECTO A LA FORMULA DE BAKI)
# normalizada en el sentido que S_00=1
############################
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



############################
# S.matrix normalizada
# Normalizada en el sentido que S^2=C
# con C_ij=\delta_{i,j^*}
############################
SMatrix_Normalized:=function(R,k,labels)
local D,S;
S:=SMatrix(R,k,labels);
D:=Sqrt(S[1]*S[1]);
return D^(-1)*S;
end;



########################################################
########################################################
## Qinta Parte:
## Dato modular como un registro
########################################################
########################################################


############################
# Calcula dato modular de A_n en nivel principal k
#  En otras palabras el dato modular de SU(n+1)_k
############################
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









