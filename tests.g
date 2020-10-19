
###############################
# Calcula la Multiplicative and additive Central charges
# el output es un par de c la aditiva y la cmultiplicativa
###############################
CentralCharge:=function(R,k)
local c,n,x;
n:=Size(SimpleSystem(R));
x:=k+n+1; # nivel 
c:= k*(n+1)*x^(-1); # additive central charge
return [c,E(12*DenominatorRat(c))^NumeratorRat(c)];
end;


###############################
# Transpuesta conjugado
###############################
DagaMatrix:=function(U)
return TransposedMat( ComplexConjugate(U) );
end;

###############################
# comprueba si U es unitaria
###############################

IsUnitaryMatrix:=function(U)
local v,n;
n:=Size(U);
v:=false;
    if U*DagaMatrix(U)=IdentityMat(n) then 
        v:=true;
    fi;
return v;
end;



###############################
# Suma de Gauss de un dato modular
#  S y T es el dato modular
# signo debe ser 1 o -1
###############################
Gauss_Sum:=function(S,T,signo)
local k,s,x;
s:=Size(S);
k:=0;
for x in [1..s] do
    k:=k+(S[1][x])^2*T[x][x]^(signo);
od;
return k;
end;



###############################
# revisa si el par es un dato modular para un dato cualquiera
###############################
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

###############################
# Test: Revisa si el dato obtenido con ModularData_TipoA(n,k)
# es dato modular
###############################
Test_modular:=function(n,k)
local MD,test;
MD:=ModularData_TipoA(n,k);
test:=IsModularData(MD.Smatrix,MD.Tmatrix);
return test;
end;


