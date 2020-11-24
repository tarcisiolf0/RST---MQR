clear
clc
%% Planta
Ta=0.1;                         %Tempo de amostra
num=1;                          %Numerador Continuo
den=[0.1 1.1 1];                %Denominador Continuo
gp=tf(num,den);                 %Função de transferencia
ftz=c2d(gp,Ta, 'zoh');          %Planta Discreta

[num,den] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(num,den, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Função de Transferência em TD

A = den;
B = num;

Na = 2;
Nb = 1;
d = 1;

Nr = Na - 1;
Ns = Nb + d -1;

%Condições do polinômio P(z^-1)
%0.25 <= w0Ta <= 1.5 ; 0.7 <= zeta <= 1

Ts=1;                     %Tempo de estabelecimento desejado malha fechada
ep=0.6;                            %Epsilon (Coeficiente de amortecimento)
wn=4/(ep*Ts);                      %Frequencia natural do sistema
Mp = exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 

%Função de Transferência em Malha Fechada desejada
z = ep;
[numd,dend]=ord2(wn,z);

gpd = tf(numd,dend);
ftzd= c2d(gpd,Ta, 'zoh');           %Planta Discreta
[Bmd, Amd] = tfdata(ftzd, 'v');
sysd = filt(Bmd,Amd, Ta);
tfd = set(sysd, 'variable', 'z^-1'); %Função de Transferência em TD

%Polos conjugados
t1 = -2*exp(-ep*wn*Ta);
t2 = sqrt(1-ep^2);
t3 = wn*Ta;
t4 = cos(wn*Ta*sqrt(1-ep^2));

%Coeficientes do polinomio desejado de acordo 
%com as especificações de desempenho
p1=-2*exp(-ep*wn*Ta)*cos(wn*Ta*sqrt(1-ep^2));
p2=exp(-2*ep*wn*Ta);

%Coeficientes do polinomio desejado 
%com um polo insignificante igual a zero
Am=[1 p1 p2 0];

%Resolução de Equações metodo matricial
% M2=[1 0 B(2); A(2) B(2) B(3); A(3) B(3) 0];
% q2=[Am(2)-A(2); Am(3)-A(3);0];
% X2=inv(M2)*q2;
% R2=[1 X2(1)];
% S2=[X2(3) X2(2)];


M = [1      0       0       0;
    A(2)    1       B(2)    0;
    A(3)    A(2)    B(3)    B(2);
    0       A(3)    0       B(3)];

p = [1; p1; p2; 0];

X=inv(M)*p;

%Polinomio R e S
S = [X(1) X(2)];
R = [X(3) X(4)];

%Polinomio T
T=sum(Am)/sum(B);

tc = 100;

denhcl1 = conv(A,S);
denhcl2 = conv(B,R);
numHCL = (T*B);
denHCL = (denhcl1+denhcl2);