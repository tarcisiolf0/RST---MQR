clear
clc
%% Planta
T_a=0.1;                          %Tempo de amostra
num=1;                          %Numerador Continuo
den=[0.1 1.1 1];                %Denominador Continuo
gp=tf(num,den);                 %Função de transferencia
ftz=c2d(gp,T_a, 'zoh');           %Planta Discreta

[num,den] = tfdata(ftz, 'v')         %num e den discreto
sys = filt(num,den, T_a);
A = den;
B = num;

% degA = polynomialDegree(polyA);
% degB = polynomialDegree(polyB);

Ts=1;                      %Tempo de estabelecimento desejado malha fechada
ep=0.6;                     %Epsilon (Coeficiente de amortecimento)
wn=4/(ep*Ts);               %Frequencia natural do sistema
Mp = exp((-ep*pi)/sqrt(1 - ep^2));% Overshoot 

%Polos conjugados
t1 = -2*exp(-ep*wn*T_a);
t2 = sqrt(1-ep^2);
t3 = wn*T_a;
t4 = cos(wn*T_a*sqrt(1-ep^2));


p1=-2*exp(-ep*wn*T_a)*cos(wn*T_a*sqrt(1-ep^2));
p2=exp(-2*ep*wn*T_a);

% Polinomio deseado con un polo insignificante igual a zero
Am=[1 p1 p2 0];
 
%Resolução de Equações metodo matricial
%M=[1 0 B(2); A(2) B(2) B(3); A(3) B(3) 0];

%q=[Am(2)-A(2); Am(3)-A(3);0;0];


M = [1 0 0 0;
    A(2) 1 B(2) 0;
    A(3) A(2) B(3) B(2);
    0 A(3) 0 B(3)];

p = [1; p1; p2; 0];

X=inv(M)*p;

Na = 2;
Nb = 1;
d = 1;

Nr = Na - 1;
Ns = Nb +d -1;

%R=[1 X(1)];             %Polinomio R
%S=[X(3) X(2)];            %Polinomio S

R = [X(1) X(2)];
S = [X(3) X(4)];
%Polinomio T
T=sum(Am)/sum(B);

