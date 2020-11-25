clc;
close all;
clear;
%% Modelagem Nise (Sistema Massa Mola 4ª Ordem)
% K1 = K2 = K3 = 1 N/m
% fv1 = fv2 = fv3 = 1 Ns/m
% M1 = M2 = 1Kg = 1 Ns^2/m
%img2 = imread('nise_two_mass_figure.png');
%figure('Name', 'Sistema Massa Mola Nise')
%imshow(img2);

% Função de Transferênica Gs = X2(s)/F(s)
syms s
M1=1; M2=1;
fv1=1; fv2=1; fv3=1;
K1=1; K2=1; K3=1;

% delta = [((M1*s^2)+(fv1 + fv3)*s + (K1+K2)) -(fv3*s + K2);
%          -(fv3*s + K2)  ((M2*s^2)+(fv2 + fv3)*s + (K2+K3))];  
     
a = [M1 (fv1 + fv3) (K1+K2)];
b = [0 -fv3 -K2];
c = [0 -fv3 -K2];
d = [M2 (fv2 + fv3) (K2+K3)];

delta1 = conv(a, d);
delta2 = -conv(c, b);
d = delta1 + delta2;

num4 = [fv3 K2];
den4 = d;

% k = 3;
gs4 = tf(num4, den4);
%% Aproximando o Sistema para um modelo de 2ª Ordem (Modelagen Nise)

% Mp = 24.2, ts = 8.55s, tr = 1.47
% Fator de Amortecimento - zeta
% zeta = -log(Mp)/sqrt(pi^2 + log(Mp)^2);
% tp = pi/(wn * sqrt(1 - zeta^2))
mp = 0.242;
ts = 8.55;
tp = 3.62;

zeta = -log(mp)/sqrt(pi^2 + log(mp)^2);
wn = pi/(tp * sqrt(1 - zeta^2));
ts = 4/wn*zeta;

num2 = wn^2;
den2 = [1 2*wn*zeta wn^2];
k = 1/3;

gs2 = tf(k*num2, den2);

figure('Name','Resposta ao Degrau - Sistema Massa Mola Nise');
t=0:0.01:15;
step(gs4, t);
hold on
step(gs2, t);

xlabel('Tempo');
ylabel('Amplitude');
title('Comparação das 2 Repostas');
legend('Sistema Original', 'Sistema de 2ª Ordem Aproximado');
grid on
hold off

%% Controlador RST

Ta = 0.01;
ftz=c2d(gs2,Ta, 'zoh');          %Planta Discreta

[num,den] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(num,den, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Função de Transferência em TD

A = den;
B = num;
Na = 2;
Nb = 1;
d = 1;

Np = Na + Nb + d -1;
Nr = Na - 1;
Ns = Nb +d -1;

%Condições do polinômio P(z^-1)
%0.25 <= w0Ta <= 1.5 ; 0.7 <= zeta <= 1

Ts = 1;                     %Tempo de estabelecimento desejado malha fechada
ep = 0.7;                            %Epsilon (Coeficiente de amortecimento)
%tp=3;                             % Tempo de pico desejado
%wn=pi/(tp * sqrt(1 - zeta^2));    % Frequencia Natural
wn = 4/(ep*Ts);
Mp = exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 

%Função de Transferência em Malha Fechada desejada
z = ep;
[numd,dend]= ord2(wn,z);

gpd = tf(numd,dend);
ftzd = c2d(gpd,Ta, 'zoh');           %Planta Discreta
[Bmd, Amd] = tfdata(ftzd, 'v');
sysd = filt(Bmd,Amd, Ta);
tfd = set(sysd, 'variable', 'z^-1'); %Função de Transferência em TD

% 1ª-Coeficientes do polinomio desejado de acordo 
%com as especificações de desempenho
p1=-2*exp(-ep*wn*Ta)*cos(wn*Ta*sqrt(1-ep^2));
p2=exp(-2*ep*wn*Ta);

%Coeficientes do polinomio desejado 
Am=[1 p1 p2 0];
p = [1;p1;p2;0];
M = [1      0       0       0;
    A(2)    1       B(2)    0;
    A(3)    A(2)    B(3)    B(2);
    0       A(3)    0       B(3)];

X=inv(M)*p;

%Polinomio R e S
S = [X(1) X(2)];
R = [X(3) X(4)];

%Polinomio T 2ª Equação Diofantina Entrada Degrau
% Xdeg = 1 - z^-1
X1 = [1 -1 0 0];

M1 = [X1(1) 0     B(1)    0
      X1(2) X1(1) B(2)   B(1)
      X1(3) X1(2) B(3)   B(2)
      X1(4) X1(3) 0      B(3)];
  
X3 = inv(M1)*p;

L = [X3(1) X3(2)];
T = [X3(3) X3(4)];

tc = 55;
w0 = 1;
Tas = 1;

% denhcl1 = conv(A,S);
% denhcl2 = conv(B,R);
% numHCL = (X3(3)*B);
% denHCL = (denhcl1+denhcl2);