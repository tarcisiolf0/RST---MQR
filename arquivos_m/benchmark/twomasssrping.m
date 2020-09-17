clc;
close all;
clear;
%% Modelagem Nise (Sistema Massa Mola 4ª Ordem)
% K1 = K2 = K3 = 1 N/m
% fv1 = fv2 = fv3 = 1 Ns/m
% M1 = M2 = 1Kg = 1 Ns^2/m
img2 = imread('nise_two_mass_figure.png');
figure('Name', 'Sistema Massa Mola Nise')
imshow(img2);

% Função de Transferênica Gs = X2(s)/F(s)
syms s
M1=1; M2=1;
fv1=1; fv2=1; fv3=1;
K1=1; K2=1; K3=1;

delta = [((M1*s^2)+(fv1 + fv3)*s + (K1+K2)) -(fv3*s + K2);
         -(fv3*s + K2)  ((M2*s^2)+(fv2 + fv3)*s + (K2+K3))];  
     
d = det(delta);

num4 = [fv3 K2];
den4 = [1 4 7 6 3];

% k = 3;
Gs4 = tf(num4, den4);
%% Aproximando o Sistema para um modelo de 2ª Ordem (Modelagen Nise)

% Mp = 24.2, ts = 8.55s, tr = 1.47
% Fator de Amortecimento - zeta
% zeta = -log(Mp)/sqrt(pi^2 + log(Mp)^2);
% tp = pi/(wn * sqrt(1 - zeta^2))
Mp = 0.242;
ts = 8.55;
tp = 3.62;

zeta = -log(Mp)/sqrt(pi^2 + log(Mp)^2);
wn = pi/(tp * sqrt(1 - zeta^2));
ts = 4/wn*zeta;

num2 = wn^2;
den2 = [1 2*wn*zeta wn^2];
k = 1/3;

Gs2 = tf(k*num2, den2);

figure('Name','Resposta ao Degrau - Sistema Massa Mola Nise');
t=0:0.01:15;
step(Gs4, t);
hold on
step(Gs2, t);

xlabel('Tempo');
ylabel('Amplitude');
title('Comparação das 2 Repostas');
legend('Sistema Original', 'Sistema de 2ª Ordem Aproximado');
grid on
hold off

%% Controlador RST

Ta = 0.1;
ftz=c2d(Gs2,Ta, 'zoh');          %Planta Discreta

[num,den] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(num,den, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Função de Transferência em TD

A = den;
B = num;

% Sistema com 1 Massa
% num = 1;
% den = [1 1 1];
% Gs2 = tf(num, den);
% Ta = 0.1;
% ftz=c2d(Gs2,Ta, 'zoh');          %Planta Discreta

[num,den] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(num,den, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Função de Transferência em TD

A = den;
B = num;
Na = 2;
Nb = 2;
d = 1;

Nr = Na - 1;
Ns = Nb +d -1;

%Condições do polinômio P(z^-1)
%0.25 <= w0Ts <= 1.5 ; 0.7 <= zeta <= 1

Ts=50;                     %Tempo de estabelecimento desejado malha fechada
ep=0.7;                            %Epsilon (Coeficiente de amortecimento)
%tp=3;                             % Tempo de pico desejado
%wn=pi/(tp * sqrt(1 - zeta^2));    % Frequencia Natural
wn = 4/(ep*Ts);
Mp=exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 

%Função de Transferência em Malha Fechada desejada
z = ep;
[numd,dend]=ord2(wn,z);

gpd = tf(numd,dend);
ftzd = c2d(gpd,Ta, 'zoh');           %Planta Discreta
[numdd, dendd] = tfdata(ftzd, 'v');
sysd = filt(numdd,dendd, Ta);
tfd = set(sysd, 'variable', 'z^-1'); %Função de Transferência em TD

%Coeficientes do polinomio desejado de acordo 
%com as especificações de desempenho
p1=-2*exp(-ep*wn*Ta)*cos(wn*Ta*sqrt(1-ep^2));
p2=exp(-2*ep*wn*Ta);

%Coeficientes do polinomio desejado 
%com um polo insignificante igual a zero
Am=[1 p1 p2 0];

M = [1 0 0 0;
    A(2) 1 B(2) 0;
    A(3) A(2) B(3) B(2);
    0 A(3) 0 B(3)];

p = [1; p1; p2; 0];

X=inv(M)*p;

%Polinomio R e S
R = [X(1) X(2)];
S = [X(3) X(4)];

%Polinomio T
T=sum(Am)/sum(B);