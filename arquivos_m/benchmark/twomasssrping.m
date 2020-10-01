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

Ta = 0.1;
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
ep = 0.9;                            %Epsilon (Coeficiente de amortecimento)
%tp=3;                             % Tempo de pico desejado
%wn=pi/(tp * sqrt(1 - zeta^2));    % Frequencia Natural
wn = 4/(ep*Ts);
Mp = exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 

%Função de Transferência em Malha Fechada desejada
z = ep;
[numd,dend]= ord2(wn,z);

gpd = tf(numd,dend);
% figure
% step(gpd);
ftzd = c2d(gpd,Ta, 'zoh');           %Planta Discreta
[numdd, dendd] = tfdata(ftzd, 'v');
sysd = filt(numdd,dendd, Ta);
tfd = set(sysd, 'variable', 'z^-1'); %Função de Transferência em TD

% 1ª-Coeficientes do polinomio desejado de acordo 
%com as especificações de desempenho
p1=-2*exp(-ep*wn*Ta)*cos(wn*Ta*sqrt(1-ep^2));
p2=exp(-2*ep*wn*Ta);

% 2ª-Coeficiente do polinomio desejado de acordo  
% com as especificações de desempenho
p3=-exp(-ep*wn*Ta);
% p3 é a parte real dos polos complexos conjugados
% P(z) = (1 + p1z^-1 + p2z^-2) (1 + p3z^-1)

% 3ª-Polo 5 vezes a parte real dos polos complexos conjugados
% raizesdesejadas = roots(dend);
% polo1 = raizesdesejadas(1);
% polo2 = raizesdesejadas(2);
% polo3 = 5*real(polo1);
% RPMF = poly([polo1 polo2 polo3]);
% k1 = 6.68;
% PMF = tf(k1,RPMF);
% PMFD = c2d(PMF, Ta, 'zoh');
% [NUMPMFD, DENPMFD] = tfdata(PMFD, 'v');
% 
% figure
% t1=0:0.01:5;
% step(gpd,t1);
% hold on
% step(PMF,t1);
% xlabel('Tempo');
% ylabel('Amplitude');
% title('Comparação das 2 Repostas');
% legend('Sistema 2ª Ordem Desejado', 'Sistema de 3ª Ordem');
% grid on
% hold off

% 4 Multiplicação do polo auxiliar pelos polos desejados
% p3 = -1,3333; Multiplicando por 5, p3 = -6,66667
nump3 = 1;
denp3 = [1 13.333];

tf3 = tf(nump3, denp3);
tfd3 = c2d(tf3, Ta, 'zoh');
[numdp3, dendp3] = tfdata(tfd3, 'v');
tfd3z = set(tfd3, 'variable', 'z^-1');
p3t = dendp3(2);
%Coeficientes do polinomio desejado 
% Am=[1 p1 p2 0];
% Am=[1 (p1+p3) (p2+p1*p3) (p2*p3)];
% Am = [1 DENPMFD(2) DENPMFD(3) DENPMFD(4)];
Am = [1 (p3t+p1) (p2+(p1*p3t)) (p2*p3t)];
M = [1 0 0 0;
    A(2) 1 B(2) 0;
    A(3) A(2) B(3) B(2);
    0 A(3) 0 B(3)];

% p = [1; p1; p2; 0];
% p = [1; (p1+p3); (p2+p1*p3); (p2*p3)];
% p = [1; DENPMFD(2); DENPMFD(3); DENPMFD(4)];
p = [1; (p3t+p1); (p2+(p1*p3t)); (p2*p3t)];
X=inv(M)*p;

%Polinomio R e S
R = [X(1) X(2)];
S = [X(3) X(4)];

% S = [X(3) (X(4) - X(3)) -X(4)];
%Polinomio T
T=sum(Am)/sum(B);

denhcl1 = conv(A,S);
denhcl2 = conv(B,R);
numHCL = (T*B);
denHCL = (denhcl1+denhcl2);

tc = 1;