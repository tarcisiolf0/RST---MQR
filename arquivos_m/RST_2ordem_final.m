clc
close all
clear
%% Planta
Ta=0.1;                         %Tempo de amostra
num=1;                          %Numerador Continuo
den=[0.1 1.1 1];                %Denominador Continuo
gp=tf(num,den);                 %Função de transferencia
ftz=c2d(gp,Ta, 'zoh');          %Planta Discreta

[numd,dend] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(numd,dend, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Função de Transferência em TD

A = dend;
B = numd;

Na = 2;
Nb = 1;
d = 1;

Nr = Na - 1;
Ns = Nb + d -1;

%Condições do polinômio P(z^-1)
%0.25 <= w0Ta <= 1.5 ; 0.7 <= zeta <= 1

Ts=1;                              %Tempo de estabelecimento desejado malha fechada
ep=0.7;                            %Epsilon (Coeficiente de amortecimento)
wn=4/(ep*Ts);                      %Frequencia natural do sistema
Mp = exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 

initial_theta = [-0.4 -0.4 0.005 0.10];
w0 = 1;
Tas = Ta;

theta_switch_time = 10*Ta;