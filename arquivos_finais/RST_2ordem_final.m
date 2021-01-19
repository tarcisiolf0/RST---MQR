clc
close all
clear
%% Planta
Ta=0.1;                         %Tempo de amostra
num=1;                          %Numerador Continuo
den=[0.1 1.1 1];                %Denominador Continuo
gp=tf(num,den);                 %Funcao de transferencia
ftz=c2d(gp,Ta, 'zoh');          %Planta Discreta

[numd,dend] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(numd,dend, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Funcao de Transferencia em TD

A = dend;
B = numd;

Na = 2;
Nb = 1;
d = 1;

Nr = Na - 1;
Ns = Nb + d -1;

%Condicoes do polinomio P(z^-1)
%0.25 <= w0Ta <= 1.5 ; 0.7 <= zeta <= 1

Ts=1;                              %Tempo de estabelecimento desejado malha fechada
ep=0.7;                            %Epsilon (Coeficiente de amortecimento)
wn=4/(ep*Ts);                      %Frequencia natural do sistema
Mp = exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 

% Funcao de transferencia em malha fechada desejada
z = ep;
[numd,dend]=ord2(wn,z);
gpd = tf(numd,dend);
ftzd= c2d(gpd,Ta, 'zoh');           %Planta Discreta
[Bmd, Amd] = tfdata(ftzd, 'v');
sysd = filt(Bmd,Amd, Ta);
tfd = set(sysd, 'variable', 'z^-1'); %Funcao de Transferencia em TD

initial_theta = [-0.4 -0.4 0.005 0.10];
w0 = 1;
Tas = Ta;
theta_switch_time = 10*Ta;

t= 0:.01:7;
y = step(gp,t);
% step(gp);
plot(t,y, 'LineWidth', 2)
title('Resposta ao Degrau')
xlabel('Tempo(s)')
ylabel('Amplitude')
legend('y(t)')
