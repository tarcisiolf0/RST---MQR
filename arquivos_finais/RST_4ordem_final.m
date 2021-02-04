clc;
close all;
clear;
%% Modelagem Nise (Sistema Massa Mola 4� Ordem)

% Fun��o de Transfer�nica Gs = X2(s)/F(s)
syms s
M1=1; M2=4;
fv1=0.1; fv2=0.4; fv3=1; 
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
gs4 = tf(num4, den4);

% Planta Discreta
% 1/10 da maior constante de tempo
% 10% do tempo de acomoda��o
% T_95/Ta = 5 ... 15

% Ta = 1.4992;
% Ta = 5.92;
% Ta = 3.94;

% Ta = 1.4992;
Ta = 1;
% Ta = 4.88;
ftz=c2d(gs4,Ta, 'zoh');          %Planta Discreta
[numd,dend] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(numd,dend, Ta);
sys_d = set(sys, 'variable', 'z^-1'); %Fun��o de Transfer�ncia em TD

% step(gs4)
% hold on
% step(ftz)

%initial_theta = [numd(2) numd(3) numd(4) numd(5) dend(2) dend(3) dend(4) dend(5)];

initial_theta = [-1.8 1.58 -0.66 0.23 0.03 0.09 -0.01 -0.009];
% initial_theta = [-0.200 0.0200 0.3000 0.1700 0.0200 0.3000 0.1000 -0.1200];
% initial_theta = [-1.3 -0.02 0.07 0.2 -0.001 -6e-05 0.003 -7e-06];
% initial_theta = [-0.45 0.14 0.45 0.07 0.1 0.27 0.9 -0.05];
w0 = 1/10;
Tas = Ta;

theta_switch_time = 10*Ta;

A = dend;
B = numd;
Na = 4;
Nb = 3;
d = 1;

Np = Na + Nb + d -1;
Nr = Na - 1;
Ns = Nb +d -1;
 
% Condi��es do polin�mio P(z^-1)
% 0.25 <= w0Ta <= 1.5 ; 0.7 <= zeta <= 1

Ts = 10;                     %Tempo de estabelecimento desejado malha fechada
ep = 0.7;                            %Epsilon (Coeficiente de amortecimento)
wn = 4/(ep*Ts);
Mp = exp((-ep*pi)/sqrt(1 - ep^2)); %Overshoot 
Tao = (2*pi)/(25*wn*(1-2*ep^2)+sqrt(4*ep^4 -4*ep^2 + 2)); % Tempo de amostragem Ogata

%% Plot
close all
t= 0:.01:90;
y = step(gs4,t);
plot(t,y, 'LineWidth', 2)
set(gca,'FontSize',15)
grid on
title('Resposta ao Degrau', 'FontSize', 30)
xlabel('Tempo(s)', 'FontSize', 25)
ylabel('Amplitude', 'FontSize', 25)
legend('Saída do Sistema', 'FontSize', 20)
figure 

% Malha fechada
sys = feedback(gs4,1);
y2 = step(sys,t);
plot(t,y2, 'LineWidth', 2)
set(gca,'FontSize',15)
grid on
title('Resposta ao Degrau', 'FontSize', 30)
xlabel('Tempo(s)','FontSize', 25)
ylabel('Amplitude', 'FontSize', 25)
legend('Saída do Sistema', 'FontSize', 20)
