clear;
clc;
close all;
%% 1ª Ordem
num = 1; den = [10 1];  % Numerador e Denominador Contínuo
Gs = tf(num,den);       % Função de Transferência
Gz = c2d(Gs, 1, 'zoh'); % Planta Discretizada

[B, A] = tfdata(Gz, 'v');   % B = Numerador e A = Denominador
nA = 1; 
nB = 0;
d = 1;
Ta = 1;
Tal = den(1); 
Ctt = 1;

% Calculando a ordem dos polinômios R e S
nR = nB + d - 1;
nS = nA - 1;

% Calculando o parâmetro p1
p1 = exp(-Ta/Tal);
PMF = [1 -p1];

% Calculando os polinômios R e S por método matricial
% PMF(z) = A*R + z^(-d)*B*S

M = [1 0 ; A(2) B(2)];
q = [1 ; PMF(2)];
X = inv(M)*q;

R = X(1);
S = X(2:end);

% Calculando o Polinômio T
T = sum(PMF)/sum(B);