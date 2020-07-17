clear;
clc;
close all;
%% 1� Ordem
num = 1; den = [10 1];  % Numerador e Denominador Cont�nuo
Gs = tf(num,den);       % Fun��o de Transfer�ncia
Gz = c2d(Gs, 1, 'zoh'); % Planta Discretizada

[B, A] = tfdata(Gz, 'v');   % B = Numerador e A = Denominador
nA = 1; 
nB = 0;
d = 1;
Ta = 1;
Tal = den(1); 
Ctt = 1;

% Calculando a ordem dos polin�mios R e S
nR = nB + d - 1;
nS = nA - 1;

% Calculando o par�metro p1
p1 = exp(-Ta/Tal);
PMF = [1 -p1];

% Calculando os polin�mios R e S por m�todo matricial
% PMF(z) = A*R + z^(-d)*B*S

M = [1 0 ; A(2) B(2)];
q = [1 ; PMF(2)];
X = inv(M)*q;

R = X(1);
S = X(2:end);

% Calculando o Polin�mio T
T = sum(PMF)/sum(B);