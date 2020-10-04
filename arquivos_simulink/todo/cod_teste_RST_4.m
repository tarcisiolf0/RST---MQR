clc;
close all;
clear;
%% Modelagem Nise (Sistema Massa Mola 4ª Ordem)
% K1 = K2 = K3 = 1 N/m
% fv1 = fv2 = fv3 = 1 Ns/m
% M1 = M2 = 1Kg = 1 Ns^2/m
% img2 = imread('nise_two_mass_figure.png');
% figure('Name', 'Sistema Massa Mola Nise')
% imshow(img2);

% Função de Transferênica Gs = X2(s)/F(s)
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

%% Controlador RST

Ta = 0.1;
ftz=c2d(gs4,Ta, 'zoh');          %Planta Discreta

[num,den] = tfdata(ftz, 'v');         %num e den discreto
sys = filt(num,den, Ta);
sysz = set(sys, 'variable', 'z^-1'); %Função de Transferência em TD

A = den;
B = num;
Nas = size(A);
Nbs = size(B);
d = 1;
Na = Nas(2) - 1;
Nb = Nbs(2) - 1 - d;

Np = Na + Nb + d -1;
Nr = Na - 1;
Ns = Nb +d -1;


%Condições do polinômio P(z^-1)
%0.25 <= w0Ta <= 1.5 ; 0.7 <= zeta <= 1

Ts = 1;                     %Tempo de estabelecimento desejado malha fechada
ep = 0.9;                            %Epsilon (Coeficiente de amortecimento)
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

polosdesejados = roots(dend);
p3 = 5 * real(polosdesejados(1));
p4 = 10 * real(polosdesejados(1));

denp3 = [1 -p3];
denp4 = [1 -p4];

tfcp3 = tf(1, denp3);
tfcp4 = tf(1, denp4);

ftzp3 = c2d(tfcp3,Ta, 'zoh');           %Planta Discreta
[Bp3, Ap3] = tfdata(ftzp3, 'v');
sysp3 = filt(Bp3,Ap3, Ta);
tfdp3 = set(sysp3, 'variable', 'z^-1'); %Função de Transferência em TD

ftzp4 = c2d(tfcp4,Ta, 'zoh');           %Planta Discreta
[Bp4, Ap4] = tfdata(ftzp4, 'v');
sysp4 = filt(Bp4,Ap4, Ta);
tfdp4 = set(sysp4, 'variable', 'z^-1'); %Função de Transferência em TD

% 1ª-Coeficientes do polinomio desejado de acordo 
%com as especificações de desempenho
p1=-2*exp(-ep*wn*Ta)*cos(wn*Ta*sqrt(1-ep^2));
p2=exp(-2*ep*wn*Ta);
p3=Ap3(2);
p4=Ap4(2);
%Coeficientes do polinomio desejado 
Am=[1 (p1+p3+p4) (p2+(p3*p4)+(p1*p3)+(p1*p4)) ((p1*p3*p4)+(p2*p3)+(p2*p4)) (p2*p3*p4)];
p = [1; (p1+p3+p4); (p2+(p3*p4)+(p1*p3)+(p1*p4)); ((p1*p3*p4)+(p2*p3)+(p2*p4)); (p2*p3*p4);0;0;0];
M = [1   0    0    0    0    0    0    0;
    A(2) 1    0    0    0    0    0    0;
    A(3) A(2) 1    0    B(2) 0    0    0;
    A(4) A(3) A(2) 1    B(3) B(2) 0    0;
    A(5) A(4) A(3) A(2) B(4) B(3) B(2) 0;    
    0    A(5) A(4) A(3) B(5) B(4) B(3) B(2);
    
    0    0   A(5)  A(4) 0    B(5) B(4) B(3);
    0    0   0     A(5) 0    0    B(5) B(4);
    
    ];
  
X=inv(M)*p;
% %Polinomio R e S
S = [X(1) X(2) X(3) X(4)];
R = [X(5) X(6) X(7) X(8)];

T=sum(Am)/sum(B);

tc=1;