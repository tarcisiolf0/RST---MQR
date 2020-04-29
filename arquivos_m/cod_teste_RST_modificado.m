clear
clc
%% Planta
T_a=0.1;                          %Tempo de amostra
num=1;                          %Numerador Continuo
den=[0.1 1.1 1];                %Denominador Continuo
gp=tf(num,den);                 %Funçao de transferencia
ftz=c2d(gp,T_a, 'zoh');           %Planta Discreta

[B,A] = tfdata(ftz, 'v')         %num e den discreto em função de z
%%
%transformar os coeficientes em função de Z para em função de z^-1
degr_A = size(pole(ftz), 1);
degr_B = size(zero(ftz), 1);

if degr_A > degr_B
    degr_dif = degr_A - degr_B;
    B = [B(1+degr_dif:end) zeros(degr_dif)];
elseif degr_B > degr_A
    degr_dif = degr_B - degr_A;
    A = [A(1+degr_dif:end) zeros(degr_dif)];
end

%%
gpz=filt(B,A,T_a,'name','Planta'); %planta na forma z⁻1


Ts=1;                      %Tempo de estabelecimento desejado malha fechada
ep=0.6;                     %Epsilon (Coeficiente de amortecimento)
wn=4/(ep*Ts);               %Frequencia natural do sistema

%%
%----------------------------------%
%Polos conjugados
t1 = -2*exp(-ep*wn*T_a);
t2 = sqrt(1-ep^2);
t3 = wn*T_a;
t4 = cos(wn*T_a*sqrt(1-ep^2));


p1=-2*exp(-ep*wn*T_a)*cos(wn*T_a*sqrt(1-ep^2));
p2=exp(-2*ep*wn*T_a);

% Polinomio deseado con un polo insignificante igual a zero
Am=[1 p1 p2 0];
 
%Resolução de Equaçõees metodo matricial
M=[1 0 B(2); A(2) B(2) B(3); A(3) B(3) 0];
q=[Am(2)-A(2); Am(3)-A(3);Am(4)];
X=inv(M)*q;
 
R=[1 X(1)];             %Polinomio R
S=[X(3) X(2)];            %Polinomio S
%Polinomio T
T=sum(Am)/sum(B);