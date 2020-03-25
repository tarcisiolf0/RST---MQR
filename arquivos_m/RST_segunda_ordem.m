clear
close all
clc
%%
nit = 400; umin = 0; umax = 4.9; 
u(1:10) = 0; erro(1:10) = 0; y(1:10) = 0;
% Referência
r(1:100) = 1; r(101:200) = 4; r(201:300) = 2; r(301:nit) = 3;
%% Planta
T=0.1;                          %Tempo de amostra
num=1;                          %Numerador Continuo
den=[0.1 1.1 1];                %Denominador Continuo
gp=tf(num,den);                 %Função de transferencia
ftz=c2d(gp,T);                  %Planta Discreta
[B,A]=tfdata(ftz,'v');          %Divide Numerador em B e denominador em A
 
gpz=filt(B,A,T,'name','Planta');      %Funcao de transferenca discreta
 
%% PMF
%Atribuição de Polos
Ts=1;                      %Tempo de estabelecimento desejado malha fechada
ep=0.6;                     %Epsilon (Coeficiente de amortecimento)
wn=4/(ep*Ts);               %Frequencia natural do sistema
 
%Polos conjugados
p1=-2*exp(-ep*wn*T)*cos(wn*T*sqrt(1-ep^2));
p2=exp(-2*ep*wn*T);
 
% Polinomio deseado con un polo insignificante igual a zero
Am=[1 p1 p2 0];
 
%Resolução de Equações metodo matricial
M=[1 0 B(2); A(2) B(2) B(3); A(3) B(3) 0];
q=[Am(2)-A(2); Am(3)-A(3);Am(4)];
X=inv(M)*q;
 
R=[1 X(1)];             %Polinomio R
S=[X(3) X(2)];            %Polinomio S
%Polinomio T
Tz=sum(Am)/sum(B);
%%
for k = 4:nit
   
    y(k)=B(2)*u(k-1)+B(3)*u(k-2)-A(2)*y(k-1)-A(3)*y(k-2);
    u(k)=Tz*r(k)-S(1)*y(k)-S(2)*y(k-1)+R(2)*u(k-1);
    %     Saturação
    if u(k) >= umax 
     u(k) = umax;
     elseif u(k) <= umin
     u(k) = umin;
    end
end
 
t = 0:T:(nit-1)*T;
 
figure
subplot(2,1,1),plot(t,r,'--k',t,y,'-r','Linewidth',2.5)
xlabel('Tempo (s)');
ylabel('Velocidade');
legend('y_r','y','Location','SouthEast')
grid on;
hold
subplot(2,1,2),plot(t,u,'--b','Linewidth',3)
xlabel('Tempo (s)');
ylabel('Controle (volts)');
legend('u')
grid on;