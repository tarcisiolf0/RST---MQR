clear
clc
close all
%Condições para realizar o ciclo de controle
nit = 400; umin = 0; umax = 4.9 ; ts = 1;
u(1:10) = 0; erro(1:10) = 0; y(1:10) = 0;
% Referência
r(1:100) = 1; r(101:200) = 4; r(201:300) = 2; r(301:nit) = 3;
%% Planta
num=1;                          %Numerador Continuo
den=[10 1];                     %Denominador Continuo
gp=tf(num,den);                 %Funcion de transferencia
ftz=c2d(gp,ts);                 %Planta Discreta
[B,A]=tfdata(ftz,'v');          %Divide Numerador em B e denominador em A
%PMF (Polo asignado)
p1=exp(-ts/den(1));
pmf=[1 -p1];
%%
%Solução das equações por metodo matricial
M=[A(1) B(1); A(2) B(2)];
q=[pmf(1); pmf(2)];
X= M\q;%inv(M)*q
R=X(1);                 %Polinomio R
S=X(2:end)';            %Polinomio S
%Polinomio T
T=sum(pmf)/sum(B);
%%
%Loop de Controle
for k = 4:nit
y(k)=B(2)*u(k-1)-A(2)*y(k-1);       %Saida do sistema
u(k)=(T*r(k)-S*y(k))/R;             %Lei de controle
%     Saturacção para a lei de controle min=0 e max=5
if u(k) >= umax 
    u(k) = umax;
elseif u(k) <= umin
    u(k) = umin;
end
end
%%
%Codigo para plotar a resposta transitoria e a ação de controle
t = 0:ts:(nit-1)*ts;
figure
subplot(2,1,1),plot(t,r,'--k',t,y,'-r','Linewidth',2.5)
xlabel('Tempo (s)');
ylabel('Velocidade');
legend('W_t','Y_t','Location','SouthEast')
grid on;
hold
subplot(2,1,2),plot(t,u,'--b','Linewidth',3)
xlabel('Tempo (s)');
ylabel('Controle (volts)');
legend('U_t')
grid on;