% Mediu elastic 1D
clear; close all; clc;
% Selectori:
TR=1; % dinamica in timp real - 1 / dinamica cu timp controlat - 0
TIP=1; % unda longitudinala - 1 / unda transversala - 0
ok=0; % parametru pentru realizarea unei singura semi-oscilatii
% Parametrii fizici ai sistemului:
P=100; % numarul de corpuri
m=0.75*ones(1,P); % kg; masele corpurilor primei jumatati
m(1,(P/2)+1:P)=5; % kg; masele corpurilor celei de-a doua jumatati
k=25*ones(1,P+1); % N/m; constantele elastice ale resorturilor primei jumatati
k(1,round((P+1)/2):(P+1))=50; % N/m; constantele elastice ale resorturilor 
%celei de-a doua jumatati
a=0.5; % m; lungimile resorturilor nedeformate
Tmax=2*pi*sqrt(mean(m)/mean(k)); % timp caracteristic sistemului
ti=0; tf=20*Tmax; N=100000; t=linspace(ti,tf,N); dt=t(2)-t(1); % timp discret
eta=zeros(N,P); % prealocare deplasari (N momente de timp; P pozitii)
% Deplasari initiale si valori de start:
eta0=zeros(1,P); % m; mediu neperturbat initial
eta(1,:)=eta0; eta(2,:)=eta0; % pasul 1 si pasul 2 temporale
% Conditii la capete (frontiera):
A=2*a; % amplitudinea perturbatiei
OM=pi/2; % pulsatia perturbatiei
etas=A*sin(OM*t); etad=0; % functii de timp (capat dreapta fixat)
for i=2:N-1 % ciclul temporal
    for j=2:P-1 % ciclul spatial
    % Recurenta de ordinul II (vezi Curs 4):
    eta(i+1,j)=2*eta(i,j)-eta(i-1,j)+...
        dt^2/m(j)*(k(j)*(eta(i,j-1)-eta(i,j))+k(j+1)*(eta(i,j+1)-eta(i,j)));
    end % ciclul spatial
    % Relatie particulara la j=1 (primul oscilator)
    eta(i+1,1)=2*eta(i,1)-eta(i-1,1)+...
        dt^2/m(1)*(k(1)*(etas(i)-eta(i,1))+k(2)*(eta(i,2)-eta(i,1)));
    % Relatie particulara la j=P (ultimul oscilator)
  
    if ok==0 && A-etas(i)<=1e-6
      ok=1;
    end
    if ok==1 && etas(i)<=1e-6
      etas=zeros(N,P);
      ok=2;
    end
      eta(i+1,P)=2*eta(i,P)-eta(i-1,P)+...
        dt^2/m(P)*(k(P)*(eta(i,P-1)-eta(i,P))+k(P+1)*(etad-eta(i,P)));
end % ciclul temporal
x=a:a:P*a; % coordonatele de echilibru ale corpurilor
figure(1); % Simularea dinamica a undelor
set(1,'Position',[50,100,1200,300]); % redimensionare fereastra
tic; simt=0; % porneste cronometrul si initializeaza timpul simularii
while simt<=tf % ciclul grafic
    hold off; % sterge plot precedent
    index=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t de simt
    if TIP==1
      plot((x(:)+eta(index,:)')*[1,1],[-A,A],'-b','MarkerSize',20); hold on;
      title('Unda longitudinala');
    else
      plot(x(:),eta(index,:),'.b','MarkerSize',20); hold on;
      ylabel('eta / m');
      title('Unda transversala');
    end
    xlabel('x / m');
    axis([x(1)-a x(P)+a,-2*A,2*A]);
    if TR==1 % 1 - in timp real
        simt=toc; % actualizeaza timpul simularii cu ceasul sistemului
        text(0.8*x(P),1.5*A,['t = ',num2str(round(t(index))),' s']);
    else
        simt=simt+5e-2; % incrementeaza timpul simularii
        text(0.8*x(P),1.5*A,['t = ',num2str(round(t(index)*10)),' ds']);
    end
    pause(1e-6)
end

    
    






