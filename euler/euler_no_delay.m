clear all;
close all;
clc;

alphaH = 0.0233;
miuH = 0.0233;
betaH = 0.0206;
gammaH = 0.5;
k = 0.2;
miuO = 0.0384;
miuC = 0.0115;
alphaC = 0.0115;
betaC = 0.962;
gammaC = 0.5;

Sh0 = 0.96;
Ih0 = 0.04;
Rh0 = 0;
Sc0 = 0.84;
Ic0 = 0.16;
Rc0 = 0;
O0 = 0.2;
t0 = 0;
tEnd = 200;
h = 0.1;
N = (tEnd-t0)/h;

%Solusi awal
T = (t0:h:tEnd);
Sh = zeros(1,(N+1));
Sh(1)=Sh0;
Ih = zeros(1,(N+1));
Ih(1)=Ih0;
Rh = zeros(1,(N+1));
Rh(1)=Rh0;
O = zeros(1,(N+1));
O(1)=O0;
Sc = zeros(1,(N+1));
Sc(1)=Sc0;
Ic = zeros(1,(N+1));
Ic(1)=Ic0;
Rc = zeros(1,(N+1));
Rc(1)=Rc0;

%Metode Euler
for i = 1:N;
    Fi=alphaH-miuH*Sh(i)-betaH*Sh(i)*O(i);
    Sh(i+1) = Sh(i)+h*Fi;
    Gi =betaH*Sh(i)*O(i)-(miuH+gammaH)*Ih(i);
    Ih(i+1) = Ih(i)+h*Gi;
    Ri = gammaH*Ih(i)-miuH*Rh(i);
    Rh(i+1) = Rh(i)+h*Ri;
    Hi=k*Ic(i)-miuO*O(i);
    O(i+1) = O(i)+h*Hi;
    Ji=alphaC-miuC*Sc(i)-betaC*Sc(i)*O(i);
    Sc(i+1) = Sc(i)+h*Ji;
    Ki =betaC*Sc(i)*O(i)-(miuC+gammaC)*Ic(i);
    Ic(i+1) = Ic(i)+h*Ki;
    Ni = gammaC*Ic(i)-miuC*Rc(i);
    Rc(i+1) = Rc(i)+h*Ni;
end

%plot
plot (T,Sc,'-k',T,Ic,'--k',T,Rc,':k','LineWidth',1.5);
xlabel ('Waktu'), ylabel ('Jumlah Populasi');
grid on;
legend('Kucing Susceptible', 'Kucing Infected', 'Kucing Recovered');
title('Dinamika Populasi Kucing');