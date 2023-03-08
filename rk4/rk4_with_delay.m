clear all;
close all;
clc;

dt=0.1;
tau=5;
t=0:dt:500;
Tt=(-tau:dt:0);
n=length(t);
Sh0 = 0.96;
Ih0 = 0.04;
Rh0 = 0;
Sc0 = 0.84;
Ic0 = 0.16;
Rc0 = 0;
O0 = 0.2;

%parameter
miuH = 0.0233;
alphaH = 0.0233;
betaH = 0.0206;
gammaH = 0.3;
k = 0.2;
miuO = 0.0384;
miuC = 0.0115;
alphaC = 0.0115;
betaC = 0.962;
gammaC = 0.2;

%initial condition
a=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) alphaH-miuH*Sh-betaH*Sh*O;
b=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) betaH*Sh*O-(miuH+gammaH)*Ih;
c=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) gammaH*Ih-miuH*Rh;
d=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) k*Ik-miuO*O;
e=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) alphaC-betaC*Sc*O-miuC*Sc;
f=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) betaC*Sc*O-(miuC+gammaC)*Ic;
g=@(Sh,Ih,Rh,O,Sc,Ic,Ik,Rc,t) gammaC*Ic-miuC*Rc;
Sh = zeros(1,n);
Sh(1)=Sh0;
Ih = zeros(1,n);
Ih(1)=Ih0;
Rh = zeros(1,n);
Rh(1)=Rh0;
O = zeros(1,n);
O(1)=O0;
Sc = zeros(1,n);
Sc(1)=Sc0;
Ic = zeros(1,n);
Ic(1)=Ic0;
Rc = zeros(1,n);
Rc(1)=Rc0;
Ik=Ic0+Tt-Tt;

for x=1:(n-1)
    g1=dt*a(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    j1=dt*b(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    k1=dt*c(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    l1=dt*d(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    m1=dt*e(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    o1=dt*f(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    p1=dt*g(Sh(x),Ih(x),Rh(x),O(x),Sc(x),Ic(x),Ik(x),Rc(x),t(x));
    
    g2=dt*a(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    j2=dt*b(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    k2=dt*c(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    l2=dt*d(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    m2=dt*e(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    o2=dt*f(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    p2=dt*g(Sh(x)+g1/2,Ih(x)+j1/2,Rh(x)+k1/2,O(x)+l1/2,Sc(x)+m1/2,Ic(x)+o1/2,Ik(x)+o1/2,Rc(x)+p1/2,t(x)+dt/2);
    
    g3=dt*a(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    j3=dt*b(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    k3=dt*c(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    l3=dt*d(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    m3=dt*e(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    o3=dt*f(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    p3=dt*g(Sh(x)+g2/2,Ih(x)+j2/2,Rh(x)+k2/2,O(x)+l2/2,Sc(x)+m2/2,Ic(x)+o2/2,Ik(x)+o2/2,Rc(x)+p2/2,t(x)+dt/2);
    
    g4=dt*a(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    j4=dt*b(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    k4=dt*c(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    l4=dt*d(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    m4=dt*e(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    o4=dt*f(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    p4=dt*g(Sh(x)+g3,Ih(x)+j3,Rh(x)+k3,O(x)+l3,Sc(x)+m3,Ic(x)+o3,Ik(x)+o3,Rc(x)+p3,t(x)+dt);
    
    Sh(x+1)=Sh(x)+(g1+2*g2+2*g3+g4)/6;
    Ih(x+1)=Ih(x)+(j1+2*j2+2*j3+j4)/6;
    Rh(x+1)=Rh(x)+(k1+2*k2+2*k3+k4)/6;
    O(x+1)=O(x)+(l1+2*l2+2*l3+k4)/6;
    Sc(x+1)=Sc(x)+(m1+2*m2+2*m3+m4)/6;
    Ic(x+1)=Ic(x)+(o1+2*o2+2*o3+o4)/6;
    Rc(x+1)=Rc(x)+(p1+2*p2+2*p3+p4)/6;
    
    Ik((tau/dt)+1+x)=Ic(x);
end

plot(t,Sc,'-k', t, Ic, '--k',t, Rc, ':k', 'LineWidth', 1.5);
hold on
title('Dinamika Populasi Kucing')
legend('Sc', 'Ic', 'Rc')
xlabel('Waktu (per Hari)')
ylabel('Jumlah Populasi')
grid on