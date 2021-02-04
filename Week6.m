%% Prepare working space
clear
close all
clc

%1 Current Point Source Potential Field
%Parameters
rhom=300;    %ohm-cm  
I=1;         %mA
x=0:50;      %um
y=0:50;      %um

%1.1 Calculating potential field across slice
for i=1:length(x)
    for j=1:length(y)
Vp(i,j)=(rhom/4*pi)*(I/sqrt(10^2+(x(i)-25)^2+(y(j)-25)^2));
end
end

%Plot
figure(1); 
imagesc(x,y,Vp);
colormap(jet);
c = colorbar;

xlabel('x (\mum)');
ylabel('y (\mum)');
c.Label.String = 'V (V)';
print -depsc2 Plot1.eps

%1.2 Potential fields, electrical fields and activation functions along axon
for i=1:length(x)
%I=1uA
Va(i)=(rhom/4*pi)*(I/sqrt(10^2+(x(i)-25)^2));
Ea(i)=(rhom/4*pi)*(-I*(x(i)-25)/(10^2+(x(i)-25)^2)^1.5);
Aa(i)=(rhom/4*pi)*((-I/(10^2+(x(i)-25)^2)^1.5))*(-3*(x(i)-25)^2/(10^2+(x(i)-25)^2)+1);
%I=-1uA
Van(i)=(rhom/4*pi)*(-I/sqrt(10^2+(x(i)-25)^2));
Ean(i)=(rhom/4*pi)*(I*(x(i)-25)/(10^2+(x(i)-25)^2)^1.5);
Aan(i)=(rhom/4*pi)*((I/(10^2+(x(i)-25)^2)^1.5))*(-3*(x(i)-25)^2/(10^2+(x(i)-25)^2)+1);
end

figure(2); 
subplot(2,1,1);
plot(x,Va);
title('I = 1 mA');
xlabel('x (\mum)');
ylabel('V (V)');
subplot(2,1,2);
plot(x,Van);
title('I = -1 mA');
xlabel('x (\mum)');
ylabel('V (V)');
print -depsc2 Plot2.eps

figure(3) ;
subplot(2,1,1);
plot(x,Ea);
title('I = 1 mA');
xlabel('x (\mum)');
ylabel('E (V/m)');
subplot(2,1,2);
plot(x,Ean);
title('I = -1 mA');
xlabel('x (\mum)');
ylabel('E (V/m)');
print -depsc2 Plot3.eps

figure(4);
subplot(2,1,1);
plot(x,Aa);
title('I = 1 mA');
xlabel('x (\mum)');
ylabel('A (V/m^2)');
subplot(2,1,2);
plot(x,Aan);
title('I = -1 mA');
xlabel('x (\mum)');
ylabel('A (V/m^2)');
print -depsc2 Plot4.eps

%External currents for axon stimulation
t=30;
ts=5;
tx=1;
dt=0.01;
tv=0:dt:t;
n=1:100;
I1 = zeros(1,t/dt+1);
I1(1,ts/dt:(ts+tx)/dt)=-0.25;
I2 = zeros(1,t/dt+1);
I2(1,ts/dt:(ts+tx)/dt)=-1;
I3 = zeros(1,t/dt+1);
I3(1,ts/dt:(ts+tx)/dt)=-0.5;
I3(1,(ts+tx)/dt+1:(ts+2*tx)/dt)=0.5;
I4 = zeros(1,t/dt+1);
I4(1,ts/dt:(ts+tx)/dt)=-2;
I4(1,(ts+tx)/dt+1:(ts+2*tx)/dt)=2;
I5 = zeros(1,t/dt+1);
I5(1,ts/dt:(ts+tx)/dt)=0.25;
I6 = zeros(1,t/dt+1);
I6(1,ts/dt:(ts+tx)/dt)=5;

V=multicompartment_model(I1);
figure(5);
imagesc(tv,n,V);
colormap(jet);
c=colorbar;
xlabel('t (ms)');
ylabel('Compartment Nr.');
c.Label.String = 'V (mV)';
title('Mono-phasic pulse -0.25 mA');
set(gca,'FontSize',12)
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot5.eps
V=multicompartment_model(I2); 
figure(6);
imagesc(tv,n,V);
colormap(jet);
c=colorbar;
xlabel('t (ms)');
ylabel('Compartment Nr.');
c.Label.String = 'V (mV)';
title('Mono-phasic pulse -1 mA');
set(gca,'FontSize',12)
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot6.eps
V=multicompartment_model(I3); 
figure(7);
title('Bi-phasic pulse 0.5 mA');
imagesc(tv,n,V);
colormap(jet);
c=colorbar;
xlabel('t (ms)');
ylabel('Compartment Nr.');
c.Label.String = 'V (mV)';
title('Bi-phasic pulse 0.5 mA');
set(gca,'FontSize',12)
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot7.eps

V=multicompartment_model(I4);
figure(8);

imagesc(tv,n,V);
colormap(jet);
c=colorbar;
xlabel('t (ms)');
ylabel('Compartment Nr.');
c.Label.String = 'V (mV)';
title('Bi-phasic pulse 2 mA');
set(gca,'FontSize',12)
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot8.eps

V=multicompartment_model(I5); 
figure(9);

imagesc(tv,n,V);
colormap(jet);
c=colorbar;
xlabel('t (ms)');
ylabel('Compartment Nr.');
c.Label.String = 'V (mV)';
title('Mono-phasic pulse 0.25 mA');
set(gca,'FontSize',12)
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot9.eps
V=multicompartment_model(I6); 
figure(10);

imagesc(tv,n,V);
colormap(jet);
c=colorbar;
xlabel('t (ms)');
ylabel('Compartment Nr.');
c.Label.String = 'V (mV)';
title('Mono-phasic pulse 5 mA');
set(gca,'FontSize',12)
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot10.eps