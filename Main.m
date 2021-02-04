%% Prepare working space
clear
close all
clc

%% Stimulation currents
nc=100;

%Time in ms
t=100;
tx=5;
dt=0.025;

% Stimulation currents in uA/cm2
Istim1 = zeros(nc,t/dt+1);
Istim1(1,1:(tx/dt))=10;

Istim2 = zeros(nc,t/dt+1);
Istim2(20,1:(tx/dt))=10;
Istim2(80,1:(tx/dt))=10;

%% Function calls and plots
%Stimulation at 1st compartment
V = multicompartment_model(Istim1);


figure(1)                               
imagesc(0:.025:100, 1:100, V)
set(gca,'FontSize',12)
c = colorbar;
ax = gca;
ax.YDir = 'normal';
c.Label.String = 'V (mV)';
xlabel('t (ms)')
ylabel('Compartment Nr.')
title('\fontsize{10}Action potential propagation for stimulation at compartment 1')
print -depsc2 Plot1.eps

%Stimulation at 20th and 80th compartments
V = multicompartment_model(Istim2);

figure(2)                              
imagesc(0:.025:100, 1:100, V)
set(gca,'FontSize',12)
c = colorbar;
ax = gca;
ax.YDir = 'normal';
c.Label.String = 'V (mV)';
xlabel('t (ms)')
ylabel('Compartment Nr.')
title('\fontsize{10}Action potential propagation for stimulation at compartments 20 and 80')
print -depsc2 Plot2.eps