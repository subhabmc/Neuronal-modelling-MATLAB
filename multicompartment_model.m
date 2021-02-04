function V = multicompartment_model(I) 

%% Initialize variables
nc=100;
t=30;
dt=0.01;

% Conductances in mS/cm^2
gNa = 120;
gK = 36;
gL = 0.3;

% Resting potentials in mV
ENa = 115;
EK = -12;
EL = 10.6;

% Neuronal circuit parameters
Cm = 1;
T = 6.3;
Ra=(0.01)*(0.5e-4)/(pi*(1.5e-4)^2);

V = zeros(nc, t/dt+1);
m = V;
n = V;
h = V;
i_hh = V;

%Matrix C
C=zeros(nc,nc);
for i=1:nc
    for j=1:nc
    if abs(i-j)==1
        C(i,j)=1;
    elseif i==j
        C(i,j)=-2;
    end
    end
end
C(1,1)=-1;
C(nc,nc)=-1;

 A=eye(nc,nc)-(dt/(Cm*Ra))*C;
 
%% Equations

% Ionic currents:
iNa = @(V, m, h) gNa*m.^3.*h.*(V - ENa);
iK = @(V, n) gK*n.^4.*(V - EK);
iL = @(V) gL*(V - EL);

% Temperature correction (T in °C):
k = @(T) 3^(.1*(T - 6.3));

% Rate equations (V in mV):
am = @(V) ((2.5 - .1*V))./(exp(2.5 - .1*V) - 1);
an = @(V) (.1 - .01*V)./(exp(1 - .1*V) - 1);
ah = @(V) .07*exp(-V/20);
bm = @(V) 1*4.*exp(-V/18);
bn = @(V) 1*.125*exp(-V/80);
bh = @(V) 1./(exp(3 - .1*V) + 1);

% Initialization
m(:, 1) = am(V(1))/(am(V(1))+bm(V(1)));    % Initial m-value
n(:, 1) = an(V(1))/(an(V(1))+bn(V(1)));    % Initial n-value
h(:, 1) = ah(V(1))/(ah(V(1))+bh(V(1)));    % Initial h-value

rhom=300;
Ve=zeros(nc,t/dt+1);
x=Ve;
x(1,:)=0.5;
for i=1:100 
Ve(i,:)=(rhom/4*pi)*(I(1,:)./sqrt(10^2+(x(i,:)-25).^2));
if i==100
    break;
end
x(i+1,:)=x(i,:)+0.5;
end

%% HH --> Multi-compartment

% Compartment
for i = 1:(t/dt)
    
    % Exponential Euler method for calculating the Gating equations
    Ax = -k(T)*(am(V(:, i)) + bm(V(:, i)));
    Bx = k(T)*am(V(:, i));
    m(:, i+1) = m(:, i).*exp(Ax*dt) + (Bx/Ax)*(exp(Ax*dt)-1);
    
    Ax = -k(T)*(an(V(:, i)) + bn(V(:, i)));
    Bx = k(T)*an(V(:, i));
    n(:, i+1) = n(:, i).*exp(Ax*dt) + (Bx/Ax)*(exp(Ax*dt)-1);   
    
    Ax = -k(T)*(ah(V(:, i)) + bh(V(:, i)));
    Bx = k(T)*ah(V(:, i));
    h(:, i+1) = h(:, i).*exp(Ax*dt) + (Bx/Ax)*(exp(Ax*dt)-1);
    
    % Calculation of the i_ion current with the m, n and h parameters
    i_hh(:, i+1) = iNa(V(:, i), m(:, i+1), h(:, i+1)) + iK(V(:, i), n(:, i+1)) + iL(V(:, i));
    
    % Implicit Euler Method to calculate the voltage
    b = V(:, i) + (dt/Cm) * (-i_hh(:, i+1))+(dt/Cm*Ra)*C*Ve(:,i+1);
    
    V(:, i+1) = A\b;
    V(100,i+1)=V(99,i+1);
    V(1,i+1)=V(2,i+1);
end


end