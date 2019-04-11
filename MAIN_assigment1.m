%% MAIN assignment 1
clear all;

seed = 33;
rng(seed)



% Another seed was used in the report. I know that this is the point of 
% having a seed in the first place, but I have been unable to find it.
% It is late now and have to hand in in a moment, but the qualitative finds
% are still the same. Just saw now that all figures must print without 
% altering the code at all. That is the reason for the ugly, ugly placement
% of brownian motion all over the place. Soryy about that. 

%% Brownian motion

% Time domain, coarseness of grid and number of realizations.  
T=10;
N=10^4;

% the dimensionality of the brownian motion. 
nB=1;

% Time step
dt = T/N;

% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
dB = sqrt(dt) .* randn(nB,N);
B = [ zeros(nB,1) , cumsum(dB,2)];
Tw = 0:dt:T;
%% Stochastic Integral

% Initialization of the solutions. 
nB=3;
x0 = zeros(nB,1);
nx = size(x0,1);
X = zeros(nx,N+1);


sigma = 2;

% Euler-Maruyama, with no drift term. 
for k=1:N                                       
    g = sigma.*sqrt(1+B(:,k).^2);  
    X(:,k+1) = X(:,k)+g.*dB(:,k);                   
end

pvar = (sigma^2).*(Tw+(1/2).*Tw.^2);
mvar = -(sigma^2.*(Tw+(1/2).*Tw.^2));

mean = zeros(1,length(Tw));
pstd = sqrt(pvar);
mstd = -sqrt(pvar);
%% Plots for 3 realizations with mean and +- std. 
figure(1)
plot(Tw,X(1,:))
hold on 
plot(Tw,X(2,:))
plot(Tw,X(3,:))
plot(Tw,mean)
plot(Tw,pstd)
plot(Tw,mstd)
hold off
xlabel("t")
ylabel("I_t")

%% Histograms for the distribution at t=10


% Time domain, coarseness of grid and number of realizations.  
T=10;
N=10^4;

% the dimensionality of the brownian motion. 
nB=1000;

% Time step
dt = T/N;

% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
dB = sqrt(dt) .* randn(nB,N);
B = [ zeros(nB,1) , cumsum(dB,2)];
Tw = 0:dt:T;

x0 = zeros(nB,1);
nx = size(x0,1);
X = zeros(nx,N+1);


sigma = 2;

% Euler-Maruyama, with no drift term. 
for k=1:N                                       
    g = sigma.*sqrt(1+B(:,k).^2);  
    X(:,k+1) = X(:,k)+g.*dB(:,k);                   
end

pvar = (sigma^2).*(Tw+(1/2).*Tw.^2);
mvar = -(sigma^2.*(Tw+(1/2).*Tw.^2));

mean = zeros(1,length(Tw));
pstd = sqrt(pvar);
mstd = -sqrt(pvar);
figure(2)
std = sqrt(4*(T+(T^2)/2));
hist(X(:,length(X)),75)
xlabel("I_{10}")
ylabel("count")
line([std std], [0 140])
line([-std -std], [0 140])
%% 3b

% Time domain, coarseness of grid and number of realizations.  
T=10;
N=10^4;

% the dimensionality of the brownian motion. 
nB=1;

% Time step
dt = T/N;

% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
dB = sqrt(dt) .* randn(nB,N);
B = [ zeros(nB,1) , cumsum(dB,2)];
Tw = 0:dt:T;
sigma = 100;
lambda = 100;
s = 2;

% Pre-allocation and IC's. 
x0 = [0;0];
nx = size(x0,1);
X1 = zeros(nx,N+1); 
X1(:,1) = x0;    

% The Euler Maruyama Method, with time and white noise steps. 
for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = [-lambda*X1(1,k);-sin(X1(2,k))+s.*X1(1,k).*cos(X1(2,k))];     
    g = [sigma;0.0];       
    X1(:,k+1) = X1(:,k)+f.*dt+g.*dB;                   
end

figure(3)
plot(Tw,X1(1,:))
xlabel("t")
ylabel("X_t")
figure(4)
plot(Tw,X1(2,:))
xlabel("t")
ylabel("Y_t")

% no need to include? 
figure(5)
plot(X1(1,1:500),X1(2,1:500))
xlabel("X_t")
ylabel("Y_t")



%% 3c
figure(6)
plot(Tw,B)
xlabel("t")
ylabel("B_t")

% very rough approximation of integral
intX=cumsum(X1(1,:));

figure(100)
plot(Tw,intX/1000)
xlabel("t")
ylabel("int_0^t X_s ds")
%%
diffB = intX/1000 - B;

plot(diffB)
%% 3d


sigma = 100;
lambda = 100;
s = 2;

% Pre-allocation and IC's. 
x0 = 0;
nx = size(x0,1);
X2 = zeros(nx,N+1); 
X2(:,1) = x0;    

% The Euler Maruyama Method, with time and white noise steps. 
for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = -sin(X2(k));     
    g = s.*cos(X2(k));       
    X2(:,k+1) = X2(:,k)+f.*dt+g.*dB;    
end

figure(8)
plot(Tw,X2(1,:))
xlabel("t")
ylabel("Y_t^I")

figure(11)
plot(Tw,X1(2,:))
hold on
plot(Tw,X2(1,:))
xlabel("t")
ylabel("Y")
legend("Y_t","Y_t^I")
%%

X1mX2 = X1(2,:)-X2(1,:);

figure(12)
plot(Tw,X1mX2)
xlabel("t")
ylabel("Y_t-Y_t^I")
%% 3e

x0 = 0;
nx = size(x0,1);
X3 = zeros(nx,N+1); 
X3(:,1) = x0;   

for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = -sin(X3(k));     
    g = s*cos(X3(k));
    % Euler-Maruyama predictor step. 
    X3(:,k+1) = X3(:,k)+f.*dt+g.*dB;
    % Modification of drift and diffusion. 
    f_tilde = (1/2).*(-sin(X3(k))-sin(X3(k+1)));
    g_tilde = (1/2).*(s.*cos(X3(k))+s.*cos(X3(k+1)));
    X3(:,k+1) = X3(:,k)+f_tilde.*dt+g_tilde.*dB;
end

figure(9)
plot(Tw,X3(1,:))
xlabel("t")
ylabel("Y_t^S")


X1mX3=X1(2,:)-X3(1,:);

figure(10)
plot(Tw,X1(2,:))
hold on
plot(Tw,X3(1,:))
xlabel("t")
ylabel("Y")
legend("Y_t","Y_t^S")

figure(15)
plot(Tw,X1mX3)
xlabel("t")
ylabel("Y_t-Y_t^S")



%% 3f Numerical verification of the conversion between stratonovich and Ito


% Time domain, coarseness of grid and number of realizations.  
T=2;
N=10^3;
nB=1;

% Time step
dt = T/N;

% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
dB = sqrt(dt) * randn(nB,N);
B = [ zeros(nB,1) , cumsum(dB,2)];
Tw = 0:dt:T;
s = 2;

x0 = 0;
nx = size(x0,1);
X4 = zeros(nx,N+1); 
X4(:,1) = x0;   

for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = X4(k)*sqrt(1-X4(k).^2);     
    g = s*(1-X4(k).^2);
    % Euler-Maruyama predictor step. 
    X4(:,k+1) = X4(:,k)+f*dt+g*dB;
    % Modification of drift and diffusion. 
    f_tilde = (1/2)*(-X4(k)*sqrt(1-X4(k).^2)-X4(k+1)*sqrt(1-X4(k+1).^2));
    g_tilde = (1/2)*(s*(1-X4(k).^2)+s*(1-X4(k+1).^2));
    X4(:,k+1) = X4(:,k)+f_tilde*dt+g_tilde*dB;
end


% Pre-allocation and IC's. 
x0 = 0;
nx = size(x0,1);
X5 = zeros(nx,N+1); 
X5(:,1) = x0;    

% The Euler Maruyama Method, with time and white noise steps. 
for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = -X5(k)*sqrt(1-X5(k).^2)-s.^2*(X5(k)-X5(k).^3);     
    g = s*(1-X5(k).^2);       
    X5(:,k+1) = X5(:,k)+f*dt+g*dB;    
end

figure(10)
plot(Tw,X4(1,:))
hold on
plot(Tw,X5(1,:))
legend("Stratonovich","Ito")
title("Stratonovich vs Ito")
xlabel("t")
ylabel("Z_t")

