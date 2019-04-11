%% MAIN assignment 2

clear all;

seed = 43;
rng(seed)
tic
%% Predator-Prey
tic
% Time domain, coarseness of grid and number of realizations.  
T=1000;
N=1e5;

% the dimensionality of the brownian motion. 
nB=2;

% number of instances of nb-dimensional BM. 

% Time step
dt = T/N;

% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
dB = sqrt(dt) .* randn(nB,N);
B = [ zeros(nB,1) , cumsum(dB,2)];
Tw = 0:dt:T;

%Euler Maruyama

r = 1;
K = 1;
beta = -10;
epsilon = 0.1;
mu = 0.05;
sigma_N = 0.1;
sigma_P = 0.1;

x0 = ones(nB,1).*1e-4;
nx = size(x0,1);
X = zeros(nx,N+1);
X(:,1)=x0;
fb = @(x,y) [r*x*(1-x/K)-beta*x*y;epsilon*beta*x*y-mu*y];
gb = @(x,y) [sigma_N*x;sigma_P*y];
for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = fb(X(1,k),X(2,k));     
    g = gb(X(1,k),X(2,k));       
    X(:,k+1) = X(:,k)+f.*dt+g.*dB; 
end
toc
figure(12)
plot(Tw,X(1,:))
xlabel('t')
ylabel('N_t')
figure(13)
plot(Tw,X(2,:))
xlabel('t')
ylabel('P_t')
figure(14)
plot(X(1,:),X(2,:))
xlabel('N_t')
ylabel('P_t')

mean_N=mean(X(1,1:end));
mean_P=mean(X(2,1:end));

var_N=var(X(1,1:end));
var_P=var(X(2,1:end));

cov_NP=cov(X(1,1:end),X(2,1:end));
%% Lamperti

tic
% Stochastic Heun for the Ito equation, (include the transformation
% invariance argument)
%T=1000;
%N=10e5;

% the dimensionality of the brownian motion. 
%nB=2;

% number of instances of nb-dimensional BM. 

% Time step
%dt = T/N;
%Tw = 0:dt:T;
% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
%dB = sqrt(dt) .* randn(nB,N);
%B = [ zeros(nB,1) , cumsum(dB,2)];
%Tw = 0:dt:T;

%r = 1;
%K = 1;
%beta = 10;
%epsilon = 0.1;
%mu = 0.05;
%sigma_N = 0.1;
%sigma_P = 0.1;

% For the Lamperti transformed system 
x20 =log((ones(nB,1)).*1e-4).*sigma_N^-1;
nx = size(x20,1);
Y = zeros(nx,N+1);
Y(:,1)=x20;

fa = @(x,y) [((r*x*(1-x/K))-beta*x*y)/(sigma_N*x)-(1/2)*sigma_N;(epsilon*beta*x*y-mu*y)/(sigma_P*y)-(1/2)*sigma_P];

for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);
    f = fa(exp(sigma_N*Y(1,k)),exp(sigma_N*Y(2,k)));
    g = [1;1];
    % Euler-Maruyama step as predictor. 
    Y(:,k+1) = Y(:,k)+f.*dt+g.*dB;
    % Corrector. 
    f_tilde = (1/2).*(fa(exp(sigma_N*Y(1,k)),exp(sigma_N*Y(2,k)))+fa(exp(sigma_N*Y(1,k+1)),exp(sigma_N*Y(2,k+1))));
    g_tilde = (1/2).*([1;1]+[1;1]);
    Y(:,k+1) = Y(:,k)+f_tilde.*dt+g_tilde.*dB;
end
toc
figure(50)
plot(Tw,Y(1,:))
xlabel('t')
ylabel('Y_t')
figure(60)
plot(Tw,Y(2,:))
xlabel('t')
ylabel('Y_t')


%%
% runtime analysis done with paper and pencil. Can be reproduced by tictoc
A=[1e5,1e6,1e7];
B=[0.82,6.76,65.21];
C=[1.046,9.93,96.89];
figure(58)
plot(A,B)
xlabel('N')
ylabel('runtime')
hold on
plot(A,C)
legend('Euler-Maruyama','Heun with Lamperti')



%% We back-transform the solutions from the Heun method used on the lamperti transformed system



Y2(1,:)=exp(sigma_N.*Y(1,:));
Y2(2,:)=exp(sigma_P.*Y(2,:));


figure(1000)
plot(X(1,:))
xlabel('t')
ylabel('N_t')
hold on
plot(Y2(1,:))
legend('X_t','Z_t')
figure(2000)
plot(X(2,:))
xlabel('t')
ylabel('P_t')
hold on
plot(Y2(2,:))
legend('X_t','Z_t')


%% Sensitivity, NOT USED

% sensitivities for the zero solution. 

% Time domain, coarseness of grid and number of realizations.  
T=100;
N=1e5;

% the dimensionality of the brownian motion. 
nB=2;

% number of instances of nb-dimensional BM. 

% Time step
dt = T/N;

% Simulation of brownian motion as cumulative sum of the increments.
% Time axis included. 
dB = sqrt(dt) .* randn(nB,N);
B = [ zeros(nB,1) , cumsum(dB,2)];
Tw = 0:dt:T;

%Euler Maruyama

r = 1;
K = 1;
beta = 5;
epsilon = 0.1;
mu = 0.05;
sigma_N = 0.1;
sigma_P = 0.1;
n = 2;
x0 = ones(nB,1);
nx = size(x0,1);
X = zeros(nx,N+1);
X(:,1)=x0;
fb = @(x,y) [(r-((2*r*n)/K))*x-beta*n*y;(epsilon*beta*n-mu)*y];
gb = @(x,y) [sigma_N;sigma_P];
for k=1:N           
    dt = Tw(k+1)-Tw(k);                
    dB = B(:,k+1)-B(:,k);           
    f = fb(X(1,k),X(2,k));     
    g = gb(X(1,k),X(2,k));       
    X(:,k+1) = X(:,k)+f.*dt+g.*dB;                   
end


figure(120)
plot(Tw,X(1,:))
figure(130)
plot(Tw,X(2,:))

figure(140)
plot(X(1,:),X(2,:))


%% Transition probabilities and filter
clear;

% constants
sigma = 2;
lambda = 1;

% state grid
h = 0.2;
x = -20:h:20.1;

% Advection-diffusion form of Forward Kolmogorov
u = @(y) -lambda.*y-(sigma.^2).*y;

D = @(y) (1/2).*(sigma.^2).*(1+y.^2);

% we solve for different values of t

G = fvade(u,D,x,'r');

dt = 0.05;

P1 = expm(G*dt);

figure(10)
imagesc(x,x,P1)
colorbar

dt = 0.1;

P2 = expm(G*dt);

figure(20)
imagesc(x,x,P2)
colorbar
dt = 0.2;

P3 = expm(G*dt);

figure(30)
imagesc(x,x,P3)
colorbar

dt = 0.5;

P4 = expm(G*dt);

figure(40)
imagesc(x,x,P4)
colorbar

%%
figure(41)
plot(x(1,1:end-1),P1(100,:))
hold on
plot(x(1,1:end-1),P2(100,:))
plot(x(1,1:end-1),P3(100,:))
plot(x(1,1:end-1),P4(100,:))
legend('dt=0.05','dt=0.1','dt=0.2','dt=0.5')
xlabel('State')
ylabel('P')
figure(49)
plot(x(1,1:end-1),P1(5,:))
hold on
plot(x(1,1:end-1),P2(5,:))
plot(x(1,1:end-1),P3(5,:))
plot(x(1,1:end-1),P4(5,:))
legend('dt=0.05','dt=0.1','dt=0.2','dt=0.5')
xlabel('State')
ylabel('P')
%%

Y = importdata('Assignment2-observations.txt',' ');
Y = Y.data;
tY = Y(:,1);
Y = Y(:,2);
plot(tY,Y,'.')


%%

s = 1;

l = @ (x,y) (1./(s.*sqrt(1+x.^2).*sqrt(2.*pi))).*exp((-((y-x).^2))./(2.*s^2.*(1+x.^2)));  % state likelihood fun
xc = x(1:(end-1)) + 0.5 * diff(x);
L = l(xc,Y);

figure(13)
plot(x,l(1,x))  
hold on
plot(x,l(2,x))  
plot(x,l(3,x))
xlabel('x')
ylabel('l(x)')
legend('X_t=1','X_t=2','X_t=3')
hold off
%%

G = fvade(u,D,x,'r');

dt = 0.1;

P = expm(G*dt);

c = zeros(1,length(tY));  % normalization
phi = zeros(length(tY)+1,length(P));  % predictor
psi = zeros(length(tY),length(P));  % estimator

%%
% estimate prior probability of x as the stationary distribution
X0 = null(G');
X0 = X0./sum(X0);
phi(1,:) = X0;


for i = 1:length(tY)-1
    c(i) = sum(phi(i,:) .* L(i+1,:));
    psi(i,:) = (1./c(i)) * phi(i,:).*L(i+1,:);  % data update
    phi(i+1,:) = psi(i,:) * P;  % time update
end

figure(1)
plot(tY,sum(psi.*xc,2))
hold on
tY2 = tY:0.1:20.1;
tY2 = tY2';
plot(tY2,sum(phi.*xc,2))
plot(tY,Y,'.')
xlabel('t')
ylabel('X_t')
legend('psi','phi')

%% CAREFUL, RUNNING MORE THAN ONCE REMOVES MORE OF PHI, YOU WILL HAVE TO RUN THE WHOLE SCRIPT
phi(end,:)=[];
sum1 = abs(Y(1:length(Y(:,1)),1)-sum(phi.*xc,2));
sumPHI=sum(sum1);

sum2 = abs(Y(1:length(Y(:,1)),1)-sum(psi.*xc,2));
sumPSI=sum(sum2);




figure(2)
imagesc(psi)
figure(3)
imagesc(phi)




