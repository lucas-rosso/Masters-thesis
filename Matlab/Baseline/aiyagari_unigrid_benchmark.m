% --------------------------------------------------- %
%% The Portfolio Choice Channel of Wealth Inequality %%
% Author: Lucas Rosso, based on Ben Moll's code with non-uniform grid %
% uniform grid, Aiyagari (1994) + Portfolio Choice a la Merton (1969)
% Date: 18-02-2021 %
% Extremely Preliminar %
% --------------------------------------------------- %
%% AIYAGARI + FAT TAIL

clear all; clc; close all;

tic;

ga = 2;
rho = 0.053;

% Guvenen et al. (2019)
the = -log(0.9);
sigma2 = (0.04*2*the)/(1-exp(-2*the));
sigma = sqrt(sigma2);

z1 = 1 - sigma;
z2 = 1 + sigma;

l_z1 = log(z1); l_z2 = log(z2);

p_z = sqrt(the/(pi*sigma2*(1-exp(-2*the))))*exp(-(the*(l_z1-(l_z2*exp(-the)))^(2))/(sigma2*(1-exp(-2*the))));

z     = [round(z1,2),round(z2,2)]; % state values (row vector)
la    = [round(p_z,2),round(p_z,2)]; % Poisson intensities

w = 1;
R = 0.06;


sig2 = 0.18^2; 
sig = 0.18;

I= 10000;
phi = 1;
amin = -phi;
amax = 75;

a = linspace(amin,amax,I)';
da = a(2)-a(1);
da2 = da^2;

aa = [a,a];

zz = ones(I,1)*z;
kk = ones(I,2);

% safe asset
rbPos = 0.02;      % return on safe asset 
rbNeg = 0.08;      % 6 percent wedge as in Kaplan, Moll and Violante (AER, 2018)
r = rbPos.*(aa-kk>=0) + rbNeg.*(aa-kk<0);

maxit= 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
dV0 = zeros(I,2);
dV2 = zeros(I,2);
c = zeros(I,2);
c0 = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

%INITIAL GUESS
v0(:,1) = (w*z(1) + rbPos.*a).^(1-ga)/(1-ga)/rho;
v0(:,2) = (w*z(2) + rbPos.*a).^(1-ga)/(1-ga)/rho;

v = v0;

for n=1:maxit
    disp(n)
    V = v; 
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z + r(I,:).*amax + (R-r(I,:)).^2/(ga*sig2)*amax).^(-ga); %will never be used
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z + r(1,:).*amin).^(-ga); %state constraint boundary condition

    %second derivative: approximation only differs at amax
    dV2(2:I-1,:) = (V(3:I,:)-2*V(2:I-1,:)+V(1:I-2,:))/da2; 
    dV2(I,:) = -ga*dVb(I,:)/amax;
    
    %consumption and savings with forward difference
    cf = max(dVf,10^(-10)).^(-1/ga);
    kf = max(-(dVf./dV2).*(R-r)./sig2,0);
    kf = min(kf,aa+phi);
    ssf = w*zz + (R-r).*kf + r.*aa - cf;
    %consumption and savings with backward difference
    cb = max(dVb,10^(-10)).^(-1/ga);
    kb = max(-(dVb./dV2).*(R-r)./sig2,0);
    kb = min(kb,aa+phi);
    ssb = w*zz + (R-r).*kb + r.*aa - cb;
    %consumption and derivative of value function at steady state
    k0 = (kb + kf)/2; %could do something fancier here but the following seems to work well. And more importantly, it's fast
    c0 = w*zz + (R-r).*k0 + r.*aa;
    dV0 = max(c0,10^(-10)).^(-ga);
   
    
    % dV_upwind makes a choice of forward or backward differences based on the sign of the drift    
    If = real(ssf) > 10^(-12); %positive drift --> forward difference
    Ib = (real(ssb) < -10^(-12)).*(1-If); %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
   
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = max(dV_Upwind,10^(-10)).^(-1/ga);
    u = c.^(1-ga)/(1-ga);
    k = max(-dV_Upwind./dV2.*(R-r)./sig2,0);
    k = min(k,aa+phi);
    
   
    %CONSTRUCT MATRIX
    X = -Ib.*ssb./da + (sig2/2).*(k.^2)/da2;
    Y = - If.*ssf./da + Ib.*ssb./da - sig2.*(k.^2)/da2;
    Z = If.*ssf./da + (sig2/2).*(k.^2)/da2; 
    
    xi = -amax*(R-r).^2/(2*ga*sig2);
    X(I,:) = -min(ssb(I,:),0)/da - xi(I,:)/da;
    Y(I,:) = -max(ssf(I,:),0)/da + min(ssb(I,:),0)/da + xi(I,:)/da;
    Z(I,:) = max(ssf(I,:),0)./da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A1(I,I) = Y(I,1) + Z(I,1);
    A2(I,I) = Y(I,2) + Z(I,2);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    B = (1/Delta + rho)*speye(2*I) - A;
    
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;
 
    % check for proper transition matrix
    error = full(sum(A,2));
    check = abs(error)>1e-9;
    check_ = sum(check);
    
    if max(abs(sum(A,2)))>10^(-9)
        disp('Improper Transition Matrix')
        max(abs(sum(A,2)))
        how_many = check_ 
    break
    end

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%RECOMPUTE TRANSITION MATRIX WITH REFLECTING BARRIER AT amax
X = -min(ssb,0)./da + sig2/2.*k.^2./da2;
Y = -max(ssf,0)./da + min(ssb,0)./da - sig2.*k.^2./da2;
Z = max(ssf,0)./da + sig2/2.*k.^2./da2; 

A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
A1(I,I) = Y(I,1) + Z(I,1);
A2(I,I) = Y(I,2) + Z(I,2);
A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;

%WORK WITH RESCALED DENSITY \tilde{g} BELOW
AT = A';
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 157;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
% g_tilde = AT\b;
g_stacked = AT\b;

%rescale \tilde{g} so that it sums to 1
g_sum = g_stacked'*ones(2*I,1)*da;
g_stacked = g_stacked./g_sum;

g = reshape(g_stacked,I,2);

% check1 = g(:,1)'*da_tilde;
% check2 = g(:,2)'*da_tilde;
% check = sum(g)*da

adot = w*zz + (R-r).*k + r.*aa - c;

%% Top shares
ga_cum = sum(g(1,:))*da;
for i=2:I
    ga_cum = [ga_cum; ga_cum(i-1) + sum(g(i,:))*da];
end

b50 = (round(ga_cum,2) <= 0.50);
m40 = (round(ga_cum,2)>0.50 & round(ga_cum,2)<0.90);
p99 = (round(ga_cum,2) >= 0.99);
p90 = (round(ga_cum,2) >= 0.90);
p95 = (round(ga_cum,2) >= 0.95);

wealth = a.*sum(da*g')';
tot_wealth     = sum(wealth);
top1_share     = (wealth'*p99)/tot_wealth;
top5_share     = (wealth'*p95)/tot_wealth;
top10_share    = (wealth'*p90)/tot_wealth;
mid40_share    = (wealth'*m40)/tot_wealth;
bottom50_share = (wealth'*b50)/tot_wealth;

% top shares
wealth_stats = {'Top 1%', round(top1_share,3)*100; ...
                'Top 5%', round(top5_share,3)*100; ...
                'Top 10%', round(top10_share,3)*100; ...
                'Middle 40%', round(mid40_share,3)*100; ...
                'Bottom 50%', round(bottom50_share,3)*100 }
