% ---------------------------------------------- %
%% Two asset Portfolio Allocation Model %%
% Author: Lucas Rosso, based on Ben Moll's code with non-uniform grid %
% Effect of permanent income shock on top shares
% Date: 15-02-2021 %
% Extremely Preliminar %
% ---------------------------------------------- %
%% AIYAGARI + FAT TAIL

clear all; clc; close all;

tic;

ga = 2;
rho = 0.053;

% Guvenen et al. (2019)
the = -log(0.9);
sigma2 = (0.0324*2*the)/(1-exp(-2*the));
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

%non-uniform grid
% x = linspace(0,1,I)';
% coeff = 5; power = 10;
% xx  = x + coeff*x.^power;
% xmax = max(xx); xmin = min(xx);
% a = (amax-amin)/(xmax - xmin)*xx + amin;
% daf = ones(I,1);
% dab = ones(I,1);
% daf(1:I-1) = a(2:I)-a(1:I-1);
% dab(2:I) = a(2:I)-a(1:I-1);
% daf(I)=daf(I-1); dab(1)=dab(2);

a = linspace(amin,amax,I)';
da = a(2)-a(1);
da2 = da^2;

aa = [a,a];
% daaf = daf*ones(1,2);
% daab = dab*ones(1,2);

%objects for approximation of second derivatives
% denom = 0.5*(daaf + daab).*daab.*daaf;
% weightf = daab./denom;
% weight0 = -(daab + daaf)./denom;
% weightb = daaf./denom;


zz = ones(I,1)*z;
kk = ones(I,2);

% safe asset
rbPos = 0.02;      % return on safe asset 
rbNeg = 0.08; % 6 percent wedge as in Kaplan, Moll and Violante (AER, 2018)
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

% v0(:,1) = (w*z(1) + r(:,1).*a + a.*((R-r(:,1)).^2)./(ga*sig2)).^(1-ga)/(1-ga)/rho;
% v0(:,2) = (w*z(2) + r(:,2).*a + a.*((R-r(:,2)).^2)./(ga*sig2)).^(1-ga)/(1-ga)/rho;

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
%     dV2b(2:I-1,:) = (daab(2:I-1,:).*V(3:I,:) - (daab(2:I-1,:)+daaf(2:I-1,:)).*V(2:I-1,:) + daaf(2:I-1,:).*V(1:I-2,:))./denom(2:I-1,:);
%     dV2f(2:I-1,:) = (daab(2:I-1,:).*V(3:I,:) - (daab(2:I-1,:)+daaf(2:I-1,:)).*V(2:I-1,:) + daaf(2:I-1,:).*V(1:I-2,:))./denom(2:I-1,:);
%     dV2b(I,:) = -ga*dVb(I,:)/amax;
%     dV2f(I,:) = -ga*dVf(I,:)/amax;
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
%     X = -Ib.*ssb./daab + sig2/2.*k.^2.*weightb;
%     Y = - If.*ssf./daaf + Ib.*ssb./daab + sig2/2.*k.^2.*weight0;
%     Z = If.*ssf./daaf + sig2/2.*k.^2.*weightf; 

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
% da_tilde = 0.5*(dab + daf);
% da_tilde(1) = 0.5*daf(1); da_tilde(I) = 0.5*dab(I);
% da_stacked = [da_tilde;da_tilde];
% grid_diag = spdiags(da_stacked,0,2*I,2*I);

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

% 
% gg = grid_diag\g_tilde; %convert from \tilde{g} to g
% 
% g = [gg(1:I),gg(I+1:2*I)];
% 
% check1 = g(:,1)'*da_tilde;
% check2 = g(:,2)'*da_tilde;
% check = sum(g)*da

%CALCULATE THEORETICAL POWER LAW EXPONENT
% zeta = ga*(2*sig2*(rho - r)/(R-r)^2 -1);

adot = w*zz + (R-r).*k + r.*aa - c;

% risky_share = (R-r)./(ga*sig2);
% cslope = (rho - (1-ga)*r)./ga - (1-ga)/(2*ga)*(R-r).^2/(ga*sig2);
% sbar = (r-rho)./ga + (1+ga)/(2*ga)*(R-r).^2/(ga*sig2);

% plot(a,k,a,risky_share.*a)
% plot(a,c,a,cslope.*a)

%COMPUTE DISTRIBUTION OF x=log(a)
% for i=1:I
%     G(i,1) = sum(g(1:i,1).*da_tilde(1:i));
%     G(i,2) = sum(g(1:i,2).*da_tilde(1:i));
% end

% f = zeros(I,2);
% x = log(max(a,0));
% dx = zeros(I,1);
% for i =2:I
% dx(i) = x(i)-x(i-1);
% f(i,1) = (G(i,1)-G(i-1,1))/dx(i);
% f(i,2) = (G(i,2)-G(i-1,2))/dx(i);
% end
% f(1)=0;
% 
% xmin = log(1); xmax = log(amax);
% ga = g(1:i,1) + g(1:i,2);
% fx = f(:,1)+f(:,2);


% amax1 = 50;
% amin1 = amin-2;

% figure
% h1 = plot(a,g(:,1),'b',a,g(:,2),'r','LineWidth',2)
% set(gca,'FontSize',16)
% legend(h1,'g_1(a)','g_2(a)')
% text(-0.555,-.015,'$\underline{a}$','FontSize',16,'interpreter','latex')
% line([amin amin], [0 max(max(g))],'Color','Black','LineStyle','--')
% xlabel('Wealth, $a$','interpreter','latex')
% ylabel('Densities, $g_i(a)$','interpreter','latex')
% xlim([amin1 amax1])
% ylim([0 0.3])
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% print -depsc densities_fat.eps

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
