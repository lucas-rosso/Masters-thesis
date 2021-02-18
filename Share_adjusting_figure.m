% ---------------------------------------------- %
%% Share of households adjusting portfolio %%
% Author: Lucas Rosso %
% Date: 18-02-2021 %
% Extremely Preliminar %
% ---------------------------------------------- %

clear all; close all; clc;

load('calibration.mat') % internally calibrated parameters 
tic 

K = [0.1,optimal_params,(3:10)/10];

adj_share = zeros(length(K),1);

for k = 1:length(K)
%% Parameters and initial objetcts

UseNoAdjustmentAsGuess = 1;
AssumeConcavity = 0; % this is relevant in the part with adjustment
BilinW = 1; % when all neighboring points of (b',a') are non-adjustment points, use bilinear weights in dividing mass to be rebalanced 
            % across them in the KF algorithm

% Parameters:
s     = 2;         % risk aversion in CRRA utility
rbPos = 0.02;      % return on liquid asset (e.g. checking account)
rbNeg = 0.10;
rho   = 0.053;     % discount rate (beta = 0.96)

ka    = K(k)       % cost of adjustment (assuming the wage w is 1)

% Safe Asset
I     = 100; 
bmin  = -1;
bmax  = 25;
b     = linspace(bmin,bmax,I)'; % column vector
db    = b(2)-b(1);

% Risky asset 
J     = 100; 
amin  = 0;
amax  = 50; 
a     = linspace(amin,amax,J); % row vector
da    = a(2)-a(1);
da2 = da^2;

% Drift and Variance of the stock return
mu    = 0.06;  
s2    = 0.18^2; % Gomes and Michaelides (JF, 2005)                   

% Poisson shocks:
Nz    = 2; % number of exog. states
the = -log(0.9);
sigma2 = (0.04*2*the)/(1-exp(-2*the));
sigma = sqrt(sigma2);

z1 = 1 - sigma;
z2 = 1 + sigma;

l_z1 = log(z1); l_z2 = log(z2);

p_z = sqrt(the/(pi*sigma2*(1-exp(-2*the))))*exp(-(the*(l_z1-(l_z2*exp(-the)))^(2))/(sigma2*(1-exp(-2*the))));

z     = [round(z1,2),round(z2,2)]; % state values (row vector)
la    = [round(p_z,2),round(p_z,2)]; % Poisson intensities

% Create mesh: b in 1st dimension, a in 2nd dimension, z in 3rd dimension
bb = b*ones(1,J); % repeat b column vector J times in the 2nd dimension
aa = ones(I,1)*a; % repeat a row vector I times in the 1st dimension
zz = ones(J,1)*z; % repeat z row vector J times in the 1st dimension (will then take the 2 columns with 1st and 2nd state values separately)

bbb = zeros(I,J,Nz); aaa = zeros(I,J,Nz); zzz = zeros(I,J,Nz); % these are the meshes
bbb(:,:,1) = bb; bbb(:,:,2) = bb;
aaa(:,:,1) = aa; aaa(:,:,2) = aa;
zzz(:,:,1) = z(1);zzz(:,:,2)= z(2);
    
%allow for differential borrowing and lending rates (only matters when bmin<0)
Rb = rbPos.*(bbb>=0) + rbNeg.*(bbb<0); % mesh of returns on liquid asset to allow for different borrowing rate if borrowing is allowed

L = I*J*Nz; % total number of nodes

% Adjustment decision set-up:
Nx      = 600; % number of grid points for cash-in-hand when adjusting
NaP     = 600; % number of grid points for a' when adjusting
xmin      = amin + bmin + ka; % minimum cash-in-hand
xmax      = amax + bmax; % maximum
x     = linspace(xmin,xmax,Nx);

sumGrid_aux   = aa + bb; % mesh of cash-in-hand from (b, a) nodes
sumGrid_aux   = sumGrid_aux(:); % transform to column vector

sumGrid = sumGrid_aux;

% Iteration set-up:
Delta     = 1000; % step size (for time) in the implicit method (as opposed to db and db which are steps in the state dimensions)
maxit     = 1000;
tol       = 1e-6;
iter      = 0;
dist      = zeros(maxit,1); % hold errors during iteration

% Initial guess
V0(:,:,1) = (1-s)^(-1)*(z(1) + rbPos.*bb).^(1-s)/rho; % assuming u(c) = c^(1 - s)/(1 - s). This is similar to value of "staying put" u(z + rb*b)/rho                                                          
V0(:,:,2) = (1-s)^(-1)*(z(2) + rbPos.*bb).^(1-s)/rho;
v         = V0;  

%  tau = 15; % from two_assets_kinked.m, this is justified as: if ra>>rb, impose tax on ra*a at high a, otherwise some households accumulate infinite illiquid wealth
%            % (not needed if ra is close to - or less than - rb). It implies aDrift = ra*a*(1 - (a/(amax * 0.999))^(tau - 1))
%  tau0 = mu.*(amax*.999)^(1-tau); % in two_assets_kinked.m the multiplier is 1.33 instead of 0.999. That multiplier makes the tax much less aggressive, as with 
%                                  % 0.999 the drift of a goes down to slightly negative towards amax, while with 1.33 the drift becomes just slightly concave near amax
%  T = tau0*a.^tau;
%  aDrift = mu*a - T;
% plot(a,aDrift,a,zeros(J,1),a,ra.*a)

aDrift = mu*a;
diff  = s2*a.^2;

%% Build matrix Aswitch summarizing evolution of a % i.e. states which are exogenous in the no adjustment case: risky asset, its return and the income process
chi = diff/(2*da2); % this is analogous to X, Y, Z constructed below for b (but call chi instead of x since x is cash-in-hand). This (and below) has 
                         % similarities with p. 16/17 in Achdou et al. 2017 - Online Appendix
yy  = - aDrift/da - diff/da2;
zeta= aDrift/da + diff/(2*da2);

%This will be the upperdiagonal of the A_switch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined (since this will be the Ith upper diagonal, the first I elements are dropped)
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)]; % this is analogous to zeta(j) * ones(I, 1) i.e. repeating zeta(j) I times in a vector
end

%This will be the center diagonal of the A_switch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the A_switch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)]; % this will have length I*(J-1) because the last I elements in constructing the Ith lower diagonal will be dropped anyway
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Aaux=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

Aswitch = [Aaux - speye(I*J)*la(1), speye(I*J)*la(1);speye(I*J)*la(2),Aaux - speye(I*J)*la(2)]; % this adds the transitions between z1 and z2 (keep in mind they have
                                                                                                                % to balance along the row and sum to 0)
                                                                                                              
%% SOLVE WITHOUT ADJUSTMENT
if UseNoAdjustmentAsGuess == 1

for n=1:maxit
    V = v;
    
    % forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db; % remember b changes along the 1st dimension (column vector)
    Vbf(I,:,:,:) = (zzz(I,:,:) + Rb(I,:,:).*bmax).^(-s); %will never be used, but impose state constraint b<=bmax just in case
                                                           % this is v'b(I,j,k)_F = u'(z(k) + rb(I) * b(I))
    % backward difference
    Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
    Vbb(1,:,:) = (zzz(1,:,:) + Rb(1,:,:).*bmin).^(-s); %state constraint boundary condition: v'b(1,j,k)_B = u'(z(k) + rb(1) * b(1))

    Vbf = max(Vbf,10^(-6)); % this is to avoid dealing with too small numbers in parts where v is almost flat wrt b
    Vbb = max(Vbb,10^(-6));
    
     %consumption and savings with forward difference
     cf = Vbf.^(-1/s); % from FOC
     sf = zzz + Rb.*bbb - cf;
     %consumption and savings with backward difference
     cb = Vbb.^(-1/s);
     sb = zzz + Rb.*bbb - cb;
     %consumption and derivative of value function at steady state
     c0 = zzz + Rb.*bbb; % this is where s=0 i.e. db = 0
     
     % make a choice of forward or backward differences based on the sign of the drift (upwind scheme)
     Ib = (sb < 0); %negative drift --> backward difference. This is where we assume concavity, because we don't check sf <= 0 as well (see *** reference)
     If = (sf > 0).*(1-Ib); %positive drift --> forward difference. Note that the (1-Ib) makes sure we take only points where also sb>=0, so Ib and If are 
                               % disjoint. Concavity of v wrt b would imply that, because with concavity sf <= sb, but this enforces it in a setting
                               % where there could be some numerical instability
     I0 = (1-If-Ib); %at steady state. These have to be points where sb>=0 & sf<=0. If we didn't impose disjointness in the previous step, this could = -1 where
                        % (sb < 0 & sf > 0) and we would need to split these points (see *** reference). Instead we assume concavity i.e. these points don't arise 

     c = cf.*If + cb.*Ib + c0.*I0;
     u = c.^(1-s)/(1-s);
     
     %CONSTRUCT MATRIX AA (see p. 5 Achdou et al. 2017 - Online Appendix)
     X = -Ib.*sb/db; % lower
     Y = -If.*sf/db + Ib.*sb/db; % center
     Z = If.*sf/db; % upper    

    updiag = zeros(I*J,Nz); % notice that endogenous nodes (b, a) are listed in the 1st dimension, while exogenous state z is listed in the 2nd
    lowdiag = zeros(I*J,Nz);
    centdiag = zeros(I*J,Nz);

    % center diagonal
    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    % lower and upper diagonals
    % (* reference)
    lowdiag(1:I-1,:) = X(2:I,1,:); % exclude the 1st row i.e. b(1). This is only at a(1) (same for updiag below). Notice lowdiag(I,:)=0
    updiag(2:I,:) = Z(1:I-1,1,:); % exclude the last row i.e. b(I). Notice upperdiag(1,:)=0
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)]; % this stacks the columns of X in the right way (adding a 0 at the end, 
                                                                                       % with size [1, Nz] to account for exogenous states)
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end

AA1=spdiags(centdiag(:,1),0,I*J,I*J)+spdiags(updiag(:,1),1,I*J,I*J)+spdiags(lowdiag(:,1),-1,I*J,I*J); % notice that when filling the upper diagonal
                                                                                       % the first element of the vector is dropped. When filling the lower diagonal
                                                                                       % the last element is dropped. This is consistent with how lowdiag and updiag
                                                                                       % are constructed initially (see * reference). Note: removed unneeded ; 0]
AA2=spdiags(centdiag(:,2),0,I*J,I*J)+spdiags(updiag(:,2),1,I*J,I*J)+spdiags(lowdiag(:,2),-1,I*J,I*J);       
AA = [AA1, sparse(I*J,I*J); sparse(I*J,I*J), AA2]; % this completes the I*J matrixes to bring them to I*J*Nz
A = AA + Aswitch; % adds the transition of "exogenous states" (in the no adjustment case, risky asset evolves exogenously)

 % checking that rows in A sum to zero
error = full(sum(A,2));
check = abs(error)>1e-9;
check_ = sum(check);

if max(abs(sum(A,2)))>10^(-9)
disp('Improper Transition Matrix')
max(abs(sum(A,2)))
how_many = check_ 
break
end
 
  B = (1/Delta + rho)*speye(I*J*Nz) - A; % this and the next 3 steps implement eq. 15 p.6 Achdou et al. 2017 - Online Appendix

  u_stacked = [reshape(u(:,:,1),I*J,1);reshape(u(:,:,2),I*J,1)]; % stack into vector I*J*Nz
  V_stacked = [reshape(V(:,:,1),I*J,1);reshape(V(:,:,2),I*J,1)];

  vec = u_stacked + V_stacked/Delta;

  %IMPLICIT UPDATING
  V_stacked_12 = B\vec; %SOLVE SYSTEM OF EQUATIONS

  V(:,:,1) = reshape(V_stacked_12(1:I*J),I,J); % bring back into mesh form
  V(:,:,2) = reshape(V_stacked_12(I*J+1:I*J*Nz),I,J);

  Vchange = V - v;
  v = V;
  
  dist(n) = max(max(max(max(abs(Vchange)))));
  disp(['Initial guess, no adjustment, iter ', int2str(n), ', dist ' num2str(dist(n)) ]);
      if dist(n)<tol
         disp(['Value Function Converged, Iteration = ' int2str(n)]);
         break
      end
end
end

%% SOLVE WITH ADJUSTMENT

for n=1:maxit

    % except for Hf, Hb, H0 and the AssumeConcavity == 0, this part (up to adjustment decision for x grid) is the same as above without adjustment
    V = v;
    % forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db; % remember b changes along the 1st dimension (column vector)
    Vbf(I,:,:) = (zzz(I,:,:) + Rb(I,:,:).*bmax).^(-s); %will never be used, but impose state constraint b<=bmax just in case
                                                           % this is v'b(I,j,k)_F = u'(z(k) + rb(I) * b(I))
    % backward difference
    Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
    Vbb(1,:,:,:) = (zzz(1,:,:) + Rb(1,:,:).*bmin).^(-s); %state constraint boundary condition: v'b(1,j,k)_B = u'(z(k) + rb(1) * b(1))

    Vbf = max(Vbf,10^(-6)); % this is to avoid dealing with too small numbers in parts where v is almost flat wrt b
    Vbb = max(Vbb,10^(-6));
    
    %consumption and savings with forward difference
    cf = Vbf.^(-1/s);
    sf = zzz + Rb.*bbb - cf;
    Hf = cf.^(1-s)/(1-s) +Vbf.*sf;
    
    %consumption and savings with backward difference
    cb = Vbb.^(-1/s);
    sb = zzz + Rb.*bbb - cb;
    Hb = cb.^(1-s)/(1-s) +Vbb.*sb;
    
    %consumption and derivative of value function at steady state
    c0 = zzz + Rb.*bbb;
    %hamiltonians;
    
    % make a choice of forward or backward differences based on the sign of the drift
    if AssumeConcavity == 1 % same as above
        Ib = (sb < 0); %negative drift --> backward difference
        If = (sf > 0).*(1-Ib); %positive drift --> forward difference
        I0 = (1-If-Ib); %at steady state
    else % (*** reference)
        I0 = (1-(sf>0)) .* (1-(sb<0)); % points where sf <= 0 & sb >= 0
        Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0); % points where (sf <= 0 & sb < 0) U (sf > 0 & sb >= 0). "Unique" because we know here we can easily take sb
                                                           % over (sf <= 0 & sb < 0) and sf over (sf > 0 & sb >= 0)
        Iboth = (sb<0).*(sf>0);
        Ib = Iunique.*(sb<0) + Iboth.*(Hb>Hf); % split the points where (sb < 0 & sf > 0) based on u(cb) + Vbb * sb </> u(cf) + Vbf * sf: if > then backward. Note
                                               % that the other elements in the HJB are the same between backward and forward
        If = Iunique.*(sf>0) + Iboth.*(Hb<=Hf);
    end
    
    c = cf.*If + cb.*Ib + c0.*I0;
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = -Ib.*sb/db;
    Y = -If.*sf/db + Ib.*sb/db;
    Z = If.*sf/db;
    
    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    lowdiag(1:I-1,:) = X(2:I,1,:);
    updiag(2:I,:) = Z(1:I-1,1,:);
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end

AA1=spdiags(centdiag(:,1),0,I*J,I*J)+spdiags(updiag(:,1),1,I*J,I*J)+spdiags(lowdiag(:,1),-1,I*J,I*J);
    AA2=spdiags(centdiag(:,2),0,I*J,I*J)+spdiags(updiag(:,2),1,I*J,I*J)+spdiags(lowdiag(:,2),-1,I*J,I*J);
    AA = [AA1, sparse(I*J,I*J); sparse(I*J,I*J), AA2];
A = AA + Aswitch; 

 % checking that rows in A sum to zero
error = full(sum(A,2));
check = abs(error)>1e-9;
check_ = sum(check);

if max(abs(sum(A,2)))>10^(-9)
disp('Improper Transition Matrix')
max(abs(sum(A,2)))
how_many = check_ 
break
end

%%%%% Adjustment decision for x grid: this constructs objects on the predetermined grid for x, before interpolating at the x values corresponding to nodes (b, a)
vstarAux = zeros(Nz,Nx); % The Aux elements are those defined on the grid for x
bAdjAux  = zeros(Nz,Nx);
aAdjAux  = zeros(Nz,Nx);
vAdj = zeros(Nz,NaP); % for each x in the grid, search over a grid for a' - backing out b' = x - a' - ka (see ** reference) - to find the maximum. This vadj
                         % holds the values over such grid for a'
aPmin = amin;
G1 = griddedInterpolant(bb,aa,V(:,:,1)); % this is the interpolant for v over the nodes (b, a) at z(1)
G2 = griddedInterpolant(bb,aa,V(:,:,2));

    for i = 1:Nx
        aPmax = min(x(i) - ka - bmin,amax); %CONSTRAINT a'<=aMax, taking into account the minimum value for b
        aP = linspace(aPmin,aPmax,NaP);
        aP = max(aP,x(i)-ka-bmax); %CONSTRAINT b'<=bMax
        bP = x(i) - ka - aP; % (** reference)
        vAdj(1,:) = G1(bP',aP'); % collect values over the grid for a' (and implied grid for b', given x)
        vAdj(2,:) = G2(bP',aP');

        [VstarAux(1,i),idx1] = max(vAdj(1,:)); % find maximum and argmax, and store maximum
        [VstarAux(2,i),idx2] = max(vAdj(2,:));
        aAdjAux(1,i) = aP(idx1); aAdjAux(2,i) = aP(idx2); % store argmax
        bAdjAux(1,i) = bP(idx1); bAdjAux(2,i) = bP(idx2);
    end
    
%%%%% Interpolate for all possible values of a+b:
Vstar = [lininterp1(x',VstarAux(1,:)',sumGrid);lininterp1(x',VstarAux(2,:)',sumGrid)]; % recall that sumGrid is the vector of x=b+a over nodes (b,a)

      
% SOLVE USING LCP  
    
  B = (1/Delta + rho)*speye(I*J*Nz) - A;

    u_stacked = [reshape(u(:,:,1),I*J,1);reshape(u(:,:,2),I*J,1)];
    V_stacked = [reshape(V(:,:,1),I*J,1);reshape(V(:,:,2),I*J,1)];
    
    vec = u_stacked + V_stacked/Delta; % from B = ... to here it's again the same as in the without adjustment case
    
  q = -vec + B*Vstar;
    
    %using Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
    z0 = V_stacked-Vstar; lbnd = zeros(I*J*Nz,1); ubnd = Inf*ones(I*J*Nz,1); % set initial value z0 where we guess v(n+1)=v(n) and lower/upper bound for solution z
    z = LCP(B,q,lbnd,ubnd,z0,0); % last argument is display of iteration data
   
    LCP_error = max(abs(z.*(B*z + q))); % check the accuracy of the LCP solution (should all be 0)
    if LCP_error > 10^(-6)
        disp('LCP not solved')
        %break
    end
    
    V_stacked = z+Vstar; %calculate value function
    
    % bring back into mesh form (exactly as in no adjustment case)
    V(:,:,1) = reshape(V_stacked(1:I*J),I,J); % this part is again the same as in the without adjustment case
    V(:,:,2) = reshape(V_stacked(I*J+1:I*J*Nz),I,J);
    
  Vchange = V - v;
  v = V;

    dist(n) = max(max(max(abs(Vchange))));
    disp(['With adjustment, adjustment, iter ', int2str(n), ', dist ' num2str(dist(n)) ]);
    if dist(n)<tol
        disp(['Value Function Converged, Iteration = ' int2str(n)]);
        break
    end

end 

toc


adj = V_stacked==Vstar; %INDICATOR FOR ADJUSTING
adj = (abs(V_stacked - Vstar)<10^(-6)); % allow for some tolerance in defining solution for adjustment decision

adjRegion1 = reshape(adj(1:I*J),I,J).*ones(I,J); % multiplication by ones(I,J) transforms logicals into floats
adjRegion2 = reshape(adj(I*J+1:end),I,J).*ones(I,J);
adjRegion(:,:,1)=adjRegion1; adjRegion(:,:,2)=adjRegion2; % this reshapes the adj. decision into (I, J, Nz)

adj_share(k) = sum(sum(sum(adjRegion)))/(I*J*Nz);
end

adj_share = [100;adj_share*100];

cd 'Figures'
figure
hold on
plot(optimal_params*ones(length(K)+1,1),[0,K]*100,'k--','LineWidth',1.5);
plot([0,K],adj_share,'o-r','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','r'); grid on
set(gca, 'GridLineStyle','--')
hold off
% legend('Low Income','High Income','Location','NorthEast');
% v=vline(0.2,'k--');
% set(v,'LineWidth', 3);
xlim([0 1]);
ylim([0 100]);
xlabel('Fixed Cost ($ \kappa $)','FontSize',12,'interpreter', 'latex');
ylabel('Percentage of Households Adjusting','FontSize',12,'interpreter', 'latex');

print -dpng adj_Share.png
