% ----------------------------------------------------- %
%% The Portfolio Choice Channel of Wealth Inequality   %%
% Based on Liquid-Illiquid LCP code (from Ben Moll) %
% Code to generate and save figures
% Author: Lucas Rosso %
% Date: 18-02-2021 %
% Extremely Preliminar %
% ----------------------------------------------------- %

% load('DATA_04012021')

cd 'Figures'

%% Figure 1: Adj Region

mymap = [0 0 1
    1 0 0];

figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
surf(b,a,adjRegion1','FaceAlpha',0.6)
hold on
view([0 90])
colormap(mymap)
xlabel('Safe Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Risky Wealth (a)','FontSize',14, 'interpreter','latex')
title('Low Productivity ($z_{1}$)','FontSize',14, 'interpreter','latex')
xlim([bmin bmax])
ylim([amin amax])
shading flat
grid off
box on

subplot(1,2,2)
surf(b,a,adjRegion2','FaceAlpha',0.6)
view([0 90])
colormap(mymap)
xlabel('Safe Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Risky Wealth (a)','FontSize',14, 'interpreter','latex')
title('High Productivity ($z_{2}$)','FontSize',14, 'interpreter','latex')
xlim([bmin bmax])
ylim([amin amax])
shading flat
grid off
box on

print -dpng adjRegion.png
% print -depsc adjRegion.eps

%% Figure 2: Adjustment Targets a,b

bAdj_plot = bAdj;
bAdj_plot(adjRegion ~= 1)=NaN;
% bAdj_plot = bAdj_plot ./ bbb - 1;

aAdj_plot = aAdj;
aAdj_plot(adjRegion ~= 1)=NaN;
% aAdj_plot = aAdj_plot ./ aaa - 1;

% b
figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(2,2,1)
surf(b,a,bAdj_plot(:,:,1)')
xlabel('Safe Wealth (b)','FontSize',11, 'interpreter','latex')
title('Safe Adj. Target, Low Income','FontSize',11, 'interpreter','latex')
xlim([bmin bmax])
ylim([amin amax])

subplot(2,2,2)
surf(b,a,bAdj_plot(:,:,2)')
xlabel('Safe Wealth (b)','FontSize',11, 'interpreter','latex')
title('Safe Adj. Target, High Income','FontSize',11, 'interpreter','latex')
xlim([bmin bmax])
ylim([amin amax])

% print -dpng adjTarget_b.png
% print -depsc adjTarget_b.eps
% -----------------------------------------------

% a
subplot(2,2,3)
surf(b,a,aAdj_plot(:,:,1)')
xlabel('Safe Wealth (b)','FontSize',11, 'interpreter','latex')
ylabel('Risky Wealth (a)','FontSize',11, 'interpreter','latex')
title('Risky Adj. Target, Low Income','FontSize',11, 'interpreter','latex')
xlim([bmin bmax])
ylim([amin amax])

subplot(2,2,4)
surf(b,a,aAdj_plot(:,:,2)')
xlabel('Safe Wealth (b)','FontSize',11, 'interpreter','latex')
ylabel('Risky Wealth (a)','FontSize',11, 'interpreter','latex')
title('Risky Adj. Target, High Income','FontSize',11, 'interpreter','latex')
xlim([bmin bmax])
ylim([amin amax])

print -dpng adjTarget.png
% print -depsc adjTarget_a.eps
% -----------------------------------------------

% Figure 3: KF equation
% figure
% set(gcf,'PaperPosition',[0 0 15 5])
% subplot(1,2,1)
% surf(b,a,g(:,:,1)')
% xlabel('Safe Wealth (b)','FontSize',11, 'interpreter','latex')
% ylabel('Risky Wealth (a)','FontSize',11, 'interpreter','latex')
% title('Stationary Density, Low Income','FontSize',11, 'interpreter','latex')
% xlim([0 bmax])
% ylim([0 amax])
% 
% subplot(1,2,2)
% surf(b,a,g(:,:,2)')
% xlabel('Safe Wealth (b)','FontSize',11, 'interpreter','latex')
% ylabel('Risky Wealth (a)','FontSize',11, 'interpreter','latex')
% title('Stationary Density, High Income','FontSize',11, 'interpreter','latex')
% xlim([0 bmax])
% ylim([0 amax])
% 
% print -dpng distribution.png
% print -depsc distribution.eps
% -----------------------------------------------

% Figure: A transition matrix
figure
spy(A)
print -dpng spy_A.png
% print -depsc spy_A.eps
% -----------------------------------------------

%% Wealth Distribution

% Figure: Marginals 
auxi_g = permute(g, [2,1,3]); 

% risky
marg_a(:,1) = sum(g(:,:,1))*db; 
marg_a(:,2) = sum(g(:,:,2))*db;

%safe
marg_b(:,1) = sum(auxi_g(:,:,1))*da; 
marg_b(:,2) = sum(auxi_g(:,:,2))*da;

% figure 
% subplot(1,2,1)
% plot(b,marg_b); grid on, hold on
% legend('Low income', 'High income','location','northeast')
% xlabel('Safe Asset ($b$)','Fontsize',12,'interpreter','latex')
% xlim([0 bmax])
% 
% subplot(1,2,2)
% plot(a,marg_a); grid on, hold on
% legend('Low income', 'High income','location','northeast')
% xlabel('Risky Asset ($a$)','Fontsize',12,'interpreter','latex')
% xlim([0 bmax])
% hold off

% bar graph

% both, editing xlims
aux = (a<18)';
new_marg_a(1:sum(aux),:) = marg_a(1:sum(aux),:);
new_marg_a(sum(aux)+1,1) = sum(marg_a(:,1).*(1-aux));
new_marg_a(sum(aux)+1,2) = sum(marg_a(:,2).*(1-aux));

aux2 = (b<8);
new_marg_b(1:sum(aux2),:) = marg_b(1:sum(aux2),:);
new_marg_b(sum(aux2)+1,1) = sum(marg_b(:,1).*(1-aux2));
new_marg_b(sum(aux2)+1,2) = sum(marg_b(:,2).*(1-aux2));

% Risky Asset b(1:sum(aux2)+1)
figure 
subplot(1,2,1)
a1 = bar(a(1:sum(aux)+1),new_marg_a(:,1),'r'); grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Risky Asset ($a$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_a(a)$','Fontsize',12,'interpreter','latex')
xlim([-da/2 a(sum(aux)+1)+da/2])
text(0.7,0.25,['$\leftarrow Pr(a=0| \, y_{1})=' num2str(round(da*marg_a(1,1),2)) '$'],'FontSize',12,'interpreter','latex','Color','k');

subplot(1,2,2)
a3 = bar(a(1:sum(aux)+1),new_marg_a(:,1),'r'); grid on, hold on
a2 = bar(a(1:sum(aux)+1),new_marg_a(:,2),'b'); grid on, hold on
set(gca, 'GridLineStyle','--')
set(a3, 'visible','off')
xlabel('Risky Asset ($a$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_a(a)$','Fontsize',12,'interpreter','latex')
xlim([-da/2 a(sum(aux)+1)+da/2])
text(0.7,0.25,['$\leftarrow Pr(a=0| \, y_{2})=' num2str(round(da*marg_a(1,2),2)) '$'],'FontSize',12,'interpreter','latex','Color','k');
legend('Low income','High income','location','northeast')
hold off

print -dpng marg_risky.png

% -----------------------------------------------

% Safe Asset
figure 
subplot(1,2,1)
a1 = bar(b(1:sum(aux2)+1),new_marg_b(:,1),'r'); grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Asset ($b$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_b(b)$','Fontsize',12,'interpreter','latex')
xlim([bmin-(db/2) b(sum(aux2)+1)+db/2])

subplot(1,2,2)
a3 = bar(b(1:sum(aux2)+1),new_marg_b(:,1),'r'); grid on, hold on
a2 = bar(b(1:sum(aux2)+1),new_marg_b(:,2),'b'); grid on, hold on
set(gca, 'GridLineStyle','--')
set(a3, 'visible','off')
xlabel('Safe Asset ($b$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_b(b)$','Fontsize',12,'interpreter','latex')
xlim([bmin-(db/2) b(sum(aux2)+1)+db/2])
legend('Low income','High income','location','northeast')
hold off

print -dpng marg_safe.png

clear aux aux2 new_marg_a new_marg_b

% -----------------------------------------------

% both  
auxi_ga = permute(g, [1,3,2]);
auxi_gb = permute(g, [2,3,1]);
marginal_a = reshape(sum(sum(auxi_ga))*db,1,I);
marginal_b = reshape(sum(sum(auxi_gb))*da,1,J);

figure 
subplot(1,2,1)
a1 = bar(b,marginal_b,'r'); grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Asset ($b$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_b(b)$','Fontsize',12,'interpreter','latex')
xlim([bmin-(db/2) bmax])
% text(0.2,0.9,['$\leftarrow Pr(b=0)=' num2str(round(db*marginal_b(1),2)) '$'],'FontSize',14,'interpreter','latex','Color','k');

subplot(1,2,2)
a2 = bar(a,marginal_a,'b'); grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Risky Asset ($a$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_a(a)$','Fontsize',12,'interpreter','latex')
xlim([-da/2 amax])
text(2,0.6,['$\leftarrow Pr(a=0)=' num2str(round(da*marginal_a(1),2)) '$'],'FontSize',14,'interpreter','latex','Color','k');
hold off
print -dpng marg_both.png

% both, editing xlims
aux = (a<18);
new_marginal_a(1:sum(aux)) = marginal_a(1:sum(aux));
new_marginal_a(sum(aux)+1) = sum(marginal_a.*(1-aux));

aux2 = (b<8)';
new_marginal_b(1:sum(aux2)) = marginal_b(1:sum(aux2));
new_marginal_b(sum(aux2)+1) = sum(marginal_b.*(1-aux2));

figure 
subplot(1,2,1)
a1 = bar(b(1:sum(aux2)+1),new_marginal_b,'r'); grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Asset ($b$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_b(b)$','Fontsize',12,'interpreter','latex')
xlim([bmin-(db/2) b(sum(aux2)+1)+db/2])
% text(0.2,0.9,['$\leftarrow Pr(b=0)=' num2str(round(db*marginal_b(1),2)) '$'],'FontSize',14,'interpreter','latex','Color','k');

subplot(1,2,2)
a2 = bar(a(1:sum(aux)+1),new_marginal_a,'b'); grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Risky Asset ($a$)','Fontsize',12,'interpreter','latex')
ylabel('Marginal Dist., $g_a(a)$','Fontsize',12,'interpreter','latex')
xlim([-da/2 a(sum(aux)+1)+da/2])
text(0.8,0.6,['$\leftarrow Pr(a=0)=' num2str(round(da*marginal_a(1),2)) '$'],'FontSize',14,'interpreter','latex','Color','k');
hold off
print -dpng marg_both2.png

clear aux aux2 new_marginal_a new_marginal_b

% only risky asset, both income types
% figure
% bar(a,marginal_a,'b'); grid on, hold on
% xlabel('Risky Asset ($a$)','Fontsize',12,'interpreter','latex')
% xlim([-da/2 amax])
% text(2,0.5,['$\leftarrow Pr(a=0)=' num2str(round(da*marginal_a(1),2)) '$'],'FontSize',14,'interpreter','latex','Color','k');
% hold off

%% Kappa and the Merton Rule

a_final = aAdj.*adjRegion + aaa.*(adjRegion==0);
b_final = bAdj.*adjRegion + bbb.*(adjRegion==0);
wealth = a_final + b_final;
share_a = a_final./wealth;

share_a_aux = permute(share_a, [2,1,3]);

% figure 
% subplot(1,2,1)
% plot(a(1:(J-1)), share_a_aux(1:(J-1),2:J,1)), grid on, hold on
% set(gca, 'GridLineStyle','--')
% xlabel('Risky Wealth, $a$','FontSize',12,'Interpreter','latex')
% ylabel('Risky Share $\omega$','FontSize',12,'Interpreter','latex')
% title('Low Income ($y_{1}$)','FontSize',12,'Interpreter','latex')
% 
% subplot(1,2,2)
% plot(a(1:(J-1)), share_a_aux(1:(J-1),2:I,2)), grid on, hold on
% set(gca, 'GridLineStyle','--')
% xlabel('Risky Wealth, $a$','FontSize',12,'Interpreter','latex')
% ylabel('Risky Share $\omega$','FontSize',12,'Interpreter','latex')
% title('High Income ($y_{2}$)','FontSize',12,'Interpreter','latex')
% hold off
% print -dpng risky_share_policy.png
% -----------------------------------------------


%% Expected returns across wealth dist.
% exp_r = (mu - rbPos).*share_a + rbPos;
% 
% figure 
% subplot(1,2,1)
% plot(a, exp_r(:,:,1)), grid on, hold on
% set(gca, 'GridLineStyle','--')
% xlabel('Risky Wealth, $a$','FontSize',12,'Interpreter','latex')
% ylabel('Expected Return to Wealth','FontSize',12,'Interpreter','latex')
% title('Low Income ($y_{1}$)','FontSize',12,'Interpreter','latex')
% 
% subplot(1,2,2)
% plot(a, exp_r(:,:,2)), grid on, hold on
% set(gca, 'GridLineStyle','--')
% xlabel('Risky Wealth, $a$','FontSize',12,'Interpreter','latex')
% ylabel('Expected Return to Wealth','FontSize',12,'Interpreter','latex')
% title('High Income ($y_{2}$)','FontSize',12,'Interpreter','latex')
% hold off

%% Savings and Consumption Policy Function

% savings
figure 
subplot(1,2,1)
plot(b, ss(:,:,1)), grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Wealth, $b$','FontSize',12,'Interpreter','latex')
ylabel('Savings in Safe Asset','FontSize',12,'Interpreter','latex')
xlim([0 25])
title('Low Income ($y_{1}$)','FontSize',12,'Interpreter','latex')

subplot(1,2,2)
plot(b, ss(:,:,2)), grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Wealth, $b$','FontSize',12,'Interpreter','latex')
ylabel('Savings in Safe Asset','FontSize',12,'Interpreter','latex')
xlim([0 25])
title('High Income ($y_{2}$)','FontSize',12,'Interpreter','latex')
hold off
print -dpng savings_policy.png


% consumption
figure 
subplot(1,2,1)
plot(b, c(:,:,1)), grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Wealth, $b$','FontSize',12,'Interpreter','latex')
ylabel('Consumption','FontSize',12,'Interpreter','latex')
xlim([0 25])
title('Low Income ($y_{1}$)','FontSize',12,'Interpreter','latex')

subplot(1,2,2)
plot(b, c(:,:,2)), grid on, hold on
set(gca, 'GridLineStyle','--')
xlabel('Safe Wealth, $b$','FontSize',12,'Interpreter','latex')
ylabel('Consumption','FontSize',12,'Interpreter','latex')
xlim([0 25])
title('High Income ($y_{2}$)','FontSize',12,'Interpreter','latex')
hold off

print -dpng consumption_policy.png


%% Wealth figures and stats
wealth_vec = reshape(aa+bb,I*J,1);
g_auxiliar = [g_stacked(1:I*J),g_stacked(I*J +1:I*J*Nz)];
g_cum = [sum(g_auxiliar')]'*da*db;

WEALTH = [wealth_vec, g_cum];
WEALTH = sortrows(WEALTH,1);

gcum = WEALTH(:,2);
wealth_vec = WEALTH(:,1);

g_sum = gcum(1);
for j = 2:I*J
    g_sum = [g_sum; g_sum(j-1)+gcum(j)];
end

% g_sum = round(g_sum, 2);

% indicators

% kappa as a share of initial stock of adjusters
ka_pctg_cost = ka/(sum(sum(sum(adjRegion.*(aaa+bbb))))/sum(sum(sum(adjRegion))))

b50 = (round(g_sum,2) <= 0.50);
m40 = (round(g_sum,2)>0.50 & round(g_sum,2)<0.90);
p90 = (round(g_sum,2) >= 0.90);
p95 = (round(g_sum,2) >= 0.95);
p99 = (round(g_sum,2) >= 0.99);

% risky wealth
risky_wealth = sum(sum(sum(a_final.*g)))*da*db

% total wealth
total_wealth = wealth_vec'*gcum

% top 1%
top1 = wealth_vec'*(gcum.*p99);
top1_share = top1/total_wealth;

% top 5%
top5 = wealth_vec'*(gcum.*p95);
top5_share = top5/total_wealth;

% top 10%
top10 = wealth_vec'*(gcum.*p90);
top10_share = top10/total_wealth;

% middle 40%
mid40 = wealth_vec'*(gcum.*m40);
mid40_share = mid40/total_wealth;
% bottom 50%
bottom50 = wealth_vec'*(gcum.*b50);
bottom50_share = bottom50/total_wealth;

wealth_stats = {'Top 1%', round(top1_share,3)*100; ...
                'Top 5%', round(top5_share,3)*100; ...
                'Top 10%', round(top10_share,3)*100; ...
                'Middle 40%', round(mid40_share,3)*100; ...
                'Bottom 50%', round(bottom50_share,3)*100 }
            
%% avg risky share
risky_share_data = [0.00185908,0.02011726,0.06025168,0.13573499,0.21601109,0.28422712,0.32045861,0.38150663,0.43721997,0.52562847]

share_a = min(max(reshape(share_a, I*J,2),0),1); % restricting portfolio to (0,1)
share_a = mean(share_a')';

WEALTH = [WEALTH, share_a];
WEALTH = sortrows(WEALTH,1);

aux = (g_sum<1/20);
rshare = sum((share_a.*gcum.*aux)/sum(gcum.*aux));
for i=2:20
    aux = (g_sum>(i-1)/20 & g_sum<i/20);
    rshare(i) = sum((share_a.*gcum.*aux)/sum(gcum.*aux));
end

merton_rule = (mu-rbPos)/(s2*s) %,'LineWidth',1.5,'k--'

figure
plot((1:20)*5,merton_rule*ones(20,1),'k--','LineWidth',1.5); grid on, hold on
plot((1:20)*5,rshare,'o-r','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','r'); grid on, hold on
% plot((1:10)*10,round(risky_share_data,2),'^-b','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','b'); grid on, hold on
set(gca, 'GridLineStyle','--')
legend('Merton Rule','Model','location','northwest')
xlabel('Financial Wealth Percentile','FontSize',12, 'interpreter','latex')
ylabel('Risky Asset Share','FontSize',12, 'interpreter','latex')
xlim([0.05 100])
hold off
print -dpng risky_share_deciles.png

%% gini (in stata)
% T = table(WEALTH);
% writetable(T,'gini_data.xlsx');
