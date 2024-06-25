clear
clc
close all

outdir = fullfile('..','output');
figdir = fullfile('..','figures');

JJ = loadArray(fullfile(outdir,'JJ.txt'),[1,1]);
NA = loadArray(fullfile(outdir,'NA.txt'),[1,1]);
NP = loadArray(fullfile(outdir,'NP.txt'),[1,1]);
NS = loadArray(fullfile(outdir,'NS.txt'),[1,1]);
NS_all = NS*NS;

a = loadArray(fullfile(outdir,'a.txt'),[NA,1]);
%a_grid = loadArray(fullfile(outdir,'a_grid.txt'),[NA+1,1]);

ages   = loadArray(fullfile(outdir,'ages.txt'),[JJ,1]);
c_coh  = loadArray(fullfile(outdir,'c_coh.txt'),[JJ,1]);
ym_coh = loadArray(fullfile(outdir,'ym_coh.txt'),[JJ,1]);
yf_coh = loadArray(fullfile(outdir,'yf_coh.txt'),[JJ,1]);
l_coh  = loadArray(fullfile(outdir,'l_coh.txt'),[JJ,1]);
cv_c  = loadArray(fullfile(outdir,'cv_c.txt'),[JJ,1]);
cv_y  = loadArray(fullfile(outdir,'cv_y.txt'),[JJ,1]);
tot_earn = loadArray(fullfile(outdir,'tot_earn.txt'),[JJ,1]);
nchild = loadArray(fullfile(outdir,'nchild.txt'),[JJ,1]);
pchild = loadArray(fullfile(outdir,'pchild.txt'),[JJ,1]);
% Survival probs
psi_vec = loadArray(fullfile(outdir,'psi.txt'),[JJ+1,1]);
% Age-determ efficiency profile
eff_vec = loadArray(fullfile(outdir,'eff.txt'),[JJ,1]);

% Distribution
stat_dist = loadArray(fullfile(outdir,'phi.txt'),[NA,NS_all,NP,JJ]);
% Policy for a'(a,theta,eta,age)
aplus = loadArray(fullfile(outdir,'aplus.txt'),[NA,NS_all,NP,JJ]);

tot_income = ym_coh+yf_coh;
load('old_res.mat','tot_income_save');
load('aplus_old.mat','aplus_old');

myerr = mean(abs(tot_income_save-tot_income));
disp(myerr)
aplus_old = permute(aplus_old,[1,3,2,4]);
myerr2 = mean(abs(aplus_old(2:end,:,:,:)-aplus),"all");
disp(myerr2)

mu_a = sum(stat_dist,[2,3,4]);

i_fig = 1; % Initialize figure counter

figure
i_fig = i_fig+1;
plot([ages;ages(end)+1],psi_vec,'LineWidth',2)
title('Survival probabilities')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(ages,eff_vec,'LineWidth',2)
title('Age-determ efficiency profile')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(ages,c_coh,'LineWidth',2)
hold on
plot(ages,tot_earn,'LineWidth',2)
legend('Consumption','Total earnings')
xlabel('Age j')
ylabel('Mean')
ylim([0,5])
title('Consumption and Earnings, life-cycle')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(ages,ym_coh,'LineWidth',2)
hold on
plot(ages,yf_coh,'LineWidth',2)
hold on
plot(ages,ym_coh+yf_coh,'LineWidth',2)
legend('Male','Female','Total')
ylim([0,5])
xlabel('Age j')
ylabel('Mean')
title('Earnings, life-cycle')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(ages,l_coh,'LineWidth',2)
ylim([0,1])
xlabel('Age j')
ylabel('Mean')
title('Female LF particip., life-cycle')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(ages,cv_c,'LineWidth',2)
hold on 
plot(ages,cv_y,'LineWidth',2)
legend('Consumption','Income')
%ylim([0,1])
xlabel('Age j')
ylabel('Mean')
title('Variance of the logs, life-cycle')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(ages,nchild,'LineWidth',2)
hold on
plot(ages,pchild,'LineWidth',2)
legend('Num. of children','Childcare costs')
xlabel('Age of mother, j')
title('Fertility and childcare costs')
print(fullfile(figdir,'nchild'),'-dpng')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(a,mu_a/JJ,'LineWidth',2)
title('Distribution of assets, pdf')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(a,cumsum(mu_a/JJ),'LineWidth',2)
title('Distribution of assets, cdf')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')

figure
i_fig = i_fig+1;
plot(a,a,'--','LineWidth',2)
hold on
plot(a,aplus(:,1,1,20),'LineWidth',2)
hold on
plot(a,aplus(:,2,1,20),'LineWidth',2)
legend('45 line','Low type','High type')
title('Policy for saving')
print(fullfile(figdir,['fig',num2str(i_fig)]),'-dpng')
