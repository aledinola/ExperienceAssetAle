%% This script is called from main.m

tot_income = stat_age.ym_coh.Mean+stat_age.yf_coh.Mean+stat_age.pen_coh.Mean;

ra = find(ages==65); % First age in retirement 

%% Make figures

figure
plot(ages,Params.s_j,'LineWidth',2)
title('Survival probabilities')
print(fullfile(par.figdir,'surv_prob'),'-dpng')

figure
plot(ages,Params.eff_j,'LineWidth',2)
title('Age-determ efficiency profile')
print(fullfile(par.figdir,'age_profile'),'-dpng')

figure
plot(ages,Params.nchild_j,'LineWidth',2)
hold on
plot(ages,Params.pchild_j,'LineWidth',2)
legend('Num. of children','Childcare costs')
xlabel('Age of mother, j')
title('Fertility and childcare costs')
print(fullfile(par.figdir,'nchild'),'-dpng')
print(fullfile(par.figdir,'children'),'-dpng')


%% Replication Figure 10.3 

figure
subplot(2,2,1)
plot(ages,tot_income,'LineWidth',2)
hold on
plot(ages,stat_age.c_coh.Mean,'LineWidth',2)
hold on
plot(ages,stat_age.c_eq_coh.Mean,'--','LineWidth',2)
legend('Income','Consumption','Cons. (equivalence)','Location','southeast')
xlabel('Age j')
ylabel('Mean')
ylim([0,5])

subplot(2,2,2)
plot(ages(1:ra-1),stat_age.yf_coh.Mean(1:ra-1),'LineWidth',2)
hold on
plot(ages(1:ra-1),stat_age.ym_coh.Mean(1:ra-1),'LineWidth',2)
legend('Female','Male','Location','southeast')
ylim([0,3])
xlabel('Age j')
ylabel('Labor earnings (mean)')

subplot(2,2,3)
plot(ages(1:ra-1),stat_age.l_coh.Mean(1:ra-1),'LineWidth',2)
ylim([0,1])
xlabel('Age j')
ylabel('Female labor force particip.')

subplot(2,2,4)
plot(ages(1:ra-1),stat_age.h_coh.Mean(1:ra-1),'LineWidth',2)
hold on
plot(ages(1:ra-1),Params.eff_j(1:ra-1),'LineWidth',2)
ylim([0,2.5])
xlabel('Age j')
ylabel('Human capital')
legend('Female','Male','Location','southeast')
print(fullfile(par.figdir,'fig_10_3'),'-dpng')

%% Replication Figure 10.4

figure
plot(ages,stat_age.a_coh.Mean,'LineWidth',2)
xlabel('Age j')
ylabel('Wealth (mean)')
title('Savings in the presence of children')
print(fullfile(par.figdir,'fig_10_4'),'-dpng')
