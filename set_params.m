function [par] = set_params()

% Set flag
% Important note: for small grid, say NA<=250, better loop version
% For larger grid, say NA=600, vectorized version (on GPU) is faster
par.do_vfi_vec = 2; % 0 = on CPU, vectorized only (a',a)
                    % 1 = on GPU, vectorized over (a',a,h), low memory use
                    % 2 = on GPU, vectorized over (a',a',h,z), high memeory
par.do_refine  = 0; % SET IT TO 0, not implemented yet
par.n_ref      = 20;
discretize_method = 'Rouwenhorst'; %'Rouwenhorst','Tauchen'

par.inpdir = 'inputs'; 
par.outdir = 'output';
par.figdir = 'figures';

par.JJ = 80;            % Maximum age
par.JR = 45;            % Retirement age
par.agejshifter = 20;   % age 1 in model = age 20 in the data

%% Grid dimensions
par.NA = 601;           % No. points asset grid (201 or 601)
par.NH = 11;            % No. points human capital grid
par.NP = 2;             % No. of permanent household types
par.NS = 5;             % No. points shocks for each gender (NS*NS in total)
%NS_all = NS*NS;

%% Economic parameters
% -- Household preference parameters
par.crra = 2.0;
par.nu = 0.12;
par.beta = 0.98;
% -- Household risk process
par.sigma_theta = 0.242; %variance fixed effect
par.sigma_eps   = 0.022; %variance AR1 shock
par.rho         = 0.985;
% -- Human capital parameters
par.xi    = [0.05312, -0.00188];
par.del_h = 0.074;
% -- Prices 
par.r   = 0.04;
par.w_m = 1.00;
par.w_f = 0.75;
% -- size of the asset grid
par.a_l    = 0.0;
par.a_u    = 450;
par.a_grow = 0.05;
a_space    = 3.0;
% -- Fixed effects, distribution and grid
par.theta_dist = [0.5,0.5]';
theta_grid_log = [-sqrt(par.sigma_theta),sqrt(par.sigma_theta)]';
par.theta_grid = exp(theta_grid_log);

% -- Use Rouwenhorst to discretize the stochastic process for (eta_m,eta_f),
%    under the assumption that eta_m,eta_f are INDEPENDENT
opt.parallel = 0;
switch discretize_method
    case 'Rouwenhorst'
        [eta_m_log,par.pi_m]=discretizeAR1_Rouwenhorst(0.0,par.rho,sqrt(par.sigma_eps),par.NS,opt);
        [eta_f_log,par.pi_f]=discretizeAR1_Rouwenhorst(0.0,par.rho,sqrt(par.sigma_eps),par.NS,opt);
    case 'Tauchen'
        %[z_grid,P]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q, tauchenoptions)
        Tauchen_q = 2.0;
        [eta_m_log,par.pi_m]=discretizeAR1_Tauchen(0.0,par.rho,sqrt(par.sigma_eps),par.NS,Tauchen_q,opt);
        [eta_f_log,par.pi_f]=discretizeAR1_Tauchen(0.0,par.rho,sqrt(par.sigma_eps),par.NS,Tauchen_q,opt);
end
par.eta_m = exp(eta_m_log);
par.eta_f = exp(eta_f_log);

% -- Grid for assets. 
par.a = make_grid(par.a_l,par.a_u,par.NA,a_space,1);

% -- Import age-related stuff
par.pen     = loadArray(fullfile(par.inpdir,'pen.txt'),[par.JJ,1]);
% Labor efficiency units
par.eff     = loadArray(fullfile(par.inpdir,'eff.txt'),[par.JJ,1]);
% Survival probabilities
par.psi     = loadArray(fullfile(par.inpdir,'psi.txt'),[par.JJ+1,1]);
% No. of children per household
par.nchild  = loadArray(fullfile(par.inpdir,'nchild.txt'),[par.JJ,1]);
%par.nchild  = zeros(par.JJ,1);
% Childcare costs
par.pchild  = loadArray(fullfile(par.inpdir,'pchild.txt'),[par.JJ,1]);
%par.pchild  = zeros(par.JJ,1);

% -- Grid for human capital (in logs)
par.h_l = 0;
% Compute max human capital by age (assume female always works)
h_max_age = zeros(par.JJ,1);
h_max_age(1) = 0;
for j = 2:par.JR-1
    h_max_age(j) = h_max_age(j-1) + par.xi(1) + par.xi(2)*(j-1);
end
h_max_age(par.JR:par.JJ) = h_max_age(par.JR-1);
h_u = max(h_max_age);
h_space = 2.0;
par.h_f_grid = make_grid(par.h_l,h_u,par.NH,h_space,1);
par.h_max_age = h_max_age;

% Initial condition for (eta_m,eta_f) in the distribution
par.is_initial_m = round((par.NS+1)/2);
par.is_initial_f = round((par.NS+1)/2);

end %end function "set_params"