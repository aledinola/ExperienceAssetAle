%% Wage risk and spousal insurance
% Author: Alessandro Di Nola
% NOTE: It seems that VFI toolkit has a small mistake in the Distribution
% code which leads to some differences in the life-cycle profiles.
% do_toolkit = 0 is the recommended value.

clear
clc
close all
format long g
addpath('tools')
% toolkit on laptop
%addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab'))
% toolkit on desktop PC
addpath(genpath('C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab'))

%% Set parameters

[par] = set_params();

do_toolkit = 0;

%% Declare grids and probabilities for toolkit
n_d    = 2;
d_grid = [0,1]'; % Female labor supply
a_grid = [par.a;par.h_f_grid]; % Grids for assets and female human capital
n_a    = [par.NA,par.NH];  % Number of grid points for asset, and female human capital
n_z    = [par.NS,par.NS,par.NP]; % we have eta_m, eta_f, theta
z_grid = [par.eta_m;par.eta_f;par.theta_grid];
pi_theta = eye(par.NP); % transition matrix for perm type is I(2)
% Kronecker in reverse order theta*eta_f*eta_m
pi_z = kron(kron(pi_theta,par.pi_f),par.pi_m);
if size(pi_z,1)~=prod(n_z)
    error('Incorrect number of elements in pi_z')
end
% Set initial conditions for shocks (eta_m,eta_f,theta)
eta_m_init = par.is_initial_m;
eta_f_init = par.is_initial_f;
theta_init = (1:par.NP)';
% We consider all fixed effects theta

% Combine all initial conditions for (eta_m,eta_f,theta) in "z":
z_init = sub2ind([par.NS,par.NS,par.NP],eta_m_init*ones(par.NP,1),eta_f_init*ones(par.NP,1),theta_init);

%% Set parameters for the toolkit
% Set demographic parameters
Params.J    = par.JJ; 
Params.Jr   = par.JR;   
N_j         = Params.J;     
Params.agej = 1:1:Params.J;
Params.agejshifter=par.agejshifter; % Makes keeping track of actual age easy in terms of model age

% Set age-dependent parameters. The toolkit recognizes automatically that
% these are age-dependent
Params.eff_j    = par.eff;
Params.pen_j    = par.pen;
Params.s_j      = par.psi(1:Params.J); % Conditional survival probabilities 
Params.pchild_j = par.pchild; % Cost of child-care
Params.nchild_j = par.nchild; % Number of children

% Set all other parameters
Params.beta = par.beta;
Params.w_m  = par.w_m;
Params.w_f  = par.w_f;
Params.r    = par.r;
Params.crra = par.crra;
Params.nu   = par.nu;
Params.xi_1 = par.xi(1);
Params.xi_2 = par.xi(2);
Params.del_h = par.del_h;
Params.h_l   = par.h_l;

%% Now, create the return function 

% TOOLKIT: (d,a',a,z) where d is the vector of decision variables, z is
% the vector of exogenous shocks and a is the vector of endogenous states.
% IN THIS MODEL: (l_f,a',a,h_f,eta_m,eta_f,theta)
DiscountFactorParamNames={'beta','s_j'};

ReturnFn = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,agej,Jr) ...
    f_ReturnFn(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,agej,Jr);

%% Create the function for human capital accumulation
% aprimeFn gives the value of h_fprime
vfoptions.aprimeFn=@(l_f,h_f,agej,xi_1,xi_2,del_h,h_l) f_HC_accum(l_f,h_f,agej,xi_1,xi_2,del_h,h_l);

%% Now solve the value function iteration problem

vfoptions.experienceasset=1; % Using an experience asset
% Note: by default, assumes it is the last d variable that controls the
% evolution of the experience asset (and that the last a variable is
% the experience asset).
vfoptions.verbose          = 1;
vfoptions.verboseparams    = 1;
vfoptions.lowmemory        = 1;

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_z]); % Put no households anywhere on grid
% All agents start with zero assets, one unit of human capital and the 
% median value of each AR1 shock and the mass for permanent type
for ip=1:par.NP
    jequaloneDist(1,1,eta_m_init,eta_f_init,ip)=par.theta_dist(ip); 
end

%% Options for the 'stationary distribution' of households

% Start with a mass of one at initial age, use the conditional survival
% probabilities sj to calculate the mass of those who survive to next
% period, repeat. Once done for all ages, normalize to one
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:Params.J
    Params.mewj(jj)=Params.s_j(jj-1)*Params.mewj(jj-1);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age

%% Set options for distribution and life-cycle profiles
simoptions=struct(); % Use the default options
simoptions.verbose = 1;
% We also need to tell simoptions about the experience asset
simoptions.experienceasset=1;
simoptions.d_grid   = d_grid; % because we are using an experience asset
simoptions.a_grid   = a_grid; % because we are using an experience asset
simoptions.aprimeFn = vfoptions.aprimeFn;

%http://discourse.vfitoolkit.com/t/how-to-cut-runtimes-for-life-cycle-profiles-or-allstats-when-you-dont-need-all-the-stats-simoptions-whichstats/253
simoptions.whichstats=zeros(1,7);
simoptions.whichstats(1) = 1;
simoptions.whichstats(3) = 1;

%% FnsToEvaluate are how we say what we want to graph the life-cycles of
% Like with return function, we have to include (d,aprime,a,z) as first
% inputs, then just any relevant parameters.
% (h,aprime,a,z)==(l_f,aprime,a,eta_m,eta_f,theta)

% l_coh is fraction of time worked, or the share of women who work
FnsToEvaluate.l_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) l_f; 
FnsToEvaluate.a_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) a;
FnsToEvaluate.h_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) exp(h_f);
FnsToEvaluate.ym_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j) w_m*eff_j*theta*eta_m;
FnsToEvaluate.yf_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_f) w_f*exp(h_f)*theta*eta_f*l_f;
FnsToEvaluate.pen_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,pen_j) pen_j;
FnsToEvaluate.y_coh = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f) w_m*eff_j*theta*eta_m+w_f*exp(h_f)*theta*eta_f*l_f;
FnsToEvaluate.c_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr) f_ReturnFn_cons(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr);
FnsToEvaluate.c_eq_coh = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr) f_ReturnFn_cons_equiv(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr);

%% Call function to do VFI, distribution and profiles

if do_toolkit==1
[V,Policy,StationaryDist,stat_age,vf_time,time_distrib,age_time] = ...
    solve_model(par,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,...
    DiscountFactorParamNames,vfoptions,jequaloneDist,AgeWeightsParamNames,...
    simoptions,FnsToEvaluate);
else
[V,Policy,StationaryDist,stat_age,vf_time,time_distrib,age_time] = ...
    solve_model_ale(par,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,...
    DiscountFactorParamNames,vfoptions,jequaloneDist,AgeWeightsParamNames,...
    simoptions,FnsToEvaluate);

end

% figure
% plot(1:Params.Jr-1,stat_age.l_coh(1:Params.Jr-1))

%% Calculate additional moments
%mu_a = sum(StationaryDist,[2,3,4,5,6]);

ages = Params.agejshifter+(1:Params.J);

fprintf('Time (in secs) for VFI: %f \n',vf_time)
fprintf('Time (in secs) for Distrib.: %f \n',time_distrib)
fprintf('Time (in secs) for Age profiles.: %f \n',age_time)

%% Call script to plot life-cycle moments

if do_toolkit==1
    make_plots_toolkit
else
    make_plots 
end
