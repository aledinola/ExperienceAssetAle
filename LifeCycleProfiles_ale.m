function [stat_age] = LifeCycleProfiles_ale(pi_z,Mu,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions)
% There are 4 STATE VARIABLES
% Distribution(a,h,z,j)
%   Policy(1:2,a,h,z,j)
% where the exogenous state z = (eta_m,eta_f,theta)
% j denotes age and is always the last state variable.

verbose = simoptions.verbose;

% Convert (a,h,eta_m,eta_f,theta,j) into (a,h,z,j)
Mu = gather(Mu);
Policy = gather(Policy);
Mu = reshape(Mu,[n_a(1),n_a(2),prod(n_z),N_j]);
Policy = reshape(Policy,[2,n_a(1),n_a(2),prod(n_z),N_j]);

[n_a,n_h,a_grid,h_grid,z_grid,n_z] = prepare_grids(n_a,a_grid,z_grid,n_z,pi_z);

Pol_d_ind      = squeeze(Policy(1,:,:,:,:)); % indexes (a,h,z,j)
Pol_aprime_ind = squeeze(Policy(2,:,:,:,:)); % indexes (a,h,z,j)

Pol_d  = d_grid(Pol_d_ind);      % values (a,h,z,j)
Pol_ap = a_grid(Pol_aprime_ind); % values (a,h,z,j)

if ~isequal(size(Mu),[n_a,n_h,n_z,N_j])
    error('Size of stationary distribution is not correct')
end
if ~isequal(size(Pol_d),[n_a,n_h,n_z,N_j])
    error('Size of Pol_d is not correct')
end
if ~isequal(size(Pol_ap),[n_a,n_h,n_z,N_j])
    error('Size of Pol_ap is not correct')
end
if any(Pol_d<0,"all") || any(Pol_d>1,"all")
    error('Pol_d is not in [0,1]')
end

%% Grid dimensions

% Set parameters that do not depend on age
r    = Params.r;
w_m  = Params.w_m;
w_f  = Params.w_f;
%crra = Params.crra;
%nu   = Params.nu;
Jr   = Params.Jr;
%xi_1  = Params.xi_1;
%xi_2  = Params.xi_2;
%del_h = Params.del_h;
%h_l   = Params.h_l;

%% Compute conditional averages over the life-cycle

% l_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta)

MuVec = reshape(Mu,[n_a*n_h*n_z,N_j]); 

l_coh = zeros(N_j,1);
a_coh = zeros(N_j,1);
h_coh = zeros(N_j,1);

y_coh = zeros(N_j,1);
ym_coh = zeros(N_j,1);
yf_coh = zeros(N_j,1);
pen_coh = zeros(N_j,1);
c_coh = zeros(N_j,1);
c_eq_coh = zeros(N_j,1);

for j = 1:N_j

    MuVec_j = MuVec(:,j);
    MuVec_j = MuVec_j/sum(MuVec_j);

    if verbose==1
        fprintf('Age = %d out of %d \n',j,N_j)
    end
    % Set age-dependent parameters
    eff_j     = Params.eff_j(j);
    pchild_j  = Params.pchild_j(j);
    pen_j     = Params.pen_j(j);
    nchild_j  = Params.nchild_j(j);
    age_j     = Params.agej(j);
    %s_j       = Params.s_j(j);

%     FnsToEvaluate.l_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) l_f; 
% FnsToEvaluate.a_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) a;
% FnsToEvaluate.h_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) exp(h_f);
% FnsToEvaluate.ym_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j) w_m*eff_j*theta*eta_m;
% FnsToEvaluate.yf_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_f) w_f*exp(h_f)*theta*eta_f*l_f;
% FnsToEvaluate.pen_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,pen_j) pen_j;
% FnsToEvaluate.y_coh = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f) w_m*eff_j*theta*eta_m+w_f*exp(h_f)*theta*eta_f*l_f;
% FnsToEvaluate.c_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr) f_ReturnFn_cons(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr);
% FnsToEvaluate.c_eq_coh = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr) f_ReturnFn_cons_equiv(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr);
% 
    Values_l_coh = zeros(n_a,n_h,n_z);
    Values_a_coh = zeros(n_a,n_h,n_z);
    Values_h_coh = zeros(n_a,n_h,n_z);

    Values_y_coh = zeros(n_a,n_h,n_z);
    Values_ym_coh = zeros(n_a,n_h,n_z);
    Values_yf_coh = zeros(n_a,n_h,n_z);
    Values_pen_coh = zeros(n_a,n_h,n_z);
    Values_c_coh = zeros(n_a,n_h,n_z);
    Values_c_eq_coh = zeros(n_a,n_h,n_z);

    for z_c=1:n_z
        eta_m = z_grid(z_c,1);
        eta_f = z_grid(z_c,2);
        theta = z_grid(z_c,3);
        for h_c=1:n_h
            h_f = h_grid(h_c);
            %for a_c=1:n_a
                a = a_grid;
                l_f = Pol_d(:,h_c,z_c,j);
                aprime = Pol_ap(:,h_c,z_c,j);
                Values_l_coh(:,h_c,z_c)=FnsToEvaluate.l_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta);
                Values_a_coh(:,h_c,z_c)=FnsToEvaluate.a_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta);
                Values_h_coh(:,h_c,z_c)=FnsToEvaluate.h_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta);
                Values_ym_coh(:,h_c,z_c)=FnsToEvaluate.ym_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j);
                Values_yf_coh(:,h_c,z_c)=FnsToEvaluate.yf_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_f);

                Values_y_coh(:,h_c,z_c)=FnsToEvaluate.y_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f);
                Values_pen_coh(:,h_c,z_c)=FnsToEvaluate.pen_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,pen_j);
                Values_c_coh(:,h_c,z_c)=FnsToEvaluate.c_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,age_j,Jr);
                Values_c_eq_coh(:,h_c,z_c)=FnsToEvaluate.c_eq_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,age_j,Jr);
            %end
        end
    end

    l_coh(j)    = compute_mean(Values_l_coh,MuVec_j,n_a,n_h,n_z);
    a_coh(j)    = compute_mean(Values_a_coh,MuVec_j,n_a,n_h,n_z);
    h_coh(j)    = compute_mean(Values_h_coh,MuVec_j,n_a,n_h,n_z);
    y_coh(j)    = compute_mean(Values_y_coh,MuVec_j,n_a,n_h,n_z);
    ym_coh(j)   = compute_mean(Values_ym_coh,MuVec_j,n_a,n_h,n_z);
    yf_coh(j)   = compute_mean(Values_yf_coh,MuVec_j,n_a,n_h,n_z);
    pen_coh(j)  = compute_mean(Values_pen_coh,MuVec_j,n_a,n_h,n_z);
    c_coh(j)    = compute_mean(Values_c_coh,MuVec_j,n_a,n_h,n_z);
    c_eq_coh(j) = compute_mean(Values_c_eq_coh,MuVec_j,n_a,n_h,n_z);

end %end j

stat_age.l_coh   = l_coh;
stat_age.a_coh   = a_coh;
stat_age.h_coh   = h_coh;
stat_age.ym_coh  = ym_coh;
stat_age.yf_coh  = yf_coh;
stat_age.y_coh   = y_coh;
stat_age.pen_coh = pen_coh;
stat_age.c_coh    = c_coh;
stat_age.c_eq_coh = c_eq_coh;


end %end function

function Mean_j = compute_mean(Values,MuVec_j,n_a,n_h,n_z)

ValuesVec = reshape(Values,[n_a*n_h*n_z,1]);

Mean_j = sum(ValuesVec.*MuVec_j);

end