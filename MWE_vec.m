clear,clc,close all
% check 13.0085
% STATE VARIABLES
% V(a,h,z,j)
% a: asset holdings
% h: human capital (female)
% z: labor productivity shocks (eta_m, eta_f are shocks, theta is permanent)
% j: age (from 1 to N_j)
% CHOICE VARIABLES
% d:  Female labor supply (only extensive margin: either 0 or 1)
% a': Next-period assets
% h' is implied by (d,a) and consumption is implied by the budget
% constraint
% DYNAMIC PROGRAMMING PROBLEM
% V(a,h,z,j) = max_{d,a'} F(d,a',a,h,z,j)+beta*s_j*E[V(a',h',z',j+1)|z]
% subject to
% h'=G(d,h), law of motion for human capital

verbose = 1;

%% Define grids and grid sizes
N_j = 80;
n_a = 201;
n_h = 11;
n_z = 50;
n_d = 2;
a_grid = linspace(0,450,n_a)';
h_grid = linspace(0,0.72,n_h)';
z_grid = linspace(0.9,1.1,n_z)';
d_grid = [0,1]';

z_grid = repmat(z_grid,[1,3]);

pi_z = rand(n_z,n_z);
pi_z = pi_z./sum(pi_z,2);

aprime_val = a_grid;  %(a',1)
a_val      = a_grid'; %(1,a)
h_grid3    = shiftdim(h_grid,-2);

%% Set parameters that do not depend on age
beta = 0.98;
r    = 0.04;
w_m  = 1;
w_f  = 0.75;
crra = 2;
nu   = 0.12;
Jr   = 45;
xi_1  = 0.05312;
xi_2  = -0.00188;
del_h = 0.074;
h_l   = 0;

p.eff_j     = ones(N_j,1);
p.pchild_j  = ones(N_j,1);
p.pen_j     = ones(N_j,1);
p.nchild_j  = ones(N_j,1);
p.s_j       = ones(N_j,1);
p.age_j     = (1:1:N_j)';

% Initialize output arrays
V  = zeros(n_a,n_h,n_z,N_j);
Policy = zeros(2,n_a,n_h,n_z,N_j);

tic

%% Solve problem in the last period

% Set age-dependent parameters
eff_j     = p.eff_j(N_j);
pchild_j  = p.pchild_j(N_j);
pen_j     = p.pen_j(N_j);
nchild_j  = p.nchild_j(N_j);
age_j     = p.age_j(N_j);
s_j       = p.s_j(N_j);

% V(a,h,z,N_j) = max_{d,a'}
V_d = zeros(n_a,n_h,n_z,n_d);
Pol_aprime_d = zeros(n_a,n_h,n_z,n_d);

for d_c = 1:n_d
    d_val = d_grid(d_c);
    for z_c = 1:n_z
        eta_m_val = z_grid(z_c,1);
        eta_f_val = z_grid(z_c,2);
        theta_val = z_grid(z_c,3);
        for h_c = 1:n_h
            h_val = h_grid(h_c);
            % RetMat is (a',a)
            RetMat = ReturnFn(d_val,aprime_val,a_val,h_val,eta_m_val,eta_f_val,theta_val,...
                w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,age_j,Jr);
            [max_val,max_ind] = max(RetMat,[],1);
            V_d(:,h_c,z_c,d_c) = max_val;
            Pol_aprime_d(:,h_c,z_c,d_c) = max_ind;
        end %end h
    end %end z
end %end d

[V(:,:,:,N_j),d_max] = max(V_d,[],4);

Policy(1,:,:,:,N_j) = d_max; % Optimal d
for z_c=1:n_z
    for h_c=1:n_h
        for a_c = 1:n_a
            d_star = d_max(a_c,h_c,z_c);
            Policy(2,a_c,h_c,z_c,N_j) = Pol_aprime_d(a_c,h_c,z_c,d_star); % Optimal a'
        end
    end
end

%% Backward iteration over age

for j = N_j-1:-1:1

    if verbose==1; fprintf('Age %d out of %d \n',j,N_j); end

    V_next = V(:,:,:,j+1); %V(a',h',z')

    % Set age-dependent parameters
    eff_j     = p.eff_j(j);
    pchild_j  = p.pchild_j(j);
    pen_j     = p.pen_j(j);
    nchild_j  = p.nchild_j(j);
    age_j     = p.age_j(j);
    s_j       = p.s_j(j);

    for z_c = 1:n_z

        eta_m_val = z_grid(z_c,1);
        eta_f_val = z_grid(z_c,2);
        theta_val = z_grid(z_c,3);

        % Compute EV(a',h'), given z
        % EV = zeros(n_a,n_h);
        z_prob = pi_z(z_c,:)';
        % for zprime_c = 1:n_z
        %     EV = EV+V_next(:,:,zprime_c)*z_prob(zprime_c);
        % end %end z'
        EV = V_next.*shiftdim(z_prob,-2); %V(a',h',z')*Prob(1,1,z')
        EV = sum(EV,3); %EV(a',h')

        for d_c=1:n_d
            d_val = d_grid(d_c);
            % Ret_mat is (a',a)
            Ret_mat = ReturnFn(d_val,aprime_val,a_val,h_grid3,eta_m_val,eta_f_val,theta_val,...
                w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,age_j,Jr);
            h_val      = h_grid;
            hprime_val = f_HC_accum(d_val,h_val,age_j,xi_1,xi_2,del_h,h_l); %(h,1)
            % ind_l and weight_l have size (n_h,1)
            [ind_l,weight_l] = find_loc_vec(h_grid,hprime_val); %(h,1)
            % EV_interp is (a',h), size (n_a,n_h)
            EV_interp = EV(:,ind_l).*weight_l'+EV(:,ind_l+1).*(1-weight_l');
            % entireRHS is (a',a,h) since (a',a,h) + (a',1,h)
            entireRHS = Ret_mat+beta*s_j*permute(EV_interp,[1,3,2]);
            [max_val,max_ind] = max(entireRHS,[],1);
            % max_val and max_ind are (1,a)
            V_d(:,:,z_c,d_c) = max_val;  %best V given d
            Pol_aprime_d(:,:,z_c,d_c) = max_ind; %best a' given d
        end %end d
    end %end z

    [V(:,:,:,j),d_max] = max(V_d,[],4);

    Policy(1,:,:,:,j) = d_max; % Optimal d
    for z_c=1:n_z
        for h_c=1:n_h
            for a_c = 1:n_a
                d_star = d_max(a_c,h_c,z_c);
                Policy(2,a_c,h_c,z_c,j) = Pol_aprime_d(a_c,h_c,z_c,d_star); % Optimal a'
            end
        end
    end

end %end j

toc

mean(Policy,'all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = ReturnFn(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,...
    pchild_j,pen_j,r,nchild_j,crra,nu,agej,Jr)

% Calculate earnings (incl. child care costs) of men and women
y_m = w_m*eff_j*theta*eta_m;
y_f = w_f*l_f*(exp(h_f)*theta*eta_f - pchild_j);
% l_f can be either 0 or 1
% calculate available resources
cash = (1+r)*a + pen_j*(agej>=Jr) + (y_m + y_f)*(agej<Jr);
cons = cash-aprime;

F = (cons/(sqrt(2+nchild_j))).^(1-crra)/(1-crra) - nu*l_f;
F(cons<=0) = -inf;

%F = -inf(size(cons)); %(a',a)
%pos = cons>0;         %(a',a)
%F(pos) = (cons(pos)/(sqrt(2+nchild_j))).^(1-crra)/(1-crra) - nu*l_f;%(a',a)

end %end function "f_ReturnFn"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h_f_prime] = f_HC_accum(l_f,h_f,age_j,xi_1,xi_2,del_h,h_l)
% l_f: d variable that affects h'
% h_f: current-period value of h

h_f_prime = h_f + (xi_1 + xi_2*age_j)*l_f - del_h*(1-l_f);
h_f_prime = max(h_f_prime, h_l);

%h_f_prime = h_f;
%h_f_prime = max(h_f_prime, h_l);


end %end function f_HC_accum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jl,omega] = find_loc_vec(x_grid,xi)
%-------------------------------------------------------------------------%
% DESCRIPTION
% Find jl s.t. x_grid(jl)<=xi<x_grid(jl+1)
% for jl=1,..,N-1
% omega is the weight on x_grid(jl) so that
% omega*x_grid(jl)+(1-omega)*x_grid(jl+1)=xi
% INPUTS
% x_grid must be a strictly increasing column vector (nx,1)
% xi can be a N-dim array with dim (s1,s2,..)
% OUTPUTS
% jl: Left point, same size as xi
% omega: weight on the left point, same size as xi
% NOTES
% Matlab recommneds to replace "histc" with "discretize".
%-------------------------------------------------------------------------%
nx = size(x_grid,1);

% For each 'xi', get the position of the 'x' element bounding it on the left [p x 1]
[~,jl] = histc(xi,x_grid); %#ok<HISTC>
jl(xi<=x_grid(1)) = 1;
jl(xi>=x_grid(nx)) = nx-1;

%Weight on x_grid(j)
omega = (x_grid(jl+1)-xi)./(x_grid(jl+1)-x_grid(jl));
omega = max(min(omega,1),0);

end %end function "find_loc_vec"
