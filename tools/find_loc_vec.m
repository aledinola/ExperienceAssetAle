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
