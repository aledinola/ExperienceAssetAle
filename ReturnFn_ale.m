function F = ReturnFn_ale(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,...
    pchild_j,pen_j,r,nchild_j,crra,nu,agej,Jr)

% Calculate earnings (incl. child care costs) of men and women
y_m = w_m*eff_j*theta*eta_m;
y_f = w_f*l_f*(exp(h_f)*theta*eta_f - pchild_j);
% l_f can be either 0 or 1
% calculate available resources
cash = (1+r)*a + pen_j*(agej>=Jr) + (y_m + y_f)*(agej<Jr);
cons = cash-aprime;

%pos    = cons>0;
F = (cons/(sqrt(2+nchild_j))).^(1-crra)/(1-crra) - nu*l_f;
F(cons<=0) = -inf;

end %end function "f_ReturnFn"