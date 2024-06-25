function cons_equiv = f_ReturnFn_cons_equiv(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,...
    pchild_j,pen_j,r,nchild_j,agej,Jr)

% Calculate earnings (incl. child care costs) of men and women
y_m = w_m*eff_j*theta*eta_m;
y_f = w_f*l_f*(exp(h_f)*theta*eta_f - pchild_j);
% l_f can be either 0 or 1
% calculate available resources
cash = (1+r)*a + pen_j*(agej>=Jr) + (y_m + y_f)*(agej<Jr);
cons = cash-aprime;
cons_equiv = cons/(sqrt(2+nchild_j));

end %end function "f_ReturnFn"