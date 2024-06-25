function [h_f_prime] = f_HC_accum(l_f,h_f,age_j,xi_1,xi_2,del_h,h_l)
% l_f: d variable that affects h'
% h_f: current-period value of h

h_f_prime = h_f + (xi_1 + xi_2*age_j)*l_f - del_h*(1-l_f);
h_f_prime = max(h_f_prime, h_l);

%h_f_prime = h_f;
%h_f_prime = max(h_f_prime, h_l);


end %end function