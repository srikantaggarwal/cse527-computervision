function [coeff] = getCubic(t1, t2, c)
    A = [t1^3, t1^2, t1, 1; 
         t2^3, t2^2, t2, 1;
         3*t1^2, 2*t1, 1, 0;
         3*t2^2, 2*t2, 1, 0;];
    b = [c;0;0;0];
    coeff = pinv(A)*b;
%    coeff = A\b;
end