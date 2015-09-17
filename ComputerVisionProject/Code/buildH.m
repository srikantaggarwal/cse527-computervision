function H = buildH(lambda, gaama, n)
    H = zeros(n,n);
    for i=1:n
        if i > 1
            H(i-1,i) = -2*gaama*lambda;
            H(i,i-1) = -2*gaama*lambda;
        end
        H(i,i) = gaama + 4*gaama*lambda;
    end
    H = 2*H;
end