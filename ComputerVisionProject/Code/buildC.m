function C = buildC(gaama, optModels)
    n = length(optModels);
    C = zeros(n,1);
    for i=1:n
        c0 = optModels(i).c;
        C(i) = -2*gaama*c0;
    end
end