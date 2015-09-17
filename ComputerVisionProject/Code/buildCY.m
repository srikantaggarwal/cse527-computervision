function C = buildCY(gaama, optModelsY)
    n = length(optModelsY);
    C = zeros(n,1);
    for i=1:n
        c0 = optModelsY(i).c;
        C(i) = -2*gaama*c0;
    end
end