function C = buildT2(gaama, optModelsX)
    n = length(optModelsX);
    C = zeros(n,1);
    for i=1:n
        t2 = optModelsX(i).t2;
        C(i) = -2*gaama*(t2.x);
    end
end