function C = buildT1(gaama, optModelsX)
    n = length(optModelsX);
    C = zeros(n,1);
    for i=1:n
        t1 = optModelsX(i).t1;
        C(i) = -2*gaama*(t1.x);
    end
end