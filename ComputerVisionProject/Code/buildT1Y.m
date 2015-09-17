function C = buildT1Y(gaama, optModelsY)
    n = length(optModelsY);
    C = zeros(n,1);
    for i=1:n
        t1 = optModelsY(i).t1;
        C(i) = -2*gaama*(t1.y);
    end
end