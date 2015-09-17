function C = buildT2Y(gaama, optModelsY)
    n = length(optModelsY);
    C = zeros(n,1);
    for i=1:n
        t2 = optModelsY(i).t2;
        C(i) = -2*gaama*(t2.y);
    end
end