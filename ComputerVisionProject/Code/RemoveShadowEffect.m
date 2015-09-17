function [outGradX, outGradY] = RemoveShadowEffect(ModelsX, ModelsY, gradX, gradY, I)
%solve for horizontal lines
outGradX = gradX;
outGradY = gradY;
[ny, nx] = size(I);
numPoints = length(ModelsX);
tmpX = zeros(ny, nx);
tmpY = zeros(ny, nx);
for i = 1:numPoints
    t1 = ModelsX(i).t1;
    t2 = ModelsX(i).t2;
    c = ModelsX(i).c;
    
    if I(t1.y, t1.x) > I(t2.y, t1.x)
        tmp = t2;
        t2 = t1;
        t1 = tmp;
    end
    
    coeff = getCubic(t1.x, t2.x, c);
    for j = t1.x:t2.x
        C_prime = 3*coeff(1)*j^2 + 2*coeff(2)*j + coeff(3);
        outGradX(t1.y, j) = -abs(abs(gradX(t1.y, j)) - abs(C_prime));
        tmpX(t1.y, j) = 1;
    end
    for j = t2.x:t1.x
        C_prime = 3*coeff(1)*j^2 + 2*coeff(2)*j + coeff(3);
        outGradX(t1.y, j) = abs(abs(gradX(t1.y, j)) - abs(C_prime));
        tmpX(t1.y, j) = 1;
    end

end

numPoints = length(ModelsY);
for i = 1:numPoints
    t1 = ModelsY(i).t1;
    t2 = ModelsY(i).t2;
    c = ModelsY(i).c;
    
    if I(t1.y, t1.x) > I(t2.y, t2.x)
        tmp = t2;
        t2 = t1;
        t1 = tmp;
    end
    
    coeff = getCubic(t1.y, t2.y, c);
    for j = t1.y:t2.y
        C_prime = 3*coeff(1)*j^2 + 2*coeff(2)*j + coeff(3);
        outGradY(j, t1.x) = -abs(abs(gradY(j, t1.x)) - abs(C_prime));
        tmpY(j, t1.x) = 1;
    end
    for j = t2.y:t1.y
        C_prime = 3*coeff(1)*j^2 + 2*coeff(2)*j + coeff(3);
        outGradY(j, t1.x) = abs(abs(gradY(j, t1.x)) - abs(C_prime));
        tmpY(j, t1.x) = 1;
    end
end

% figure, imshow(abs(outGradX),[]);
% figure, imshow(abs(gradX),[]);
%figure, imshow(tmpX,[]);
%figure, imshow(tmpY,[]);
%figure, imshow(I);
%NSImg = poisson_solver_function(outGradX, outGradY, I);
%NSImg = exp(NSImg-1);
%figure, imshow(NSImg,[]);
end