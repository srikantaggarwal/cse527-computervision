%% Input image into logarithmic image
function outputImageText = RemoveShadowRGB(InputImage, MaskImage)
%InputImage = im2double(imread('inputImages/DSC_0394.jpg'));
%MaskImage = imread('gtMask/DSC_0394.png');
Radius = 8;
InputImage = im2double(InputImage);
meanval = zeros(1,3);
outGradX = InputImage;
outGradY = InputImage;
outGradXtext = InputImage;
outGradYtext = InputImage;
r = Radius;
for comp=1:3
    I = InputImage(:,:,comp);
    
    %% Getting boundary points of image
    BW = im2bw(MaskImage);
    dim = size(BW);
    col = round(dim(2)/2)-90;
    row = min(find(BW(:, col)));
    boundary = bwtraceboundary(BW, [row, col], 'N');
    % imshow(inputImage);
    % hold on;
    % plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 3);

    %% Calculating optimal illumination model
    % Taking log of whole image
    % [ny, nx] = size(I);
    R = log(I+1);
    meanval(1, comp) = mean(I(:));
    %Gx = gradX(:,:,comp);
    %Gy = gradY(:,:,comp);
    [Gx, Gy] = imgradientxy(R, 'CentralDifference');
    % figure, imshow(R);
    models = struct('t1', {}, 't2', {}, 'c', {});
    optModelsX = struct('t1', {}, 't2', {}, 'c', {});
    optModelsY = struct('t1', {}, 't2', {}, 'c', {});
    % Let Gx be the gradient of original image along x direction
    % Let Gy be the gradient of original image along y direction
    % [Gx, Gy] = imgradientxy(R, 'Sobel');
    
    numPoints = length(boundary);
    t0 = struct('x', {}, 'y', {});
    for i = 1:numPoints
        t0(i).x = boundary(i, 2);
        t0(i).y = boundary(i, 1);
    end

    %solve for horizontal lines
    % if comp==1
    for i = 1:numPoints
        gprod = 1;
        g_t1 = t0(i);
        g_t2 = t0(i);
        g_c = 0;
        for j = 1:r
            for k = 1:r
                t1 = t0(i);
                t2 = t0(i);
                t1.x = t1.x - j;
                t2.x = t2.x + k;

                if t1.x < 1
                    t1.x = 1;
                end
                if t2.x < 1
                    t2.x = 1;
                end
                if t1.x > size(R,2)
                    t1.x = size(R,2);
                end
                if t2.x > size(R,2)
                    t2.x = size(R,2);
                end
                if R(t2.y, t2.x) < R(t1.y, t1.x)
                    tmp = t2;
                    t2 = t1;
                    t1 = tmp;
                end
                c = R(t1.y, t1.x) - R(t2.y, t2.x);
                coeff = getCubic(t1.x, t2.x, c);
                points = zeros(2*r,1);
                for l = -r:r
                    t = t0(i);
                    t.x = t.x + l;
                    if t.x < 1
                        t.x = 1;
                    end
                    if t.x > size(R,2)
                        t.x = size(R,2);
                    end
                    tmp = t.x;
                    C_prime = 3*coeff(1)*tmp^2 + 2*coeff(2)*tmp + coeff(3);
                    G_hat = Gx(t.y, t.x) - C_prime;
                    points(l+r+1) = G_hat;
                end
                mu = mean(points);
                covar = var(points);
                product = -1;
                for l = -r:r
                    product = product * phi(points(l+r+1), mu, covar);
                end

                if product <= gprod
                    gprod = product;
                    g_t1 = t1;
                    g_t2 = t2;
                    g_c = c;
                end
            end
        end
        optModelsX(i).t1 = g_t1;
        optModelsX(i).t2 = g_t2;
        optModelsX(i).c = g_c;
    end

    %solve for horizontal lines
    for i = 1:numPoints
        gprod = 1;
        g_t1 = t0(i);
        g_t2 = t0(i);
        g_c = 0;
        for j = 1:r
            for k = 1:r
                t1 = t0(i);
                t2 = t0(i);
                t1.y = t1.y - j;
                t2.y = t2.y + k;

                if t1.y < 1
                    t1.y = 1;
                end
                if t2.y < 1
                    t2.y = 1;
                end
                if t1.y > size(R,1)
                    t1.y = size(R,1);
                end
                if t2.y > size(R,1)
                    t2.y = size(R,1);
                end
                if R(t2.y, t2.x) < R(t1.y, t1.x)
                    tmp = t2;
                    t2 = t1;
                    t1 = tmp;
                end
                c = R(t1.y, t1.x) - R(t2.y, t2.x);
                coeff = getCubic(t1.y, t2.y, c);
                points = zeros(2*r,1);
                for l = -r:r
                    t = t0(i);
                    t.y = t.y + l;
                    if t.y < 1
                        t.y = 1;
                    end
                    if t.y > size(R,1)
                        t.y = size(R,1);
                    end
                    tmp = t.y;
                    C_prime = 3*coeff(1)*tmp^2 + 2*coeff(2)*tmp + coeff(3);
                    G_hat = Gy(t.y, t.x) - C_prime;
                    points(l+r+1) = G_hat;
                end
                mu = mean(points);
                covar = var(points);
                product = -1;
                for l = -r:r
                    product = product * phi(points(l+r+1), mu, covar);
                end

                if product <= gprod
                    gprod = product;
                    g_t1 = t1;
                    g_t2 = t2;
                    g_c = c;
                end
            end
        end
        optModelsY(i).t1 = g_t1;
        optModelsY(i).t2 = g_t2;
        optModelsY(i).c = g_c;
    end
    % end
    %{
    figure, imshow(I);
    for i = 1:length(optModelsY)
        t1 = optModelsY(i).t1;
        t2 = optModelsY(i).t2;
        hold on;
        plot([t1.x, t2.x], [t1.y, t2.y], 'g', 'LineWidth', 3);
    end

    figure, imshow(I);
    for i = 1:length(optModelsX)
        t1 = optModelsX(i).t1;
        t2 = optModelsX(i).t2;
        hold on;
        plot([t1.x, t2.x], [t1.y, t2.y], 'g', 'LineWidth', 3);
    end
    %}
    lambda = 10;
    gaama = 0.9;

    H1 = buildH(lambda, gaama, numPoints);
    f1 = buildC(gaama, optModelsX);
    H2 = buildH(lambda, 1-gaama, numPoints);
    f2 = buildT1(1-gaama, optModelsX);
    H3 = buildH(lambda, 1-gaama, numPoints);
    f3 = buildT2(1-gaama,optModelsX);
    Z = zeros(numPoints, numPoints);
    H = [H1, Z, Z; Z, H2, Z; Z, Z, H3];
    f = [f1; f2; f3];

    A1 = eye(numPoints);
    b1 = zeros(numPoints,1);
    b2 = zeros(numPoints,1);
    b3 = zeros(numPoints,1);
    b4 = zeros(numPoints,1);
    b5 = zeros(numPoints,1);

    for i = 1:numPoints
        % b2(i) = t0(i).x;
        b2(i) = min(t0(i).x + r, size(I, 2));
        b3(i) = min(t0(i).x + r, size(I, 2));
        b4(i) = -1 * max(1, t0(i).x - r);
        % b5(i) = -1 * min(t0(i).x + 1, size(I, 2));
        b5(i) = -1 * max(1, t0(i).x - r);
    end

    b = [b1;b2;b3;b4;b5];
    pad = zeros(numPoints);
    A = [A1,pad,pad;
         pad,A1,pad;
         pad,pad,A1;
         pad,-A1,pad;
         pad,pad,-A1];

    vals = quadprog(H, f, A, b);
    cx = vals(1:numPoints);
    tx1 = vals((numPoints+1):(2*numPoints));
    tx2 = vals((2*numPoints)+1:end);
    tx1 = round(tx1);
    tx2 = round(tx2);
    modelsX = optModelsX;
    for i = 1:numPoints
        modelsX(i).c = cx(i);
        tmp = t0(i);
        tmp.x = tx1(i);
        modelsX(i).t1 = tmp;
        tmp.x = tx2(i);
        modelsX(i).t2 = tmp;
    end
    %{
    figure, imshow(I);
    for i = 1:length(modelsX)
        t1 = modelsX(i).t1;
        t2 = modelsX(i).t2;
        hold on;
        plot([t1.x, t2.x], [t1.y, t2.y], 'r', 'LineWidth', 3);
    end
    %}
    lambda = 10;
    gaama = 0.9;

    H1 = buildH(lambda, gaama, numPoints);
    f1 = buildC(gaama, optModelsY);
    H2 = buildH(lambda, 1-gaama, numPoints);
    f2 = buildT1Y(1-gaama, optModelsY);
    H3 = buildH(lambda, 1-gaama, numPoints);
    f3 = buildT2Y(1-gaama,optModelsY);
    Z = zeros(numPoints, numPoints);
    H = [H1, Z, Z; Z, H2, Z; Z, Z, H3];
    f = [f1; f2; f3];

    A1 = eye(numPoints);
    b1 = zeros(numPoints,1);
    b2 = zeros(numPoints,1);
    b3 = zeros(numPoints,1);
    b4 = zeros(numPoints,1);
    b5 = zeros(numPoints,1);

    for i = 1:numPoints
        b2(i) = min(t0(i).y + r, size(I, 1));
        % b2(i) = t0(i).y;
        b3(i) = min(t0(i).y + r, size(I, 1));
        b4(i) = -1 * max(1, t0(i).y - r);
        % b5(i) = -1 * min(t0(i).y + 1, size(I, 1));
        b5(i) = -1 * max(1, t0(i).y - r);
    end

    b = [b1;b2;b3;b4;b5];
    pad = zeros(numPoints);
    A = [A1,pad,pad;
         pad,A1,pad;
         pad,pad,A1;
         pad,-A1,pad;
         pad,pad,-A1];

    vals = quadprog(H, f, A, b);
    cy = vals(1:numPoints);
    ty1 = vals((numPoints+1):(2*numPoints));
    ty2 = vals((2*numPoints)+1:end);
    ty1 = round(ty1);
    ty2 = round(ty2);
    modelsY = optModelsY;
    for i = 1:numPoints
        modelsY(i).c = cy(i);
        tmp = t0(i);
        tmp.y = ty1(i);
        modelsY(i).t1 = tmp;
        tmp.y = ty2(i);
        modelsY(i).t2 = tmp;
    end
    
    %{
    figure, imshow(I);
    for i = 1:length(modelsY)
        t1 = modelsY(i).t1;
        t2 = modelsY(i).t2;
        hold on;
        plot([t1.x, t2.x], [t1.y, t2.y], 'r', 'LineWidth', 3);
    end
    %}
    %[img, OGx, OGy] = test2(modelsX, modelsY, Gx, Gy, R);
    % cancel the shadow effect in gradient field
    [OGx, OGy] = RemoveShadowEffect(modelsX, modelsY, Gx, Gy, R);
    extMask = extendMask(BW, boundary, modelsX, modelsY);

    gtcX = calcGradientST(OGx, boundary, extMask, r, 1);
    gtcY = calcGradientST(OGy, boundary, extMask, r, 0);
    %figure, imshow(Gx);
    %figure, imshow(Gy);
    %figure, imshow(OGx);
    %figure, imshow(OGy);
    %figure, imshow(gtcX);
    %figure, imshow(gtcY);
    
%     G = sqrt(Gx.*Gx + Gy.*Gy);
%     OG = sqrt(OGx.*OGx + OGy.*OGy);
%     Gtc = sqrt(gtcX.*gtcX + gtcY.*gtcY);
%     x = 1:size(I,2);
%     yp = modelsX(1).t1.y;
%     disp(yp);
%     figure, plot(x, I(yp, :), 'r');
%     figure, plot(x, G(yp, :), 'g');
%     figure,plot(x, OG(yp, :), 'b');
%     figure,plot(x, Gtc(yp, :), 'y');
    %plot(1:size(I,2), G(end-50, :));
    %plot(1:size(I,2), OG(end-50, :));
    %plot(1:size(I,2), Gtc(end-50, :));
    OGx(:,end) = 0;
    OGy(end,:) = 0;
    outGradX(:,:,comp) = OGx;
    outGradY(:,:,comp) = OGy;
    % BW = BW + addGrad;
    % make gradient texture consistent
    % should I process gradients separately or as single magnitude?
    % extMask=BW;
    
    
    
    %gtcX = calcGradientST(OGx, boundary, BW, r);
    %gtcY = calcGradientST(OGy, boundary, BW, r);
    gtcX(:,end) = 0;
    gtcY(end,:) = 0;
    outGradXtext(:,:,comp) = gtcX;
    outGradYtext(:,:,comp) = gtcY;
    
    %reconstruct
    % RecImg = poisson_solver_function(OGx, OGy, I);
    % TextRecImg = poisson_solver_function(gtcX, gtcY, I);
%     recMean = mean(R(:));
%     textRecMean = mean(R(:));
%     poissonOn = 0;
%     OGx(:,end) = 0;
%     OGx(end,:) = 0;
%     gtcX(:,end) = 0;
%     gtcY(end,:) = 0;
%     
%     RecImg = ImageRecH(OGx,OGy,recMean,poissonOn);
%     TextRecImg = ImageRecH(gtcX,gtcY,textRecMean,poissonOn);
%     
%     RecImg = exp(RecImg-1);
%     TextRecImg = exp(TextRecImg-1);
%     %initial image
%     figure, imshow(R, []);
%     %reconstructed
%     figure, imshow(RecImg,[]);
%     figure, imshow(TextRecImg, []);
%    outputImage = RecImg;
%    outputImageText = TextRecImg;
end

 poissonOn = 1;
 outputImage = ImageRecH(outGradX, outGradY, meanval, poissonOn);
 outputImageText = ImageRecH(outGradXtext, outGradYtext, meanval, poissonOn);

figure; imshow(InputImage, []);
%figure; imshow(outGradX, []);
%figure; imshow(outGradY, []);
%figure; imshow(outGradXtext, []);
%figure; imshow(outGradYtext, []);
figure; imshow(outputImage, []);
figure; imshow(outputImageText, []);
%imshow(img);
end