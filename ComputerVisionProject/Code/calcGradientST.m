%% Transform gradient in shadow area to be texture-consistent
% G: gradient of shadow-free image
% bp: boundary points of shadow
% mask: mask image which contains 1 where shadow is, 0 otherwise
% r: radius
function Result = calcGradientST(G, bp, mask, r, horiz, modelsX, modelsY)
    n = size(bp, 1); % number of boundary points
    t0y=bp(1:n,1); 
    t0x=bp(1:n,2);
    szY = size(G, 1);
    szX = size(G, 2);
    Result = G;
    %% get the gradient only of umbra+penumbra, zeros outside
    Gsh = G.*mask;    
    
    muS = mean(Gsh(:));% sampling mean    
    stdS = std(Gsh(:), 1); % sampling std deviation
    %% find grad.values in umbra and lit area (along horizontal and vert. lines)
    umbArr=zeros(n,1);
    litArr=zeros(n,1);
    extentSz = 3;
    for i=1:n % for each boundary point       
        if horiz==1        
            %along horizontal line  
            xua=[];%array for umbra area extent
            xla=[];%array for lit area extent
            offsetx1 =max(t0x(i)-r-extentSz, 1);
            % suppose X1 is for umbra, check later
            for counter=1:extentSz
                extX1= offsetx1 + counter;
                xua(counter)=extX1;            
            end
            % suppose extX2 is for lit area
            offsetx2 = min(t0x(i)+r, szX);
            for counter=1:extentSz
                extX2= min(offsetx2 + counter-1, szX);
                xla(counter)=extX2;
            end
            
            % check last points in arr, umbra or lit area (the same is true for
            % other points in arr).
            if (mask(t0y(i), offsetx1+1) > mask(t0y(i), offsetx2))%x1 is lit area
                tmpa=xua;
                xua = xla;
                xla = tmpa;                
            elseif(mask(t0y(i), offsetx1+1) > mask(t0y(i), offsetx2))
                disp('equal!!!!')
            end     
        else
            yua=[];%size of extentSz
            yla=[];
            offsety1 = max(t0y(i)-r-extentSz, 1);
            % suppose Y1 is for umbra, check later
            for counter=1:extentSz
                extY1= offsety1 + counter;
                yua(counter)=extY1;            
            end
            % suppose extX2 is for lit area
            offsety2 = min(t0y(i)+r, szY);
            for counter=1:extentSz
                extY2=min(offsety2 + counter-1, szY);
                yla(counter)=extY2;
            end
            % check if need to swap lit area and umbra array
            if (mask(offsety1+1, t0x(i)) > mask(offsety2, t0x(i))) % y2 is umbra
                tmpa=yua;
                yua=yla;
                yla=tmpa;
            end
        end        
        
        % sample the gradient value in the umbra and lit area
        sumU=0;
        sumL = 0;
        if horiz==1
            for counter=1:extentSz
                sumU=sumU+G(t0y(i),xua(counter));
                sumL=sumL+G(t0y(i),xla(counter));
            end
        else
            for counter=1:extentSz
                sumU=sumU+G(yua(counter),t0x(i));
                sumL=sumL+G(yla(counter),t0x(i));
            end
        end
        
        umbArr(i)=sumU/extentSz;
        litArr(i)=sumL/extentSz;        
    end    
    
    % calc statistics
    muUmbra = mean(umbArr);
    varUmbra = var(umbArr);
    muLa = mean(litArr);
    varLa = var(litArr);    
    
    %% results
    
    muSe = muUmbra - muLa;    
    varSe = varUmbra - varLa;
    
    disp('shadow effect params:');
    [muUmbra, muLa, varUmbra, varLa]
    
    varS = stdS^2;    
    if (varS > varSe)
        muT = muS - muSe;
        stdT = sqrt(varS - varSe);        
    else
        muT = muSe - muS;
        stdT = sqrt(varSe - varS);
        disp('varS < varSe');        
        [varS, varSe, muS, muSe, horiz]
    end
    
    for i = 1:size(Gsh, 1)
        for j = 1:size(Gsh, 2)
            if (Gsh(i, j) == 0)
                continue;
            end;
            
            Result(i,j) = muT + (Gsh(i,j) - muS)*stdT/stdS;
        end
    end    
end