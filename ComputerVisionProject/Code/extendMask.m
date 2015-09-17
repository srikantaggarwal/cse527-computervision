function em = extendMask(mask,bp, modelsX, modelsY)
    n = size(bp,1);
    em=mask;
    szY=size(mask,1);
    szX=size(mask,2);
    r=5;
    arr=size(n,1);
    for i=1:n 
        bx=bp(i,2);
        by=bp(i,1);
        t2=modelsX(i).t2;
        t2y=modelsY(i).t2;
        arr(i)=bx-t2.x;
        % search when mask appears
        for ind=-r:r            
            x=min(szX-1,max(1, bx+ind));% getting x safely
            if (mask(by,x) == 0 && mask(by,x+1)==1)  
                for xx=t2.x:x
                    em(max(1,by-1), xx)=1;
                    em(by, xx)=1;
                    em(min(szY,by+1), xx)=1;
                end                
            elseif (mask(by, x)==1 && mask(by, x+1) == 0)
                for xx=x+1:t2.x
                    em(max(1,by-1), xx)=1;
                    em(by, xx)=1;
                    em(min(szY,by+1), xx)=1;
                end
            end
        end
        
        for ind=-r:r            
            y=min(szY-1,max(1, by+ind));% getting y safely
            if (mask(y,bx) == 0 && mask(y+1,bx)==1)                  
                for yy=t2y.y:y
                    em(yy, max(bx-1,1)) = 1;
                    em(yy, bx) = 1;
                    em(yy, min(bx+1,szX)) = 1;
                end
            elseif (mask(y, bx)==1 && mask(y+1, bx) == 0)
                for yy=y+1:t2y.y
                    em(yy, max(bx-1,1)) = 1;
                    em(yy, bx) = 1;
                    em(yy, min(bx+1,szX)) = 1;
                end
            end
        end
    end
    %figure(1),imshow(mask,[]);
    %figure(2),imshow(em, []);
end