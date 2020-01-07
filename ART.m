%% Reference 
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%%
function x  = ART(A,AT,b,x,lambda,n,bpos,bfig)

if (nargin < 8)
    bfig	= false;
end

if (nargin < 7)
    bpos	= true;
end

if (nargin < 6)
    n       = 1e2;
end

ATA	= AT(A(ones(size(x), 'single')));

for i = 1:n
    
    x  	= x + lambda*AT(b - A(x))./ATA;
    
    if (bpos)
        x(x < 0) = 0;
    end
    
    if (bfig)
        figure(1); colormap gray;
        imagesc(x);
        axis image off;
         title([num2str(i) ' / ' num2str(n)]);
        drawnow();
    end
end

x   = gather(x);

end