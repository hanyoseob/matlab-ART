%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

%%
function x  = ART(A,AT,b,x,lambda,niter,bpos)

if (nargin < 7)
    bpos	= true;
end

if (nargin < 6)
    niter       = 1e2;
end

ATA	= AT(A(ones(size(x), 'single')));

for i = 1:niter
    
    x  	= x + lambda*AT(b - A(x))./ATA;
    
    if (bpos)
        x(x < 0) = 0;
    end

    figure(1); colormap gray;
    imagesc(x);
    axis image off;
    title(num2str([i, niter], '%d / %d'));
    drawnow();
end

x   = gather(x);

end