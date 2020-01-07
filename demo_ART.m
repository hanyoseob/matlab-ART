%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA 

%%
clear ;
close all ;
home ;

%% GPU Processing
% If there is GPU device on your machine, 
% then isgpu is true. Otherwise, it is false.
bgpu    = true;
bfig    = true;

%% SYSTEM SETTING
N       = 512;
VIEW    = 360;
THETA   = linspace(0, 180, VIEW + 1);   THETA(end) = [];

A       = @(x) radon(x, THETA);
AT      = @(y) iradon(y, THETA, 'none', N)/(pi/(2*length(THETA)));
AINV    = @(y) iradon(y, THETA, N);

%% DATA GENERATION
load('XCAT512.mat');
x       = imresize(double(XCAT512), [N, N]);
p       = A(x);
x_full  = AINV(p);

%% LOW-DOSE SINOGRAM GENERATION
i0              = 1e4;
pn              = exp(-p);
pn              = i0.*pn;
pn              = poissrnd(pn);
pn              = max(-log(max(pn,1)./i0),0);
pn(isnan(pn))   = 0;
pn(isinf(pn))   = 0;

y               = pn;
% y       = p;

%% Algebraic Reconstruction Technique (ART) INITIALIZATION
x_low   = AINV(y);
x0      = zeros(size(x));
lambda  = 1e0;
niter   = 2e2;
bpos    = true;

%% Run Algebraic Reconstruction Technique (ART)
if bgpu
    y   = gpuArray(y);
    x0  = gpuArray(x0);
end

x_art   = ART(A, AT, y, x0, lambda, niter, bpos, bfig);

%% CALCUATE QUANTIFICATION FACTOR 
x_low       = max(x_low, 0);
x_art       = max(x_art, 0);
nor         = max(x(:));

mse_x_low   = immse(x_low./nor, x./nor);
mse_x_art   = immse(x_art./nor, x./nor);

psnr_x_low 	= psnr(x_low./nor, x./nor);
psnr_x_art 	= psnr(x_art./nor, x./nor);

ssim_x_low  = ssim(x_low./nor, x./nor);
ssim_x_art  = ssim(x_art./nor, x./nor);


%% DISPLAY
wndImg  = [0, 0.03];

figure('name', 'Algebraic Reconstruction Technique (ART)');
colormap(gray(256));

suptitle('Algebraic Reconstruction Technique (ART)');
subplot(221);   imagesc(x,     	wndImg); 	axis image off;     title('ground truth');
subplot(222);   imagesc(x_full, wndImg);   	axis image off;     title(['full-dose_{FBP, view : ', num2str(VIEW) '}']);
subplot(223);   imagesc(x_low,  wndImg);   	axis image off;     title({['low-dose_{FBP, view : ', num2str(VIEW) '}'], ['MSE : ' num2str(mse_x_low, '%.4e')], ['PSNR : ' num2str(psnr_x_low, '%.4f')], ['SSIM : ' num2str(ssim_x_low, '%.4f')]});
subplot(224);   imagesc(x_art,  wndImg);  	axis image off;     title({['recon_{ART, view : ', num2str(VIEW) '}'], ['MSE : ' num2str(mse_x_art, '%.4e')], ['PSNR : ' num2str(psnr_x_art, '%.4f')], ['SSIM : ' num2str(ssim_x_art, '%.4f')]});

