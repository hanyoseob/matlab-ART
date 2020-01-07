%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

%%
clear ;

isgpu   = true;

%%  SYSTEM SETTING
N       = 512;
ANG     = 180;
VIEW    = 360;
THETA   = linspace(0, ANG, VIEW + 1);   THETA(end) = [];

A       = @(x) radon(x, THETA);
AT      = @(y) iradon(y, THETA, 'none', N)/(pi/(2*length(THETA)));
AINV    = @(y) iradon(y, THETA, N);

%% DATA GENERATION
load('XCAT512.mat');
x       = imresize(double(XCAT512), [N, N]);
p       = A(x);
x_full  = AINV(p);

%% LOW-DOSE SINOGRAM GENERATION
i0      = 5e4;
pn      = exp(-p);
pn      = i0.*pn;
pn      = poissrnd(pn);
pn      = max(-log(max(pn,1)./i0),0);

y       = pn;

%% Algebraic Reconstruction Technique (ART) INITIALIZATION
x_low   = AINV(y);
x0      = zeros(size(x));
lambda  = 1e0;
niter   = 2e2;
bpos    = true;

%% RUN Algebraic Reconstruction Technique (ART)
if isgpu
    y = gpuArray(y);
    x0 = gpuArray(x0);
end

x_art   = ART(A, AT, y, x0, lambda, niter, bpos);

%% CALCULATE QUANTIFICATION FACTOR
x_low       = max(x_low, 0);
x_art       = max(x_art, 0);
nor         = max(x(:));

mse_x_low   = immse(x_low./nor, x./nor);
mse_x_art   = immse(x_art./nor, x./nor);

psnr_x_low  = psnr(x_low./nor, x./nor);
psnr_x_art  = psnr(x_art./nor, x./nor);

ssim_x_low  = ssim(x_low./nor, x./nor);
ssim_x_art  = ssim(x_art./nor, x./nor);


%% DISPLAY
wndImg  = [0, 0.03];
wndPrj  = [0, 6];

figure(1);
colormap(gray(256));

% suptitle('Algebraic Reconstruction Technique');
subplot(241);   imagesc(x,      wndImg);    axis image off;     title(['ground truth']);

subplot(242);   imagesc(x_full, wndImg);    axis image off;     title(['full-dose_{view : ', num2str(VIEW) '}']);
subplot(246);   imagesc(p, wndPrj);         xlabel(['Angle ( \Delta\theta : ' num2str(ANG/VIEW) ' \circ )']);   ylabel('Detector'); title(['full-dose_{view : ', num2str(VIEW) '}']);

subplot(243);   imagesc(x_low,  wndImg);    axis image off;     title({['low-dose_{view : ', num2str(VIEW) '}' ' using I0 = ' num2str(i0, '%.2e') ], ['MSE : ' num2str(mse_x_low, '%.4e')], ['PSNR : ' num2str(psnr_x_low, '%.4f')], ['SSIM : ' num2str(ssim_x_low, '%.4f')]});
subplot(247);   imagesc(y,  wndPrj);        xlabel(['Angle ( \Delta\theta : ' num2str(ANG/VIEW) ' \circ )']);   ylabel('Detector'); title(['low-dose_{view : ', num2str(VIEW) '} using I0 = ' num2str(i0, '%.2e')]);

subplot(244);   imagesc(x_art,  wndImg);    axis image off;     title({['recon_{art}'], ['MSE : ' num2str(mse_x_art, '%.4e')], ['PSNR : ' num2str(psnr_x_art, '%.4f')], ['SSIM : ' num2str(ssim_x_art, '%.4f')]});
subplot(248);   imagesc(y - p);             xlabel(['Angle ( \Delta\theta : ' num2str(ANG/VIEW) ' \circ )']);   ylabel('Detector'); title('full-dose - low-dose');


