function wb = wishart_blockdiag(X,L,varargin)

% Wishart based detection of reflection symmetry in polSAR covmat data
% or test for diagonal-only in 3x3 or 2x2 covmat data
%
% wb = wishart_blockdiag(X,L);
%
% ShhShh* ShhShv* ShhSvv*     k    a   rho
% ShvShh* ShvShv* ShvSvv*  =  a*   ksi b
% SvvShh* SvvShv* SvvSvv*     rho* b*  zeta
%
% X is the (row-wise) upper-right elements of
% the covariance matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%   ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
% = [k Re(a) Im(a) Re(rho) Im(rho) ksi Re(b) Im(b) zeta]
%
% Input
%   X	   - covariance matrix input
%   L	   - equivalent number of looks (ENL) for X
%   nvars  - optional
%            if nvars==3: test for diagonal only in 3x3 matrix;
%            else nvars must be 4 or 9 depending on X
%   Ptest  - optional
%            logical, output P-values (and rho, omega2 and f), default: true
%
% Output: struct wb with fields
%   lnQ     - logarithm of test statistic Q
%   L       - equivalent number of looks, ENL (from input)
%   rho     - aux variable
%   omega2  - aux variable
%   f       - number of degrees of freedom
%   P       - probability of finding smaller -2 \rho \ln Q
%
% References
%
% Paul Connetable, Knut Conradsen, Allan A. Nielsen and Henning Skriver (2022).
% Test statistics for reflection symmetry, applications to quad-polarimetric
% SAR data for detection of man-made structures.
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote
% Sensing 15, 2877-2890. Open Access DOI:10.1109/JSTARS.2022.3162670.
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2025).
% A test statistic for block-diagonal covariance matrix structure in
% polSAR data. IEEE Geoscience and Remote Sensing Letters 22, 4012305.
% DOI:10.1109/LGRS.2025.3605978. 

% (c) Copyright 2021
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 6 Mar 2021
% 18 Oct 2023 /AA: complex lnQ and lnQ>0 fixed

if nargin<2 || nargin>4
    help wishart_blockdiag
    error('wishart_blockdiag: wrong input');
end
if nargout>1
    help wishart_blockdiag
    error('wishart_blockdiag: wrong output');
end

if ischar(X), X = freadenvisar(X); end

if ndims(X)~=3
    error('wishart_blockdiag: X must be 3-dim')
end
pX = size(X);
nvars = pX(3);
if nargin>2
    nvars = varargin{1};
end
Ptest = true;
if nargin>3
    Ptest = varargin{2};
end

wb.L = L;

if nvars==9
    % fprintf('wishart_blockdiag: quad/full pol, is matrix reflection symmetric?\n');
    wb.lnQ = L*(log(wishart_det(X))-log(wishart_det(X(:,:,[1 4 5 9])))-log(X(:,:,6)));
    wb.rho = 1 - 3/2/L;
    wb.omega2 = 5/12/L^2/wb.rho^2;
    p1 = 2;
    p2 = 1;
    p3 = 0;
elseif nvars==4
    if not(pX(3)==4)
        error('*** wishart_blockdiag: test 2x2 for diagonal only but matrix is not 2x2 ***')
    end
    % fprintf('wishart_blockdiag: dual pol, is matrix diagonal only?\n');
    wb.lnQ = L*(log(wishart_det(X))-log(X(:,:,1))-log(X(:,:,4)));
    wb.rho = 1 - 1/L;
    wb.omega2 = 1/12/L^2/wb.rho^2;
    p1 = 1;
    p2 = 1;
    p3 = 0;
elseif or(nvars==3,nvars==19) % '19' because '9' is already used
    if not(pX(3)==9)
        error('*** wishart_blockdiag: test 3x3 for diagonal only but matrix is not 3x3 ***')
    end
    % fprintf('wishart_blockdiag: quad/full pol, is matrix diagonal only?\n');
    wb.lnQ = L*(log(wishart_det(X))-log(X(:,:,1))-log(X(:,:,6))-log(X(:,:,9)));
    wb.rho = 1 - 4/3/L;
    wb.omega2 = 7/12/L^2/wb.rho^2;
    p1 = 1;
    p2 = 1;
    p3 = 1;
else
    error('wishart_blockdiag: input must be 9 or 4 (real) bands');
end
clear X

if ~isreal(wb.lnQ)
    warning('*** some lnQ values are complex, only real part is used ***')
    fprintf('max(imag(lnQ(:))) = %g\n', max(imag(wb.lnQ(:))));
    wb.lnQ = real(wb.lnQ);
end
if max(wb.lnQ(:))>0
    warning('*** lnQ > 0, bad ***')
    fprintf('max(lnQ(:)) = %g\n', max(wb.lnQ(:)));
    wb.lnQ(wb.lnQ>0) = 0;
end

p = p1 + p2 + p3;
% psump2 = p^2 - (p1^2 + p2^2 + p3^2);
% psump3 = p^3 - (p1^3 + p2^3 + p3^3);
% wb.f = psump2;
wb.f = p^2 - (p1^2 + p2^2 + p3^2);
% wb.rho_ = 1 - psump3/psump2/3/L;
if Ptest
    % psump4 = p^4 - (p1^4 + p2^4 + p3^4);
    % % wb.omega2_ = (-psump3^2/psump2/36 + psump4/24)/(L^2*wb.rho^2)^2; % typo
    % wb.omega2_ = (-psump3^2/psump2/36 + psump4/24)/L^2/wb.rho^2;
    wb.P = (1 - wb.omega2)*gammainc(-wb.rho*wb.lnQ,0.5*wb.f  ) ...
              + wb.omega2 *gammainc(-wb.rho*wb.lnQ,0.5*wb.f+2);
    % wb.rho = 1;
    % wb.omega2 = 0;
    % wb.P = gammainc(-wb.lnQ,0.5*wb.f);
end

if false

%% Generated data
rng('default')
Sigma = wc2mat([4 0 0 1 1 2 0 0 3]) % refection symmetry
chol(Sigma)
L = 13;
% or
cd Foulum; L64 = flipud(freadenvisar('fl064_l/fl064_l','b')); L = 13; cd ..
[simimage,coord,Sigma] = wishart_rand_patch(L,L64,1,1,1)

simimageS = wishart_rand(L,1024,1024,Sigma,1); % reflection symmetry
wb_diag = wishart_blockdiag(simimageS,L,3,true) % diag-only test
wc_Phisto(wb_diag,nan,'diag-only test') 
wb_refl = wishart_blockdiag(simimageS,L) % reflection symmetry test
wc_Phisto(wb_refl,nan,'reflection symmetry test')
% is ShhSvv* zero?
wb_ShhSvv = wishart_blockdiag(simimageS(:,:,[1 4 5 9]),L) % is ShhSvv* zero?
wc_Phisto(wb_ShhSvv,nan,'is ShhSvv* zero?')

simimage1 = wishart_rand(L,1024,1024,eye(3),1); % diag-only
wb1_diag = wishart_blockdiag(simimage1,L,3,true) % diag-only test
wc_Phisto(wb1_diag,nan,'diag-only test') 
wb1_refl = wishart_blockdiag(simimage1,L) % reflection symmetry test
wc_Phisto(wb1_refl,nan,'reflection symmetry test')
% is ShhSvv* zero?
wb1_ShhSvv = wishart_blockdiag(simimage1(:,:,[1 4 5 9]),L) % is ShhSvv* zero?
wc_Phisto(wb1_ShhSvv,nan,'is ShhSvv* zero?')

%% F-SAR X-band
cd FSAR/Vejers
load Xband_Vejers.mat
cd ../..
SAR = cat(3,HHHH,real(HHHV),imag(HHHV),real(HHVV),imag(HHVV),HVHV,real(HVVV),imag(HVVV),VVVV);
SAR = double(SAR(:,1601:1601+4800,:)); % used in paper
% SAR = double(SAR(:,1601:1601+4800-1,:));
% % SAR = mult_sqrt2(SAR);
% clear HHHH HHHV HHVV HVHV HVVV VVVV
% % enlHHHH = ENLfilter(SAR(:,:,1));
% % enlHVHV = ENLfilter(SAR(:,:,6));
% % enlVVVV = ENLfilter(SAR(:,:,9));
% enl = ENLfilter(SAR(:,:,[1 6 9]));
% figure, histogram(enl), grid on, grid minor
L = 3;
L = 4;
L = 20; % says Paul

% cd Foulum; SAR = flipud(freadenvisar('fl064_l/fl064_l','b')); L = 13; cd ..
% cd radarsat2bonn; SAR = imread('20090525.tif'); L = 12; cd ..
% cd radarsat2bonn; SAR = double(imread('20090525.tif')); L = 12; cd ..

% imshowsar(complex(SAR(:,:,4),SAR(:,:,5)))

wb = wishart_blockdiag(SAR,L)
[nr,nc,nb] = size(SAR);
SARrefl = SAR;
SARrefl(:,:,[2 3 7 8]) = zeros(nr,nc,4);
[lnQ,omega2,rho,P,f] = wishart_change(SAR,SARrefl,L);
% wc = wcRjl(cat(4,SAR,SARrefl),L,0)

% ihplot = false;
ihplot = true;
mv = 0.005; % mv = 0.02; % subplot_tight vertical margin 
mh = 0.005; % mh = 0.01; % subplot_tight horizontal margin

pval = 0.001;
hf = figure;
h(1) = subplot_tight(3,3,1,[mv,mh]);
% % imshowsarrgb(SAR)
% % imshowsarrgbprctile(SAR)
% imshowsarpauli(SAR)
vz = van_zyl(SAR)
imshowrgbprctile(vz.Lambda(:,:,[2 3 1]));
% imshowrgbprctile(10*log10(vz.Lambda(:,:,[2 3 1])));
h(2) = subplot_tight(3,3,2,[mv,mh]);
aux = -2*wb.rho*wb.lnQ;
imshow(aux,[0,prctile(aux(:),99.9)])
h(3) = subplot_tight(3,3,3,[mv,mh]);
auz = -2*rho*lnQ;
imshow(auz,[0,prctile(auz(:),99.9)])
h(4) = subplot_tight(3,3,5,[mv,mh]);
imshow(wb.P,[0 1])
h(5) = subplot_tight(3,3,6,[mv,mh]);
imshow(P,[0 1])
h(6) = subplot_tight(3,3,8,[mv,mh]);
imshow(wb.P>1-pval)
h(7) = subplot_tight(3,3,9,[mv,mh]);
imshow(P>1-pval)
if ihplot
    rho_cplx = rho_llrr(SAR); % this is rho_LLRR
    [nrows,ncols] = size(rho_cplx);
    amp = abs(rho_cplx);%.^2;
    phas = angle(rho_cplx);
    h(8) = subplot_tight(3,3,4,[mv,mh]);
    % intensity-hue/hue-value plot of amplitude and phase
    imshow(hsv2rgb(cat(3,(phas+pi)/(2*pi),ones(nrows,ncols),amp)))
%     colormap(gca,hsv(256)) % not needed for hsv2rgb but for colorbar
%     cb = colorbar('SouthOutside');
%     set(cb,'XTick',(0:64:256)/256)
%     set(cb,'XTickLabel',{'-180','-90','0','90','180'})
%     drawnow
    rho_cplx = conj(rho_cplx); % this is rho_RRLL, not rho_LLRR
    [nrows,ncols] = size(rho_cplx);
%     amp = abs(rho_cplx);
    phas = angle(rho_cplx);
    h(9) = subplot_tight(3,3,7,[mv,mh]);
    % intensity-hue/hue-value plot of amplitude and phase
    imshow(hsv2rgb(cat(3,(phas+pi)/(2*pi),ones(nrows,ncols),amp)))
%     colormap(gca,hsv(256)) % not needed for hsv2rgb but for colorbar
%     cb = colorbar('SouthOutside');
%     set(cb,'XTick',(0:64:256)/256)
%     set(cb,'XTickLabel',{'-180','-90','0','90','180'})
%     drawnow
end
linkaxes(h,'xy')

% idx = randsample(nr*nc,max(nr,nc));
% % idx = 1:(nr*nc);
maxx = max(aux(:));
maxz = max(auz(:));
auxx = max(maxx,maxz);
figure
plot(aux(:),auz(:),'.')
title('Test statistic: $-2 \rho \ln Q$','FontSize',16,...
    'Interpreter','latex')
hold on
plot([0 auxx],[0 auxx],'LineWidth',1.5)
axis equal square
grid on
xlabel('''Correct'' test statistic, df=4')
ylabel('''Naive'' test statistic, df=9')

figure
plot(wb.P(:),P(:),'.')
title('$P\{$not refl. symm.$\} = P\{<-2 \rho \ln Q\}$','FontSize',16,...
    'Interpreter','latex')
hold on
plot([0 1],[0 1],'LineWidth',1.5)
axis equal square
grid on
xlabel('''Correct'' test statistic, df=4')
ylabel('''Naive'' test statistic, df=9')

% figure(hf.Number)
% subplot_tight(3,3,4,7*[mv,mh]);
% plot(aux(:),auz(:),'.'), grid on
% subplot_tight(3,3,7,7*[mv,mh]);
% plot(wb.P(:),P(:),'.'), axis equal square, grid on

end % false

return

%% Dependence on ENL, p, rho, omega_2

init = false;

if init
    % L = [3 4 5 10 20 30 40 50];
    L = [3:10 15 20 30 40 50];
    L0 = L; % for resetting below
    p = 3; % 2 3
    % k is #time points 
    k = 1; %52; % 5 10 20 50 100
    nrep = 10; % number of repetitions
    nr = 100;
    nc = 1000;
end


%% generate data
tic
% SAR_L1 = struct('rep',nan(nr,nc,1,k));
SAR_L2 = struct('rep',nan(nr,nc,2*2,k));
SAR_L  = struct('rep',nan(nr,nc,3*3,k));
for ix=1:nrep
    % SAR_L1(ix).rep = nan(nr,nc,1,k);
    SAR_L2(ix).rep = nan(nr,nc,2*2,k);
    for ixL=1:numel(L)
        SAR_L(ix,ixL).rep = nan(nr,nc,3*3,k);
    end
end
toc
% rng('default')
tic
for ix=1:nrep
    % SAR_L1(ix).rep = wishart_rand_silent(1,nr,nc,eye(1),k); %,1);
    SAR_L2(ix).rep = wishart_rand_silent(2,nr,nc,eye(2),k); %,1);
    for ixL=1:numel(L)
        SAR_L(ix,ixL).rep = wishart_rand_silent(L(ixL),nr,nc,eye(3),k); %,1);
    end
end
toc

%% 3x3 matrices
tic
entrpy_wb_r  = nan(nrep,numel(L));
entrpy_wb_r1 = nan(nrep,numel(L));
for ixr=1:nrep
    for ixL=1:numel(L)
        wb = wishart_blockdiag(SAR_L(ixr,ixL).rep,L(ixL)); % reflection symmetry test
        pP = wc_entropy(wb,nan);
        entrpy_wb_r(ixr,ixL) = pP(3); % Shannon entropy
        wb = wishart_blockdiag_r1(SAR_L(ixr,ixL).rep,L(ixL)); % reflection symmetry test
        pP = wc_entropy(wb,nan);
        entrpy_wb_r1(ixr,ixL) = pP(3); % Shannon entropy
    end
end
ms3x3r = [mean(entrpy_wb_r);  std(entrpy_wb_r)
          mean(entrpy_wb_r1); std(entrpy_wb_r1)];
toc

%% 2x2 matrices
tic
entrpy_wb_d  = nan(nrep,numel(L)+1);
entrpy_wb_d1 = nan(nrep,numel(L)+1);
for ixr=1:nrep
    ixL = 1; % start with L=2 from SAR_L2
    wb = wishart_blockdiag(SAR_L2(ixr).rep,2); % diag only test
    pP = wc_entropy(wb,nan);
    entrpy_wb_d(ixr,ixL) = pP(3); % Shannon entropy
    wb = wishart_blockdiag_r1(SAR_L2(ixr).rep,2); % diag only test
    pP = wc_entropy(wb,nan);
    entrpy_wb_d1(ixr,ixL) = pP(3); % Shannon entropy
    for ixL=1:numel(L)
        wb = wishart_blockdiag(SAR_L(ixr,ixL).rep(:,:,6:9),L(ixL)); % diag only test
        pP = wc_entropy(wb,nan);
        entrpy_wb_d(ixr,ixL+1) = pP(3); % Shannon entropy
        wb = wishart_blockdiag_r1(SAR_L(ixr,ixL).rep(:,:,6:9),L(ixL)); % diag only test
        pP = wc_entropy(wb,nan);
        entrpy_wb_d1(ixr,ixL+1) = pP(3); % Shannon entropy
    end
end
ms2x2d = [mean(entrpy_wb_d);  std(entrpy_wb_d)
          mean(entrpy_wb_d1); std(entrpy_wb_d1)];
toc

%% plots

ax = nan(4,1);

mv = 0.08;
mh = 0.08;
lw = 2;

%% plot

figure('Name','3x3 matrices')

ax(1) = subplot_tight(2,1,1,[mv,mh]);
errorbar(L,ms3x3r(1:2:end,:),ms3x3r(2:2:end,:),'LineWidth',lw)
grid on, grid minor
xlabel('ENL')
ylabel('Entropy')
legend('3x3 refl \rho, \omega_2','3x3 refl','Location','southeast')
% title("Block diagonal")
xlim([0 51])
% ylim([0.75 1.01])

ax(2) = subplot_tight(2,1,2,[mv,mh]);
errorbar(L,ms3x3r(1:2:end,:),ms3x3r(2:2:end,:),'LineWidth',lw)
grid on, grid minor
xlabel('ENL')
ylabel('Entropy')
legend('3x3 refl \rho, \omega_2','3x3 refl','Location','southeast')
% title("Block diagonal")
xlim([0 51])
ylim([0.9989 1])


figure('Name','2x2 matrices')

L = [2 L0];

ax(3) = subplot_tight(2,1,1,[mv,mh]);
errorbar(L,ms2x2d(1:2:end,:),ms2x2d(2:2:end,:),'LineWidth',lw)
grid on, grid minor
xlabel('ENL')
ylabel('Entropy')
legend('2x2 diag \rho, \omega_2','2x2 diag','Location','southeast')
% title("Block diagonal")
xlim([0 51])
% ylim([0.9 1.003])

ax(4) = subplot_tight(2,1,2,[mv,mh]);
errorbar(L,ms2x2d(1:2:end,:),ms2x2d(2:2:end,:),'LineWidth',lw)
grid on, grid minor
xlabel('ENL')
ylabel('Entropy')
legend('2x2 diag \rho, \omega_2','2x2 diag','Location','southeast')
% title("Block diagonal")
xlim([0 51])
ylim([0.996 1])

L = L0;
