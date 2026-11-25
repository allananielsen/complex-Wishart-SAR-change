function [simimage,coord,Sigma] = wishart_rand_patch(L,SAR,np,varargin)

% Generate observations from polSAR imagery W(ENL,SAR) (by click or coord)
%
% [simimage,coord,Sigma] = wishart_rand_patch(L,SAR,np,varargin)
%
% Choose np pixels in SAR image bmo ginput(np) or from coord as seed Sigma
% for generation of upper right (UR) part of complex Wishart distributed
% covariance matrices
%
% Based on work by Knut Conradsen and Morton J. Canty
%
% Input
%   L       - ENL, equivalent number of looks, must be an integer  
%   SAR     - 3-D, nrows x ncols polSAR scene with UR 2x2 or 3x3 cov. mat.
%   np      - number of patches
%   nrows   - #rows generated per patch (optional, default is 256)
%   ncols   - #cols generated per patch (optional, default is nrows)
%   coord   - if present we do not use ginput but the pixel coordinates in
%             coord with one row per patch (this overwrites the value of np)
%   if nargin>6 plot 1-D histograms of generated data
%
% Remember, that here we generate covariance matrices.
% To obtain complex Wishart distributed matrices we must therefore multiply
% by ENL:
%
% C = wishart_rand_patch(L,SAR,4); % polSAR covariance matrix representation
% W = L*C; % complex Wishart distributed
%
% or
%
% np = 9;
% rng('default')
% [s1,coord,Sigma] = wishart_rand_patch(L,SAR,np,128,128);
% s2 = wishart_rand_patch(L,SAR,np,128,128,coord);
% s3 = wishart_rand_patch(L,SAR,np,128,128,coord);
% s4 = wishart_rand_patch(L,SAR,np,128,128,coord);
% s14 = cat(4,s1,s2,s3,s4);
% clear s1 s2 s3 s4
% vzx = van_zylx(s14)
% figure, imshowrgbprctile(vzx(1).eigvals,[2 3 1],[1 99],100)
% wcfull = wcRjl(s14,L)
% [nr,nc,~,~] = size(s14);
% pQRj = wc_lnRtable(wcfull,true(nr,nc));
% [pP,pE] = wc_Ptable(wcfull,true(nr,nc));
% matrix2latex(pE,6)
%
% See wishart_rand, wishart_rand_scene
%
% References
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "Generation of Sample Complex Wishart Distributed Matrices and Change
% Detection in Polarimetric SAR Data".
% International Conference on Digital Image and Signal Processing,
% DISP 2019, Article no. 031, St Hugh's College, Oxford University,
% England, 29-30 April 2019.
% https://doi.org/10.48550/arXiv.2611.22943.
% https://www2.imm.dtu.dk/pubdb/pubs/7124-full.html
%
% K. Conradsen, H. Skriver and A. A. Nielsen (2024).
% “Omnibus change detection in block diagonal covariance matrix polSAR data
% illustrated with simulated and Sentinel-1 data,”
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing
% pp. 16359–16376,
% https://doi.org/10.1109/JSTARS.2024.3453442
%
% A. A. Nielsen, H. Skriver and K. Conradsen (2025).
% "Multifrequency omnibus change detection in covariance matrix polSAR data,"
% IEEE Geoscience and Remote Sensing Letters 22, 4005505,
% https://doi.org/10.1109/LGRS.2025.3541861
%
% A. A. Nielsen, H. Skriver and K. Conradsen (2025).
% "A test statistic for block-diagonal covariance matrix structure in polSAR data",
% IEEE Geoscience and Remote Sensing Letters 22, 4012305,
% https://doi.org/10.1109/LGRS.2025.3605978
%
% A. A. Nielsen, H. Skriver and K. Conradsen (2026).
% "Generation of sample complex Wishart distributed matrices and tests for
% equality in polarimetric SAR data",
% in prep.


% (c) Copyright 2018-2025
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 23-24 Oct 2018
% 17 Jul 2024
% 8 Aug 2024 /AA: k added
% 23 Dec 2024 /AA: adapted to generate patches
% 26 Jan 2025 /AA: divide by L

% s = rng;
fprintf("\nYou might want to use rng('default')\n");
fprintf("\n***********************************************************\n");
fprintf("*** wishart_rand_patch generates covariance matrix data ***\n")
fprintf("*** multiply by L to get Wishart distributed data       ***\n");
fprintf("***********************************************************\n");

if nargin<3
    error('wishart_rand_patch(L,SAR,k), SAR must be 3-D nrows x ncols with 2x2 or 3x3 cov. mat.')
end
if nargin>7
    error('wishart_rand_patch(L,SAR,k,1), SAR must be 3-D nrows x ncols with 2x2 or 3x3 cov. mat.')
end

if ismatrix(SAR) | ndims(SAR)>4 % SAR must be 3-D or 4-D
% if ndims(SAR)<3 | ndims(SAR)>4 % SAR must be 3-D or 4-D
    error('wishart_rand_patch: SAR must be 3-D nrows x ncols with 2x2 or 3x3 cov. mat.')
end
if ndims(SAR)==4 % SAR is 4-D
    % [~,~,~,nt] = size(SAR);
    fprintf('wishart_rand_scene: SAR is 4-D, we use first ''time point'' only\n')
    SAR = squeeze(SAR(:,:,:,1));
end
[nrows,ncols,pp] = size(SAR); % SAR is 3-D
if ~(pp==4 || pp==9)
    error('wishart_rand_scene: cov. mat. must be 2x2 or 3x3')
end
p = sqrt(pp);

nrows = 256;%128;
if nargin>3
    nrows = varargin{1};
end
ncols = nrows;
if nargin>4
    ncols = varargin{2};
end
nobs = nrows*ncols;

%% ginput or coord
if nargin>5 % coord
    coord = varargin{3};
    [np,~] = size(coord); % overwrite value of np in input
else % ginput
    ii = 1;
    if pp==4
        % aux = cat(3,S1(:,:,4,ii),S1(:,:,1,ii),S1(:,:,4,ii)./S1(:,:,1,ii));
        % aux = cat(3,S1(:,:,4,ii),S1(:,:,1,ii),sqrt(S1(:,:,2,ii).^2+S1(:,:,3,ii).^2));
        aux = cat(3,SAR(:,:,4,ii),SAR(:,:,1,ii),abs(complex(SAR(:,:,2,ii),SAR(:,:,3,ii))));
    else % pp=9
        aux = cat(3,SAR(:,:,1,ii),SAR(:,:,6,ii),SAR(:,:,9,ii));
    end
    % figure, imshowrgbprctile(aux)
    % figure, imshowrgbprctile(sqrt(aux))
    figure, imshowrgbprctile(10*log10(aux));
    fprintf('\n*** click on %d representative pixels... ***\n\n',np)
    coord = int32(ginput(np));
    fprintf('*** ginput done... ***\n\n')
end
Sigma = nan(p,p,np);
for ir=1:np
    % Sigma(:,:,ir) = wc2mat(SAR(int32(coord(ir,2)),int32(coord(ir,1)),:));
    Sigma(:,:,ir) = wc2mat(SAR(coord(ir,2),coord(ir,1),:));
end

%%
simimage = nan(nobs,pp,np);
%simimage = nan(nobs,np,pp);
for ik = 1:np
    for ix=1:nobs
        X = randn(L,p);
        Y = randn(L,p);
        W = 0.5*(X'*X + Y'*Y - 1j*(X'*Y - Y'*X));
        C = chol(Sigma(:,:,ik));
        W = C'*W*C;
        if p==3
            simimage(ix,1,ik) = real(W(1,1));
            simimage(ix,2,ik) = real(W(1,2));
            simimage(ix,3,ik) = imag(W(1,2));
            simimage(ix,4,ik) = real(W(1,3));
            simimage(ix,5,ik) = imag(W(1,3));
            simimage(ix,6,ik) = real(W(2,2));
            simimage(ix,7,ik) = real(W(2,3));
            simimage(ix,8,ik) = imag(W(2,3));
            simimage(ix,9,ik) = real(W(3,3));
        else
            simimage(ix,1,ik) = real(W(1,1));
            simimage(ix,2,ik) = real(W(1,2));
            simimage(ix,3,ik) = imag(W(1,2));
            simimage(ix,4,ik) = real(W(2,2));
        end
    end
end

simimage = reshape(simimage,nrows,ncols,pp,np);
simimage = permute(simimage,[1 2 4 3]);
simimage = reshape(simimage,nrows,ncols*np,pp);
simimage = simimage/L;

%% re-arrange output image
if mod(np,4)==0 && np>4 % 4 rows of patches
    naux = ncols*np/4;
    s1 = simimage(:,1:naux,:);
    s2 = simimage(:,naux+1:2*naux,:);
    s3 = simimage(:,2*naux+1:3*naux,:);
    s4 = simimage(:,3*naux+1:4*naux,:);
    simimage = cat(1,s1,s2,s3,s4);
    % % 4 cols of patches
    % simimage = permute(simimage,[2 1 3]);
    clear s1 s2 s3 s4
elseif mod(np,3)==0 && np>3 % 3 rows of patches
    naux = ncols*np/3;
    s1 = simimage(:,1:naux,:);
    s2 = simimage(:,naux+1:2*naux,:);
    s3 = simimage(:,2*naux+1:3*naux,:);
    simimage = cat(1,s1,s2,s3);
    clear s1 s2 s3
elseif mod(np,2)==0 && np>2% np is even, 2 rows of patches
    naux = ncols*np/2;
    s1 = simimage(:,1:naux,:);
    s2 = simimage(:,naux+1:2*naux,:);
    simimage = cat(1,s1,s2);
    clear s1 s2
end

% simimage = reshape(simimage,nrows,ncols,np,pp);

% if k==1
%     simimage = squeeze(simimage);
% end

% rng(s);

%%
if nargin>6
    figure('Name','Generated matrix data')
    for ix=1:p*p
        ax(ix) = subplot_tight(p,p,ix,[0.06,0.06]);
        h = histogram(simimage(:,:,ix));
        if nrows*ncols>100*100+1 % or some other 'high' number
            h.EdgeColor = 'none';
        end
        grid on
    end
    % linkaxes(ax,'xy')
    if p==3
        linkaxes([ax(1) ax(6) ax(9)],'xy') % real elements
        % xlim([900 1350])
        linkaxes([ax(2) ax(3) ax(4) ax(5) ax(7) ax(8)],'xy') % complex
        % xlim([900 1350])
    else
        linkaxes([ax(1) ax(4)],'xy') % real elements
        % xlim([900 1350])
        linkaxes([ax(2) ax(3)],'xy') % complex
        % xlim([900 1350])
    end
end

return

%% Simulated VV from Sentinel-1

cd S1_117
L = 6;
L = 9; % says Henning Skriver
F = dir('*_.hdr');
S1 = freadenvisar1(F(1).name(1:10),'n');
cd ..
S1 = permute(S1,[2 1 3]);
SAR = S1(681:1280,1:2000,:);
np = 4; % number of patches/Sigmas
nr = 256;
nc = 256;
rng('default')
[s1,coord,Sigma] = wishart_rand_patch(L,SAR,np,nr,nc);
aux = 10*log10(cat(3,s1(:,:,4),s1(:,:,1),abs(complex(s1(:,:,2),s1(:,:,3)))));
figure, imshowrgbprctile(aux);
s2 = wishart_rand_patch(L,SAR,np,nr,nc,coord);
s3 = wishart_rand_patch(L,SAR,np,nr,nc,coord);
s4 = wishart_rand_patch(L,SAR,np,nr,nc,coord);
s5 = wishart_rand_patch(L,SAR,np,nr,nc,coord);
s6 = wishart_rand_patch(L,SAR,np,nr,nc,coord);
s16 = cat(4,s1,s2,s3,s4,s5,s6);
clear s1 s2 s3 s4 s5 s6
wc = wcRjl(s16(:,:,4,:),L,false) % VV only
pQR = wc_lnRtable(wc,NaN,'Generated, VV only')
pE = wc_Ptable(wc,NaN,'Generated, VV only')

%% Constructed Sigma

% single bounce
a = complex(2,0);
a2 = abs(a).^2;
S1 = [a2 0 0 real(a) imag(a) 0 0 0 1];
%S1 = S1 + 0.001*rand(1,9);
% double bounce
b = complex(-2,0);
b2 = abs(b).^2;
S2 = [b2 0 0 real(b) imag(b) 0 0 0 1];
%S2 = S2 + 0.001*rand(1,9);
% volume scattering (vegetation)
S3 = [0 0 0 0 0 2 0 0 0];
%S3 = S3 + 0.001*rand(1,9);

% mostly single
c = complex(1,1);
c2 = abs(c)^2;
S4 = [c2 0 0 real(c) imag(c) 0 0 0 1];
%S4 = S4 + 0.001*rand(1,9);

% mostly double
d = complex(-1,1);
d2 = abs(d).^2;
S5 = [d2 0 0 real(d) imag(d) 0 0 0 1];
%S5 = S5 + 0.001*rand(1,9);

% single + volume
S6 = S1 + S3;

% double + volume
S7 = S2 + S3;

% 
e = complex(0.1,2);
% e = complex(-0.1,2);
e2 = abs(e).^2;
S8 = [e2 0 0 real(e) imag(e) 0 0 0 1];

S9 = S8 + S3;


simimage = reshape(cat(1,S1,S2,S3,S4,S5,S6,S7,S8,S9),1,9,9);
vz = van_zyl(simimage)
vz.Lambda(1,1:9,:)
vzx = van_zylx(simimage)
vzx.eigvals(1,1:9,:)
vzf = van_zyl_full(simimage)
vzf.eigvals(1,1:9,:)

figure, imshowrgbprctile(vz.Lambda,[2 3 1],[1 99],10000) % magn 100x = 10000%
figure, imshowrgbprctile(vzx.eigvals,[2 3 1],[1 99],10000)
figure, imshowrgbprctile(real(vzf.eigvals),[2 3 1],[1 99],10000)

figure, imshowrgb(vz.Lambda,[2 3 1],2,10000)
figure, imshowrgb(vzx.eigvals,[2 3 1],2,10000)
figure, imshowrgb(real(vzf.eigvals),[2 3 1],2,10000)

