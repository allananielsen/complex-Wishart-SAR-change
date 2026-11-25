function vz = van_zylx(SAR)

% van_zylx calculates van Zyl's (Cloude's) decomposition of covariance data
%   reflection symmetry is assumed/enforced
%
% SAR input images SAR are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g., for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%  ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   SAR     - 3-D or 4-D array of 3x3 complex covariance matrix
%             (polarimetric SAR image in covariance matrix form,
%             reading order is alphabetical quad/full pol
%             (as read by freadenvisar)
%               [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*)
%                        ShvShv*                 Re(ShvSvv*) Im(ShvSvv*)
%                                                SvvSvv*]
%             or azumuthal symmetry
%               [ShhShh*                         Re(ShhSvv*) Im(ShhSvv*)
%                        ShvShv*
%                                                SvvSvv*]
%
%             If quad/full pol data are input, reflection symmetry is
%             enforced
%
% Output
%   vz - struct with fields
%           eigvals (unlike in polsar_eigvals, eigenvalues are not sorted)
%           Lambda
%           entropy (Shannon entropy)
%           alpha
%           beta
%
% Eigenvalues are not sorted (unlike in polsar_eigvals), but Lambda1 and
% Lambda2 (and alpha and beta) values are swapped if Re(alpha)<0 so that
% vz.Lambda(:,:,1) is associated with single (surface, odd) bounce scattering (blue)
% vz.Lambda(:,:,2) is associated with double (even) bounce scattering (red)
% vz.Lambda(:,:,3) is the remainder associated with volume scattering (green)
%
% See also cloude_pottier and polsar_entropy
%
% References
%
% [1] Jakob J. van Zyl (1992).
%     "Application of Cloude's target decomposition theorem to polarimetric
%     imaging radar data".
%     Proc. SPIE 1748, Radar Polarimetry, 184-191 (12 February 1993).
%     DOI:10.1117/12.140615
%
% [2] Jakob J. van Zyl, Motofumi Arii and Yunjin Kim (2011).
%     "Model-Based Decomposition of Polarimetric SAR Covariance Matrices
%     Constrained for Nonnegative Eigenvalues".
%     IEEE Transactions on Geoscience and Remote Sensing, 49(9), 3452-3459,
%     September 2011.
%     DOI:10.1109/TGRS.2011.2128325
%
% [3] Allan A. Nielsen (2019).
%     "Fast matrix based computation of eigenvalues and the Loewner order
%     in PolSAR data".
%     IEEE Geoscience and Remote Sensing Letters 17(10), 1727-1731, 2020.
%     DOI:10.1109/LGRS.2019.2952202

% (c) Copyright 2021
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 27 April-5 June 2021

warning('*** van_zylx is (experi)mental - use van_zyl ***')

%%
pX = size(SAR);
nrows = pX(1);
ncols = pX(2);
if  ndims(SAR)==3
    ntims = 1;
elseif ndims(SAR)==4
    ntims = pX(4);
else
    error('van_zylx: wrong dimensionality of input images');
end
nvars = pX(3);

if nvars==9 % full pol
    p = 3;
    fprintf('van_zylx: quad/full pol input, reflection symmetry enforced\n')
%     warning('van_zylx: quad/full pol input, reflection symmetry enforced')
elseif nvars==5 % full pol, reflection symmetry
    p = 3;
%     warning('van_zylx: reflection symmetry')
else
    error('van_zylx: input must be quad/full pol or reflection symmetry')
end

% full pol
if nvars==9
    SAR = SAR(:,:,[1 4 5 6 9],:);
end

for it=1:ntims
    %%
    if class(SAR)=='double'
        vz(it).eigvals = nan(nrows,ncols,p,'double');
        vz(it).Lambda = nan(nrows,ncols,p,'double');
    else
        vz(it).eigvals = nan(nrows,ncols,p,'single');
        vz(it).Lambda = nan(nrows,ncols,p,'single');
    end

    %%
    % reflection symmetry
    X = SAR(:,:,:,it);
    disc = sqrt((X(:,:,1)-X(:,:,5)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    vz(it).eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,5) + disc);
    vz(it).eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,5) - disc);
    vz(it).eigvals(:,:,3) = X(:,:,4); % eigvals(:,:,3) = Lambda(:,:,3)
    B = X(:,:,1)+X(:,:,4)+X(:,:,5); % trace
    p1 = vz(it).eigvals(:,:,1)./B;
    p2 = vz(it).eigvals(:,:,2)./B;
    p3 = vz(it).eigvals(:,:,3)./B;
    clear B
    vz(it).entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    clear p1 p2 p3
%     auxp = X(:,:,5)-X(:,:,1)+disc; % van Zyl
%     auxm = X(:,:,5)-X(:,:,1)-disc; % van Zyl
    auxp = X(:,:,5)-X(:,:,1)+disc; % AA
    auxm = X(:,:,5)-X(:,:,1)-disc; % AA
    clear disc
    tworho = 2*complex(X(:,:,2),X(:,:,3));
    idx = X(:,:,2)<0; % double bounce
    clear X
    
    %%
    tworho2 = tworho.*conj(tworho);
    auxp2 = auxp.^2;
    vz(it).Lambda(:,:,1) = real(vz(it).eigvals(:,:,1).*auxp2./(auxp2+tworho2));
    clear auxp2
    auxm2 = auxm.^2;
    vz(it).Lambda(:,:,2) = real(vz(it).eigvals(:,:,2).*auxm2./(auxm2+tworho2));
    clear tworho2
    vz(it).Lambda(:,:,3) = real(vz(it).eigvals(:,:,3)); % waste of memory
    clear auxm2
    vz(it).alpha = tworho./auxp;
    clear auxp
    vz(it).beta = tworho./auxm; 
    clear auxm tworho

    %% check must be zero - and it is! :-)
    % alpha2 = vz(it).alpha.*conj(vz(it).alpha);
    % beta2 = vz(it).beta.*conj(vz(it).beta);
    % [nrows,ncols,~] = size(X);
    % vz(it).check = X -...
    % repmat(vz(it).Lambda(:,:,1),1,1,5).*cat(3,alpha2,real(vz(it).alpha),imag(vz(it).alpha),zeros(nrows,ncols),ones(nrows,ncols)) -...
    % repmat(vz(it).Lambda(:,:,2),1,1,5).*cat(3,beta2,real(vz(it).beta),imag(vz(it).beta),zeros(nrows,ncols),ones(nrows,ncols)) -...
    % repmat(vz(it).Lambda(:,:,3),1,1,5).*cat(3,zeros(nrows,ncols),zeros(nrows,ncols),zeros(nrows,ncols),ones(nrows,ncols),zeros(nrows,ncols));
    % aux = vz(it).check;
    % % [min(vz(it).check(:)) max(vz(it).check(:))]

    %% if real(vz(it).alpha)<0: Lambda1 and alpha associated with even (double) bounce
    % after changes: Lambda1 and alpha always associated with odd (single) bounce, surface
% 	idx = real(tworho)<0; % double bounce
% %     idx = real(vz(it).alpha)<0;
%     % vz(it).idx = idx;
    aux = vz(it).alpha(idx);
    vz(it).alpha(idx) = vz(it).beta(idx);
    vz(it).beta(idx) = aux;
    clear aux
    Lamb1 = vz(it).Lambda(:,:,1);
    Lamb2 = vz(it).Lambda(:,:,2);
    aux = Lamb1(idx);
    Lamb1(idx) = Lamb2(idx);
    Lamb2(idx) = aux;
    clear aux
    vz(it).Lambda(:,:,1) = Lamb1;
    clear Lamb1
    vz(it).Lambda(:,:,2) = Lamb2;
    clear Lamb2

    %% if real(rho)<0: eigvals 1 and 2 both associated with odd (single) and even (double) bounce
    % after changes: eigvals 1 always associated with odd  (single) bounce, surface
    % after changes: eigvals 2 always associated with even (double) bounce
    Lamb1 = vz(it).eigvals(:,:,1);
    Lamb2 = vz(it).eigvals(:,:,2);
    aux = Lamb1(idx);
    Lamb1(idx) = Lamb2(idx);
    Lamb2(idx) = aux;
    clear aux idx
    vz(it).eigvals(:,:,1) = Lamb1;
    clear Lamb1
    vz(it).eigvals(:,:,2) = Lamb2;
    clear Lamb2

    %%
    vz(it).comment1 = 'Lambda1 and alpha are associated with single (odd) bounce (surface scattering)';
    vz(it).comment2 = 'Lambda2 and beta are associated with double (even) bounce';
    vz(it).comment3 = 'Lambda3 is associated with the remainder (typically volume scattering)';
    vz(it).comment4 = 'Low entropy is associated with a dominant scattering mechanism, high entropy with randomness';
    vz(it).comment5 = 'eigvals(:,:,3) = Lambda(:,:,3)';
    vz(it).comment6 = 'There is a typo in ref [1] eqs (16) and (17): the numerators under the squareroot sign should be squared';
    vz(it).comment7 = 'Visualization: RGB of eigvals (or Lambda) 2, 3, 1';
%     vz(it).comment = 'beta = -1./conj(alpha) is not stored';
    
    %% check must be zero - and it is! :-)
    % alpha2 = vz(it).alpha.*conj(vz(it).alpha);
    % beta2 = vz(it).beta.*conj(vz(it).beta);
    % [nrows,ncols,~] = size(X);
    % vz(it).check = X -...
    % repmat(vz(it).Lambda(:,:,1),1,1,5).*cat(3,alpha2,real(vz(it).alpha),imag(vz(it).alpha),zeros(nrows,ncols),ones(nrows,ncols)) -...
    % repmat(vz(it).Lambda(:,:,2),1,1,5).*cat(3,beta2,real(vz(it).beta),imag(vz(it).beta),zeros(nrows,ncols),ones(nrows,ncols)) -...
    % repmat(vz(it).Lambda(:,:,3),1,1,5).*cat(3,zeros(nrows,ncols),zeros(nrows,ncols),zeros(nrows,ncols),ones(nrows,ncols),zeros(nrows,ncols));
    % % [min(vz(it).check(:)) max(vz(it).check(:))]

    %% R: double bounce, G: volume scattering, B: single bounce
    % figure, imshowrgb(10*log10(vz(it).Lambda(:,:,[2 3 1])))
    % figure, imshowrgbprctile(10*log10(vz(it).Lambda(:,:,[2 3 1])))
    % figure, imshowrgbprctile(vz(it).Lambda(:,:,[2 3 1]))

%%
end % for it=1:ntims

%% Visualization
while false
figure('Name','van Zyl/Cloude')
beta = -vz.alpha./abs(vz.alpha).^2;
rr = vz.Lambda(:,:,2).*(abs(beta).^2+1);
clear beta
gg = vz.Lambda(:,:,3);
bb = vz.Lambda(:,:,1).*(abs(vz.alpha).^2+1);
imshowrgbprctile(cat(3,rr,gg,bb))
% imshowrgbprctile(10*log10(cat(3,rr,gg,bb)))
end

%%
while false

%% examples from van Zyl (1992)
X = nan(3,3,5);
% ocean surface
X(1,1,1) = 1.0000;
X(1,1,2) = 1.8298;
X(1,1,3) = 0.2557;
X(1,1,4) = 0.0952;
X(1,1,5) = 3.7877;
% clearcut area, P-band
X(2,1,1) = 1.0000;
X(2,1,2) = 0.0928;
X(2,1,3) = 0.0582;
X(2,1,4) = 0.5261;
X(2,1,5) = 0.5642;
% clearcut area, L-band
X(2,2,1) = 1.0000;
X(2,2,2) = 0.2324;
X(2,2,3) = 0.1057;
X(2,2,4) = 0.5308;
X(2,2,5) = 0.7580;
% clearcut area, C-band
X(2,3,1) = 1.0000;
X(2,3,2) = 0.3559;
X(2,3,3) = 0.0440;
X(2,3,4) = 0.4083;
X(2,3,5) = 0.7159;
% forested area, P-band
X(3,1,1) = 1.0000;
X(3,1,2) = 0.2803;
X(3,1,3) = 0.0167;
X(3,1,4) = 0.3301;
X(3,1,5) = 0.6529;
% forested area, L-band
X(3,2,1) = 1.0000;
X(3,2,2) = 0.3950;
X(3,2,3) = 0.0582;
X(3,2,4) = 0.3485;
X(3,2,5) = 0.7122;
% forested area, C-band
X(3,3,1) = 1.0000;
X(3,3,2) = 0.3669;
X(3,3,3) = 0.0016;
X(3,3,4) = 0.2416;
X(3,3,5) = 0.4685;
% [nrows,ncols,~] = size(X);
nrows = 1;
ncols = 1;
vzx = van_zylx(X)

aa = cart2polc(vzx.alpha);
bb = cart2polc(vzx.beta);
aa-bb % phase is pi

alpha2 = vzx.alpha(1,1).*conj(vzx.alpha(1,1));
% beta = vzx.beta;
beta = -vzx.alpha./abs(vzx.alpha).^2;
beta2 = beta(1,1).*conj(beta(1,1));
mata = repmat(vzx.Lambda(1,1,1),1,1,5).*cat(3,alpha2,real(vzx.alpha(1,1)),imag(vzx.alpha(1,1)),zeros(nrows,ncols),ones(nrows,ncols));
matb = repmat(vzx.Lambda(1,1,2),1,1,5).*cat(3,beta2,real(beta(1,1)),imag(beta(1,1)),zeros(nrows,ncols),ones(nrows,ncols));
% mata = repmat(vzx.Lambda(:,:,1),1,1,5).*cat(3,alpha2,real(vzx.alpha),imag(vzx.alpha),zeros(nrows,ncols),ones(nrows,ncols));
% matb = repmat(vzx.Lambda(:,:,2),1,1,5).*cat(3,beta2,real(vzx.beta),imag(vzx.beta),zeros(nrows,ncols),ones(nrows,ncols));
squeeze(mata)
squeeze(matb)
[0; 0; vzx.Lambda(1,1,3); 0; 0]
[vzx.Lambda(1,1,1)*(alpha2+1) vzx.Lambda(1,1,2)*(beta2+1) vzx.Lambda(1,1,3)]
squeeze(vzx.alpha(1,1,:))'

% analytical eigenvectors for X(ir,ic,:)
ir = 2, ic=1,
xi   = X(ir,ic,1);
rho2 = X(ir,ic,2)^2 + X(ir,ic,3)^2;
eta  = X(ir,ic,4);
zeta = X(ir,ic,5);
Delta = (zeta-xi)^2 + 4*rho2;
auxp = zeta - xi + sqrt(Delta);
auxm = zeta - xi - sqrt(Delta);
e1 = sqrt((auxp^2)/(auxp^2+4*rho2)) * [vzx.alpha(ir,ic);0;1]
e2 = sqrt((auxm^2)/(auxm^2+4*rho2)) * [vzx.beta(ir,ic);0;1]
e3 = [0;1;0]
squeeze(vzx.eigvals(ir,ic,:)).'
% numerical eigenvectors for X(ir,ic,:)
[V,D] = eig(wc2mat(X(ir,ic,:)))
% Xmat = wc2mat(X)

%%
X = nan(1,3,9);
% single bounce
rho = complex(1,0);
X(1,1,1) = rho*conj(rho);
X(1,1,2) = 0;
X(1,1,3) = 0;
X(1,1,4) = real(rho);
X(1,1,5) = imag(rho);
X(1,1,6) = 0;
X(1,1,7) = 0;
X(1,1,8) = 0;
X(1,1,9) = 1;
% double bounce
rho = complex(-1,0);
X(1,2,1) = rho*conj(rho);
X(1,2,2) = 0;
X(1,2,3) = 0;
X(1,2,4) = real(rho);
X(1,2,5) = imag(rho);
X(1,2,6) = 0;
X(1,2,7) = 0;
X(1,2,8) = 0;
X(1,2,9) = 1;
% double bounce
X(1,3,1) = 0;
X(1,3,2) = 0;
X(1,3,3) = 0;
X(1,3,4) = 0;
X(1,3,5) = 0;
X(1,3,6) = 2;
X(1,3,7) = 0;
X(1,3,8) = 0;
X(1,3,9) = 0;
vzX = van_zylx(X)
vzX.eigvals
%vzX.Lambda
imtool(vzX.eigvals(:,:,[2 3 1])) % magnify and zoom


%% for target detection: plot without logarithms
cd FSAR/Vejers
load Xband_Vejers.mat
cd ../..
% set(gca,'DataAspectRatio',[0.29 0.25 1])
imshowsar(HHHV,'HHHV');
imshowsar(HVVV,'HVVV');
imshowsar(HHVV,'HHVV');
[nr,nl] = size(HHHH);
figure('Name','RG power HHHV, HVVV')
imshowrgbprctile(cat(3,abs(HHHV).^2,abs(HVVV).^2,zeros(nr,nl)));
figure('Name','RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,abs(HHHV).^2,abs(HVVV).^2,abs(HHVV).^2));
SAR = cat(3,HHHH,real(HHHV),imag(HHHV),real(HHVV),imag(HHVV),HVHV,real(HVVV),imag(HVVV),VVVV);
% SAR = mult_sqrt2(SAR);
% clear HHHH HHHV HHVV HVHV HVVV VVVV
% clear nr nl
vzx = van_zylx(SAR)
vz = van_zyl(double(SAR(:,1601:1601+4800,:)))
vzq = polsar_eigvals_rgb(double(SAR(:,1601:1601+4800,:)))
figure('Name','van Zyl/Cloude decomposition')
imshowrgbprctile(vzx.Lambda(:,:,[2 3 1]));
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'DataAspectRatio',[0.30 0.25 1])
% figure('Name','van Zyl/Cloude decomposition')
% imshowrgbprctile(vz.Lambda,[2 3 1],[0.5 99.5]);
% figure('Name','van Zyl/Cloude decomposition [10*log10]')
% imshowrgb(10*log10(vz.Lambda(:,:,[2 3 1])));

% Residual from azim model (one obs only)
row = 1; col = 1;
SAR1 = double([SAR(row,col,1) complex(SAR(row,col,2),SAR(row,col,3)) complex(SAR(row,col,4),SAR(row,col,5));
       complex(SAR(row,col,2),-SAR(row,col,3)) SAR(row,col,6) complex(SAR(row,col,7),SAR(row,col,8))
       complex(SAR(row,col,4),-SAR(row,col,5)) complex(SAR(row,col,7),-SAR(row,col,8)) SAR(row,col,9)]);
fprintf('det of covmat = %g\n',det(SAR1))
res = SAR1 - ...
      (vzx.Lambda(row,col,1)*[abs(vzx.alpha(row,col)).^2 0 vzx.alpha(row,col);0 0 0;conj(vzx.alpha(row,col)) 0 1] + ...
       vzx.Lambda(row,col,2)*[abs(vzx.beta(row,col)).^2  0 vzx.beta(row,col) ;0 0 0;conj(vzx.beta(row,col)) 0 1] + ...
       vzx.Lambda(row,col,3)*[0 0 0;0 1 0;0 0 0])
fprintf('det of residual (azim) = %g\n',det(res))
% [V,D] = eig(SAR1);
% res1 = SAR1 - ...
%       (vzx.eigvals(row,col,1)*[abs(vzx.alpha(row,col)).^2 0 vzx.alpha(row,col);0 0 0;conj(vzx.alpha(row,col)) 0 1] + ...
%        vzx.eigvals(row,col,2)*[abs(vzx.beta(row,col)).^2  0 vzx.beta(row,col) ;0 0 0;conj(vzx.beta(row,col)) 0 1] + ...
%        vzx.eigvals(row,col,3)*[0 0 0;0 1 0;0 0 0])
% fprintf('det of residual (full pol) = %g\n',det(res1))

% Another visualization
figure('Name','van Zyl/Cloude')
beta = -vzx.alpha./abs(vzx.alpha).^2;
rr = vzx.Lambda(:,:,2).*(abs(beta).^2+1);
clear beta
gg = vzx.Lambda(:,:,3);
bb = vzx.Lambda(:,:,1).*(abs(vzx.alpha).^2+1);
imshowrgbprctile(cat(3,rr,gg,bb))
% imshowrgbprctile(10*log10(cat(3,rr,gg,bb)))

rrggbb = cat(3,rr,gg,bb);
rrggbb = 10*log10(cat(3,rr,gg,bb));
clear rr gg bb
minL = prctile(rrggbb(:),1);
maxL = prctile(rrggbb(:),99);
imshow((rrggbb-minL)/(maxL-minL))
clear rrggbb

[rho,rho0,theta,tau] = rho_llrr(SAR,true);

%% Useless: IH-plot of entropy and alpha
% V = vzx.entropy;
% amin = prctile(V(:),1);
% amax = prctile(V(:),99);
% V = -(V-amin)/(amax-amin);
% % H = (phas-min(phas(:)))/(max(phas(:)-min(phas(:))));
% % H = (phas-min(phas(:)))/12.0;
% H = abs(vzx.alpha);%.^2;
% minH = min(H(:));
% H = (H-minH)/(max(H(:))-minH);
% figure('Name','IH-plot entropy/abs(alpha)')
% imshow(hsv2rgb(cat(3,H,ones(nr,nl),V)))
% colormap(hsv(256))
% cb = colorbar('SouthOutside');
% set(cb, 'XTick', ((-pi:(pi/2):pi)+pi)/2/pi)
% set(cb, 'XTickLabel', {'-180','-90','0','90','180'})
% xlabel(cb, 'degrees')
% drawnow

%% % to Knut
% HHHH = HHHH(1376:1376+99,5480:5480+99);
% HHHV = HHHV(1376:1376+99,5480:5480+99);
% HHVV = HHVV(1376:1376+99,5480:5480+99);
% HVHV = HVHV(1376:1376+99,5480:5480+99);
% HVVV = HVVV(1376:1376+99,5480:5480+99);
% VVVV = VVVV(1376:1376+99,5480:5480+99);
% Amat = [HHHH(:),real(HHHV(:)),imag(HHHV(:)),real(HHVV(:)),imag(HHVV(:)),HVHV(:),real(HVVV(:)),imag(HVVV(:)),VVVV(:)];
% writematrix(Amat,'Vejers_misc.txt');

%%
cd radarsat2bonn/
F = dir('*.tif');
p = length(F);
nr = 1000;
nc = 1000;
SAR = nan(nr,nc,9,p);
time = nan(p,1);
h = waitbar(0,'Reading data...','name','Read SAR data');
for ii=1:p
    waitbar(ii/p,h);
    SAR(:,:,:,ii) = imread(F(ii).name);
    year = str2num(F(ii).name(1:4));
    month = str2num(F(ii).name(5:6));
    day = str2num(F(ii).name(7:8));
    datevec = [year,month,day,0,0,0];
    time(ii) = datenum(datestr(datevec));
end
close(h)
cd ..
vz = van_zylx(coh2cov(SAR))
nr = nr - 50;
for ix=1:p
    img(ix).rgb = 10*log10(vz(ix).Lambda(51:end,:,[2 3 1]));
%     img(ix).rgb = vz(ix).Lambda(51:end,:,[2 3 1]);
end
% same stretch for R, for G and for B
for ie=1:3
    ev = [];
    for ix=1:p
        ev = cat(3,ev,img(ix).rgb(:,:,ie));
    end
    mm = mean(reshape(ev,nr*nc,p));
    ss = std(reshape(ev,nr*nc,p));
    for ix=1:p
        img(ix).rgb(:,:,ie) = (img(ix).rgb(:,:,ie)-mm(ix))/ss(ix);
    end
end
clear ev
figure('Name','van Zyl/Cloude decomposition, same stretch')
nstd = 1.5;
for ix=1:p
    h(ix) = subplot_tight(3,4,ix,[0.005,0.005]);
    imshow(uint8(255*((img(ix).rgb+nstd)/(2*nstd))))
end
linkaxes(h,'xy')
% figure('name','van Zyl/Cloude decomposition')
% for ix=1:p
%     h(ix) = subplot_tight(3,4,ix,[0.005,0.005]);
% %     imshowrgb(10*log10(vz(ix).Lambda(51:end,:,[2 3 1])))
%     imshowrgbprctile(10*log10(vz(ix).Lambda(51:end,:,[2 3 1])),[1 2 3],[1 1 1;99 99 99])
% end
% linkaxes(h,'xy')

%%
cd Foulum
L62 = flipud(freadenvisar('fl062_l/fl062_l','b'));
L63 = flipud(freadenvisar('fl063_l/fl063_l','b'));
L64 = flipud(freadenvisar('fl064_l/fl064_l','b'));
L65 = flipud(freadenvisar('fl065_l/fl065_l','b'));
L68 = flipud(freadenvisar('fl068_l/fl068_l','b'));
L74 = flipud(freadenvisar('fl074_l/fl074_l','b'));
cd ..
figure('Name','L62 RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,L62(:,:,2).^2+L62(:,:,3).^2,L62(:,:,7).^2+L62(:,:,8).^2,L62(:,:,4).^2+L62(:,:,3).^5));
figure('Name','L63 RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,L63(:,:,2).^2+L63(:,:,3).^2,L63(:,:,7).^2+L63(:,:,8).^2,L63(:,:,4).^2+L63(:,:,3).^5));
figure('Name','L64 RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,L64(:,:,2).^2+L64(:,:,3).^2,L64(:,:,7).^2+L64(:,:,8).^2,L64(:,:,4).^2+L64(:,:,3).^5));
figure('Name','L65 RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,L65(:,:,2).^2+L65(:,:,3).^2,L65(:,:,7).^2+L65(:,:,8).^2,L65(:,:,4).^2+L65(:,:,3).^5));
figure('Name','L68 RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,L68(:,:,2).^2+L68(:,:,3).^2,L68(:,:,7).^2+L68(:,:,8).^2,L68(:,:,4).^2+L68(:,:,3).^5));
figure('Name','L74 RGB power HHHV, HVVV, HHVV')
imshowrgbprctile(cat(3,L74(:,:,2).^2+L74(:,:,3).^2,L74(:,:,7).^2+L74(:,:,8).^2,L74(:,:,4).^2+L74(:,:,3).^5));
SAR = cat(4,L62,L63,L64,L65,L68,L74);
clear L62 L63 L64 L65 L68 L74
[nr,nc,~,p] = size(SAR);
vz = van_zylx(SAR)
for ix=1:p
    img(ix).rgb = 10*log10(vz(ix).Lambda(:,:,[2 3 1]));
end
% same stretch for R, for G and for B
for ie=1:3
    ev = [];
    for ix=1:p
        ev = cat(3,ev,img(ix).rgb(:,:,ie));
    end
    mm = mean(reshape(ev,nr*nc,p));
    ss = std(reshape(ev,nr*nc,p));
    for ix=1:p
        img(ix).rgb(:,:,ie) = (img(ix).rgb(:,:,ie)-mm(ix))/ss(ix);
    end
end
clear ev
figure('Name','van Zyl/Cloude decomposition, same stretch')
nstd = 1.5;
for ix=1:p
    h(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshow(uint8(255*((img(ix).rgb+nstd)/(2*nstd))))
end
linkaxes(h,'xy')

figure('Name','van Zyl/Cloude decomposition')
for ix = 1:6
    h(ix) = subplot_tight(2,3,ix,[0.005,0.005]);
    imshowrgbprctile(vz(ix).Lambda(:,:,[2 3 1]))
%     imshowrgb(10*log10(vz(ix).Lambda(:,:,[2 3 1])))
end
linkaxes(h,'xy')

end % while false
