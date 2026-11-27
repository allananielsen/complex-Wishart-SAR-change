function wl = wilks_lambda_change(X,Y,L,varargin)

% Change detection in polarimetric SAR data bmo Wilks' lambda
%
% wl = wilks_lambda_change(X,Y,L);
%
% Change detection in bi-temporal polarimetric SAR images
% (multi-look covariance representation)
% by means of Wilks' lambda
%
% Perform test for equality of two complex 1x1, 2x2 or 3x3
% variance-covariance matrices representing
%
%   full polarimetry
%   full polarimetry, azimutal symmetry
%   full polarimetry, diagonal case
%   dual polarimetry
%   dual polarimetry, diagonal case
%   single polarization
%
% (to read the data from disk, see freadenvisar)
%
% ShhShh* ShhShv* ShhSvv*     k    a   rho
% ShvShh* ShvShv* ShvSvv*  =  a*   ksi b
% SvvShh* SvvShv* SvvSvv*     rho* b*  zeta
%
% Input
%   X       - covariance matrix input at first time point;
%   Y       - covariance matrix input at second time point
%             X and Y must be 3-D with either
%               9 bands (full polarimetry),
%               5 bands (full pol, azimuthal symmetry),
%               3 bands (full pol, diagonal only),
%               4 bands (dual polarimetry),
%               2 bands (dual pol, diagonal only), or
%               1 band (power data);
%             reading order is alphabetical in indexes with
%             real parts before imaginary parts (as read by
%             freadenvisar), fx. for full polarimetry:
%               hhhh, re(hhhv), im(hhhv), re(hhvv), im(hhvv),
%               hvhv, re(hvvv), im(hvvv), vvvv
%   L       - number of looks for X and Y
%   Ftest   - "true" if you want F-test (as well as chi-squared test)
%
% Output
%   struct wl with fields
%       Wilks(:,:,1) - |X|/|X+Y|
%       Wilks(:,:,2) - |Y|/|X+Y|
%       L            - number of looks for X and Y
%       df           - number of degrees of freedom for chi-squared statistic
%       chisq        - Bartlett's chi-squared statistic
%       pChisq       - p-value for Bartlett's chi-squared statistic
%       df1          - do for F-test
%       df2          -
%         F          -
%        pF          -

% (c) Copyright 2018
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, www.imm.dtu.dk/~alan
% 2-14 Sep 2018

if nargin<3
    help wilks_lambda_change
    error('wilks_lambda_change: wrong input');
end
Ftest = false;
if nargin>3
    Ftest = true;
end
if nargout>1
    help wilks_lambda_change
    error('wilks_lambda_change: wrong output');
end

if ischar(X), X = freadenvisar(X); end
if ischar(Y), Y = freadenvisar(Y); end

pX = size(X);
pY = size(Y);
if ~all(pX==pY)
    error('wilks_lambda_change: X and Y must have same dimensions');
end
nrows = pX(1);
ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    nvars = 1;
else
    error('wilks_lambda_change: wrong sizes of input images');
end

if nvars==9 % fully polarimetric data
    comment = 'full polarimetry, experimental vales for Beta-dist a and b';
    fprintf('wilks_lambda_change: %s, L = %i\n',comment,L)
    p = 3;
    % simulation based, L=3:197
    alpha = -0.75 + 0.5836*L;
    beta = 2 + 4.0825*L;
elseif nvars==5 % full polarimetry, azimuthal symmetry
    comment = '3x3, azimuthal symmetry, experimental vales for Beta-dist a and b';
    fprintf('wilks_lambda_change: %s, L = %i\n',comment,L)
    p = 3;
    alpha = -0.15523 + 0.52042*L;
    beta = 0.94751 + 3.6408*L;
%     alpha = -0.21781 + 0.56567*L;
%     beta = 0.90932 + 3.9579*L;
elseif nvars==3 % full polarimetry, diagonal only
    comment = '3x3, diagonal only, experimental vales for Beta-dist a and b';
    fprintf('wilks_lambda_change: %s, L = %i\n',comment,L)
    alpha = 0.083327 + 0.51612*L; % 1/12 + 0.51612*L;
    beta = 0.20834 + 3.6127*L;
elseif nvars==4 % dual polarization
    comment = 'dual polarization, experimental vales for Beta-dist a and b';
    fprintf('wilks_lambda_change: %s, L = %i\n',comment,L)
    p = 2;
%     alpha = L*3/4 - 0.5; % -0.48027 + 0.74996*L; % simulation based
%     beta = L*9/4 + 1/6; % 0.1597 + 2.2492*L;
    alpha = L*3/4 - 1/2; % -0.497604 + 0.74999*L; % simulation based
    beta = L*9/4; % 2.25012*L;
elseif nvars==2 % dual polarimetry, diagonal only
    comment = 'dual polarization, diagonal only';   
    fprintf('wilks_lambda_change: %s, L = %i\n',comment,L)
    p = 2;
    alpha = L*3/4; %  0.75064*L; % simulation based
    beta = L*9/4; % 2.2519*L;
elseif nvars==1 % one channel data only
    comment = 'single polarization';
    fprintf('wilks_lambda_change: %s, L = %i\n',comment,L)
    p = 1;
    alpha = L; % 0.9996*L; % simulation based
    beta = L; % 0.9996*L;
else
    error('wilks_lambda_change_change: input (9, 5, 3), (4, 2) or 1 band(s)');
end

det1 = wishart_det(X);
det2 = wishart_det(Y);
X = X+Y;
clear Y
det3 = wishart_det(X);
clear X

wl.Wilks = nan(nrows,ncols,2);
wl.Wilks(:,:,1) = det1./det3;
% wl.pseudoQ = det1.*det2./det3.^2;
% % wl.Q = (2^(2*p)*wl.pseudoQ).^L; % (2^(2*p)*det1.*det2./det3.^2).^L;
% wl.lnQ = L*(2*p*log(2)+log(det1)+log(det2)-2*log(det3));
clear det1
wl.Wilks(:,:,2) = det2./det3;
clear det2 det3

% No theory, simulation based
wl.pBeta = cdf('Beta',wl.Wilks,alpha,beta);
wl.a = alpha;
wl.b = beta;
wl.L = L;
wl.Mode = comment;
wl.pBeta_comment ='p-value pBeta is simulation based, no theory, two-sided';

while false
    
m = L;
n = L;
% wl.Lx = m;
% wl.Ly = n;
wl.df = n*p;
% Bartlett's chi-squared approximation
wl.chisq = ((p-n+1)/2-m)*log(wl.Wilks);
wl.pChisq = gammainc(0.5*wl.chisq,0.5*wl.df,'upper'); % P{> chisq} = P{<= Lambda}
% wl.pChisq = gammainc(0.5*wl.chisq,0.5*wl.df); % P{<= chisq} = P{> Lambda}
% wcRjl is not 'upper', and here we use -2*ln(Q)
wl.pQ = gammainc(-wl.lnQ,0.5*nvars); % P{<= -2*lnQ}
wl.dfQ = nvars;
wl.pChisq_comment = 'p-value pChisq is P{<= Lambda} = P{> chisq}, two-sided';
wl.pQ_comment = 'p-value pQ is for change: P{<= -2*lnQ}, one-sided';

if Ftest
    r = L; %??????????????
    m = p; %??????????????
    n = p; %??????????????
    if m^2+n^2==5 % must be checked, but it doesn't happen here
        t = 1;
        Ut = wl.Wilks;
    else
        t = sqrt((m^2*n^2-4)/(m^2+n^2-5));
        Ut = wl.Wilks.^(1/t);
    end
    v = 0.5*(2*r+n-m-1);
    wl.df1 = m*n;
    wl.df2 = v*t+1-0.5*wl.df1;
    wl.F = (1-Ut)./Ut*wl.df2/wl.df1;
    clear Ut
    wl.pF = betainc(wl.df2./(wl.df2+wl.df1*wl.F),0.5*wl.df2,0.5*wl.df1,'upper'); % P(>F)
end

end % while false

while false
    
% Radarsat-2 Bonn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % coherency matrix, T (strange values, so I divide by 1e8)
% T = T*1e-8;
% N = [1 0 1;1 0 -1;0 sqrt(2) 0]/sqrt(2);
% C = N'*T*N; % covariance matrix
% T = N*C*N'; % coherency matrix
% % det(C) - det(T) % is zero

cd radarsat2bonn
x1 = imread('20100731.tif');
x2 = imread('20100824.tif');
info = geotiffinfo('20100731.tif'); % R-2 data are poorly georectified here
% info = geotiffinfo('20100824.tif');
cd ..

x1(1:50,1:250,:) = nan(50,250,9);
x2(1:50,1:250,:) = nan(50,250,9);

% x1 = x1(51:end,:,:);
% x2 = x2(51:end,:,:);
% 
% x1 = x1(1:450,351:end,:);
% x2 = x2(1:450,351:end,:);

% lw = loewner_order(x1,x2)
% wc = wcRjl(cat(4,x1,x2),12,'full',false)
wl = wilks_lambda_change(x1,x2,12);

xx1 = 10*log10(x1(:,:,[6 9 1])) - 80; % dB? weird raw data
xx2 = 10*log10(x2(:,:,[6 9 1])) - 80;
xx1(:,:,2) = xx1(:,:,2) + 6; % add 6 to green
xx2(:,:,2) = xx2(:,:,2) + 6;
xx1 = (xx1 + 24)/18; % [-24,0] mapped to [0,1] for imshow
xx2 = (xx2 + 24)/18;

[nrows,ncols,~] = size(x1);
mv = 0.005; % mv = 0.02; % subplot_tight vertical margin 
mh = 0.005; % mh = 0.01; % subplot_tight horizontal margin 

pval = 0.99995; % pval = 0.9999;
figure
h(1,1) = subplot_tight(3,2,1,[mv,mh]);
imshow(xx1)
h(1,2) = subplot_tight(3,2,2,[mv,mh]);
imshow(xx2)
h(1,3) = subplot_tight(3,2,3,[mv,mh]);
imshow((cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols))))
h(1,4) = subplot_tight(3,2,4,[mv,mh]);
imshow(255*(cat(3,wl.pBeta(:,:,1)>pval,wl.pBeta(:,:,2)>pval,zeros(nrows,ncols))))
h(1,5) = subplot_tight(3,2,5,[mv,mh]);
imshow(255*(cat(3,wl.pBeta(:,:,2)<1-pval,wl.pBeta(:,:,1)<1-pval,zeros(nrows,ncols))))
h(1,6) = subplot_tight(3,2,6,[mv,mh]);
imshow(255*(cat(3,wl.pBeta(:,:,1)>pval & wl.pBeta(:,:,2)<1-pval,...
    wl.pBeta(:,:,2)>pval & wl.pBeta(:,:,1)<1-pval,zeros(nrows,ncols))))
linkaxes(h,'xy')

figure
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(wl.Wilks(:,:,1),'Normalization','pdf');
hold on
xx = 0:0.005:hh.BinLimits(2);
pd = makedist('Beta','a',wl.a,'b',wl.b);
% pd = makedist('Beta','a',5.0,'b',43.0);
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
histogram(wl.Wilks(:,:,2),'Normalization','pdf');
hold on
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
% clear h xx f
linkaxes(hp,'xy')

% beta for no-change region(s)
figure
data = squeeze(wl.Wilks(:,:,1));
data = data(roiNC);
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(data,'Normalization','pdf');
hold on
xx = 0:0.005:hh.BinLimits(2);
pd = makedist('Beta','a',wl.a,'b',wl.b);
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
data = squeeze(wl.Wilks(:,:,2));
data = data(roiNC);
hh = histogram(data,'Normalization','pdf');
hold on
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
% clear h xx f
linkaxes(hp,'xy')

%%%%%%%%

% Pauli RGB
xx1 = 10*log10(x1(:,:,[6 9 1])) - 80; % dB? weird raw data
xx2 = 10*log10(x2(:,:,[6 9 1])) - 80;
xx1(:,:,2) = xx1(:,:,2) + 6; % add 6 to green
xx2(:,:,2) = xx2(:,:,2) + 6;
xx1 = (xx1 + 24)/18; % [-24,0] mapped to [0,1] for imshow
xx2 = (xx2 + 24)/18;

pval = 0.9999;
figure
h(1) = subplot_tight(2,2,1,[0.02,0.01]);
imshow(xx1)
% title('VV t_1')
h(2) = subplot_tight(2,2,2,[0.02,0.01]);
imshow(xx2)
% title('VV t_2')
h(3) = subplot_tight(2,2,3,[0.02,0.01]);
imshow((cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols))))
h(4) = subplot_tight(2,2,4,[0.02,0.01]);
imshow(255*(cat(3,wl.pChisq(:,:,1)>pval,wl.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
linkaxes(h,'xy')

% follow by "File | Open" in Google Earth Pro
outim = uint8(255*cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols)));
geotiffwrite('Wilks_change_Bonn',outim,info.SpatialRef,...
    'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)
  
% Sentinel-1 Fraport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd fraportfloat
x1 = imread('sentinel1_VVVH_02_warp.tif');
x2 = imread('sentinel1_VVVH_03_warp.tif');
info = geotiffinfo('sentinel1_VVVH_02_warp.tif');
cd ..

xx1 = 10*log10(cat(3,x1,x1(:,:,1)./x1(:,:,2)));
xx2 = 10*log10(cat(3,x2,x2(:,:,1)./x2(:,:,2)));

xx1(:,:,1) = (xx1(:,:,1)+22)/26; % -22 dB to  4 dB
xx1(:,:,2) = (xx1(:,:,2)+28)/24; % -28 dB to -4 dB
xx1(:,:,3) = (xx1(:,:,3))/16;    %   0 db to 16 dB
xx2(:,:,1) = (xx2(:,:,1)+22)/26;
xx2(:,:,2) = (xx2(:,:,2)+28)/24;
xx2(:,:,3) = (xx2(:,:,3))/16;

wl = wilks_lambda_change(x1,x2,4.9);

[nrows,ncols,~] = size(x1);
mv = 0.005; % mv = 0.02; % subplot_tight vertical margin 
mh = 0.005; % mh = 0.01; % subplot_tight horizontal margin 

pval = 0.99995; % pval = 0.9999;
figure%('Name','RGB VV-VH-VH/VV')
h(1,1) = subplot_tight(3,2,1,[mv,mh]);
imshow(xx1)
% imshowrgb(xx1,[1 2 3])
h(1,2) = subplot_tight(3,2,2,[mv,mh]);
imshow(xx2)
% imshowrgb(xx2,[1 2 3])
h(1,3) = subplot_tight(3,2,3,[mv,mh]);
imshow((cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols))))
h(1,4) = subplot_tight(3,2,4,[mv,mh]);
imshow(255*(cat(3,wl.pBeta(:,:,1)>pval,wl.pBeta(:,:,2)>pval,zeros(nrows,ncols))))
%h(1,5) = subplot_tight(3,2,5,[mv,mh]);
%imshow((cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols))))
h(1,5) = subplot_tight(3,2,5,[mv,mh]);
imshow(255*(cat(3,wl.pBeta(:,:,2)<1-pval,wl.pBeta(:,:,1)<1-pval,zeros(nrows,ncols))))
h(1,6) = subplot_tight(3,2,6,[mv,mh]);
imshow(255*(cat(3,wl.pBeta(:,:,1)>pval & wl.pBeta(:,:,2)<1-pval,...
    wl.pBeta(:,:,2)>pval & wl.pBeta(:,:,1)<1-pval,zeros(nrows,ncols))))
linkaxes(h,'xy')

% follow by "File | Open" in Google Earth Pro
outim = uint8(255*(cat(3,wl.pBeta(:,:,1)>pval & wl.pBeta(:,:,2)<1-pval,...
    wl.pBeta(:,:,2)>pval & wl.pBeta(:,:,1)<1-pval,zeros(nrows,ncols))));
geotiffwrite('wl_change',outim,info.SpatialRef,...
    'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)

figure
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(wl.Wilks(:,:,1),'Normalization','pdf');
hold on
xx = 0:0.005:hh.BinLimits(2);
pd = makedist('Beta','a',wl.a,'b',wl.b);
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
histogram(wl.Wilks(:,:,2),'Normalization','pdf');
hold on
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
% clear h xx f
linkaxes(hp,'xy')

% beta for no-change region(s)
figure
data = squeeze(wl.Wilks(:,:,1));
data = data(roiNC);
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(data,'Normalization','pdf');
hold on
xx = 0:0.005:hh.BinLimits(2);
pd = makedist('Beta','a',wl.a,'b',wl.b);
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
data = squeeze(wl.Wilks(:,:,2));
data = data(roiNC);
hh = histogram(data,'Normalization','pdf');
hold on
plot(xx,pdf(pd,xx),'LineWidth',1.5,'Color','Red')
grid on
% clear h xx f
linkaxes(hp,'xy')

% Kolmogorov-Smirnov test
aa = hh.BinEdges + 0.01;
aa = aa(1:end-1);
[h,p,ks2stat] = kstest2(hh.Values,pdf(pd,aa))

%%%%%%%%

pval = 0.99995; % pval = 0.9999;
figure
h(1,1) = subplot_tight(2,2,1,[0.02,0.01]);
imshow(10*log10(x1(:,:,1)),[-24 6])
% 1,title('VV t_1')
h(1,2) = subplot_tight(2,2,2,[0.02,0.01]);
imshow(10*log10(x2(:,:,1)),[-24 6])
% title('VV t_2')
h(1,3) = subplot_tight(2,2,3,[0.02,0.01]);
imshow((cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols))))
h(1,4) = subplot_tight(2,2,4,[0.02,0.01]);
imshow(255*(cat(3,wl.pChisq(:,:,1)>pval,wl.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
linkaxes(h,'xy')

pval = 0.00005;
figure('Name','Wilks''s Lambda: P{<= Lambda_X}>0.99995 | P{<= Lambda_Y}<=0.00005 (R), vv (G)')
h(2,1) = subplot_tight(2,2,1,[0.02,0.01]);
imshow(10*log10(x1(:,:,1)),[-24 6])
h(2,2) = subplot_tight(2,2,2,[0.02,0.01]);
imshow(10*log10(x2(:,:,1)),[-24 6])
h(2,3) = subplot_tight(2,2,3,[0.02,0.01]);
imshow((cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols))))
h(2,4) = subplot_tight(2,2,4,[0.02,0.01]);
imshow(cat(3,or(wl.pChisq(:,:,1)>1-pval,wl.pChisq(:,:,2)<=pval),or(wl.pChisq(:,:,2)>1-pval,wl.pChisq(:,:,1)<=pval),zeros(nrows,ncols)))
linkaxes(h,'xy')

% % follow by "File | Open" in Google Earth Pro
% outim = uint8(255*cat(3,wl.Wilks(:,:,1),wl.Wilks(:,:,2),zeros(nrows,ncols)));
% geotiffwrite('wl_change',outim,info.SpatialRef,...
%     'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)

% EMISAR L-band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Foulum
L62 = flipud(freadenvisar('fl062_l/fl062_l','b'));
L63 = flipud(freadenvisar('fl063_l/fl063_l','b'));
L64 = flipud(freadenvisar('fl064_l/fl064_l','b'));
L65 = flipud(freadenvisar('fl065_l/fl065_l','b'));
L68 = flipud(freadenvisar('fl068_l/fl068_l','b'));
L74 = flipud(freadenvisar('fl074_l/fl074_l','b'));
cd ..

[nrows,ncols,~] = size(L65);

wl1 = wilks_lambda_change(L62,L63,13);
wl2 = wilks_lambda_change(L63,L64,13);
wl3 = wilks_lambda_change(L64,L65,13);
wl4 = wilks_lambda_change(L65,L68,13);
wl5 = wilks_lambda_change(L68,L74,13);

figure('Name','Pauli representation')
h(1) = subplot_tight(2,3,1,[0.02,0.01]);
imshowsarpauli(L62,[1 2 3],[-24 -30 -20;-3 -3 0])
h(2) = subplot_tight(2,3,2,[0.02,0.01]);
imshowsarpauli(L63,[1 2 3],[-24 -30 -20;-3 -3 0])
h(3) = subplot_tight(2,3,3,[0.02,0.01]);
imshowsarpauli(L64,[1 2 3],[-24 -30 -20;-3 -3 0])
h(4) = subplot_tight(2,3,4,[0.02,0.01]);
imshowsarpauli(L65,[1 2 3],[-24 -30 -20;-3 -3 0])
h(5) = subplot_tight(2,3,5,[0.02,0.01]);
imshowsarpauli(L68,[1 2 3],[-24 -30 -20;-3 -3 0])
h(6) = subplot_tight(2,3,6,[0.02,0.01]);
imshowsarpauli(L74,[1 2 3],[-24 -30 -20;-3 -3 0])
linkaxes(h,'xy')
clear h

% figure
% h(1) = subplot_tight(2,3,1,[0.005,0.005]);
% imshowsarpauli(L62)
% h(2) = subplot_tight(2,3,2,[0.005,0.005]);
% imshowsarpauli(L63)
% h(3) = subplot_tight(2,3,3,[0.005,0.005]);
% imshowsarpauli(L64)
% h(4) = subplot_tight(2,3,4,[0.005,0.005]);
% imshowsarpauli(L65)
% h(5) = subplot_tight(2,3,5,[0.005,0.005]);
% imshowsarpauli(L68)
% h(6) = subplot_tight(2,3,6,[0.005,0.005]);
% imshowsarpauli(L74)
% linkaxes(h,'xy')

% figure
% h(1) = subplot_tight(2,3,1,[0.02,0.01]);
% imshowrgb(10*log10(L62(:,:,[6 1 9])))
% h(2) = subplot_tight(2,3,2,[0.02,0.01]);
% imshowrgb(10*log10(L63(:,:,[6 1 9])))
% h(3) = subplot_tight(2,3,3,[0.02,0.01]);
% imshowrgb(10*log10(L64(:,:,[6 1 9])))
% h(4) = subplot_tight(2,3,4,[0.02,0.01]);
% imshowrgb(10*log10(L65(:,:,[6 1 9])))
% h(5) = subplot_tight(2,3,5,[0.02,0.01]);
% imshowrgb(10*log10(L68(:,:,[6 1 9])))
% h(6) = subplot_tight(2,3,6,[0.02,0.01]);
% imshowrgb(10*log10(L74(:,:,[6 1 9])))
% linkaxes(h,'xy')

figure('Name','Wilks''s Lambda: |X|/|X+Y| (R), |Y|/|X+Y| (G)')
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow((cat(3,wl1.Wilks(:,:,1),wl1.Wilks(:,:,2),zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow((cat(3,wl2.Wilks(:,:,1),wl2.Wilks(:,:,2),zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow((cat(3,wl3.Wilks(:,:,1),wl3.Wilks(:,:,2),zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow((cat(3,wl4.Wilks(:,:,1),wl4.Wilks(:,:,2),zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow((cat(3,wl5.Wilks(:,:,1),wl5.Wilks(:,:,2),zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')
clear h1 h2 h3 h4 h5

% X only
pval = 0.9999; % pval = 0.99995;
figure
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow(255*(cat(3,wl1.pChisq(:,:,1)>pval,wl1.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow(255*(cat(3,wl2.pChisq(:,:,1)>pval,wl2.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow(255*(cat(3,wl3.pChisq(:,:,1)>pval,wl3.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow(255*(cat(3,wl4.pChisq(:,:,1)>pval,wl4.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow(255*(cat(3,wl5.pChisq(:,:,1)>pval,wl5.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')

% % either 3x2
% pval = 0.9999;
% figure
% h1 = subplot_tight(3,2,1,[0.02,0.01]);
% imshow(255*(cat(3,wl1.pChisq(:,:,1)>pval,wl1.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h2 = subplot_tight(3,2,2,[0.02,0.01]);
% imshow(255*(cat(3,wl2.pChisq(:,:,1)>pval,wl2.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h3 = subplot_tight(3,2,3,[0.02,0.01]);
% imshow(255*(cat(3,wl3.pChisq(:,:,1)>pval,wl3.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h4 = subplot_tight(3,2,4,[0.02,0.01]);
% imshow(255*(cat(3,wl4.pChisq(:,:,1)>pval,wl4.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h5 = subplot_tight(3,2,5,[0.02,0.01]);
% imshow(255*(cat(3,wl5.pChisq(:,:,1)>pval,wl5.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% linkaxes([h1,h2,h3,h4,h5],'xy')

% or 2x3
pval = 0.9999; % pval = 0.99995;
figure('Name','Wilks''s Lambda: P{<= Lambda}>0.99995')
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow((cat(3,wl1.pChisq(:,:,1)>pval,wl1.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow((cat(3,wl2.pChisq(:,:,1)>pval,wl2.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow((cat(3,wl3.pChisq(:,:,1)>pval,wl3.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow((cat(3,wl4.pChisq(:,:,1)>pval,wl4.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow((cat(3,wl5.pChisq(:,:,1)>pval,wl5.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')

figure
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow((cat(3,wl1.pChisq(:,:,1),wl1.pChisq(:,:,2),zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow((cat(3,wl2.pChisq(:,:,1),wl2.pChisq(:,:,2),zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow((cat(3,wl3.pChisq(:,:,1),wl3.pChisq(:,:,2),zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow((cat(3,wl4.pChisq(:,:,1),wl4.pChisq(:,:,2),zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow((cat(3,wl5.pChisq(:,:,1),wl5.pChisq(:,:,2),zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')


% EMISAR C-band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Foulum
C62 = flipud(freadenvisar('fl062_c/fl062_c','b'));
% C63 = flipud(freadenvisar('fl063_c/fl063_c','b'));
C64 = flipud(freadenvisar('fl064_c/fl064_c','b'));
C65 = flipud(freadenvisar('fl065_c/fl065_c','b'));
C68 = flipud(freadenvisar('fl068_c/fl068_c','b'));
cd ..

[nrows,ncols,~] = size(C65);

wl1 = wilks_lambda_change(C62,C64,13);
wl2 = wilks_lambda_change(C64,C65,13);
wl3 = wilks_lambda_change(C65,C68,13);

figure
h(1) = subplot_tight(2,2,1,[0.02,0.01]);
imshowsarpauli(C62,[1 2 3],[-15 -18 -12;-3 -3 0])
h(2) = subplot_tight(2,2,2,[0.02,0.01]);
imshowsarpauli(C64,[1 2 3],[-15 -18 -12;-3 -3 0])
h(3) = subplot_tight(2,2,3,[0.02,0.01]);
imshowsarpauli(C65,[1 2 3],[-15 -18 -12;-3 -3 0])
h(4) = subplot_tight(2,2,4,[0.02,0.01]);
imshowsarpauli(C68,[1 2 3],[-15 -18 -12;-3 -3 0])
linkaxes(h,'xy')

figure
h(1) = subplot_tight(2,2,1,[0.005,0.005]);
imshowsarpauli(C62(513:end,513:end,:))
h(2) = subplot_tight(2,2,2,[0.005,0.005]);
imshowsarpauli(C64(513:end,513:end,:))
h(3) = subplot_tight(2,2,3,[0.005,0.005]);
imshowsarpauli(C65(513:end,513:end,:))
h(4) = subplot_tight(2,2,4,[0.005,0.005]);
imshowsarpauli(C68(513:end,513:end,:))
linkaxes(h,'xy')

% figure
% h(1) = subplot_tight(2,3,1,[0.02,0.01]);
% imshowrgb(10*log10(L62(:,:,[6 1 9])))
% h(2) = subplot_tight(2,3,2,[0.02,0.01]);
% imshowrgb(10*log10(L63(:,:,[6 1 9])))
% h(3) = subplot_tight(2,3,3,[0.02,0.01]);
% imshowrgb(10*log10(L64(:,:,[6 1 9])))
% h(4) = subplot_tight(2,3,4,[0.02,0.01]);
% imshowrgb(10*log10(L65(:,:,[6 1 9])))
% h(5) = subplot_tight(2,3,5,[0.02,0.01]);
% imshowrgb(10*log10(L68(:,:,[6 1 9])))
% h(6) = subplot_tight(2,3,6,[0.02,0.01]);
% imshowrgb(10*log10(L74(:,:,[6 1 9])))
% linkaxes(h,'xy')

% X only
pval = 0.9999; % pval = 0.99995;
figure
h1 = subplot_tight(2,2,1,[0.02,0.01]);
imshow(255*(cat(3,wl1.pChisq(:,:,1)>pval,wl1.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
h2 = subplot_tight(2,2,2,[0.02,0.01]);
imshow(255*(cat(3,wl2.pChisq(:,:,1)>pval,wl2.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
h3 = subplot_tight(2,2,3,[0.02,0.01]);
imshow(255*(cat(3,wl3.pChisq(:,:,1)>pval,wl3.pChisq(:,:,1)<=1-pval,zeros(nrows,ncols))))
linkaxes([h1,h2,h3],'xy')

% % either 3x2
% pval = 0.9999;
% figure
% h1 = subplot_tight(3,2,1,[0.02,0.01]);
% imshow(255*(cat(3,wl1.pChisq(:,:,1)>pval,wl1.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h2 = subplot_tight(3,2,2,[0.02,0.01]);
% imshow(255*(cat(3,wl2.pChisq(:,:,1)>pval,wl2.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h3 = subplot_tight(3,2,3,[0.02,0.01]);
% imshow(255*(cat(3,wl3.pChisq(:,:,1)>pval,wl3.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h4 = subplot_tight(3,2,4,[0.02,0.01]);
% imshow(255*(cat(3,wl4.pChisq(:,:,1)>pval,wl4.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% h5 = subplot_tight(3,2,5,[0.02,0.01]);
% imshow(255*(cat(3,wl5.pChisq(:,:,1)>pval,wl5.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
% linkaxes([h1,h2,h3,h4,h5],'xy')

% or 2x3
pval = 0.9999; %pval = 0.99995
figure
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow(255*(cat(3,wl1.pChisq(:,:,1)>pval,wl1.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow(255*(cat(3,wl2.pChisq(:,:,1)>pval,wl2.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow(255*(cat(3,wl3.pChisq(:,:,1)>pval,wl3.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow(255*(cat(3,wl4.pChisq(:,:,1)>pval,wl4.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow(255*(cat(3,wl5.pChisq(:,:,1)>pval,wl5.pChisq(:,:,2)>pval,zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')

figure
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow((cat(3,wl1.Wilks(:,:,1),wl1.Wilks(:,:,2),zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow((cat(3,wl2.Wilks(:,:,1),wl2.Wilks(:,:,2),zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow((cat(3,wl3.Wilks(:,:,1),wl3.Wilks(:,:,2),zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow((cat(3,wl4.Wilks(:,:,1),wl4.Wilks(:,:,2),zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow((cat(3,wl5.Wilks(:,:,1),wl5.Wilks(:,:,2),zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(wl.chisq(:,:,1),'Normalization','pdf');
hold on
xx = 0:1:hh.BinLimits(2);
% f = chi2pdf(xx,wl.df);
f = pdf('Chisquare',xx,wl.df);
plot(xx,f,'LineWidth',1.5,'Color','Red')
grid on
%title(strcat('Bartlett''s \chi^2(',num2str(wl.df),') approximation'))
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
histogram(wl.chisq(:,:,2),'Normalization','pdf');
hold on
plot(xx,f,'LineWidth',1.5,'Color','Red')
grid on
% clear h xx f
%title(strcat('Bartlett''s \chi^2(',num2str(wl.df),') approximation'))
linkaxes(hp,'xy')

% chi-squared for no-change region(s)

figure
data = squeeze(wl.chisq(:,:,1));
data = data(roiNC2);
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(data,'Normalization','pdf');
hold on
xx = 0:0.2:hh.BinLimits(2);
% f = chi2pdf(xx,wl.df);
f = pdf('Chisquare',xx,wl.df);
plot(xx,f,'LineWidth',1.5,'Color','Red')
%title(strcat('Bartlett''s \chi^2(',num2str(wl.df),') approximation'))
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
data = squeeze(wl.chisq(:,:,2));
data = data(roiNC2);
hh = histogram(data,'Normalization','pdf');
hold on
% xx = 0:0.2:hh.BinLimits(2);
% % f = chi2pdf(xx,wl.df);
% f = pdf('Chisquare',xx,wl.df);
plot(xx,f,'LineWidth',1.5,'Color','Red')
% clear h xx f
%title(strcat('Bartlett''s \chi^2(',num2str(wl.df),') approximation'))
linkaxes(hp,'xy')


% F for no-change region(s) - there is something wrong with the F-values

figure
data = squeeze(wl.F(:,:,1));
data = data(roiNC);
hp(1) = subplot_tight(2,1,1,[0.05,0.05]);
hh = histogram(data,'Normalization','pdf');
hold on
xx = 0:0.01:hh.BinLimits(2);
f = pdf('F',xx,wl.df1,wl.df2);
plot(xx,f,'LineWidth',1.5,'Color','Red')
% title(strcat('Rao''s F(',num2str(wl.df1),',',num2str(wl.df2),') approximation'))
hp(2) = subplot_tight(2,1,2,[0.05,0.05]);
data = squeeze(wl.F(:,:,2));
data = data(roiNC);
histogram(data,'Normalization','pdf');
hold on
plot(xx,f,'LineWidth',1.5,'Color','Red')
% clear h xx f
% title(strcat('Rao''s F(',num2str(wl.df1),',',num2str(wl.df2),') approximation'))
linkaxes(h,'xy')


end % false



while false
if nvars==9 % fully polarimetric data
    fprintf('wilks_lambda_change: full pol\n')    
    p = 3;
    % ShhShh* ShhShv* ShhSvv*    k     a    rho   | 1   2,3     4,5
    % ShvShh* ShvShv* ShvSvv* =  a*    ksi  b     |     6       7,8
    % SvvShh* SvvShv* SvvSvv*    rho*  b*   zeta  |             9
    % 
    % det = k*ksi*zeta+a*b*rho*+rho*a**b*-rho^2*ksi-b^2*k-zeta*a^2
    %
    % a*b*rho* = a**b**rho = 2*Re(a*b*rho*)
    det1 = X(:,:,1).*X(:,:,6).*X(:,:,9) ...
        + 2*(X(:,:,4).*(X(:,:,2).*X(:,:,7)-X(:,:,3).*X(:,:,8)) ...
           + X(:,:,5).*(X(:,:,3).*X(:,:,7)+X(:,:,2).*X(:,:,8))) ...
        - X(:,:,6).*(X(:,:,4).^2+X(:,:,5).^2) ...
        - X(:,:,1).*(X(:,:,7).^2+X(:,:,8).^2) ...
        - X(:,:,9).*(X(:,:,2).^2+X(:,:,3).^2);
    det2 = Y(:,:,1).*Y(:,:,6).*Y(:,:,9) ...
        + 2*(Y(:,:,4).*(Y(:,:,2).*Y(:,:,7)-Y(:,:,3).*Y(:,:,8)) ...
           + Y(:,:,5).*(Y(:,:,3).*Y(:,:,7)+Y(:,:,2).*Y(:,:,8))) ...
        - Y(:,:,6).*(Y(:,:,4).^2+Y(:,:,5).^2) ...
        - Y(:,:,1).*(Y(:,:,7).^2+Y(:,:,8).^2) ...
        - Y(:,:,9).*(Y(:,:,2).^2+Y(:,:,3).^2);
    X = X + Y;
    clear Y
    det3 = X(:,:,1).*X(:,:,6).*X(:,:,9) ...
        + 2*(X(:,:,4).*(X(:,:,2).*X(:,:,7)-X(:,:,3).*X(:,:,8)) ...
           + X(:,:,5).*(X(:,:,3).*X(:,:,7)+X(:,:,2).*X(:,:,8))) ...
        - X(:,:,6).*(X(:,:,4).^2+X(:,:,5).^2) ...
        - X(:,:,1).*(X(:,:,7).^2+X(:,:,8).^2) ...
        - X(:,:,9).*(X(:,:,2).^2+X(:,:,3).^2);
elseif nvars==5 % full polarimetry, azimuthal symmetry
    fprintf('wilks_lambda_change: full pol, azimuthal symmetry\n')    
    p = 3;
    % ShhShh* ShhShv* ShhSvv*    k     0    rho  |  1           2,3
    % ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     4
    % SvvShh* SvvShv* SvvSvv*    rho*  0    zeta |              5
    % 
    % det = ksi*(k*zeta-rho^2)
    det1 = X(:,:,4).*(X(:,:,1).*X(:,:,5) - X(:,:,2).^2 - X(:,:,3).^2);
    det2 = Y(:,:,4).*(Y(:,:,1).*Y(:,:,5) - Y(:,:,2).^2 - Y(:,:,3).^2);
    X = X + Y;
    clear Y
    det3 = X(:,:,4).*(X(:,:,1).*X(:,:,5) - X(:,:,2).^2 - X(:,:,3).^2);
elseif nvars==3 % full polarimetry, diagonal only
    fprintf('wilks_lambda_change: full pol, diagonal only\n')    
    p = 3;
    % ShhShh* ShhShv* ShhSvv*    k     0    0    |  1
    % ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     2
    % SvvShh* SvvShv* SvvSvv*    0     0    zeta |        3
    % 
    % det = k*ksi*zeta
    det1 = X(:,:,1).*X(:,:,2).*X(:,:,3);
    det2 = Y(:,:,1).*Y(:,:,2).*Y(:,:,3);
    X = X + Y;
    clear Y
    det3 = X(:,:,1).*X(:,:,2).*X(:,:,3);
elseif nvars==4 % dual polarimetry
    fprintf('wilks_lambda_change: dual pol\n')    
    p = 2;
    % ShhShh* ShhSvv*    k     rho   | 1   2,3
    % SvvShh* SvvSvv*    rho*  ksi   |     4
    %
    % or
    %
    % ShhShh* ShhShv*    k     a     | 1   2,3
    % ShvShh* ShvShv*    a*    zeta  |     4
    %
    % or
    %
    % ShvShv* ShvSvv*    ksi   b     | 1   2,3
    % SvvShv* SvvSvv*    b*    zeta  |     4
    %
    % det  = k*ksi-rho2;
    det1 = X(:,:,1).*X(:,:,4) - X(:,:,2).^2 - X(:,:,3).^2;
    det2 = Y(:,:,1).*Y(:,:,4) - Y(:,:,2).^2 - Y(:,:,3).^2;
    X = X + Y;
    clear Y
    det3 = X(:,:,1).*X(:,:,4) - X(:,:,2).^2 - X(:,:,3).^2;
elseif nvars==2 % dual polarimetry, diagonal only
    fprintf('wilks_lambda_change: dual pol, diagonal only\n')    
    p = 2;
    % ShhShh* ShhSvv*    k     0     |  1
    % SvvShh* SvvSvv*    0     zeta  |     2
    %
    % det  = k*zeta
    det1 = X(:,:,1).*X(:,:,2);
    det2 = Y(:,:,1).*Y(:,:,2);
    X = X + Y;
    clear Y
    det3 = X(:,:,1).*X(:,:,2);
elseif nvars==1 % one channel data only
    fprintf('wilks_lambda_change: single pol\n')    
    p = 1;
    det1 = X(:,:,1);
    det2 = Y(:,:,1);
    X = X + Y;
    clear Y
    det3 = X(:,:,1);
else
    error('wilks_lambda_change_change: input (9, 5, 3), (4, 2) or 1 band(s), real (and complex)');
end
end % while false
