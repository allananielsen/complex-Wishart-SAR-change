function [rho,phd] = rho_hhvv(X)

% rho_hhvv calculates complex correlation between SHH and SVV in polSAR data
%
% rho = rho_hhvv(X);
% or
% [rho,phd] = rho_hhvv(X);
%
%   full polarimetry
%   full polarimetry, azimutal symmetry
%   dual polarimetry
%
% (to read the data from disk, see freadenvisar)
%
% ShhShh* ShhShv* ShhSvv*     k    a   rho
% ShvShh* ShvShv* ShvSvv*  =  a*   ksi b
% SvvShh* SvvShv* SvvSvv*     rho* b*  zeta
%
% Input
%   X   - covariance matrix representation, must be 3-D with
%             9 (quad/full polarimetry) bands,
%             order is alphabetical in indexes with
%             real parts before imaginary parts (as read by
%             freadenvisar), fx. for full polarimetry:
%               hhhh, re(hhhv), im(hhhv), re(hhvv), im(hhvv),
%               hvhv, re(hvvv), im(hvvv), vvvv
% Output
%   rho    - (complex) correlation between Shh and Svv
%   phd    - phase difference between Shh and Svv

% (c) Copyright 2014-2018
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, https://people.compute.dtu.dk/alan
% 2 Oct 2014 - 22 Sep 2018

if ischar(X), X = freadenvisar(X); end

pX = size(X);
nrows = pX(1);
ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
else
    error('rho_hhvv: wrong dimensionality of input');
end

if nvars==9 % fully polarimetric data
%     X(:,:,6) = 2*X(:,:,6);
%     X(:,:,[2 3 7 8]) = sqrt(2)*X(:,:,[2 3 7 8]);
%     Y(:,:,6) = 2*Y(:,:,6);
%     Y(:,:,[2 3 7 8]) = sqrt(2)*Y(:,:,[2 3 7 8]);
    fprintf('rho_hhvv: full polarimetry\n');
    % ShhShh* ShhShv* ShhSvv*    k     a    rho   | 1   2,3     4,5
    % ShvShh* ShvShv* ShvSvv* =  a*    ksi  b     |     6       7,8
    % SvvShh* SvvShv* SvvSvv*    rho*  b*   zeta  |             9
    % 
    rho = complex(X(:,:,4),X(:,:,5))./sqrt((X(:,:,1).*X(:,:,9)));
    if nargout>1
        rho = abs(rho);
        phd = atan2(X(:,:,5),X(:,:,4));
    end
elseif nvars==5 % full polarimetry, azimuthal symmetry
%     X(:,:,4) = 2*X(:,:,4);
%     Y(:,:,4) = 2*Y(:,:,4);
    fprintf('rho_hhvv: full polarimetry, azimuthal symmetry\n');
    % ShhShh* ShhShv* ShhSvv*    k     0    rho  |  1           2,3
    % ShvShh* ShvShv* ShvSvv* =  0     ksi  0    |     4
    % SvvShh* SvvShv* SvvSvv*    rho*  0    zeta |              5
    % 
    rho = complex(X(:,:,2),X(:,:,3))./sqrt((X(:,:,1).*X(:,:,5)));
    if nargout>1
        rho = abs(rho);
        phd = atan2(X(:,:,3),X(:,:,2));
    end
elseif nvars==4 % dual polarimetry
    fprintf('rho_hhvv: dual polarimetry\n');
    % ShhShh* ShhSvv*     k     rho   | 1   2,3
    % SvvShh* SvvSvv*  =  rho*  ksi   |     4
    %
    rho = complex(X(:,:,2),X(:,:,3))./sqrt((X(:,:,1).*X(:,:,4)));
    if nargout>1
        rho = abs(rho);
        phd = atan2(X(:,:,3),X(:,:,2));
    end
else
    error('rho_hhvv: input 9, 5, or 4 bands, real (and complex)');
%     fprintf('%s\n', '*** order of input images is ***very*** important ***')
%     fprintf('%s\n', '*** hhhh (real), hhhv (complex), hhvv (complex),  ***')
%     fprintf('%s\n', '*** hvhv (real), hvvv (complex), vvvv (real)      ***')
end

return

% plot phase
figure
if isreal(phd)
    imshow(phd,[-pi pi])
else
    imshow(angle(phd),[-pi pi])
end
hsvaa = hsv(256);
% start in blue, go spectral
% hsvaa = flip([hsvaa(172:end,:);hsvaa(1:171,:)]);
colormap(gca,hsvaa)
cb = colorbar('SouthOutside');
set(cb,'XTick',-pi:pi/2:pi)
set(cb,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})

% or

figure
if isreal(phd)
    imshow(phd*180/pi,[-180 180])
else
    imshow(angle(phd)*180/pi,[-180 180])
end
hsvaa = hsv(256);
% start in blue, go spectral
% hsvaa = flip([hsvaa(172:end,:);hsvaa(1:171,:)]);
colormap(gca,hsvaa)
cb = colorbar('SouthOutside');
set(cb,'XTick',-180:90:180)
set(cb,'XTickLabel',{'-180','-90','0','90','180'})

% or plot IH version of magnitude and phase

[nrows,ncols] = size(phd);
amp = abs(phd);
% amp = phd.*conj(phd);
phas = angle(phd);
% intensity-hue/hue-value plot of amplitude and phase
figure('Name','rho_hhvv'), imshow(hsv2rgb(cat(3,(phas+pi)/(2*pi),ones(nrows,ncols),amp)))
% % hue-saturation plot of amplitude and phase
% figure, imshow(hsv2rgb(cat(3,(phas+pi)/(2*pi),amp,ones(nrows,ncols))))
colormap(gca,hsv(256)) % not needed for hsv2rgb but for colorbar
cb = colorbar('SouthOutside');
set(cb,'XTick',(0:64:256)/256)
set(cb,'XTickLabel',{'-180','-90','0','90','180'})
drawnow
