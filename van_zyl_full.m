function vzf = van_zyl_full(X)

% Calculates Cloude/van Zyl for full PolSAR + eigvals, entropy and anisotropy
%
% SAR input images X are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%  ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   X       - 3-D array of 3x3 complex covariance matrix
%             (polarimetric SAR image in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%                        ShvShv*                 Re(ShvSvv*) Im(ShvSvv*)
%                                                SvvSvv*]
%
% Output
%   vzf - struct with fields
%           eigvals
%           entropy (Shannon entropy)
%           anisotropy
%           anisotropy12
%
% See also cloude_pottier and van_zyl
%
% References
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "The Loewner Order and Direction of Detected Change in Sentinel-1 and
% Radarsat-2 Data".
% IEEE Geoscience and Remote Sensing Letters 17(2), 242-246 (2020).
% DOI:10.1109/LGRS.2019.2918636.
%
% Allan A. Nielsen (2019).
% "Fast matrix based computation of eigenvalues and the Loewner order
% in PolSAR data".
% IEEE Geoscience and Remote Sensing Letters 17(10), 1727-1731 (2020).
% DOI:10.1109/LGRS.2019.2952202.
%
% Allan A. Nielsen (2020).
% "Fast Matrix Based Computation of Eigenvalues in PolSAR Data".
% IEEE International Geoscience and Remote Sensing Symposium,
% IGARSS 2021, On-line event, Brussels, Belgium, 12-16 July 2021.
% Accepted for IEEE International Geoscience and Remote Sensing Symposium,
% IGARSS 2020, Hawaii, USA, 19-24 July 2020.
%
% Allan A. Nielsen, Paul J. Connetable, Henning Skriver and Knut Conradsen.
% "A version of the Cloude/van Zyl decomposition for full/quad pol
% synthetic aperture radar data". DTU Compute Technical Report
% DOI:10.11581/dtu.00000237

% (c) Copyright 2019-2021
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 6 Aug 2019
% 31 May 2021 /van Zyl/Cloude "look-alike" sorting of eigvals
%              for full/quad pol

pX = size(X);
nrows = pX(1);
ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    nvars = 1;
else
    error('polsar_eigvals: wrong sizes of input images');
end

if nvars==9 % full pol
    p = 3;
elseif nvars==4 % dual pol
    p = 2;
elseif nvars==2 % dual pol, diagonal only
    p = 2;
elseif nvars==3 % full pol, diagonal only
    p = 3;
elseif nvars==5 % full pol, azimuthal symmetry
    p = 3;
elseif nvars==1 % single channel
    p = 1;
else
    error('polsar_eigvals: wrong dimensionality of input images')
end

if isa(X,'double')
    vzf.eigvals = nan(nrows,ncols,p,'double');
    Lamb1 = nan(nrows,ncols,'double');
    Lamb2 = nan(nrows,ncols,'double');
    Lamb3 = nan(nrows,ncols,'double');
else
    vzf.eigvals = nan(nrows,ncols,p,'single');
    Lamb1 = nan(nrows,ncols,'single');
    Lamb2 = nan(nrows,ncols,'single');
    Lamb3 = nan(nrows,ncols,'single');
end

if nvars==4 % dual pol
    disc = sqrt((X(:,:,1)-X(:,:,4)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    vzf.eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,4) + disc);
%     vzf.eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,4) - disc);
    B = X(:,:,1)+X(:,:,4); % trace
    clear X
    vzf.eigvals(:,:,2) = vzf.eigvals(:,:,1) - disc;
    clear disc
	p1 = vzf.eigvals(:,:,1)./B;
	p2 = vzf.eigvals(:,:,2)./B;
    clear B
	vzf.entropy = (-p1.*log(p1)-p2.*log(p2))/log(2); % log base 2
    clear p1 p2
    vzf.anisotropy12 = (vzf.eigvals(:,:,1)-vzf.eigvals(:,:,2))./ ...
        (vzf.eigvals(:,:,1)+vzf.eigvals(:,:,2));
%     vzf.gini = (p1.*(1-p1)+p2.*(1-p2))*2;
% 	vzf.energy = p1.*p1+p2.*p2;
%     clear p1 p2
%     vzf.renyi = -log(vzf.energy)/log(2); % Rényi or collision entropy
elseif nvars==9% full pol
	a2 = X(:,:,2).^2+X(:,:,3).^2;
	rho2 = X(:,:,4).^2+X(:,:,5).^2;
	b2 = X(:,:,7).^2+X(:,:,8).^2;
	% minus determinant
    D = -(X(:,:,6).*(X(:,:,1).*X(:,:,9)-rho2)-X(:,:,1).*b2-X(:,:,9).*a2 + ...
		2.0*(X(:,:,4).*(X(:,:,2).*X(:,:,7)-X(:,:,3).*X(:,:,8)) + ...
             X(:,:,5).*(X(:,:,3).*X(:,:,7)+X(:,:,2).*X(:,:,8))));
	C = X(:,:,6).*(X(:,:,1)+X(:,:,9))+X(:,:,1).*X(:,:,9)-rho2-a2-b2;
    clear a2 rho2 b2
	% minus trace
    B = -(X(:,:,1)+X(:,:,6)+X(:,:,9));
    idx = X(:,:,4)<0; % real part of rho<0: double bounce
    % A = 1;
% 	Q = (2.0/27.0*(B./A).^3-1.0/3.0*B.*C./A.^2+D./A)/2.0;
	Q = (2.0/27.0*B.^3-1.0/3.0*B.*C+D)/2.0;
    clear D
%   P = (-1.0/3.0*(B./A).^2+C./A)/3.0;
    P = (-1.0/3.0*B.^2+C)/3.0;
    clear C
%     P3 = P.^3;
%     Q2P3 = Q.^2+P3;
%     if (Q2P3>0) % negation of casus irreducibilis
%         imtool(Q2P3)
%     end
%     clear Q2P3 P3
%   % casus irreducibilis: Q2P3 must be non-positive
    smp = sqrt(-P);
    theta = acos(Q./P./smp);
%     clear Q P
	vzf.eigvals(:,:,1) =  2*smp.*cos( theta    /3.0);
	vzf.eigvals(:,:,2) = -2*smp.*cos((theta+pi)/3.0);
% 	vzf.eigvals(:,:,3) = -2*smp.*cos((theta-pi)/3.0);
    clear smp
	vzf.eigvals(:,:,3) = -sum(vzf.eigvals(:,:,1:2),3);
    vzf.eigvals = vzf.eigvals - B/3.0;
% %     vzf.eigvals = sort(vzf.eigvals,3,'descend');
%     vzf.anisotropy = (vzf.eigvals(:,:,2)-vzf.eigvals(:,:,3))./ ...
%         (vzf.eigvals(:,:,2)+vzf.eigvals(:,:,3));
	p1 = -vzf.eigvals(:,:,1)./B;
	p2 = -vzf.eigvals(:,:,2)./B;
	p3 = -vzf.eigvals(:,:,3)./B;
    clear B
	vzf.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    vzf.anisotropy = (vzf.eigvals(:,:,2)-vzf.eigvals(:,:,3))./ ...
        (vzf.eigvals(:,:,2)+vzf.eigvals(:,:,3));
    vzf.anisotropy12 = (vzf.eigvals(:,:,1)-vzf.eigvals(:,:,2))./ ...
        (vzf.eigvals(:,:,1)+vzf.eigvals(:,:,2));
%     vzf.gini = (p1.*(1-p1)+p2.*(1-p2)+p3.*(1-p3))*3/2;
% 	vzf.energy = p1.*p1+p2.*p2+p3.*p3;
%     clear p1 p2 p3
%     vzf.renyi = -log(vzf.energy)/log(3); % Rényi or collision entropy
%     vzf.cos3t = theta;
% ------- azimuthal symmetry "look-alike"
    diffHVHV = abs(vzf.eigvals - repmat(X(:,:,6),1,1,3));
    clear X
    [~,idxsort] = sort(diffHVHV,3,'descend');
    clear diffHVHV
    h = waitbar(0,'Sorting eigvals...','name','Sort eigvals');
    for ic=1:ncols
        waitbar(ic/ncols,h);
        for ir=1:nrows
            Lamb1(ir,ic) = vzf.eigvals(ir,ic,idxsort(ir,ic,1));
            Lamb2(ir,ic) = vzf.eigvals(ir,ic,idxsort(ir,ic,2));
            Lamb3(ir,ic) = vzf.eigvals(ir,ic,idxsort(ir,ic,3));
        end
    end
    close(h)
    clear idxsort
    vzf.eigvals(:,:,3) = Lamb3;
    clear Lamb3
    vzf.eigvals(:,:,1:2) = sort(cat(3,Lamb1,Lamb2),3,'descend');
    Lamb1 = vzf.eigvals(:,:,1);
    Lamb2 = vzf.eigvals(:,:,2);
    % after changes: eigvals(:,:,1) always associated with odd (single) bounce, surface
    aux = Lamb1(idx);
    Lamb1(idx) = Lamb2(idx);
    Lamb2(idx) = aux;
    clear aux idx
    vzf.eigvals(:,:,1) = Lamb1;
    clear Lamb1
    vzf.eigvals(:,:,2) = Lamb2;
    clear Lamb2
% -------
elseif nvars==2 % dual pol diag only
    vzf.eigvals = sort(X,3,'descend');
    B = X(:,:,1)+X(:,:,2); % trace
    clear X
	p1 = vzf.eigvals(:,:,1)./B;
	p2 = vzf.eigvals(:,:,2)./B;
    clear B
	vzf.entropy = (-p1.*log(p1)-p2.*log(p2))/log(2); % log base 2
    clear p1 p2
    vzf.anisotropy12 = (vzf.eigvals(:,:,1)-vzf.eigvals(:,:,2))./ ...
        (vzf.eigvals(:,:,1)+vzf.eigvals(:,:,2));
%     vzf.gini = (p1.*(1-p1)+p2.*(1-p2))*2;
% 	vzf.energy = p1.*p1+p2.*p2;
%     clear p1 p2
%     vzf.renyi = -log(vzf.energy)/log(2); % Rényi or collision entropy
elseif nvars==3 % full pol diag only
    vzf.eigvals = sort(X,3,'descend');
    B = X(:,:,1)+X(:,:,2)+X(:,:,3); % trace
    clear X
%     vzf.anisotropy = (vzf.eigvals(:,:,2)-vzf.eigvals(:,:,3))./ ...
%         (vzf.eigvals(:,:,2)+vzf.eigvals(:,:,3));
	p1 = vzf.eigvals(:,:,1)./B;
	p2 = vzf.eigvals(:,:,2)./B;
	p3 = vzf.eigvals(:,:,3)./B;
    clear B
	vzf.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    clear p1 p2 p3
    vzf.anisotropy = (vzf.eigvals(:,:,2)-vzf.eigvals(:,:,3))./ ...
        (vzf.eigvals(:,:,2)+vzf.eigvals(:,:,3));
    vzf.anisotropy12 = (vzf.eigvals(:,:,1)-vzf.eigvals(:,:,2))./ ...
        (vzf.eigvals(:,:,1)+vzf.eigvals(:,:,2));
%     vzf.gini = (p1.*(1-p1)+p2.*(1-p2)+p3.*(1-p3))*3/2;
% 	vzf.energy = p1.*p1+p2.*p2+p3.*p3;
%     clear p1 p2 p3
%     vzf.renyi = -log(vzf.energy)/log(3); % Rényi or collision entropy
elseif nvars==5 % azimuthal symmetry
    disc = sqrt((X(:,:,1)-X(:,:,5)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    vzf.eigvals(:,:,1) = X(:,:,4);
    vzf.eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,5) + disc);
%     vzf.eigvals(:,:,3) = 0.5*(X(:,:,1)+X(:,:,5) - disc);
    B = X(:,:,1)+X(:,:,4)+X(:,:,5); % trace
    idx = X(:,:,2)<0; % real part of rho<0: double bounce
    clear X
    vzf.eigvals(:,:,3) = vzf.eigvals(:,:,2) - disc;
    clear disc
    vzf.eigvals = sort(vzf.eigvals,3,'descend'); % to be on the safe side
	p1 = vzf.eigvals(:,:,1)./B;
	p2 = vzf.eigvals(:,:,2)./B;
	p3 = vzf.eigvals(:,:,3)./B;
    clear B
	vzf.entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    clear p1 p2 p3
    vzf.anisotropy = (vzf.eigvals(:,:,2)-vzf.eigvals(:,:,3))./ ...
        (vzf.eigvals(:,:,2)+vzf.eigvals(:,:,3));
    vzf.anisotropy12 = (vzf.eigvals(:,:,1)-vzf.eigvals(:,:,2))./ ...
        (vzf.eigvals(:,:,1)+vzf.eigvals(:,:,2));
%     vzf.gini = (p1.*(1-p1)+p2.*(1-p2)+p3.*(1-p3))*3/2;
% 	vzf.energy = p1.*p1+p2.*p2+p3.*p3;
%     clear p1 p2 p3
%     vzf.renyi = -log(vzf.energy)/log(3); % Rényi or collision entropy
    Lamb1 = vzf.eigvals(:,:,1);
    Lamb2 = vzf.eigvals(:,:,2);
    aux = Lamb1(idx);
    Lamb1(idx) = Lamb2(idx);
    Lamb2(idx) = aux;
    clear aux idx
    vzf.eigvals(:,:,1) = Lamb1;
    clear Lamb1
    vzf.eigvals(:,:,2) = Lamb2;
    clear Lamb2
else % scalar, single channel power data
    warning('polsar_eigvals: eigvals and entropy for scalar are meaningless')
    vzf.eigvals = X;
    vzf.entropy = nan;
    vzf.anisotropy = nan;
    clear X
end

% pp = vzf.eigvals./repmat(sum(vzf.eigvals,3),1,1,p);
% vzf.entropy = -sum(pp.*log(pp),3)/log(p);
% % if p==3
% %     vzf.entropy = -sum(pp.*log(pp),3)/log(3);
% % elseif p==2
% %     vzf.entropy = -sum(pp.*log(pp),3)/log(2);
% % end

vzf.comment1 = 'eigvals1 is associated with single (odd) bounce (surface scattering)';
vzf.comment2 = 'eigvals2 is associated with double (even) bounce';
vzf.comment3 = 'eigvals3 is associated with the remainder (typically volume scattering)';
vzf.comment4 = 'Low entropy is associated with a dominant scattering mechanism, high entropy with randomness';
vzf.comment5 = 'Visualization: RGB of eigvals 2, 3, 1';

while 0

cd FSAR/Vejers
load Xband_Vejers.mat
cd ../..
SAR = cat(3,HHHH,real(HHHV),imag(HHHV),real(HHVV),imag(HHVV),HVHV,real(HVVV),imag(HVVV),VVVV);
SAR = double(SAR(:,1601:1601+4800,:)); % used in paper
% SAR = double(SAR(:,1601:1601+4800-1,:));
% SAR = double(SAR);
[nr0,nl0] = size(HHHH);
[nr,nl,~] = size(SAR);

% figure('Name','RG power HHHV, HVVV')
% imshowrgbprctile(cat(3,abs(HHHV).^2,abs(HVVV).^2,zeros(nr0,nl0)));
% figure('Name','RGB power HHHV, HVVV, HHVV')
% imshowrgbprctile(cat(3,abs(HHHV).^2,abs(HVVV).^2,abs(HHVV).^2));

% % HHHV
% aux = complex(SAR(:,:,2),SAR(:,:,3));
% % imshowsar(aux,'HHHV'), %imshowsar(HHHV,'HHHV');
% figure('Name','magnitude HHHV'), imshow(abs(aux),[])
% imtool(abs(aux),[])
% 
% % HVVV
% aux = complex(SAR(:,:,7),SAR(:,:,8));
% % imshowsar(aux,'HVVV'), %imshowsar(HVVV,'HVVV');
% figure('Name','magnitude HVVV'), imshow(abs(aux),[])
% imtool(abs(aux),[])
% 
% % HHVV
% aux = complex(SAR(:,:,4),SAR(:,:,5));
% % imshowsar(aux,'HHVV'), %imshowsar(HHVV,'HHVV');
% figure('Name','magnitude HHVV'), imshow(abs(aux),[])
% imtool(abs(aux),[])

% clear HHHH HHHV HHVV HVHV HVVV VVVV

tic, vz = van_zyl(SAR), toc
tic, vzf = van_zyl_full(SAR), toc

%% IGARSS 2025
figure('Name','van Zyl/Cloude decomposition')
imshowrgbprctile(sqrt(vzf.eigvals),[2 3 1],[0.5 99.5]);
set(gca,'DataAspectRatio',[0.29 0.25 1])
print -depsc vzfVejers
L = 20;
wb = wishart_blockdiag(double(SAR),L)
aux = log(SAR(:,:,6)); % <Shvhv Shvhv*>
[minmaxx] = prctile(aux(:),[0.5 99.5]);
auxu8 = uint8(255*((aux-minmaxx(1))/(minmaxx(2)-minmaxx(1))));
idx = wb.P>0.99;
rr = auxu8;
rr(idx) = 255;
gb = auxu8;
gb(idx) = 0;
figure('Name','wb.P>0.99')
imshow(cat(3,rr,gb,gb));
set(gca,'DataAspectRatio',[0.29 0.25 1])
print -depsc wbVejers_99
%%

d_entropy = vz.entropy - vzf.entropy;
d_hvhv = SAR(:,:,6) - vzf.eigvals(:,:,3);

% figure('Name','van Zyl/Cloude decomposition')
% imshowrgbprctile(vzx.lambda,[2 3 1],[0.5 99.5]);
% set(gca,'DataAspectRatio',[0.29 0.25 1])
% figure('Name','van Zyl/Cloude decomposition [10*log10]')
% imshowrgb(10*log10(vzx.Lambda(:,:,[2 3 1])));
% set(gca,'DataAspectRatio',[0.29 0.25 1])

% % set(gca,'DataAspectRatio',[0.29 0.25 1])
% set(gca,'DataAspectRatio',[0.30 0.25 1])

figure('Name','Cloude/van Zyl decomposition')
h(1) = subplot_tight(2,2,1,[0.005,0.005]);
imshowrgbprctile(vzx.eigvals(:,:,[2 3 1]));
set(gca,'DataAspectRatio',[0.30 0.25 1])
h(2) = subplot_tight(2,2,2,[0.005,0.005]);
imshowrgbprctile(vzf.eigvals(:,:,[2 3 1]));
set(gca,'DataAspectRatio',[0.30 0.25 1])
h(3) = subplot_tight(2,2,3,[0.005,0.005]);
imshow(d_entropy,[]);
set(gca,'DataAspectRatio',[0.30 0.25 1])
h(4) = subplot_tight(2,2,4,[0.005,0.005]);
imshow(d_hvhv,[]);
set(gca,'DataAspectRatio',[0.30 0.25 1])
linkaxes(h,'xy')

figure, histogram(d_hvhv)

% N = [1 0 1;1 0 -1;0 sqrt(2) 0]/sqrt(2); % cov2coh: T = N*C*N'
% N is orthogonal: inv(N) is N'
% Vfcoh = N'*Vfcov;
% Vfcov = N*Vfcoh;

% Big house, large entropy difference
row = 2700;
col = 1540;
row = 2670;
col = 1550;
vzf.entropy(row,col)
vzx.entropy(row,col)
aux = wc2mat(cov2coh(SAR(row,col,:)));
[Vf,Df] = eig(aux)
[Vfcoh,Dfcoh] = eig(wc2mat(cov2coh(SAR(row,col,:))))
% aux(1,2) = 0; % good for covmat, not for cohmat
% aux(2,1) = 0; % good for covmat, not for cohmat
aux(1,3) = 0;
aux(3,1) = 0;
aux(2,3) = 0;
aux(3,2) = 0
[Vx,Dx] = eig(aux)
figure
bar([diag(Df),diag(Dx),squeeze(vzx.eigvals(row,col,[1 2 3]))])
bar([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,[1 2 3]))])
grid on
legend('full','vZyl','Location','northeast')
title('Eigenvalues - Big house, large entropy difference')
% alpha
alpha = acosd(abs(Vfcoh(1,:)))
p = diag(Dfcoh)'/sum(diag(Dfcoh))
alpha_mean = sum(p.*alpha)

ksi = aux(1,1);
rho = aux(1,3);
eta = aux(2,2);
zeta = aux(3,3);
del = sqrt((ksi-zeta)^2+4*rho*rho');
aup = (zeta-ksi+del)^2;
sqp = sqrt(aup/(aup+4*rho*rho'));
alpha = 2*rho/(zeta-ksi+del);
aum = (zeta-ksi-del)^2;
sqm = sqrt(aum/(aum+4*rho*rho'));
beta = -alpha/(alpha*alpha')

norm(Vf(:,1)-Vx(:,1))
norm(Vf(:,2)-Vx(:,2))
norm(Vf(:,3)-Vx(:,3))

% Road, small entropy difference
row = 1340; % road?
col = 3780; % road?
row = 1342; % road
col = 3502; % road
row = 1262; % road
col = 3755; % road
vzf.entropy(row,col)
vzx.entropy(row,col)
aux = wc2mat(cov2coh(SAR(row,col,:)));
[Vf,Df] = eig(aux)
[Vfcoh,Dfcoh] = eig(wc2mat(cov2coh(SAR(row,col,:))))
% aux(1,2) = 0; % good for covmat, not for cohmat
% aux(2,1) = 0; % good for covmat, not for cohmat
aux(1,3) = 0;
aux(3,1) = 0;
aux(2,3) = 0;
aux(3,2) = 0
[Vx,Dx] = eig(aux)
figure
bar([diag(Df),diag(Dx),squeeze(vzx.eigvals(row,col,[1 2 3]))])
bar([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,[1 2 3]))])
grid on
legend('full','vZyl','Location','northeast')
title('Eigenvalues - Road, small entropy difference')
% alpha
alpha = acosd(abs(Vfcoh(1,:)))
p = diag(Dfcoh)'/sum(diag(Dfcoh))
alpha_mean = sum(p.*alpha)

ksi = aux(1,1);
rho = aux(1,3);
eta = aux(2,2);
zeta = aux(3,3);
del = sqrt((ksi-zeta)^2+4*rho*rho');
aup = (zeta-ksi+del)^2;
sqp = sqrt(aup/(aup+4*rho*rho'));
alpha = 2*rho/(zeta-ksi+del);
aum = (zeta-ksi-del)^2;
sqm = sqrt(aum/(aum+4*rho*rho'));
beta = -alpha/(alpha*alpha')

norm(Vf(:,1)-Vx(:,1))
norm(Vf(:,2)-Vx(:,2))
norm(Vf(:,3)-Vx(:,3))

% Wood, small entropy difference
row = 487; % wood
col = 1955; % wood
row = 1217; % road
col = 4312; % road
row = 1340; % wood
col = 3780; % wood
row = 1065; % wood in paper
col = 4140; % wood in paper
vzf.entropy(row,col)
vzx.entropy(row,col)
aux = wc2mat(cov2coh(SAR(row,col,:)));
[Vf,Df] = eig(aux)
% figure, quiver3(0,0,0,Vf(1,1),Vf(2,1),Vf(3,1)) % no good for complex
[Vfcoh,Dfcoh] = eig(wc2mat(cov2coh(SAR(row,col,:))))
% aux(1,2) = 0; % good for covmat, not for cohmat
% aux(2,1) = 0; % good for covmat, not for cohmat
aux(1,3) = 0;
aux(3,1) = 0;
aux(2,3) = 0;
aux(3,2) = 0
[Vx,Dx] = eig(aux)
figure
bar([diag(Df),diag(Dx),squeeze(vzx.eigvals(row,col,[3 2 1]))])
legend('full','azim','vZyl','Location','northwest')
bar([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,:))])
legend('full','azim','Location','northeast')
grid on
title('Eigenvalues - Wood, small entropy difference')
% title('Eigenvalues - Road, small entropy difference')
% alpha
alpha = acosd(abs(Vfcoh(1,:)))
p = diag(Dfcoh)'/sum(diag(Dfcoh))
alpha_mean = sum(p.*alpha)

ksi = aux(1,1);
rho = aux(1,3);
eta = aux(2,2);
zeta = aux(3,3);
del = sqrt((ksi-zeta)^2+4*rho*rho');
aup = (zeta-ksi+del)^2;
sqp = sqrt(aup/(aup+4*rho*rho'));
alpha = 2*rho/(zeta-ksi+del);
aum = (zeta-ksi-del)^2;
sqm = sqrt(aum/(aum+4*rho*rho'));
beta = -alpha/(alpha*alpha');

norm(Vf(:,1)-Vx(:,1))
norm(Vf(:,2)-Vx(:,2))
norm(Vf(:,3)-Vx(:,3))
% figure, quiver3(0,0,0,Vx(1,1),Vx(2,1),Vx(3,1)) % no good for complex

% Next to swimming pool, large entropy difference
row = 1296;
col = 1647;
row = 1290;
col = 1650;
vzf.eigvals(row,col,:)
vzf.entropy(row,col)
vzx.entropy(row,col)
aux = wc2mat(cov2coh(SAR(row,col,:)));
[Vf,Df] = eig(aux)
% figure, quiver3(0,0,0,Vf(1,1),Vf(2,1),Vf(3,1)) % no good for complex
[Vfcoh,Dfcoh] = eig(wc2mat(cov2coh(SAR(row,col,:))))
% aux(1,2) = 0; % good for covmat, not for cohmat
% aux(2,1) = 0; % good for covmat, not for cohmat
aux(1,3) = 0;
aux(3,1) = 0;
aux(2,3) = 0;
aux(3,2) = 0
[Vx,Dx] = eig(aux)
figure
bar([diag(Df),diag(Dx),squeeze(vzx.eigvals(row,col,[3 2 1]))])
legend('full','azim','vZyl','Location','northwest')
bar([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,:))])
legend('full','vZyl','Location','northeast')
grid on
title('Eigenvalues - West of swimming pool, large entropy difference')
% alpha
alpha = acosd(abs(Vfcoh(1,:)))
p = diag(Dfcoh)'/sum(diag(Dfcoh))
alpha_mean = sum(p.*alpha)

ksi = aux(1,1);
rho = aux(1,3);
eta = aux(2,2);
zeta = aux(3,3);
del = sqrt((ksi-zeta)^2+4*rho*rho');
aup = (zeta-ksi+del)^2;
sqp = sqrt(aup/(aup+4*rho*rho'));
alpha = 2*rho/(zeta-ksi+del);
aum = (zeta-ksi-del)^2;
sqm = sqrt(aum/(aum+4*rho*rho'));
beta = -alpha/(alpha*alpha');

norm(Vf(:,1)-Vx(:,1))
norm(Vf(:,2)-Vx(:,2))
norm(Vf(:,3)-Vx(:,3))

%%%

% Histograms for paper
figure('Name','Cloude/van Zyl decomposition')
% row = 1340; % road?
% col = 3780; % road?
% row = 1342; % road
% col = 3502; % road
row = 1262; % road
col = 3755; % road
hh(1) = subplot_tight(1,2,1,[0.06,0.06]);
bar([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,:))])
% bar(10*log10([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,:))]))
legend('full','refl','Location','northeast')
grid on
% title('Eigenvalues - Road, small entropy difference')
% row = 2700; % big house
% col = 1540; % big house
row = 2670; % big house
col = 1550; % big house
% row = 1065, col = 4140 % Ha.alpha = 45.53, vegetation
% row = 1340, col = 3780 % Ha.alpha = 29.71, Wood
hh(2) = subplot_tight(1,2,2,[0.06,0.06]);
bar([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,:))])
% bar(10*log10([squeeze(vzf.eigvals(row,col,:)),squeeze(vzx.eigvals(row,col,:))]))
legend('full','refl','Location','northeast')
grid on
% title('Eigenvalues - Big house, large entropy difference')

end % while
