function vz = van_zyl(SAR)

% van_zyl calculates van Zyl's (Cloude's) decomposition of covariance data
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
%
% Eigenvalues are not sorted (unlike in polsar_eigvals), but Lambda1 and
% Lambda2 (and alpha and beta) values are swapped if Re(alpha)<0 so that
% vz.Lambda(:,:,1) is associated with single (surface, odd) bounce scattering (blue)
% vz.Lambda(:,:,2) is associated with double (even) bounce scattering (red)
% vz.Lambda(:,:,3) is the remainder associated with volume scattering (green)
%
% See also van_zyl_full and polsar_eigvals
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
% 27 April 2021

%%
pX = size(SAR);
nrows = pX(1);
ncols = pX(2);
if  ndims(SAR)==3
    ntims = 1;
elseif ndims(SAR)==4
    ntims = pX(4);
else
    error('van_zyl: wrong dimensionality of input images');
end
nvars = pX(3);

if nvars==9 % full pol
    p = 3;
    fprintf('van_zyl: quad/full pol input, reflection symmetry enforced\n')
elseif nvars==5 % full pol, reflection symmetry
    p = 3;
else
    error('van_zyl: input must be quad/full pol or reflection symmetry')
end

% full pol
if nvars==9
    SAR = SAR(:,:,[1 4 5 6 9],:);
end

for it=1:ntims
    %%
    if class(SAR)=='double'
        eigvals = nan(nrows,ncols,p-1,'double');
        vz(it).Lambda = nan(nrows,ncols,p,'double');
    else
        eigvals = nan(nrows,ncols,p-1,'single');
        vz(it).Lambda = nan(nrows,ncols,p,'single');
    end

    %%
    % reflection symmetry
    X = SAR(:,:,:,it);
    disc = sqrt((X(:,:,1)-X(:,:,5)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    vz(it).Lambda(:,:,3) = X(:,:,4); % eigvals(:,:,3) = Lambda(:,:,3)
    eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,5) + disc);
    eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,5) - disc);
    B = X(:,:,1)+X(:,:,4)+X(:,:,5); % trace
    p1 = eigvals(:,:,1)./B;
    p2 = eigvals(:,:,2)./B;
    p3 = vz(it).Lambda(:,:,3)./B;
    clear B
    vz(it).entropy = (-p1.*log(p1)-p2.*log(p2)-p3.*log(p3))/log(3); % log base 3
    clear p1 p2 p3
    auxp = X(:,:,5)-X(:,:,1)+disc; % van Zyl
    auxm = X(:,:,5)-X(:,:,1)-disc; % van Zyl
%     auxp = X(:,:,1)-X(:,:,5)+disc; % AA
%     auxm = X(:,:,1)-X(:,:,5)-disc; % AA
    clear disc
    tworho = 2*complex(X(:,:,2),X(:,:,3));
    clear X
    tworho2 = tworho.*conj(tworho);
    auxp2 = auxp.^2;
    vz(it).Lambda(:,:,1) = real(eigvals(:,:,1).*auxp2./(auxp2+tworho2));
    clear auxp2
    auxm2 = auxm.^2;
    vz(it).Lambda(:,:,2) = real(eigvals(:,:,2).*auxm2./(auxm2+tworho2));
    clear auxm2 tworho2 eigvals
    vz(it).alpha = tworho./auxp;
    clear auxp
    beta = tworho./auxm; 
    clear auxm tworho

    %% if real(vz(it).alpha)<0: Lambda1 and alpha associated with even (double) bounce
    % after changes: Lambda1 and alpha always associated with odd (single) bounce, surface
    idx = real(vz(it).alpha)<0;
    vz(it).alpha(idx) = beta(idx);
    clear beta
    Lamb1 = vz(it).Lambda(:,:,1);
    Lamb2 = vz(it).Lambda(:,:,2);
    aux = Lamb1(idx);
    Lamb1(idx) = Lamb2(idx);
    Lamb2(idx) = aux;
    clear aux idx
    vz(it).Lambda(:,:,1) = Lamb1;
    clear Lamb1
    vz(it).Lambda(:,:,2) = Lamb2;
    clear Lamb2
    vz(it).comment1 = 'Lambda1 and alpha are associated with single (odd) bounce (surface scattering)';
    vz(it).comment2 = 'Lambda2 and beta are associated with double (even) bounce';
    vz(it).comment3 = 'Lambda3 is associated with the remainder (typically volume scattering)';
    vz(it).comment4 = 'Low entropy is associated with a dominant scattering mechanism, high entropy with randomness';
    vz(it).comment5 = 'beta = -1./conj(alpha) = -alpha./abs(alpha).^2 is not stored';
    vz(it).comment6 = 'Visualization: RGB of Lambda 2, 3, 1';
%     vz(it).comment7 = 'There is a typo in ref [1] eqs (16) and (17): the numerators under the squareroot sign should be squared';

    %% R: double bounce, G: volume scattering, B: single bounce
    % figure, imshowrgb(10*log10(vz(it).Lambda(:,:,[2 3 1])))
    % figure, imshowrgbprctile(10*log10(vz(it).Lambda(:,:,[2 3 1])))
    % figure, imshowrgbprctile(vz(it).Lambda(:,:,[2 3 1]))

end % for it=1:ntims
