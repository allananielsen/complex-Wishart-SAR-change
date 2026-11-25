function Ha = H_alpha(SAR,varargin)

% Cloude-Pottier entropy-anisotropy-alpha (coherency matrix based)
% Mean alpha is determined bmo the eigenvector-eigenvalue identity
%
% Input
%   SAR       - 3-D array of 3x3 complex covariance matrix
%               (polarimetric SAR image in covariance matrix form,
%               reading order is alphabetical quad/full pol
%               (as read by freadenvisar)
%                 [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*)
%                          ShvShv*                 Re(ShvSvv*) Im(ShvSvv*)
%                                                SvvSvv*]
%   flag_beta - calculate also mean beta angle (optional)
%               (if flag_beta is true, default is false)
%
% Output
%   Ha.H       - entropy
%   Ha.A       - anisotropy
%   Ha.alpha   - mean alpha
%   Ha.beta    - mean beta (if flag_beta is true)
%
% References
% [0] Jacobi, 1834.
% [1] S. R. Cloude and E. Pottier, “An entropy based classification scheme
%     for land applications of polarimetric SAR," IEEE Transactions on
%     Geoscience and Remote Sensing, vol. 35, no. 1, pp. 68–78, Jan. 1997,
%     https://doi.org/10.1109/36.551935.
% [2] J.-S. Lee, M. R. Grunes, T. L. Ainsworth, L.-J. Du, D. L. Schuler,
%     and S. R. Cloude, “Unsupervised classification using polarimetric
%     decomposition and the complex Wishart classifier," IEEE Transactions
%     on Geoscience and Remote Sensing, vol. 37, no. 5, pp. 2249–2258, Sep.
%     1999, https://doi.org/10.1109/36.789621.
% [3] S. Hensley,"The eigenvector-eigenvalue identity and radar polarimetry,"
%     in International Geoscience and Remote Sensing Symposium (IGARSS),
%     Brussels, Belgium, Jul. 2021, IEEE.
%     https://doi.org/10.1109/IGARSS47720.2021.9554916.
% [4] P. B. Denton, S. J. Parke, T. Tao, and X. Zhang, "Eigenvectors from
%     eigenvalues: A survey of a basic identity in linear algebra,"
%     Bulletin of the American Mathematical Sociey, vol. 59, no. 1,
%     pp. 31–58, Jan. 2022, https://doi.org/10.1090/bull/1722.
% [5] A. A. Nielsen, "The eigenvector-eigenvalue identity applied to fast
%     calculation of polSAR scattering characterization,"
%     IEEE Geoscience and Remote Sensing Letters, Apr. 2022.
%     https://doi.org/10.1109/LGRS.2022.3169994.

% (c) Copyright 2022
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 15 Feb 2022

flag_beta = false;
if nargin>1
    flag_beta = true;
end

np = size(SAR,3);
if np~=9
    error('input must be upper-right 3x3 Hermitian matrix')
end

%%
% SAR = cov2coh(SAR);
lam = polsar_eigvals_only(SAR);
p = lam.eigvals./repmat(sum(lam.eigvals,3),1,1,3);
Ha.H = -sum(p.*log(p),3)/log(3); % entropy
Ha.A = (lam.eigvals(:,:,2)-lam.eigvals(:,:,3))./ ...
       (lam.eigvals(:,:,2)+lam.eigvals(:,:,3)); % anisotropy

%% All first components
lam1 = polsar_eigvals_only(SAR(:,:,[6 7 8 9])); % first minor M1
e1sq(:,:,1) = ((lam.eigvals(:,:,1)-lam1.eigvals(:,:,1))  .* ...
               (lam.eigvals(:,:,1)-lam1.eigvals(:,:,2))) ./ ...
              ((lam.eigvals(:,:,1)- lam.eigvals(:,:,2))  .* ...
               (lam.eigvals(:,:,1)- lam.eigvals(:,:,3)));
e1sq(:,:,2) = ((lam.eigvals(:,:,2)-lam1.eigvals(:,:,1))  .* ...
               (lam.eigvals(:,:,2)-lam1.eigvals(:,:,2))) ./ ...
              ((lam.eigvals(:,:,2)- lam.eigvals(:,:,1))  .* ...
               (lam.eigvals(:,:,2)- lam.eigvals(:,:,3)));
e1sq(:,:,3) = ((lam.eigvals(:,:,3)-lam1.eigvals(:,:,1))  .* ...
               (lam.eigvals(:,:,3)-lam1.eigvals(:,:,2))) ./ ...
              ((lam.eigvals(:,:,3)- lam.eigvals(:,:,1))  .* ...
               (lam.eigvals(:,:,3)- lam.eigvals(:,:,2)));
% Ha.e1sq = e1sq; % check with e1sq from eigvecijsq

clear lam1
alpha = acosd(sqrt(e1sq)); % alpha angles
clear e1sq
Ha.alpha = sum(p.*alpha,3); % alpha_mean
clear alpha

%%
if flag_beta

    % All second and third components
    lam1 = polsar_eigvals_only(SAR(:,:,[1 4 5 9])); % second minor M2
    e2sq(:,:,1) = ((lam.eigvals(:,:,1)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,1)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,1)- lam.eigvals(:,:,2))  .* ...
                   (lam.eigvals(:,:,1)- lam.eigvals(:,:,3)));
    e2sq(:,:,2) = ((lam.eigvals(:,:,2)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,2)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)- lam.eigvals(:,:,3)));
    e2sq(:,:,3) = ((lam.eigvals(:,:,3)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,3)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)- lam.eigvals(:,:,2)));
    
    lam1 = polsar_eigvals_only(SAR(:,:,[1 2 3 6])); % third minor M3
    clear SAR
    e3sq(:,:,1) = ((lam.eigvals(:,:,1)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,1)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,1)- lam.eigvals(:,:,2))  .* ...
                   (lam.eigvals(:,:,1)- lam.eigvals(:,:,3)));
    e3sq(:,:,2) = ((lam.eigvals(:,:,2)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,2)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)- lam.eigvals(:,:,3)));
    e3sq(:,:,3) = ((lam.eigvals(:,:,3)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,3)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)- lam.eigvals(:,:,2)));
    clear lam lam1
    
    beta = atand(sqrt(e3sq./e2sq)); % beta angles
    clear e2sq e3sq
    Ha.beta = sum(p.*beta,3); % beta_mean
    clear p beta

end % if flag_beta

Ha.comment_H = 'Low entropy is associated with a single dominant scattering mechanism, high entropy with random scattering';
Ha.comment_A = 'Anisotropy is a measure of departure from azimuthal symmetry (low anisotropy is associated with azimuthal symmetry)';
Ha.comment_alpha0 = 'alpha = 0: single bounce (surface) scattering';
Ha.comment_alpha45 = 'alpha = 45: scattering from differently oriented dipoles (e.g., vegetation, dihedral targets rotated with respect to flight line)';
Ha.comment_alpha90 = 'alpha = 90: (conductive) double bounce scattering';
