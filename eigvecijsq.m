function esq = eigvecijsq(Hermit)

% Eigenvector-eigenvalue identity for 3x3 Hermitian matrices
% (for fx Cloude-Pottier mean alpha or beta angle from 3x3 coherency matrix)
%
% Input
%   Hermit - 3-D array of 3x3 Hermitian matrices (upper left)
%         (fx polarimetric SAR image in coherency/covariance matrix form,
%         reading order is alphabetical (as read by freadenvisar):
%               [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%                        ShvShv*                 Re(ShvSvv*) Im(ShvSvv*)
%                                                SvvSvv*]
%
% Output
%   esq    - struct with squared components of eigenvectors
%
% References
%
% [0] Jacobi, 1834.
% [1] K. T'. Löwner, "Über monotone Matrixfunktionen,"
%     Mathematische Zeitschrift, vol. 38, pp. 177–216, 1934,
%     https://eudml.org/doc/168495.
% [2] R. C. Thompson, "Principal submatrices of normal and Hermitian
%     matrices,”
%     Illinois Journal of Mathematics, vol. 10, no. 2, pp. 296–308, 1966,
%     https://doi.org/10.1215/ijm/1256055111.
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
% alan@dtu.dk, https://people.compute.dtu.dk/alan
% 15 Feb 2022

np = size(Hermit,3);
if np~=9
    error('input must be upper-right 3x3 Hermitian matrix')
end

esq.comment = 'struct with squared components of eigenvectors';

%%
lam = polsar_eigvals_only(Hermit);

%% All first components
lam1 = polsar_eigvals_only(Hermit(:,:,[6 7 8 9])); % first minor M1
esq.e1sq(:,:,1) = ((lam.eigvals(:,:,1)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,1)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,1)- lam.eigvals(:,:,2))  .* ...
                   (lam.eigvals(:,:,1)- lam.eigvals(:,:,3)));
esq.e1sq(:,:,2) = ((lam.eigvals(:,:,2)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,2)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)- lam.eigvals(:,:,3)));
esq.e1sq(:,:,3) = ((lam.eigvals(:,:,3)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,3)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)- lam.eigvals(:,:,2)));

%% All second components
lam1 = polsar_eigvals_only(Hermit(:,:,[1 4 5 9])); % second minor M2
esq.e2sq(:,:,1) = ((lam.eigvals(:,:,1)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,1)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,1)- lam.eigvals(:,:,2))  .* ...
                   (lam.eigvals(:,:,1)- lam.eigvals(:,:,3)));
esq.e2sq(:,:,2) = ((lam.eigvals(:,:,2)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,2)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)- lam.eigvals(:,:,3)));
esq.e2sq(:,:,3) = ((lam.eigvals(:,:,3)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,3)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)- lam.eigvals(:,:,2)));

%% All third components
lam1 = polsar_eigvals_only(Hermit(:,:,[1 2 3 6])); % third minor M3
clear Hermit
esq.e3sq(:,:,1) = ((lam.eigvals(:,:,1)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,1)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,1)- lam.eigvals(:,:,2))  .* ...
                   (lam.eigvals(:,:,1)- lam.eigvals(:,:,3)));
esq.e3sq(:,:,2) = ((lam.eigvals(:,:,2)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,2)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,2)- lam.eigvals(:,:,3)));
esq.e3sq(:,:,3) = ((lam.eigvals(:,:,3)-lam1.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)-lam1.eigvals(:,:,2))) ./ ...
                  ((lam.eigvals(:,:,3)- lam.eigvals(:,:,1))  .* ...
                   (lam.eigvals(:,:,3)- lam.eigvals(:,:,2)));

%%
esq.comment11 = 'e1sq(:,:,1) contains first  component of first  eigenvector';
esq.comment21 = 'e1sq(:,:,2) contains first  component of second eigenvector';
esq.comment31 = 'e1sq(:,:,3) contains first  component of third  eigenvector';
esq.comment12 = 'e2sq(:,:,1) contains second component of first  eigenvector';
esq.comment22 = 'e2sq(:,:,2) contains second component of second eigenvector';
esq.comment32 = 'e2sq(:,:,3) contains second component of third  eigenvector';
esq.comment13 = 'e3sq(:,:,1) contains third  component of first  eigenvector';
esq.comment23 = 'e3sq(:,:,2) contains third  component of second eigenvector';
esq.comment33 = 'e3sq(:,:,3) contains third  component of third  eigenvector';
