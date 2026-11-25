function C = coh2cov(T)

% coh2cov transforms upper triangular polarimetric coherency matrix to covariance
%
% N = [1 0 1;1 0 -1;0 sqrt(2) 0]/sqrt(2);
% C = N'*T*N; % covariance matrix
% T = N*C*N'; % coherency matrix
%
% see also cov2coh
%
% % det(C) - det(T) % is zero

% (c) Copyright 2018-2019
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, https://people.compute.dtu.dk/alan
% 17 Sep 2018
% 15-19 Aug 2019 (ndims=4 added, data must be quad/full pol)

nd = ndims(T);
if nd<3 || nd>4
    error('ndims of input must be 3 or 4');
end

C = T;
sqrt2 = sqrt(2);

if nd==4
    [~,~,p,nt] = size(T);
    if p~=9, error('coh2cov: data must be quad/full pol'), end
    for ix=1:nt
        C(:,:,1,ix) = T(:,:,1,ix) +2*T(:,:,2,ix) + T(:,:,6,ix);
        C(:,:,2,ix) = sqrt2*(T(:,:,4,ix) + T(:,:,7,ix));
        C(:,:,3,ix) = sqrt2*(T(:,:,5,ix) + T(:,:,8,ix));
        C(:,:,4,ix) = T(:,:,1,ix) - T(:,:,6,ix);
        C(:,:,5,ix) = -2*T(:,:,3,ix);
        C(:,:,6,ix) = 2*T(:,:,9,ix);
        C(:,:,7,ix) = sqrt2*(T(:,:,4,ix) - T(:,:,7,ix));
        C(:,:,8,ix) = sqrt2*(T(:,:,8,ix) - T(:,:,5,ix));
        C(:,:,9,ix) = T(:,:,1,ix) -2*T(:,:,2,ix) + T(:,:,6,ix);
    end
else
    [~,~,p] = size(T);
    if p~=9, error('cov2coh: data must be quad/full pol'), end
    C(:,:,1) = T(:,:,1) +2*T(:,:,2) + T(:,:,6);
    C(:,:,2) = sqrt2*(T(:,:,4) + T(:,:,7));
    C(:,:,3) = sqrt2*(T(:,:,5) + T(:,:,8));
    C(:,:,4) = T(:,:,1) - T(:,:,6);
    C(:,:,5) = -2*T(:,:,3);
    C(:,:,6) = 2*T(:,:,9);
    C(:,:,7) = sqrt2*(T(:,:,4) - T(:,:,7));
    C(:,:,8) = sqrt2*(T(:,:,8) - T(:,:,5));
    C(:,:,9) = T(:,:,1) -2*T(:,:,2) + T(:,:,6);
end
C = 0.5*C;
