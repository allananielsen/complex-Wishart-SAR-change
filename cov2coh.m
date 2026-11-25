function T = cov2coh(C)

% cov2coh transforms upper triangular polarimetric covariance matrix to coherency
%
% N = [1 0 1;1 0 -1;0 sqrt(2) 0]/sqrt(2);
% C = N'*T*N; % covariance matrix
% T = N*C*N'; % coherency matrix
%
% see also coh2cov, cov2circ, mult_sqrt2
%
% wishart_det(C) - wishart_det(T) is zero

% (c) Copyright 2018-2019
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, https://people.compute.dtu.dk/alan
% 17 Sep 2018
% 15-19 Aug 2019 (ndims=4 added, data must be quad/full pol)

nd = ndims(C);
if nd<3 || nd>4
    error('ndims of input must be 3 or 4');
end

T = C;
sqrt2 = sqrt(2);

if nd==4
    [~,~,p,nt] = size(T);
    if p~=9, error('cov2coh: data must be quad/full pol'), end
    for ix=1:nt
        T(:,:,1,ix) = C(:,:,1,ix) +2*C(:,:,4,ix) + C(:,:,9,ix);
        T(:,:,2,ix) = C(:,:,1,ix) - C(:,:,9,ix);
        T(:,:,3,ix) = -2*C(:,:,5,ix);
        T(:,:,4,ix) = sqrt2*(C(:,:,2,ix) + C(:,:,7,ix));
        T(:,:,5,ix) = sqrt2*(C(:,:,3,ix) - C(:,:,8,ix));
        T(:,:,6,ix) = C(:,:,1,ix) -2*C(:,:,4,ix) + C(:,:,9,ix);
        T(:,:,7,ix) = sqrt2*(C(:,:,2,ix) - C(:,:,7,ix));
        T(:,:,8,ix) = sqrt2*(C(:,:,3,ix) + C(:,:,8,ix));
        T(:,:,9,ix) = 2*C(:,:,6,ix);
    end
else
    [~,~,p] = size(T);
    if p~=9, error('cov2coh: data must be quad/full pol'), end
    T(:,:,1) = C(:,:,1) +2*C(:,:,4) + C(:,:,9);
    T(:,:,2) = C(:,:,1) - C(:,:,9);
    T(:,:,3) = -2*C(:,:,5);
    T(:,:,4) = sqrt2*(C(:,:,2) + C(:,:,7));
    T(:,:,5) = sqrt2*(C(:,:,3) - C(:,:,8));
    T(:,:,6) = C(:,:,1) -2*C(:,:,4) + C(:,:,9);
    T(:,:,7) = sqrt2*(C(:,:,2) - C(:,:,7));
    T(:,:,8) = sqrt2*(C(:,:,3) + C(:,:,8));
    T(:,:,9) = 2*C(:,:,6);
end
T = 0.5*T;
