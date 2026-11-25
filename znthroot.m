function zroot = znthroot(z,n,varargin)

% znthroot calculates n complex nth roots of (complex) numbers
%
% Input
%   z(:,:)       - complex number in cartesian coordinates (real + 1i*imag)
%   n            - n for nth root
%   varargin     - any input causes only (n-1) roots to be calculated
%                  (the n nth roots sum to zero)
%
% Output
%   zroot(:,:,1:n)     - nth roots in cartesian coordinates
% or
%   zroot(:,:,1:(n-1)) - nth roots sum to zero, do not calculate the last

% (c) Copyright 2022
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 Feb 2022

if ~ismatrix(z)
    error('*** znthroot: input must be matrix ***')
end
[nr,nc] = size(z);

if nargin==1
    error('znthroot: you must input n for nth root')
end
if nargin==2
    zroot = nan(nr,nc,n);
end
if nargin>2
    zroot = nan(nr,nc,n-1);
end

zpol = cart2polc(z); % mag in real part of zpol, arg in imag part
zpol_mag = nthroot(real(zpol),n);
% zpol_arg = imag(zpol)/n;
% figure('Name','nth complex root')
for ix=1:(n-1)
    zroot(:,:,ix) = pol2cartc(complex(zpol_mag,(imag(zpol)-(ix-1)*2*pi)/n));
%     polarplot(zroot(1,1,ix),'-x','LineWidth',2,'MarkerSize',10), hold on
end
if nargin==2
    % zroot(:,:,n) = pol2cartc(complex(zpol_mag,(imag(zpol)-(n-1)*2*pi)/n));
    zroot(:,:,n) = -sum(zroot(:,:,1:(n-1)),3);
%     polarplot(zroot(1,1,n),'-x','LineWidth',2,'MarkerSize',10)
end
