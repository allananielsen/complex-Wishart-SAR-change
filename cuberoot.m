function croot = cuberoot(z,varargin)

% cuberoot calculates three complex cube roots of (complex) numbers
%
% Input
%   z(:,:)       - complex number in cartesian coordinates (real + 1i*imag)
%   varargin     - any input causes only two cube roots to be calculated
%                  (the three cube roots sum to zero)
%
% Output
%   croot(:,:,1:3) - cuberoots in cartesian coordinates
% or
%   croot(:,:,1:2) - cuberoots sum to zero, do not calculate the third

% (c) Copyright 2022
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 Feb 2022

if ~ismatrix(z)
    error('*** cuberoot: input must be matrix ***')
end
[nr,nc] = size(z);

if nargin==1
    croot = nan(nr,nc,3);
else
    croot = nan(nr,nc,2);
end

zpol = cart2polc(z); % mag in real part, arg in imag part of zpol
zpol_mag = nthroot(real(zpol),3);
% zpol_arg = imag(zpol)/3;
croot(:,:,1) = pol2cartc(complex(zpol_mag, imag(zpol)/3));
croot(:,:,2) = pol2cartc(complex(zpol_mag,(imag(zpol)-2*pi)/3));
if nargin==1
    %  croot(:,:,3) = pol2cartc(complex(zpol_mag,(imag(zpol)-4*pi)/3));
    croot(:,:,3) = -sum(croot(:,:,1:2),3);
end
