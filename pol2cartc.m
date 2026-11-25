function zcart = pol2cartc(zpol,varargin)
%POL2CART Transform polar to Cartesian coordinates.
%   [X,Y] = POL2CART(TH,R) transforms corresponding elements of data
%   stored in polar coordinates (angle TH, radius R) to Cartesian
%   coordinates X,Y.  The arrays TH and R must the same size (or
%   either can be scalar).  TH must be in radians.
%
%   [X,Y,Z] = POL2CART(TH,R,Z) transforms corresponding elements of
%   data stored in cylindrical coordinates (angle TH, radius R, height
%   Z) to Cartesian coordinates X,Y,Z. The arrays TH, R, and Z must be
%   the same size (or any of them can be scalar).  TH must be in radians.
%
%   Class support for inputs TH,R,Z:
%      float: double, single
%
%   See also CART2SPH, CART2POL, SPH2CART.

%   L. Shure, 4-20-92.
%   Copyright 1984-2004 The MathWorks, Inc. 
%
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 Feb 2022

% Input
%   zpol - magnitude in real part, argument in imaginary part

if nargin>1 % argument in degrees
    imag_zpol = pi*imag(zpol)/180;
else
    imag_zpol = imag(zpol);
end
zcart = real(zpol) .* (cos(imag_zpol) + 1i*sin(imag_zpol));
