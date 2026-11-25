function [zpol,zpoldeg] = cart2polc(zcart)

% cart2polc transforms Cartesian data to polar coordinates
%
% zpol = cart2polc(zcart)
% [zpol,zpoldeg] = cart2polc(zcart)
% [zpol,zpoldeg] = cart2polc(complex(1,2))
% [zpol,zpoldeg] = cart2polc(wc2mat(wc_rand(5,1,1,eye(3))/5))
%
% real(zpol) is the radius/magnitude
% imag(zpol) is the argument in radians
%
% real(zpoldeg) is the radius/magnitude
% imag(zpoldeg) is the argument in degrees (0-360)

% (c) Copyright 2021
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 April 2021

th = atan2(imag(zcart),real(zcart));
r = hypot(imag(zcart),real(zcart));
zpol = complex(r,th);
if nargout==2
    th = 180*th/pi;
    zpoldeg = complex(r,th);
end
