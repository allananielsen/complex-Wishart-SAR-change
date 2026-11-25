function out = mat2wc(in)

% mat2wc arranges output from fx randwc into upper-left reals
%
% L = 10;
% mat2wc(randwc(L,3,eye(3)))

% (c) Copyright 2025
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 Feb 2025

[p,p1,nel] = size(in);
if not(p==p1)
    error('mat2wc: input matrices must be square')
end

if p==3 % full/quad pol
    out = nan(nel,9);
    out(:,1) = real(in(1,1,:));
    out(:,2) = real(in(1,2,:));
    out(:,3) = imag(in(1,2,:));
    out(:,4) = real(in(1,3,:));
    out(:,5) = imag(in(1,3,:));
    out(:,6) = real(in(2,2,:));
    out(:,7) = real(in(2,3,:));
    out(:,8) = imag(in(2,3,:));
    out(:,9) = real(in(3,3,:));
elseif p==2 % dual pol
    out = nan(nel,4);
    out(:,1) = real(in(1,1,:));
    out(:,2) = real(in(1,2,:));
    out(:,3) = imag(in(1,2,:));
    out(:,4) = real(in(2,2,:));
elseif p==1 % scalar
    out = nan(nel,1);
    out(:,1) = in(1,1,:);
else
    warning('mat2wc: p must be 3, 2 or 1')
end
