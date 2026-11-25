function range = imshowrgbprctile(img,varargin)

%
% imshowrgbprctile(img,rgb,prctl,magn)
%
% img   - 3-D (or 2-D) matrix
% rgb   - order of bands shown as R, G and B, e.g. [4 -2 1] (or one band);
%         optional, defaults to [1 2 3]
%         if <0: negate band
% prctl - percentiles for stretch in %, e.g., [2 98]
% magn  - magnification in %

% (c) Copyright 2013
% Allan Aasbjerg Nielsen
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 21 Sep 2013

if nargin < 1, error('imshowrgbprctile: wrong input specification'), end
% if ndims(img)~=3, error('imshowrgbprctile: input image must be 3-D'); end

rgb   = [1 2 3]; if nargin > 1, rgb   = varargin{1}; end
prctl = [1 99];  if nargin > 2, prctl = varargin{2}; end
magn  = 100;     if nargin > 3, magn  = varargin{3}; end
[nr,nc,nf] = size(img);
sizergb = size(rgb,2);
if nf == 1
    sizergb = 1;
end
if sizergb == 1 % only one band specified
    img = cast(img, 'single');
    sgn = sign(rgb);
    img = img(:, :, abs(rgb));
    minimg = prctile(img(:),prctl(1));
    maximg = prctile(img(:),prctl(2));
    range = [minimg maximg];
    r = (img-minimg)/(maximg-minimg);
    if sgn==-1, r = 1-r; end
    r(r<0)=0; r(r>1)=1;
    imshow(r, 'InitialMagnification', magn)
elseif sizergb == 3 % 3-band image
    img = cast(reshape(img, nc*nr, nf), 'single');
    sgn = sign(rgb);
    img = img(:, abs(rgb));
    minimg = prctile(img,prctl(1));
    maximg = prctile(img,prctl(2));
    range = [minimg; maximg];
%     fprintf('range for RGB %g %g %g\n',minimg)
%     fprintf('              %g %g %g\n',maximg)
    img = reshape(img, nr, nc, 3);
%     nf = 3;
    r = (img(:,:,1)-minimg(1))/(maximg(1)-minimg(1));
    g = (img(:,:,2)-minimg(2))/(maximg(2)-minimg(2));
    b = (img(:,:,3)-minimg(3))/(maximg(3)-minimg(3));
    if sgn(1)==-1, r=1-r; end
    if sgn(2)==-1, g=1-g; end
    if sgn(3)==-1, b=1-b; end
    r(r<0)=0; r(r>1)=1;
    g(g<0)=0; g(g>1)=1;
    b(b<0)=0; b(b>1)=1;
    imshow(cat(3,r,g,b), 'InitialMagnification', magn)
else
    error('wrong specification of rgb');
end
