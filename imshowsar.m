function range = imshowsar(img,varargin)

% imshowsar shows intensity SAR data, e.g., ShhShh* (real) or ShhSvv*
%   (real,imaginary); displays intensity (and phase if relevant)
%
% imshowsar(img,range)
%
% Input
%   img     - real image or (real,imaginary) image for display
%   titl    - title for plot    
%   range   - range for intensity, defaults to 1- and 99-percentile in dB
%
% imshowsar(img,[-36 0]) stretches intensity between -36 dB and 0 dB

% (c) Copyright 2012-2021
% Allan Aasbjerg Nielsen
% alan@dtu.dk, people.compute.dtu.dk/alan
% 6 Jan 2012 - 21 Sep 2013
% 24 Sep 2019 IH added
% 12 Mar 2021 percentiles added
% 26 Apr 2021 title added

realimage = false;

% range = [-36 -6]; % dB for magnitude of input
amp = 10*log10(abs(img));
titl = 'SAR';
if nargin>1
    titl = varargin{1};
end
range(1) = prctile(amp(:),1);
range(2) = prctile(amp(:),99);
if nargin>2
    range = varargin{2};
end

p = size(img);
if length(p)==3
	img = complex(img(:,:,1), img(:,:,2));
    amp = 10*log10(abs(img));
    range(1) = prctile(amp(:),1);
    range(2) = prctile(amp(:),99);
    if nargin>2
        range = varargin{2};
    end
elseif length(p)==2
    if isreal(img)
        realimage = true;
    end
end
    
figure('Name',titl)
imshow(amp, range)
cb = colorbar('SouthOutside');
xlabel(cb, 'dB')
drawnow
if realimage
    return
end
figure('Name',titl)
imshow(180*angle(img)/pi, [-180 180])
% colormap(hsv)
colormap(flip(circshift(hsv,85,1)))
cb = colorbar('SouthOutside');
set(cb, 'XTick', -180:90:180)
set(cb, 'XTickLabel', {'-180','-90','0','90','180'})
xlabel(cb, 'degrees')
drawnow
% imshow(angle(img), [-pi pi])
% colormap(hsv)
% cb = colorbar('SouthOutside')
% set(cb, 'XTick', -pi:pi/2:pi)
% set(cb, 'XTickLabel', {'-pi','-pi/2','0','pi/2','pi'})

if ~realimage
%     V = abs(img);
%     amp = 10*log10(abs(img));
%     amin = prctile(amp(isfinite(amp(:))),1);
%     amax = prctile(amp(isfinite(amp(:))),99);
    amin = range(1);
    amax = range(2);
    V = (amp-amin)/(amax-amin);
    % H = (phas-min(phas(:)))/(max(phas(:)-min(phas(:))));
    % H = (phas-min(phas(:)))/12.0;
    H = (angle(img)+pi)/2/pi;
    figure('Name',titl)
    imshow(hsv2rgbAA(cat(3,H,ones(p(1),p(2)),V)))
%     imshow(hsv2rgb(cat(3,H,ones(p(1),p(2)),V)))
    % colormap(hsv)
%     colormap(hsv(256))
    colormap(flip(circshift(hsv,85,1)))
    cb = colorbar('SouthOutside');
    set(cb, 'XTick', ((-pi:(pi/2):pi)+pi)/2/pi)
    set(cb, 'XTickLabel', {'-180','-90','0','90','180'})
    xlabel(cb, 'degrees')
    drawnow
end
