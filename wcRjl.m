function wc = wcRjl(X,L,varargin)

% Wishart based change detection in time series of k polarimetric SAR images
%
% SAR input images are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%   ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% wc = wcRjl(cat(4,L62,L63,L64,L65,L68,L74),13);
% wc = wcRjl('dir',13);
%
% wcRjl in the first example runs wc_Rj for the combinations
%   cat(4,L62,L63,L64,L65,L68,L74)  - result in wc(1)
%   cat(4,L63,L64,L65,L68,L74)      - result in wc(2)
%   cat(4,L64,L65,L68,L74)          - ...
%   cat(4,L65,L68,L74)
%   cat(4,L68,L74)
% and outputs the result to structure wc as wc(1), wc(2), ..., wc(5)
%
% In the second example wcRjl reads all files associated with files named
% '*.hdr' in directory 'dir' (for example '.'),
% concatenates them as in the first example
% followed by the same analysis (file names must be such that the
% concatenation happens in the right temporal order)
%
% Input
%   X       - directory name (in which case we read all files associated
%             with files named '*.hdr')
%             or
%             series of k 1x1, 2x2 or 3x3 complex covariance matrices
%             (polarimetric SAR images in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%   L        - k (or one) number of looks for X
%   (pol     - polarization, 'full', ..., - obsolete, retained for compatibility)
%   1/0      - do not output sumMat and logDet (optional, true or 1 is default)
%   out4byte - output in 4-byte (and not 8-byte) precision (optional, false or 0 is default)
%
% Output
%   wc      - structure with
%             lnR     - logarithm of test quantities Q and Rj
%                       ln(Q) is lnR(:,:,1), ln(R_2) is lnR(:,:,2) etc.
%             f       - number of degrees of freedom for Q and Rj test
%                       statistics
%             rho     - aux variable
%             omega2  - aux variable
%             P       - probabilities of finding smaller value of
%                       -2*rho(i)*lnR(:,:,i), i=1,...,k
%             sumMat  - updated sum of matrices
%             logDet  - logarithm of determinants of sumMat
%             L       - equivalent number of looks, ENL (from input)
%
% References
%
% Knut Conradsen, Allan Aasbjerg Nielsen and Henning Skriver (2016).
% "Determining the Points of Change in Time Series of Polarimetric SAR Data".
% IEEE Transactions on Geoscience and Remote Sensing 54(5), 3007-3024.
% https://www2.imm.dtu.dk/pubdb/p.php?6825
% https://orbit.dtu.dk/en/publications/determining-the-points-of-change-in-time-series-of-polarimetric-sar-data(4255ddb4-1fd9-4d66-9aad-e8049a976d75).html
% DOI:10.1109/TGRS.2015.2510160
%
% Allan A. Nielsen, Knut Conradsen, Henning Skriver and Morton J. Canty (2017).
% "Visualization of and software for omnibus test based change detected in
% a time series of polarimetric SAR data".
% Canadian Journal of Remote Sensing 43(6), 582-592.
% https://www2.imm.dtu.dk/pubdb/p.php?id=7027
% DOI:10.1080/07038992.2017.1394182
%
% Allan A. Nielsen, Morton J. Canty, Henning Skriver and Knut Conradsen (2019).
% "Cloud based spatio-temporal analysis of change in sequences of Sentinel
% images." ESA Conference on Big Data from Space (BiDS), 105-108.
% Munich, Germany, 19-21 February 2019.
% https://www2.imm.dtu.dk/pubdb/p.php?id=7129
% DOI:10.2760/848593

% (c) Copyright 2014-2023
% Allan Aasbjerg Nielsen, PhD
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 29 Oct 2014 - 11 Aug 2017 - 11 Apr 2019
% 10 Mar 2021 /AA: wc.L added
% 22 Nov 2021 /AA: now uses wc_Rj (which outputs a struct) and not wishart_change_Rj
%                  now handles two 3x3 or two 2x2 diagonal blocks
%                  cleanup of input outputMat and out4byte

% output 8-byte wc with ln Q/Rj, and sumMat/logDet depending on options
outputMat = true;
% to save space, output 4-byte wc with no ln Q/Rj and no sumMat/logDet
out4byte = false;

if nargin<2 || nargin>5
    help wcRjl
    error('wcRjl: wrong input');
end
if nargout>1
    help wcRjl
    error('wcRjl: wrong output');
end
if nargin>2 % output sumMat and logDet?
    if ischar(varargin{1}) % polarization (obsolete, retained for compatibility)
        if nargin>3
            outputMat = varargin{2};
        end
        if nargin>4
            out4byte = varargin{3}; % output 4 byte data
        end
    else
        outputMat = varargin{1};
        if nargin>3
            out4byte = varargin{2}; % output 4 byte data
        end
    end
end
if out4byte
    outputMat = false;
end

% if ischar(X) % read all files in directory 'X' associated with file names '*.hdr'
if ischar (X) && isfolder(X) % read all files in directory 'X' associated with file names '*.hdr'
    directory = X;
    clear X
    F = dir(strcat(directory,'/*.hdr'));
    fname = strcat(directory,'/',F(1).name(1:end-4));
    [~,aux]  = freadenvisar(fname,'b','true'); % read hdr only
    if aux(3)==6, dd = 9; end
    if aux(3)==4, dd = 5; end
    if aux(3)==3, dd = 4; end
    X = nan(aux(1),aux(2),dd,length(F));
    for ii=1:length(F)
        fname = strcat(directory,'/',F(ii).name(1:end-4));
        fprintf('wcRjl: file %s\n',fname)
        X(:,:,:,ii) = freadenvisar(fname,'b');
    end
    clear F ii fname directory Y
end

pX = size(X);
nv = pX(3);
if nv==18
    fprintf('*** wcRjl: input is two 3x3 diagonal blocks in 6x6 matrix ***\n')
elseif nv==8
    fprintf('*** wcRjl: input is two 2x2 diagonal blocks in 4x4 matrix ***\n')
end
nsets = pX(4);
if nsets<2
    error('wcRjl: one set only?')
end
% k = nsets;

wcnew = struct('lnR',[],'f',nan,'rho',nan,'omega2',nan,'P',[],...
            'sumMat',[],'logDet',[],'L',nan);
for i=1:(nsets-1)
    wc(i) = wcnew;
end
clear wcnew

if out4byte % to save space, output 4-byte, no ln Q/Rj, no sumMat/logDet
    for i=1:(nsets-1)
        wci = wc_Rj(X(:,:,:,i:nsets),L);
        % wc(i).lnR = single(wci.lnR);
        wc(i).f = single(wci.f);
        wc(i).rho = single(wci.rho);
        wc(i).omega2 = single(wci.omega2);
        wc(i).P = single(wci.P);
        wc(i).L = single(wci.L);
        clear wci
    end
else % output 8-byte
    for i=1:(nsets-1)
            wc(i) = wc_Rj(X(:,:,:,i:nsets),L);
    end
end

if outputMat
    if nv==18 || nv==8
        wcnew = struct('lnR',[],'f',nan,'rho',nan,'omega2',nan,'P',[],...
            'sumMat',L*X(:,:,:,nsets),...
            'logDet',log(wishart_det6(L*X(:,:,:,nsets))),'L',nan);
    else
        wcnew = struct('lnR',[],'f',nan,'rho',nan,'omega2',nan,'P',[],...
            'sumMat',L*X(:,:,:,nsets),...
            'logDet',log(wishart_det(L*X(:,:,:,nsets))),'L',nan);
    end
    wc = [wc,wcnew];
    clear wcnew 
else
    wc = rmfield(wc,{'sumMat','logDet'});
end

if out4byte
    wc = rmfield(wc,'lnR');
end
