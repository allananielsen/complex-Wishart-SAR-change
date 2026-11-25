function [coord,Sigma,pE,pQ,SARsim,wc,vz] = polsar_find_coord(SAR,varargin)

% polsar_find_coord finds col and row numbers of extreme pixels in polSAR data
%
% coord = polsar_find_coord(SAR)
% [coord,Sigma] = polsar_find_coord(SAR,13);
% [coord,Sigma] = polsar_find_coord(SAR,13);
%
% [coord,Sigma,pE,pQ,SARsim,wc,vz] = polsar_find_coord(vzSAR.eigvals,L,SAR);
%
% Input
%   SAR    - 
%   L      - (optional) ENL (if nargin>1 generate data and do change analysis over six
%            simulated time points with wishart_rand_patch and wcRjl)
%   SARori - (optional)
%            (1) set to anything, fx "1" for generation of 6 instances plus
%                change detection in a no change situation
%            (2) meant to go with input "SAR" as any 3-variable input,
%                fx eigvals from vz_zylx
%                if nargin>2 and if input "SAR" has 3 variables, input
%                original SAR scene here and generate patches from
%                that scene from coord found in the 3-variable input
%   1      - (optional) set to anything, fx 1 
%            (1) set to anything, fx "1" for generation of 6 instances plus
%                change detection in a no change situation if "SARori"is
%                specified
%
% Output
%   coord  - row and col numbers of interesting pixels, input to
%            wishart_rand_patch
%   Sigma  - covmats for simulated data, only if nargin>1
%   pE     - entropies for P from wcRjl, only if nargin>1
%   pQ     - -2*rho*lnQ from wcRjl, only if nargin>1
%   SARsim - the generated data, only if nargin>1
%   wc     - struct from wcRjl, only if nargin>1
%   vz     - van Zyl decomposition of (first realisation of) SARsim, only if nargin>1

% (c) Copyright 2025
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, http://people.compute.dtu.dk/alan
% 17-22 Jan 2025

L = NaN;
if nargin>1
    L = varargin{1};
end
Sigma = NaN;
pE = NaN;
pQ = NaN;
SARsim = NaN;
wc = NaN;
vz = NaN;

[nr,nc,nv] = size(SAR);
if nv==1
    fprintf('polsar_find_coord: input is %d x %d and has %d element\n',nr,nc,nv)
else
    fprintf('polsar_find_coord: input is %d x %d and has %d elements\n',nr,nc,nv)
end

ix = 1;

if nv==9

    coord = nan(9,2);
    % or
    % coord = nan(12,2);
    % or
    % coord = nan(3,2);

    span = SAR(:,:,1) + SAR(:,:,6) + SAR(:,:,9);
 
    if size(coord,1)==9 || size(coord,1)==12

        % max (relative) Shhhh
        aux = SAR(:,:,1) ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        
        % max (relative) Shvhv, i.e., (Pauli) volume scattering
        aux = SAR(:,:,6) ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        
        % max (relative) Svvvv
        aux = SAR(:,:,9) ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        
        % max det
        aux = wishart_det(SAR);
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        % min det
        idx = min(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        
        %%
        % % max tr
        % aux = SAR(:,:,1) + SAR(:,:,6) + SAR(:,:,9);
        % idx = max(aux(:));
        % [row(ix),col(ix)] = find(aux==idx);
        % ix = ix + 1;
        % % min tr
        % idx = min(aux(:));
        % [row(ix),col(ix)] = find(aux==idx);
        % ix = ix + 1;
        % 
        % % max ShhSvv* (volume, vegetation scattering)
        % aux1 = complex(SAR(:,:,4),SAR(:,:,5));
        % aux = aux1.*conj(aux1);
        % clear aux1
        % idx = max(aux(:));
        % [row(ix),col(ix)] = find(aux==idx);
        % ix = ix + 1;
        % % min ShhSvv*
        % idx = min(aux(:));
        % [row(ix),col(ix)] = find(aux==idx);
        % ix = ix + 1;
    
        % max (relative) |ShhSvv*|
        aux1 = complex(SAR(:,:,4),SAR(:,:,5));
        aux = abs(aux1) ./ span;
        % aux = aux1.*conj(aux1);
        clear aux1
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
    
        % % max norm of complex correlation
        % aux1 = complex(SAR(:,:,4),SAR(:,:,5));
        % aux = abs(aux1 ./ sqrt(SAR(:,:,1)) ./ sqrt(SAR(:,:,9)));
        % clear aux1
        % idx = max(aux(:));
        % [row(ix),col(ix)] = find(aux==idx);
        % ix = ix + 1;
    
        aux = (SAR(:,:,1) + SAR(:,:,6)) ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        aux = (SAR(:,:,1) + SAR(:,:,9)) ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
        aux = (SAR(:,:,6) + SAR(:,:,9)) ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;

    end

    if size(coord,1)==12 || size(coord,1)==3

        % max (relative) |Shh + Svv*|^2, i.e., Pauli single bounce (surface)
        aux = SAR(:,:,1) + SAR(:,:,9) + 2*SAR(:,:,4);
        aux = aux ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
    
        % max (relative) |Shh - Svv*|^2 i.e., Pauli double bounce
        aux = SAR(:,:,1) + SAR(:,:,9) - 2*SAR(:,:,4);
        aux = aux ./ span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;
    
        % max span (trace)
        aux = span;
        idx = max(aux(:));
        [row(ix),col(ix)] = find(aux==idx);
        ix = ix + 1;

    end

elseif nv==4

    coord = nan(9,2);
    
    trc = SAR(:,:,1) + SAR(:,:,4);

    % max (relative) Shvhv
    aux = SAR(:,:,1) ./ trc;
    idx = max(aux(:));
    [row(1),col(1)] = find(aux==idx);
    
    % max (relative) Svvvv
    aux = SAR(:,:,4) ./ trc;
    idx = max(aux(:));
    [row(2),col(2)] = find(aux==idx);
    
    % max det
    aux = wishart_det(SAR);
    idx = max(aux(:));
    [row(3),col(3)] = find(aux==idx);
    % min det
    idx = min(aux(:));
    [row(4),col(4)] = find(aux==idx);

    % max tr
    aux = trc;
    idx = max(aux(:));
    [row(5),col(5)] = find(aux==idx);
    % min tr
    idx = min(aux(:));
    [row(6),col(6)] = find(aux==idx);
    
    % max ShvSvv*
    aux1 = complex(SAR(:,:,2),SAR(:,:,3));
    aux = aux1.*conj(aux1);
    idx = max(aux(:));
    [row(7),col(7)] = find(aux==idx);
    
    % max (relative) ShvSvv*
    aux = aux ./ trc;
    idx = max(aux(:));
    [row(8),col(8)] = find(aux==idx);
    % min (relative) ShvSvv*
    idx = min(aux(:));
    [row(9),col(9)] = find(aux==idx);

elseif nv==3 % meant for van Zyl decomposition eigvals as input

    coord = nan(6,2);
    
    % max single bounce
    aux = SAR(:,:,1);
    idx = max(aux(:));
    [row(ix),col(ix)] = find(aux==idx);
    ix = ix + 1;
    % min double bounce
    idx = min(aux(:));
    [row(ix),col(ix)] = find(aux==idx);
    ix = ix + 1;

    % max double bounce
    aux = SAR(:,:,2);
    idx = max(aux(:));
    [row(ix),col(ix)] = find(aux==idx);
    ix = ix + 1;
    % min double bounce
    idx = min(aux(:));
    [row(ix),col(ix)] = find(aux==idx);
    ix = ix + 1;
    
    % max volume scattering or rotated dipoles
    aux = SAR(:,:,3);
    idx = max(aux(:));
    [row(ix),col(ix)] = find(aux==idx);
    ix = ix + 1;
    % min double bounce
    idx = min(aux(:));
    [row(ix),col(ix)] = find(aux==idx);
    ix = ix + 1;

    coord = [row' col'];

elseif nv==2

    coord = nan(4,2);

    % max Shvhv
    aux = SAR(:,:,1);
    idx = max(aux(:));
    [row(1),col(1)] = find(aux==idx);
    
    trc = SAR(:,:,1) + SAR(:,:,2);

    % max (relative) Shvhv
    aux = aux ./ trc;
    idx = max(aux(:));
    [row(2),col(2)] = find(aux==idx);
    
    % max Svvvv
    aux = SAR(:,:,2);
    idx = max(aux(:));
    [row(3),col(3)] = find(aux==idx);
    
    % max (relative) Svvvv
    aux = aux ./ trc;
    idx = max(aux(:));
    [row(4),col(4)] = find(aux==idx);

elseif nv==1

    coord = nan(2,2);

    % max
    idx = max(SAR(:));
    [row(1),col(1)] = find(SAR==idx);
    
    % min
    idx = min(SAR(:));
    [row(2),col(2)] = find(SAR==idx);

else
    error('input data must be 9, 4, 2 or 1 variable(s)')
end

clear aux trc

% coord = [col' row']; % no, use fliplr(coord) in call to wishart_rand_patch
coord = [row' col']; % yes, that's right

if nargin<2
    return
end

%% ENL is input so we simulate data and perform change detection

nrsim = 256;
ncsim = nrsim;
rng('default');

% below fliplr(coord) because ginput used in wishart_rand_patch
% flips rows and cols, so we must flip back again

if nv==9
    [simim,~,Sigma] = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    figure('Name','RGB dB 1/6/9')
    imshowrgbprctile(10*log10(cat(3,simim(:,:,1),simim(:,:,6),simim(:,:,9))));
    vz = van_zylx(simim);
    figure('Name','RGB van Zyl dB eigvals 2/3/1')
    imshowrgbprctile(10*log10(vz.eigvals),[2 3 1]);
    figure('Name','RGB van Zyl dB Lambdas 2/3/1')
    imshowrgbprctile(10*log10(vz.Lambda),[2 3 1]);
    if nargin<3
        return
    end
    s2 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s3 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s4 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s5 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s6 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    SARsim = cat(4,simim,s2,s3,s4,s5,s6);
    clear s2 s3 s4 s5 s6
    wc = wcRjl(SARsim,L);
    pQ = wc_lnRtable(wc,NaN)
    pE = wc_Ptable(wc,NaN)
    % matrix2latex(pE,6)
elseif nv==4
    [simim,~,Sigma] = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    aux = complex(simim(:,:,2),simim(:,:,3));
    figure('Name','RGB dB 4/1/(2,3)')
    imshowrgbprctile(10*log10(cat(3,simim(:,:,4),simim(:,:,1),aux.*conj(aux))));
    clear aux
    if nargin<3
        return
    end
    s2 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s3 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s4 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s5 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    s6 = wishart_rand_patch(L,SAR,size(coord,1),nrsim,ncsim,fliplr(coord));
    SARsim = cat(4,simim,s2,s3,s4,s5,s6);
    clear s2 s3 s4 s5 s6
    wc = wcRjl(SARsim,L);
    pQ = wc_lnRtable(wc,NaN)
    pE = wc_Ptable(wc,NaN)
    % matrix2latex(pE,6)
elseif nv==3 % run this on original SARori scene and not on 3-variable input,
             % fx van Zyl decomp
    if nargin>2
        SARori = varargin{2};
        [nrori,ncori,nvori] = size(SARori);
        if nr~=nrori || nc~=ncori
            error('SAR and SARori must have same nrows and ncols')
        end
        if nvori~=9
            error('SARori must have 9 variables')
        end
    end
    % only min eigvals (or whatever other 3-variable input)
    fprintf('only min eigvals from van Zyl decomposition (or other 3-variable input)\n')
    for ix = 2:2:6 % ix = 1:2:6
        wc2mat(SARori(coord(ix,1),coord(ix,2),:))
    end
    [simim,~,Sigma] = ...
        wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(2:2:end,:)));
    figure('Name','RGB dB 1/6/9 (min)')
    imshowrgbprctile(10*log10(cat(3,simim(:,:,1),simim(:,:,6),simim(:,:,9))));
    vz = van_zylx(simim);
    figure('Name','RGB van Zyl dB eigvals 2/3/1 (min)')
    imshowrgbprctile(10*log10(vz.eigvals),[2 3 1]);
    figure('Name','RGB van Zyl dB Lambdas 2/3/1 (min)')
    imshowrgbprctile(10*log10(vz.Lambda),[2 3 1]);
    % only max eigvals (or whatever other 3-variable input)
    fprintf('only max eigvals from van Zyl decomposition (or other 3-variable input)\n')
    for ix = 1:2:6 % ix = 1:2:6
        wc2mat(SARori(coord(ix,1),coord(ix,2),:))
    end
    [simim,~,Sigma] = ...
        wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(1:2:end,:)));
    figure('Name','RGB dB 1/6/9 (max)')
    imshowrgbprctile(10*log10(cat(3,simim(:,:,1),simim(:,:,6),simim(:,:,9))));
    vz = van_zylx(simim);
    figure('Name','RGB van Zyl dB eigvals 2/3/1 (max)')
    imshowrgbprctile(10*log10(vz.eigvals),[2 3 1]);
    figure('Name','RGB van Zyl dB Lambdas 2/3/1 (max)')
    imshowrgbprctile(10*log10(vz.Lambda),[2 3 1]);
    % % both max and min eigvals (or whatever other 3-variable input)
    % fprintf('both max and min eigvals from van Zyl decomposition (or other 3-variable input)\n')
    % for ix = 1:6
    %     wc2mat(SARori(coord(ix,1),coord(ix,2),:))
    % end
    % [simim,~,Sigma] = ...
    %     wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord));
    % figure('Name','RGB dB 1/6/9 (max and min)')
    % imshowrgbprctile(10*log10(cat(3,simim(:,:,1),simim(:,:,6),simim(:,:,9))));
    % vz = van_zylx(simim);
    % figure('Name','RGB van Zyl dB eigvals 2/3/1 (max and min)')
    % imshowrgbprctile(10*log10(vz.eigvals),[2 3 1]);
    % figure('Name','RGB van Zyl dB Lambdas 2/3/1 (max and min)')
    % imshowrgbprctile(10*log10(vz.Lambda),[2 3 1]);
    fprintf('\nrows and cols for max eigvals in coord(1:2:end,:)\n')
    fprintf('rows and cols for min eigvals in coord(2:2:end,:)\n\n')
    if nargin<4
        return
    end
    s2 = wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(1:2:end,:)));
    s3 = wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(1:2:end,:)));
    s4 = wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(1:2:end,:)));
    s5 = wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(1:2:end,:)));
    s6 = wishart_rand_patch(L,SARori,size(coord,1),nrsim,ncsim,fliplr(coord(1:2:end,:)));
    SARsim = cat(4,simim,s2,s3,s4,s5,s6);
    clear s2 s3 s4 s5 s6
    wc = wcRjl(SARsim,L);
    pQ = wc_lnRtable(wc,NaN)
    pE = wc_Ptable(wc,NaN)
    % matrix2latex(pE,6)
elseif nv==2 % generate 2x2 covmat and extract diagonal
    SAR4 = cat(3,SAR(:,:,1),zeros(nr,nc),zeros(nr,nc),SAR(:,:,2));
    [simim4,~,Sigma] = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    simim = simim4(:,:,[1 4]); clear simim4
    aux = simim(:,:,1) ./ simim(:,:,2);
    figure('Name','RGB dB 2/1/(1/2)')
    imshowrgbprctile(10*log10(cat(3,simim(:,:,2),simim(:,:,1),aux)));
    clear aux
    if nargin<3
        return
    end
    s24 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s2 = s24(:,:,[1 4]); clear s24
    s34 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s3 = s34(:,:,[1 4]); clear s34
    s44 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s4 = s44(:,:,[1 4]); clear s44
    s54 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s5 = s54(:,:,[1 4]); clear s54
    s64 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s6 = s64(:,:,[1 4]); clear s64
    SARsim = cat(4,simim,s2,s3,s4,s5,s6);
    clear s2 s3 s4 s5 s6
    wc = wcRjl(SARsim,L);
    pQ = wc_lnRtable(wc,NaN)
    pE = wc_Ptable(wc,NaN)
    % matrix2latex(pE,6)
elseif nv==1 % generate 2x2 covmat and extract first diagonal element
    SAR4 = cat(3,SAR,zeros(nr,nc),zeros(nr,nc),SAR);
    [simim4,~,Sigma4] = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    simim = simim4(:,:,1); clear simim4
    Sigma = Sigma4(1,1,:); clear Sigma4
    figure, imshow(simim,[])
    if nargin<3
        return
    end
    s24 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s2 = s24(:,:,1); clear s24
    s34 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s3 = s34(:,:,1); clear s34
    s44 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s4 = s44(:,:,1); clear s44
    s54 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s5 = s54(:,:,1); clear s54
    s64 = wishart_rand_patch(L,SAR4,size(coord,1),nrsim,ncsim,fliplr(coord));
    s6 = s64(:,:,1); clear s64
    SARsim = cat(4,simim,s2,s3,s4,s5,s6);
    clear s2 s3 s4 s5 s6
    wc = wcRjl(SARsim,L);
    pQ = wc_lnRtable(wc,NaN)
    pE = wc_Ptable(wc,NaN)
    % matrix2latex(pE,6)
end

while false

    % Sentinel-1 inspired data generation for Tech Rep

    nt = 60;
    L = 5;
    SAR = wishart_rand(L,100,1000,eye(2),nt); %,1);
    wcd = wc_Q(SAR(:,:,[1 4],:),L)
    wc = wc_Q(SAR,L)
    wcd1 = wc_Q(SAR(:,:,[1 4],:),L,1)
    wc1 = wc_Q(SAR,L,1)

    mv = 0.075;
    mh = 0.075;
    
    figure
    ax(1) = subplot_tight(2,2,1,[mv,mh]);
    wc_lnRhisto1(wc1,nan)
    ax(2) = subplot_tight(2,2,2,[mv,mh]);
    wc_lnRhisto1(wc,nan)
    ax(3) = subplot_tight(2,2,3,[mv,mh]);
    wc_Phisto1(wc1,nan)
    ax(4) = subplot_tight(2,2,4,[mv,mh]);
    wc_Phisto1(wc,nan)
    
    figure
    ax(1) = subplot_tight(2,2,1,[mv,mh]);
    wc_lnRhisto1(wcd1,nan)
    ax(2) = subplot_tight(2,2,2,[mv,mh]);
    wc_lnRhisto1(wcd,nan)
    ax(3) = subplot_tight(2,2,3,[mv,mh]);
    wc_Phisto1(wcd1,nan)
    ax(4) = subplot_tight(2,2,4,[mv,mh]);
    wc_Phisto1(wcd,nan)

    %%
    % test: varying nt, ENL, nt, p
    
    nt = 10; % 5 10 20 50 100
    L = 5; % 5 10 20 50
    p = 3; % 2 3

    nr = 100;
    nc = 100;
    
    SAR = wishart_rand(L,nr,nc,eye(p),nt); %,1);
    if p==2
        subset = [1 4];
    elseif p==3
        subset = [1 6 9];
    else
        error('p must be 2 or 3')
    end
    wcd = wc_Q(SAR(:,:,subset,:),L);
    wc = wc_Q(SAR,L);
    wcd1 = wc_Q(SAR(:,:,subset,:),L,1);
    wc1 = wc_Q(SAR,L,1);

    wc_Phisto(wcd,nan,'wcd ');
    wc_Phisto(wc,nan,'wcf ');
    wc_Phisto(wcd1,nan,'wcd1');
    wc_Phisto(wc1,nan,'wcf1');

end