function out = wc2mat(in)

% wc2mat arranges output from fx wc_rand in matrix form
%
% L = 10;
% wc2mat(wc_rand(L,1,1,eye(2))/L)
% [zpol,zpoldeg] = cart2polc(wc2mat(wc_rand(L,1,1,eye(3))/L))

% (c) Copyright 2021-2025
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 28 Apr 2021
%  3 Jan 2022 /AA: reflection symmetry added
%  8 Dec 2024 /AA: in as vector (ismatrix) added
% 28 Jan 2025 /AA: p=1 added

% if ismatrix(in)
%     [nr,nc] = size(in);
%     if ~or(nr==1,nc==1)
%         error('wc2mat: with vector input, nr or nc must be 1')
%     end
%     nel = max(nr,nc);
%     nr = 1;
%     nc = 1;
% else
    [nr,nc,nel] = size(in);
% end

if and(nr==1,nc==1)
    if nel==9 % full/quad pol
        out = nan(3,3);
        out(1,1) = in(1);
        out(1,2) = complex(in(2),in(3));
        out(1,3) = complex(in(4),in(5));
        out(2,2) = in(6);
        out(2,3) = complex(in(7),in(8));
        out(3,3) = in(9);
        out(2,1) = conj(out(1,2));
        out(3,1) = conj(out(1,3));
        out(3,2) = conj(out(2,3));
    elseif nel==5 % reflection symmetry
        out = nan(3,3);
        out(1,1) = in(1);
        out(1,3) = complex(in(2),in(3));
        out(1,2) = 0;
        out(2,2) = in(4);
        out(2,3) = 0;
        out(3,3) = in(5);
        out(2,1) = conj(out(1,2));
        out(3,1) = conj(out(1,3));
        out(3,2) = conj(out(2,3));
    elseif nel==4 % dual pol
        out = nan(2,2);
        out(1,1) = in(1);
        out(1,2) = complex(in(2),in(3));
        out(2,2) = in(4);
        out(2,1) = conj(out(1,2));
    elseif nel==1 % scalar
        out = nan;
        out(1,1) = in(1);
    else
        warning('wc2mat: numel(in) must be 9, 5, 4 or 1')
    end
else
    if nel==9
        out = nan(3,3,nr,nc);
        out(1,1,:,:) = in(:,:,1);
        out(1,2,:,:) = complex(in(:,:,2),in(:,:,3));
        out(1,3,:,:) = complex(in(:,:,4),in(:,:,5));
        out(2,2,:,:) = in(:,:,6);
        out(2,3,:,:) = complex(in(:,:,7),in(:,:,8));
        out(3,3,:,:) = in(:,:,9);
        out(2,1,:,:) = conj(out(1,2,:,:));
        out(3,1,:,:) = conj(out(1,3,:,:));
        out(3,2,:,:) = conj(out(2,3,:,:));
    elseif nel==5
        out = nan(3,3,nr,nc);
        out(1,1,:,:) = in(:,:,1);
        out(1,3,:,:) = complex(in(:,:,2),in(:,:,3));
        out(1,2,:,:) = 0;
        out(2,2,:,:) = in(:,:,4);
        out(2,3,:,:) = 0;
        out(3,3,:,:) = in(:,:,5);
        out(2,1,:,:) = conj(out(1,2,:,:));
        out(3,1,:,:) = conj(out(1,3,:,:));
        out(3,2,:,:) = conj(out(2,3,:,:));
    elseif nel==4
        out = nan(2,2,nr,nc);
        out(1,1,:,:) = in(:,:,1);
        out(1,2,:,:) = complex(in(:,:,2),in(:,:,3));
        out(2,2,:,:) = in(:,:,4);
        out(2,1,:,:) = conj(out(1,2,:,:));
    elseif nel==1
        out = nan(1,1,nr,nc);
        out(1,1,:,:) = in(:,:,1);
    else
        warning('wc2mat: numel(in) must be 9, 5, 4 or 1')
    end
end    
    