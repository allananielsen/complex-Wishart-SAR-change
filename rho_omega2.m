function [rho,omega2,df,check] = rho_omega2(nv,nsets,ENL)

% calculate rho, omega2 and df for wcRjl (polSAR change detection)
%
% [rho,omega2,df] = rho_omega2(nv,k,ENL);
% [rho,omega2,df,check] = rho_omega2(nv,nsets,ENL)
%
% Input
%   nv      - number of elements in covmat
%               1 single channel power data
%               2 dual pol diag only
%               3 full pol diag only
%               4 dual pol full matrix
%               5 full pol reflection symmetry
%               9 full pol
%   k       - number of time points
%   ENL     - equivalent number of looks

% (c) Copyright 2023
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 16 Oct- 10 Nov 2023
    
if nargin<1 %| strcmp(varargin{2},"help")
    help rho_omega2
    return
end

k = nsets;
n = ENL;

if nv==1
    p = 1;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1); 
elseif nv==2
    p1 = 1;
    p2 = 1;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
elseif nv==3
    p1 = 1;
    p2 = 1;
    p3 = 1;
    p = p1+p2+p3;
    f = p1^2+p2^2+p3^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1) + p3*(2*p3^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1) + p3^2*(p3^2-1);
elseif nv==4
    p = 2;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1); 
elseif nv==5
    p1 = 1;
    p2 = 2;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
elseif nv==9
    p = 3;
    f = p^2;
    sump = p*(2*p^2-1);
    sump2 = p^2*(p^2-1); 
elseif nv==14 || nv==8 % this is for 2 blocks with 2x2 full matrices ('4' is taken hence '14')
    p1 = 2;
    p2 = 2;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
elseif nv==16 || nv==18 % this is for 2 blocks with 3x3 full matrices 
    p1 = 3;
    p2 = 3;
    p = p1+p2;
    f = p1^2+p2^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1);
elseif nv==19 % this is for 3 blocks with 3x3 full matrices ('9' is taken hence '19')
    p1 = 3;
    p2 = 3;
    p3 = 3;
    p = p1+p2+p3;
    f = p1^2+p2^2+p3^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1) + p3*(2*p3^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1) + p3^2*(p3^2-1);
elseif nv==26 % this is for 3 blocks with 2x2 full matrices ('16' is taken hence '26')
    p1 = 2;
    p2 = 2;
    p3 = 2;
    p = p1+p2+p3;
    f = p1^2+p2^2+p3^2;
    sump = p1*(2*p1^2-1) + p2*(2*p2^2-1) + p3*(2*p3^2-1);
    sump2 = p1^2*(p1^2-1) + p2^2*(p2^2-1) + p3^2*(p3^2-1);
else
    error('*** rho_omega2: wrong nv ***')
end

if nsets==2; nsets = 1; end

df = nan(nsets,1);
df(1) = (k-1)*f;
for j=2:nsets
    df(j) = f;
end

rho = nan(nsets,1);
omega2 = nan(nsets,1);
rho(1) = 1 - (k/n-1/(n*k))/(6*df(1))*sump;
omega2(1) = -df(1)/4*(1-1/rho(1))^2 + ...
    (k/n^2-1/n^2/k^2)*sump2/(24*rho(1)^2);
for j=2:nsets
    rho(j) = 1 - (1+1/(j*(j-1)))/(6*df(j)*n)*sump;
    omega2(j) = -f/4*(1-1/rho(j))^2 + ...
        (1+(2*j-1)/(j^2*(j-1)^2))*sump2/(24*n^2*rho(j)^2);
end

if nargout>3
    % Check
    check = [rho,omega2,df];
    if nv==2
        df(1) = 2*(k-1);
        rho(1) = 1 - 1/(3*df(1))*(k/n-1/(n*k));
        omega2(1) = -df(1)/4*(1-1/rho(1))^2;
        for j=2:nsets
            df(j) = 2;
            rho(j) = 1 - (1+1/(j*(j-1)))/(6*n);
            omega2(j) = -df(j)/4*(1-1/rho(j))^2;
        end
    elseif nv==3
        df(1) = 3*(k-1);
        rho(1) = 1 - 1/(2*df(1))*(k/n-1/(n*k));
        omega2(1) = -df(1)/4*(1-1/rho(1))^2;
        for j=2:nsets
            df(j) = 3;
            rho(j) = 1 - (1+1/(j*(j-1)))/(6*n);
            omega2(j) = -df(j)/4*(1-1/rho(j))^2;
        end
    end
    check = check - [rho,omega2,df];
    if any(check,'all')
        warning('*** rho_omega2: something is wrong ***')
    end
end

return

nL = 30;
p_rof = nan(nL,3);
for L=1:nL
    [aarho,aaomega2,aaf] = rho_omega2(18,4,L);
    % [aarho,aaomega2,aaf] = rho_omega2(4,12,L);
    p_rof(L,1) = aarho(1);
    p_rof(L,2) = aaomega2(1);
    p_rof(L,3) = aaf(1);
    % title(titl(L),"FontSize",16,"FontWeight","bold")
    titl(1) = "\rho    ";
    titl(2) = "\omega_2";
end

figure
for ix=1:2
    ax(ix) = subplot_tight(2,1,ix,[0.1,0.05]);
    % plot(1:nL,p_rof(:,ix),"LineWidth",2)
    plot(2:nL,p_rof(2:end,ix),"LineWidth",2)
    grid on
    title(titl(ix),"FontSize",16,"FontWeight","bold")
end
% linkaxes(ax,'xy')

titl(3) = "df";
figure
for ix=1:3
    ax(ix) = subplot_tight(3,1,ix,[0.1,0.05]);
    % plot(1:nL,p_rof(:,ix),"LineWidth",2)
    plot(2:nL,p_rof(2:end,ix),"LineWidth",2)
    grid on
    title(ax(ix),titl(ix),"FontSize",16)
end
% linkaxes(ax,'xy')

%%%%%%%%%%%

while false
    % Check #1 - values for testing whether one matrix is block diagonal
    if nv==4 % dual pol: testing whether matrix is diagonal 
        p1 = 1;
        p2 = 1;
    elseif nv==9 % full pol: testing whether matrix is reflection symmetric
        p1 = 1;
        p2 = 2;
    end
    p = p1 + p2;
    sump2 = p^2 - (p1^2+p2^2);
    sump3 = p^3 - (p1^3+p2^3);
    sump4 = p^4 - (p1^4+p2^4);
    df(1) = sump2;
    rho(1) = 1 - 1/(3*n)*sump3/sump2;
    % omega2(1) = ( - sump3^2/(36*sump2) + sump4/24 )/(n^2*rho(1)^2)^2;
    omega2(1) = ( - sump3^2/(36*sump2) + sump4/24 )/(n^2*rho(1)^2);
end
