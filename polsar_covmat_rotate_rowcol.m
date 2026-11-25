function SARrot = polsar_covmat_rotate_rowcol(SAR,theta,varargin)

% Rotate full/quad polSAR data (slow for-loop implementation)
%
% Input
%   SAR   - covariance matrix representation, must be 3-D with
%             nine (9) (quad/full polarimetry) bands,
%             order is alphabetical in indexes with
%             real parts before imaginary parts (as read by
%             freadenvisar), for full polarimetry:
%               hhhh, re(hhhv), im(hhhv), re(hhvv), im(hhvv),
%               hvhv, re(hvvv), im(hvvv), vvvv
%   theta - rotation angle in degrees (theta from rho_llrr)
%   bias  - constant bias angle in degrees - optional, defaults to 0
%
%   Covariance matrix representation
%
%   ShhShh* ShhShv* ShhSvv*     k    a   rho
%   ShvShh* ShvShv* ShvSvv*  =  a*   ksi b
%   SvvShh* SvvShv* SvvSvv*     rho* b*  zeta%
%
% See also rho_llrr, rho_rrll
%
% References
%
% J.-S. Lee and D. L. Schuler and T. L. Ainsworth and E. Krogager and
% D. Kasilingam and W.-M. Boerner
% "On the estimation of radar polarization orientation shifts induced by
% terrain slopes"
% IEEE Transactions on Geoscience and Remote Sensing 40(1), 30-41, Jan. 2002
% doi:10.1109/36.981347
%
% T. L. Ainsworth, D. L. Schuler and J.-S. Lee
% "Polarimetric SAR characterization of man-made structures in urban areas
% using normalized circular-pol correlation coefficients"
% Remote Sensing of Environment, 112, 2876-2885, 2008
% doi:10.1016/j.rse.2008.02.005
%
% M. Sugimoto, K. Ouchi and Y. Nakamura
% "Four-Component Scattering Power Decomposition Algorithm with Rotation
% of Covariance Matrix Using ALOS-PALSAR Polarimetric Data",
% Remote Sensing 2012, 4, 2199-2209
% doi:10.3390/rs4082199
%
% Allan A. Nielsen
% "Rotation of Polarimetric SAR Covariance Matrices",
% IEEE International Geoscience and Remote Sensing Symposium, IGARSS 2023,
% Pasadena, California, USA, 16-21 July 2023.
% https://orbit.dtu.dk/en/publications/rotation-of-polarimetric-sar-covariance-matrices

% (c) Copyright 2022
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 4 Mar 2022

if nargin>2
    bias = varargin{1};
    theta = theta + bias;
end

pX = size(SAR);
nrows = pX(1);
ncols = pX(2);
if  ndims(SAR)==3
    nvars = pX(3);
elseif ismatrix(SAR)
    nvars = 1;
else
    error('polsar_covmat_rotate: wrong size of input images');
end

if nvars~=9 % full/quad pol
    error('polsar_covmat_rotate: input must be full/quad pol')
end

s2 = sqrt(2);
SARrot = nan(size(SAR));
h = waitbar(0,'Slow implementation with for-loops over rows and cols...', ...
    'name','PolSAR rotate');
for ic=1:ncols
    waitbar(ic/ncols,h);
    for ir=1:nrows
        c2t = cosd(2*theta(ir,ic));
        s2t = sind(2*theta(ir,ic));
        U = 0.5*[  1+c2t  s2*s2t  1-c2t
                 -s2*s2t   2*c2t s2*s2t
                   1-c2t -s2*s2t  1+c2t];
%         if ~isreal(U), fprintf('U complex, ir=%i, ic=%i\n',ir,ic), end
        aux = U*wc2mat(SAR(ir,ic,:))*U.';
        SARrot(ir,ic,1) = real(aux(1,1));
        SARrot(ir,ic,2) = real(aux(1,2));
        SARrot(ir,ic,3) = imag(aux(1,2));
        SARrot(ir,ic,4) = real(aux(1,3));
        SARrot(ir,ic,5) = imag(aux(1,3));
        SARrot(ir,ic,6) = real(aux(2,2));
        SARrot(ir,ic,7) = real(aux(2,3));
        SARrot(ir,ic,8) = imag(aux(2,3));
        SARrot(ir,ic,9) = real(aux(3,3));
    end
end
close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while false

% C = sym('C',3);
syms C 3
% assumptions(C)
% U = sym('U',3,'real');
syms U 3 real
% assumptions(U)

% c2v = cos(2*v);
% s2v = sin(2*v);
% U = 0.5*[  1+c2t  s2*s2t  1-c2t
%          -s2*s2t   2*c2t s2*s2t
%            1-c2t -s2*s2t  1+c2t];

% v = sym('v','real')
syms v real
U(1,1) = 1+cos(2*v);
U(2,1) = -sqrt(2)*sin(2*v);
U(3,1) = 1-cos(2*v);
U(1,2) = sqrt(2)*sin(2*v);
U(2,2) = 2*cos(2*v);
U(3,2) = -sqrt(2)*sin(2*v);
U(1,3) = 1-cos(2*v);
U(2,3) = sqrt(2)*sin(2*v);
U(3,3) = 1+cos(2*v);
U % = 0.5*U

Cr = U*C*U.'/4;
Cr22 = Cr(2,2)
Cr22es = simplify(expand(Cr22))
grad22 = gradient(Cr(2,2),v);
grad22s = simplify(grad22)
% grad22c = collect(grad22)
Hess22 = gradient(grad22s,v)
% eqn = df_dv == 0
% % fplot([lhs(eqn) rhs(eqn)],[-pi pi]) % nonsense
S1 = solve(grad22s == 0, v, 'ReturnConditions', true)
simplify(S1.v)

% Cr22es =
% C1_1/4 - C1_3/4 + C2_2/2 - C3_1/4 + C3_3/4 ...
%     - (C1_1*cos(4*v))/4 + (C1_3*cos(4*v))/4 + (C2_2*cos(4*v))/2 + (C3_1*cos(4*v))/4 - (C3_3*cos(4*v))/4 ...
%     - (2^(1/2)*C1_2*sin(4*v))/4 - (2^(1/2)*C2_1*sin(4*v))/4 + (2^(1/2)*C2_3*sin(4*v))/4 + (2^(1/2)*C3_2*sin(4*v))/4
% = 0.25*( C1_1 - C1_3 + 2*C2_2 - C3_1 + C3_3) ...
%             + 0.25*cos(4*v)*( - C1_1 + C1_3 + 2*C2_2 + C3_1 - C3_3 ) ...
%     + sqrt(2)*0.25*sin(4*v)*( - C1_2 - C2_1 + C2_3 + C3_2 )
% = 0.25*( C1_1 - 2*real(C1_3) + 2*C2_2 + C3_3) ...
%             + 0.25*cos(4*v)*( - C1_1 + 2*real(C1_3) + 2*C2_2 - C3_3 ) ...
%     - sqrt(2)*0.5*sin(4*v)*( real(C1_2) - real(C2_3) )

% grad22s =
% C1_1*sin(4*v) - C1_3*sin(4*v) - 2*C2_2*sin(4*v) - C3_1*sin(4*v) + C3_3*sin(4*v) ...
%   - 2^(1/2)*C1_2*cos(4*v) - 2^(1/2)*C2_1*cos(4*v) + 2^(1/2)*C2_3*cos(4*v) + 2^(1/2)*C3_2*cos(4*v)
% = sin(4*v)*( C1_1 - 2*C2_2 + C3_3 - 2*real(C1_3) )
% - cos(4*v)*2*sqrt(2)*( real(C1_2) - real(C2_3) ) 

% % Hess22 = 4*C1_1*cos(4*v) - 4*C1_3*cos(4*v) - 8*C2_2*cos(4*v) - 4*C3_1*cos(4*v) + 4*C3_3*cos(4*v) ...
% % + 4*2^(1/2)*C1_2*sin(4*v) + 4*2^(1/2)*C2_1*sin(4*v) - 4*2^(1/2)*C2_3*sin(4*v) - 4*2^(1/2)*C3_2*sin(4*v)
% % Hess22 = 4*cos(4*v)*( C1_1 - C1_3 - 2*C2_2 - C3_1 + C3_3 )  ...
% %     + 4*sqrt(2)*sin(4*v)*( C1_2 + C2_1 - C2_3 - C3_2);
% % Hess22 = 4*cos(4*v)*( C1_1 - 2*real(C1_3) - 2*C2_2 + C3_3 )  ...
% %     + 4*sqrt(2)*sin(4*v)*( 2*real(C1_2) - 2*real(C2_3) );
% Hess22 = 4*cosd(4*w)*( SAR(row,col,1) - 2*SAR(row,col,4) - 2*SAR(row,col,6) + SAR(row,col,9) )  ...
%     + 8*sqrt(2)*sind(4*w)*( SAR(row,col,2) - SAR(row,col,7) );

%%%%%%%%%%

syms S 2
syms R 2
syms v real
R(1,1) =  cos(v);
R(2,1) =  sin(v);
R(1,2) = -sin(v);
R(2,2) =  cos(v)

syms Sr 2
% Shh = S(1,1)
% Shv = S(1,2)
% % Svh = S(2,1)
% Svv = S(2,2)
Sr = R.'*S*R

syms C 3
C(1,1) =         S(1,1)*S(1,1)';
C(1,2) = sqrt(2)*S(1,1)*S(1,2)';
C(1,3) =         S(1,1)*S(2,2)';
C(2,2) =       2*S(1,2)*S(1,2)';
C(2,3) = sqrt(2)*S(1,2)*S(2,2)';
C(3,3) =         S(2,2)*S(2,2)';
C(2,1) = C(1,2)';
C(3,1) = C(1,3)';
C(3,2) = C(2,3)'

syms Cr1 3
Cr1(1,1) =         Sr(1,1)*Sr(1,1)';
Cr1(1,2) = sqrt(2)*Sr(1,1)*Sr(1,2)';
Cr1(1,3) =         Sr(1,1)*Sr(2,2)';
Cr1(2,2) =       2*Sr(1,2)*Sr(1,2)';
Cr1(2,3) = sqrt(2)*Sr(1,2)*Sr(2,2)';
Cr1(3,3) =         Sr(2,2)*Sr(2,2)';
Cr1(2,1) = Cr1(1,2)';
Cr1(3,1) = Cr1(1,3)';
Cr1(3,2) = Cr1(2,3)'


%%%%%%%%%%

syms a b c d z
S2 = solve(a*z^2 + b*z + c == 0, z)
S3 = solve(a*z^3 + b*z^2 + c*z + d == 0, z)
 
end % while false
