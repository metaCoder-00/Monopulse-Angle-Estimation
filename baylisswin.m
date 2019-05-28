function w = baylisswin(varargin)
%BAYLISS   Bayliss window.
%
% BAYLISSWIN(N) returns an N-point Bayliss window in a column vector.
%
% BAYLISSWIN(N,NBAR) returns an N-point Bayliss window with NBAR nearly
% constant-level sidelobes adjacent to the mainlobe. NBAR must be an
% integer greater than or equal to one.
%
% BAYLISSWIN(N,NBAR,SLL) returns an N-point Bayliss window with SLL maximum
% sidelobe level in dB relative to the mainlobe peak. SLL must be a
% negative value, e.g., -30 dB.
%
% NBAR should satisfy NBAR >= 2*A^2+0.5, where A is equal to
% acosh(10^(-SLL/20))/pi, otherwise the sidelobe level specified is not
% guaranteed. If NBAR is not specified it defaults to 4. SLL defaults to
% -30 dB if not specified.
%
% % EXAMPLE
% %  This example generates a 64-point Bayliss window with 4 sidelobes
% %  adjacent to the mainlobe that are nearly constant-level, and a peak
% %  sidelobe level of -35 dB relative to the mainlobe peak.
%
%   w = baylisswin(64,5,-35);
%   wvtool(w);
%

%
%   References:
%     [1] Carrara, Walter G., Ronald M. Majewski, and Ron S. Goodman,
%         Spotlight Synthetic Aperture Radar: Signal Processing Algorithms,
%         Artech House, October 1, 1995.
%     [2] Brookner, Eli, Practical Phased Array Antenna Systems,
%         Artech House, Inc., 1991, pg. 2-51.

%#codegen

% Validate input and set default values.

narginchk(1, 3);  
[N,NBAR,SLL] = validateinputs(varargin{:});

C = [0.30387530, -0.05042922, -0.00027989, -0.00000343, -0.00000002;...
     0.98583020, -0.03338850,  0.00014064,  0.00000190,  0.00000001;...
     2.00337487, -0.01141548,  0.00041590,  0.00000373,  0.00000001;...
     3.00636321, -0.00683394,  0.00029281,  0.00000161,  0.00000000;...
     4.00518423, -0.00501795,  0.00021735,  0.00000088,  0.00000000];

% mu = [ 
%         0.5860670;
%         1.6970509;
%         2.7171939;
%         3.7261370;
%         4.7312271;
%         5.7345205;
%         6.7368281;
%         7.7385356;
%         8.7398505;
%         9.7408945;
%         10.7417435;
%         11.7424475;
%         12.7430408;
%         13.7435477;
%         14.7439856;
%         15.7443679;
%         16.7447044;
%         17.7450030;
%         18.7452697;
%         19.7455093;
% ];

theta = 0;

A = C(1,:)*SLL.^(0:4).';
global xi_T
xi_T = [C(2,:)*SLL.^(0:4).', C(3,:)*SLL.^(0:4).', ...
        C(4,:)*SLL.^(0:4).', C(5,:)*SLL.^(0:4).'
];

% Bayliss pulse widening (dilation) factor.
% if NBAR > length(mu)
%     sp2 = (1 + (3/4)/NBAR)^2;
% else
%     sp2 = mu(NBAR)^2/(A^2 + NBAR^2);
% end
sp2 = (NBAR+0.5)^2/(A^2+NBAR^2);

Fm = zeros(NBAR,1);
summation = 0;
k = (0:N-1)';
xi = (k-0.5*N+0.5)/N;
for m = (0:NBAR-1),
    Fm(m+1) = calculateFm(m,sp2,A,NBAR);
    summation = Fm(m+1)*sin(2*pi*(m+0.5)*xi)+summation;
end
w = exp(-1j*pi*(0:N-1)'*sind(theta)).*summation;

%-------------------------------------------------------------------
function Fm = calculateFm(m,sp2,A,NBAR)
% Calculate the cosine weights.

Z = zeros(NBAR, 1);
p = [0:m-1, m+1:NBAR-1]'; % p~=m
global xi_T

if NBAR <= 4
    Z(1:NBAR) = xi_T(1:NBAR);
else
    Z(1:4) = xi_T(1:4);
    Z(5:NBAR) = sqrt(A^2+(5:NBAR).^2);
end

Num = prod((1 - ((m+0.5)^2/sp2)./(Z(1:NBAR-1).^2)));
Den = prod((1 - (m+0.5)^2./(p+0.5).^2));

Fm = ((-1)^m*(m+0.5)^2.*Num)./(2j.*Den);

%-------------------------------------------------------------------
function [N,NBAR,SLL] = validateinputs(varargin)
% Cast to enforce Precision Rules
N = signal.internal.sigcasttofloat(varargin{1},'double','taylorwin',...
  'N','allownumeric');
% Validate order
cond = (N ~= floor(N));
if cond
    if coder.target('MATLAB')
        N = round(N);
        warning(message('signal:baylisswin:WindowLengthMustBeInteger'));
    else
        coder.internal.errorIf(cond,'signal:baylisswin:WindowLengthMustBeInteger');
    end
end

cond = (N<0);
if  cond
    coder.internal.errorIf(cond,'signal:baylisswin:WindowLengthMustBePositive');
end

% Validate NBAR
if nargin < 2,
    NBAR = 4;
else
    % Cast to enforce Precision Rules
    NBAR = signal.internal.sigcasttofloat(varargin{2},'double',...
      'baylisswin','NBAR','allownumeric');
    if NBAR<=0 || NBAR~=floor(NBAR),
        error(message('signal:baylisswin:NBARMustBeInteger'));
    end
end

% Validate SLL
if nargin < 3,
    SLL = -30;
else
    % Cast to enforce Precision Rules
    SLL = signal.internal.sigcasttofloat(varargin{3},'double','taylorwin',...
      'SLL','allownumeric');    
    if SLL>0,
        error(message('signal:taylorwin:SLLMustBeNegative'));
    end
end

% [EOF]
