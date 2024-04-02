function pwr = bandpower(varargin)
%BANDPOWER Band power
%   P = BANDPOWER(X) computes the average power in the input signal vector,
%   X.  If X is a matrix, then BANDPOWER computes the average power in each
%   column independently.
%
%   P = BANDPOWER(X, Fs, FREQRANGE) specifies FREQRANGE as a two-element
%   vector of real values, specifying the two frequencies between which you
%   want to measure the power.  Fs is the sampling rate of the input
%   signal.  The power is estimated by applying a Hamming window and using
%   a periodogram of the same length as the input vector.  If the input
%   vector X contains N samples, then FREQRANGE must be contained within
%   the interval:
%      [              0, Fs/2          ] if X is real and N is even
%      [              0, Fs*(N-1)/(2*N)] if X is real and N is odd
%      [-Fs*(N-2)/(2*N), Fs/2          ] if X is complex and N is even
%      [-Fs*(N-1)/(2*N), Fs*(N-1)/(2*N)] if X is complex and N is odd
%
%   P = BANDPOWER(Pxx, F, 'psd') computes the average power via a rectangle
%   approximation of the integral of the Power Spectral Density (PSD)
%   estimate, given in the vector Pxx over the frequencies specified in
%   vector F.
%
%   P = BANDPOWER(Pxx, F, FREQRANGE, 'psd') specifies FREQRANGE as a
%   two-element vector of real values, specifying the two frequencies
%   between which you want to measure the power.  F is a vector containing
%   the frequencies that correspond to the estimates given in Pxx.
%
%   NOTE: If the frequency range values don't exactly match the frequency
%   values stored in F then the next closest value is used.
%
%   % Example
%   %   Examine the power present within 50 kHz of the carrier of a 1 kHz
%   %   sinusoid FM modulated at 250 kHz with a modulation index of 0.1.
%
%   Fs = 1.0e6; Fc = 250e3; Fd = 1e3; m = 0.1;
%   x = vco(m*sin(2*pi*Fd*(1:1000)/Fs),Fc,Fs);
%   periodogram(x,[],length(x),Fs);
%   mod_power_dB = 10*log10(bandpower(x, Fs, Fc + [-50e3 50e3]))
%
%   See also POWERBW OBW PERIODOGRAM PWELCH MEANFREQ MEDFREQ PLOMB.

%   Copyright 2012-2021 The MathWorks, Inc.
%#codegen

narginchk(1,4);

inputArgs = cell(size(varargin));
[inputArgs{:}] = convertStringsToChars(varargin{:});

idx = 0; % 'psd' flag is absent

for i = coder.unroll(1:nargin)
    if ischar(inputArgs{i})
        % input flag must be a compile time constant for code generation
        coder.internal.assert(coder.internal.isConst(inputArgs{i}),'signal:bandpower:PsdAsConst');
        % 'psd' is the only supported input flag
        coder.internal.assert(strcmpi(inputArgs{i},'psd'),'signal:bandpower:IncorrectStringPsd');
        idx = i; % 'psd' flag is present at index 'idx'
    end
end

if idx == 0
    % Time-domain input
    pwr = timedomainbandpower(inputArgs{:});
else
    % Power spectral density input
    pwr = psdbandpower(inputArgs{1:idx-1});
end

end

function pwr = timedomainbandpower(x,fs,freqrange)

% invalid function call: bandpower(x,fs)
coder.internal.assert(nargin==1 || nargin==3,'signal:bandpower:FreqRangeMissing');

% perform column vector conversion on 'x' before checking 2d matrix
if isrow(x)
    colX = x(:);
else
    colX = x;
end
validateattributes(colX,{'single','double'},{'2d','finite'},...
    'bandpower','X',1);

if nargin==1
    % full range specified. real cast required for codegen
    pwr = real(mean(conj(colX).*colX,1));
else
    validateattributes(fs,{'numeric'},{'scalar','finite','real','positive'}, ...
        'bandpower','Fs',2);
    
    % Cast to enforce precision rules
    fsScalar = double(fs(1));
    % Compute periodogram using a hamming window with the same length as
    % input
    n = size(colX,1);
    if isreal(colX)
        [Pxx,F] = periodogram(colX,hamming(n),n,fsScalar);
    else
        [Pxx,F] = periodogram(colX,hamming(n),n,fsScalar,'centered');
    end
    % Return the bandpower
    pwr = psdbandpower(Pxx,F,freqrange);
end
end

function pwr = psdbandpower(Pxx,W,freqrange)

% Function call cannot be bandpower(pxx,'psd')
coder.internal.errorIf(nargin < 2,'signal:bandpower:FreqVectorMissing')

% perform column vector conversion on 'Pxx' before checking 2d matrix
if isrow(Pxx)
    colPxx = Pxx(:);
else
    colPxx = Pxx;
end
validateattributes(colPxx,{'single','double'},{'2d','finite','real'}, ...
    'bandpower','Pxx',1);

validateattributes(W,{'numeric'},{'vector','finite','real','increasing',...
    'numel',size(colPxx,1)},'bandpower','F',2);

% Cast W to enforce Precision rules and force a column vector
colW = double(W(:));

if nargin < 3
    freqrange = [colW(1) colW(end)];
    freqrangespecified = false;
else
    validateattributes(freqrange,{'numeric'},{'vector','finite','real',...
        'numel',2,'increasing'},'bandpower','FREQRANGE',3);
    coder.internal.errorIf((freqrange(1)<colW(1) || freqrange(2)>colW(end)),...
        'signal:bandpower:FreqRangeOutOfBounds');
    freqrangespecified = true;
end

% Find indices of freq range requested.
idx = find(colW <= freqrange(1),1,'last');
idx1 = idx(1,1);
idx = find(colW >= freqrange(2),1,'first');
idx2 = idx(1,1);

% Determine the width of the rectangle used to approximate the integral.
W_diff = diff(colW);
if freqrangespecified
    lastRectWidth = 0;  % Don't include last point of PSD data.
    width = [W_diff; lastRectWidth];
else
    % There are two cases when spectrum is twosided, CenterDC or not.
    % In both cases, the frequency samples does not cover the entire
    % 2*pi (or Fs) region due to the periodicity.  Therefore, the
    % missing freq range has to be compensated in the integral.  The
    % missing freq range can be calculated as the difference between
    % 2*pi (or Fs) and the actual frequency vector span.  For example,
    % considering 1024 points over 2*pi, then frequency vector will be
    % [0 2*pi*(1-1/1024)], i.e., the missing freq range is 2*pi/1024.
    %
    % When CenterDC is true, if the number of points is even, the
    % Nyquist point (Fs/2) is exact, therefore, the missing range is at
    % the left side, i.e., the beginning of the vector.  If the number
    % of points is odd, then the missing freq range is at both ends.
    % However, due to the symmetry of the real signal spectrum, it can
    % still be considered as if it is missing at the beginning of the
    % vector.  Even when the spectrum is asymmetric, since the
    % approximation of the integral is close when NFFT is large,
    % putting it in the beginning of the vector is still ok.
    %
    % When CenterDC is false, the missing range is always at the end of
    % the frequency vector since the frequency always starts at 0.
    
    % assuming a relatively uniform interval
    missingWidth = (colW(end) - colW(1)) / (numel(colW) - 1);
    
    % if CenterDC was not specified, the first frequency point will
    % be 0 (DC).
    centerDC = colW(1) ~= 0;
    if centerDC
        width = [missingWidth; W_diff];
    else
        width = [W_diff; missingWidth];
    end
end

% Sum the average power over the range of interest.
pwr = width(idx1:idx2)'*colPxx(idx1:idx2,:);

end