function [phase_avg, profile] = calculateShimPhase(shim, expParams)

% Residual B0 field represented using Cartesian spherical harmonics
% following the convention  https://doi.org/10.1002/cmr.1820020303 (The Ancient and Honourable Art of Shimming)

% The resulting spatial frequency distribution is converted into a
% complex, time-dependent phase factor:
%
%   phase_avg(t) = sum_k w_k * exp(1i * 2*pi * df_k * t)
%
% The phase factor can subsequently be applied to a simulated FID:
%
%   FID = FID .* phase_avg;
%
% This function does not modify experimental data or performs fitting.

  % Z    = z
  % Z2   = 2*z^2 - (x^2 + y^2)
  % Z3   = z * [2*z^2 - 3*(x^2 + y^2)]
  % Z4   = 8*z^2*[z^2 - 3*(x^2 + y^2)] + 3*(x^2 + y^2)^2
  % 
  % X    = x
  % Y    = y
  % XZ   = x*z
  % YZ   = y*z
  % XY   = x*y
  % X2Y2 = x^2 - y^2

  % NOTES
%   - Spatial coordinates are normalized to a cylinder of radius 1
%     and axial extent [-1, 1].
%   - Harmonic coefficients represent frequency deviations in Hz.
%   - Equal-volume deterministic sampling ensures reproducibility.
%   - The mean frequency offset is removed to preserve spectral
%     referencing.
%   - Spatial frequency bins reduce the cost of phase calculation

%% Validate acquisition parameters

validateattributes(expParams.dwell, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});

validateattributes(expParams.TDeff, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'});

if mod(expParams.TDeff, 2) ~= 0
    error('expParams.TDeff must contain an even number of real points.');
end

nTime = expParams.TDeff / 2;
t = (0:nTime-1) * expParams.dwell;

%% Initialize outputs

phase_avg = ones(1, nTime);

profile = struct( ...
    'freqHz', 0, ...
    'weights', 1, ...
    'rmsHz', 0, ...
    'numSpatialPoints', 0);

%% Check whether shim effects are enabled

if nargin < 1 || isempty(shim)
    return;
end

if ~isstruct(shim) || ~isscalar(shim)
    error('shim must be a scalar structure.');
end

if ~isfield(shim, 'Enabled') || isempty(shim.Enabled)
    return;
end

validateattributes(shim.Enabled, {'logical', 'numeric'}, ...
    {'scalar', 'real', 'finite'});

if ~ismember(double(shim.Enabled), [0, 1])
    error('shim.Enabled must be true or false.');
end

if ~logical(shim.Enabled)
    return;
end

%% Read harmonic coefficients

coefficientNames = { ...
    'Z1', 'Z2', 'Z3', 'Z4', ...
    'X', 'Y', 'XZ', 'YZ', 'XY', 'X2Y2'};

coefficients = zeros(1, numel(coefficientNames));

for k = 1:numel(coefficientNames)

    name = coefficientNames{k};

    if isfield(shim, name)

        value = shim.(name);

        validateattributes(value, {'numeric'}, ...
            {'scalar', 'real', 'finite'});

        coefficients(k) = double(value);

    end
end

% All coefficients equal zero: no field inhomogeneity
if all(coefficients == 0)
    return;
end

Z1   = coefficients(1);
Z2   = coefficients(2);
Z3   = coefficients(3);
Z4   = coefficients(4);
X    = coefficients(5);
Y    = coefficients(6);
XZ   = coefficients(7);
YZ   = coefficients(8);
XY   = coefficients(9);
X2Y2 = coefficients(10);

%% Generate deterministic cylindrical sample grid

% Number of spatial positions in each dimension
nRadial  = 12;
nAngular = 32;
nAxial   = 32;

% Equal-area radial sampling
rGrid = sqrt(((1:nRadial) - 0.5) / nRadial);

% Uniform angular sampling
thetaGrid = 2*pi*((1:nAngular) - 0.5) / nAngular;

% Uniform axial sampling
zGrid = -1 + 2*((1:nAxial) - 0.5) / nAxial;

% Three-dimensional sampling coordinates
[r, theta, z] = ndgrid(rGrid, thetaGrid, zGrid);

x = r .* cos(theta);
y = r .* sin(theta);

rho2 = x.^2 + y.^2;

%% Construct Cartesian spherical-harmonic basis

% First-order harmonics
H_Z1 = z;
H_X  = x;
H_Y  = y;

% Second-order harmonics
H_Z2   = 2*z.^2 - rho2;
H_XZ   = x .* z;
H_YZ   = y .* z;
H_XY   = x .* y;
H_X2Y2 = x.^2 - y.^2;

% Third-order axial harmonic
H_Z3 = z .* (2*z.^2 - 3*rho2);

% Fourth-order axial harmonic
H_Z4 = 8*z.^2 .* (z.^2 - 3*rho2) + 3*rho2.^2;

%% Calculate residual B0 frequency distribution

% Frequency deviations in Hz
df = ...
      Z1   * H_Z1 ...
    + Z2   * H_Z2 ...
    + Z3   * H_Z3 ...
    + Z4   * H_Z4 ...
    + X    * H_X ...
    + Y    * H_Y ...
    + XZ   * H_XZ ...
    + YZ   * H_YZ ...
    + XY   * H_XY ...
    + X2Y2 * H_X2Y2;

df = df(:);

% Remove the mean frequency offset.
% This preserves the overall spectral reference.
df = df - mean(df);

%% Compress the spatial frequency distribution

nBins = 512;

[counts, ~, binIdx] = histcounts(df, nBins);

counts = counts(:);

% Mean actual frequency of positions assigned to each bin
binFreq = accumarray(binIdx, df, [nBins, 1], @mean, 0);

% Discard empty frequency bins
keep = counts > 0;

freqHz = binFreq(keep);
weights = counts(keep);

% Normalize distribution weights
weights = weights / sum(weights);

%% Calculate spatially averaged phase evolution

% Process the time axis in chunks to limit memory consumption.
chunkSize = 2048;

for startIdx = 1:chunkSize:nTime

    endIdx = min(startIdx + chunkSize - 1, nTime);

    idx = startIdx:endIdx;

    phaseMatrix = exp(1i * 2*pi * (freqHz * t(idx)));

    phase_avg(idx) = weights.' * phaseMatrix;

end

%% Return frequency distribution for diagnostics / GUI preview

profile.freqHz = freqHz;
profile.weights = weights;
profile.rmsHz = sqrt(mean(df.^2));
profile.numSpatialPoints = numel(df);

end