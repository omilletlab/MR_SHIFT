function result = fitAsymmetricTail_1p3( ...
    params, expParams, AddLB, experimentalReference, fitMask, options)
% Fit a single-spin reference using T2, GL, CS shift and an asymmetric tail.
%
% Inputs
%   params                Previously fitted single-spin MR-SHIFT parameters.
%   expParams             Acquisition and processing parameters.
%   AddLB                 Existing additional line-broadening setting.
%   experimentalReference Experimental spectrum with other fits subtracted.
%   fitMask               Logical mask of points used for fitting.
%   options               Optional fields:
%       Direction              +1 or -1 (default -1).
%       T2StartFactor          Starting T2 / original T2 (default 1.2).
%       T2BoundsFactors        Bounds / original T2 (default [0.5 4]).
%       GLStart                Initial GL in native MR-SHIFT units.
%       GLBounds               GL bounds (default [0 2]).
%       CSShiftStartHz         Initial reference chemical-shift offset (default 0 Hz).
%       CSShiftBoundsHz        Allowed reference shift (default [-0.5 0.5] Hz).
%       EtaStart               Initial tail fraction (default 0.1).
%       EtaBounds              Tail-fraction bounds (default [0 0.6]).
%       LambdaStartHz          Initial tail scale in Hz (default 0.5).
%       LambdaBoundsHz         Tail-scale bounds in Hz (default [0.02 3]).
%       UseNoiseInsensitiveLoss Use squared epsilon-insensitive loss.
%       NoiseSigma             Experimental noise standard deviation.
%       ResidualNoiseSigma     Threshold multiplier (default 2).
%
% T2, GL, reference CS shift, eta and lambda are nonlinear parameters.
% Amplitude is fitted linearly or with a one-dimensional noise-insensitive
% objective. The reference CS shift is a local fitting parameter; no change
% is made to params.w_CS_J or to the saved global tail modulation.
% The no-tail comparator optimises T2, GL and the same reference CS shift.
% Requires calculateAsymmetricTailPhase_1p3.m and FastFIDgeneration_1p3.m.

    if nargin < 6 || isempty(options)
        options = struct();
    end

    defaults = struct( ...
        'Direction', -1, ...
        'T2StartFactor', 1.2, ...
        'T2BoundsFactors', [0.5 8], ...
        'GLStart', [], ...
        'GLBounds', [0 2], ...
        'CSShiftStartHz', 0, ...
        'CSShiftBoundsHz', [-0.5 0.5], ...
        'EtaStart', 0.1, ...
        'EtaBounds', [0 0.9], ...
        'LambdaStartHz', 0.5, ...
        'LambdaBoundsHz', [0.02 3], ...
        'UseNoiseInsensitiveLoss', false, ...
        'NoiseSigma', NaN, ...
        'ResidualNoiseSigma', 2);

    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(options, name) || isempty(options.(name))
            options.(name) = defaults.(name);
        end
    end

    if ~isscalar(options.Direction) || ~ismember(options.Direction, [-1, 1])
        error('Direction must be -1 or +1.');
    end
    validateattributes(options.T2StartFactor, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.EtaStart, {'numeric'}, ...
        {'scalar', 'real', 'finite'});
    validateattributes(options.LambdaStartHz, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.CSShiftStartHz, {'numeric'}, ...
        {'scalar', 'real', 'finite'});
    validateattributes(options.ResidualNoiseSigma, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});

    if ~isscalar(params.w_CS_J) || ~isscalar(params.T2) || ...
            ~isscalar(params.GL) || ~isscalar(params.Amp) || ...
            ~isfinite(params.T2) || params.T2 <= 0 || ...
            ~isfinite(params.GL) || params.GL < 0 || ...
            ~isfinite(params.Amp) || params.Amp <= 0
        error('Select a previously fitted single-spin reference.');
    end

    bounds = [ ...
        params.T2 * options.T2BoundsFactors(:).'; ...
        options.GLBounds(:).'; ...
        options.EtaBounds(:).'; ...
        options.LambdaBoundsHz(:).'; ...
        options.CSShiftBoundsHz(:).'];
    if any(~isfinite(bounds(:))) || any(bounds(:, 1) >= bounds(:, 2)) || ...
            any(bounds(1, :) <= 0) || bounds(2, 1) < 0 || ...
            bounds(3, 1) < 0 || bounds(3, 2) > 1 || ...
            any(bounds(4, :) <= 0)
        error('Invalid T2, GL, eta, lambda or CS-shift bounds.');
    end

    if isfield(expParams, 'SimulationShimPhase') && ...
            ~isempty(expParams.SimulationShimPhase)
        error('Disable the existing lineshape correction before fitting.');
    end
    if ~isfield(params, 'Skeleton') || isempty(params.Skeleton)
        [~, params.Skeleton] = NSpin_Simulation_1p3( ...
            params, expParams, AddLB);
    end
    if ~isscalar(params.Skeleton) || numel(params.Skeleton.freqAll) ~= 1
        error('The reference must have a single simulated transition.');
    end

    target = real(double(experimentalReference(:)));
    fitMask = logical(fitMask(:));
    if numel(target) ~= expParams.pts(1) || ...
            numel(fitMask) ~= numel(target) || nnz(fitMask) < 10 || ...
            any(~isfinite(target(fitMask)))
        error('Invalid experimental reference or fitting window.');
    end

    if logical(options.UseNoiseInsensitiveLoss)
        validateattributes(options.NoiseSigma, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'});
        delta = options.ResidualNoiseSigma * options.NoiseSigma;
    else
        delta = 0;
    end

    validateattributes(expParams.BF1, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});

    y = target(fitMask);
    t = (0:expParams.TDeff/2-1) * expParams.dwell;
    originalSpectrum = real(params.Amp * ...
        FastFIDgeneration_1p3(params, expParams, AddLB));
    originalSpectrum = originalSpectrum(:);

    simplexOptions = optimset('Display', 'off', ...
        'MaxIter', 500, 'MaxFunEvals', 850, ...
        'TolX', 1e-4, 'TolFun', 1e-6);
    scalarOptions = optimset('Display', 'off', 'TolX', 1e-5);

    if isempty(options.GLStart)
        glStart = max(params.GL, 0.3);
    else
        validateattributes(options.GLStart, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'nonnegative'});
        glStart = options.GLStart;
    end

    t2Start = params.T2 * options.T2StartFactor;
    start = [t2Start, glStart, options.EtaStart, ...
        options.LambdaStartHz, options.CSShiftStartHz];
    start = clampToBounds(start, bounds);

    [zFit, tailLoss, tailExitFlag] = fminsearch( ...
        @lossWithTail, encode(start, bounds), simplexOptions);
    candidate = decode(zFit, bounds);
    [tailLoss, tailAmp, tailSpectrum] = evaluate( ...
        candidate(1), candidate(2), candidate(3), ...
        candidate(4), candidate(5));

    noTailBounds = bounds([1 2 5], :);
    noTailStart = clampToBounds(start([1 2 5]), noTailBounds);
    [zNoTail, noTailLoss, noTailExitFlag] = fminsearch( ...
        @lossWithoutTail, encode(noTailStart, noTailBounds), simplexOptions);
    noTail = decode(zNoTail, noTailBounds);
    [noTailLoss, noTailAmp, noTailSpectrum] = evaluate( ...
        noTail(1), noTail(2), 0, options.LambdaStartHz, noTail(3));

    if noTailLoss <= tailLoss + 1e-9 * max(1, noTailLoss)
        best = [noTail(1), noTail(2), 0, ...
            options.LambdaStartHz, noTail(3)];
        bestLoss = noTailLoss;
        bestAmp = noTailAmp;
        bestSpectrum = noTailSpectrum;
        bestFlag = noTailExitFlag;
    else
        best = candidate;
        bestLoss = tailLoss;
        bestAmp = tailAmp;
        bestSpectrum = tailSpectrum;
        bestFlag = tailExitFlag;
    end

    result.T2Original = params.T2;
    result.T2Optimized = best(1);
    result.GLOriginal = params.GL;
    result.GLOptimized = best(2);
    result.eta = best(3);
    result.lambdaHz = best(4);
    result.direction = options.Direction;
    result.CSShiftHz = best(5);
    result.CSShiftNoTailHz = noTail(3);
    result.CSShiftPPM = best(5) / expParams.BF1;
    result.CSOriginalPPM = params.w_CS_J;
    result.CSOptimizedPPM = params.w_CS_J + result.CSShiftPPM;
    result.amplitudeFactor = bestAmp;
    result.amplitudeOptimized = params.Amp * bestAmp;
    result.originalSpectrum = originalSpectrum;
    result.fittedSpectrum = bestSpectrum;
    result.residualBefore = target - originalSpectrum;
    result.residualAfter = target - bestSpectrum;
    result.rmseBefore = sqrt(mean(result.residualBefore(fitMask).^2));
    result.rmseAfter = sqrt(mean(result.residualAfter(fitMask).^2));
    result.rmseNoTail = sqrt(mean((y - noTailSpectrum(fitMask)).^2));
    result.T2NoTail = noTail(1);
    result.GLNoTail = noTail(2);
    result.improvementPercent = 100 * (1 - ...
        result.rmseAfter / max(result.rmseBefore, eps));
    result.improvementVsNoTailPercent = 100 * (1 - ...
        result.rmseAfter / max(result.rmseNoTail, eps));
    result.objective = bestLoss;
    result.objectiveNoTail = noTailLoss;
    result.noiseThreshold = delta;
    result.optimizerExitFlag = bestFlag;
    result.tailExitFlag = tailExitFlag;
    result.noTailExitFlag = noTailExitFlag;
    result.modulation = calculateAsymmetricTailPhase_1p3( ...
        t, best(3), best(4), options.Direction);

    function p = clampToBounds(p, localBounds)
        lo = localBounds(:, 1).';
        hi = localBounds(:, 2).';
        p = max(lo + 1e-4*(hi-lo), ...
            min(hi - 1e-4*(hi-lo), p));
    end

    function z = encode(p, localBounds)
        low = localBounds(:, 1).';
        high = localBounds(:, 2).';
        fraction = (p - low) ./ (high - low);
        fraction = min(max(fraction, 1e-4), 1 - 1e-4);
        z = log(fraction ./ (1 - fraction));
    end

    function p = decode(z, localBounds)
        low = localBounds(:, 1).';
        high = localBounds(:, 2).';
        p = low + (high - low) ./ (1 + exp(-z));
    end

    function loss = lossWithTail(z)
        p = decode(z, bounds);
        loss = evaluate(p(1), p(2), p(3), p(4), p(5));
    end

    function loss = lossWithoutTail(z)
        p = decode(z, noTailBounds);
        loss = evaluate(p(1), p(2), 0, options.LambdaStartHz, p(3));
    end

    function [loss, ampFactor, spectrum] = evaluate(T2, GL, eta, lambdaHz, csShiftHz)
        trialParams = params;
        trialParams.T2 = T2;
        trialParams.GL = GL;
        trialExpParams = expParams;
        trialExpParams.SimulationShimPhase = ...
            calculateAsymmetricTailPhase_1p3( ...
                t, eta, lambdaHz, options.Direction) .* ...
            exp(1i * 2*pi * csShiftHz * t);

        template = real(params.Amp * ...
            FastFIDgeneration_1p3(trialParams, trialExpParams, AddLB));
        template = template(:);
        x = template(fitMask);
        denominator = x.' * x;
        if ~isfinite(denominator) || denominator <= realmin
            loss = inf;
            ampFactor = 0;
            spectrum = zeros(size(template));
            return;
        end

        ampFactor = max(0, (x.' * y) / denominator);
        if delta > 0
            ampLimit = max([4 * ampFactor, ...
                2 * norm(y) / sqrt(denominator), 1e-3]);
            ampFactor = fminbnd( ...
                @(a) sum(max(abs(y - a*x) - delta, 0).^2), ...
                0, ampLimit, scalarOptions);
        end

        spectrum = ampFactor * template;
        residual = y - spectrum(fitMask);
        loss = sum(max(abs(residual) - delta, 0).^2);
    end
end
