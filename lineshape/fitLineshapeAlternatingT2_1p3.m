function result = fitLineshapeAlternatingT2_1p3( ...
    params, expParams, AddLB, ppmGrid, experimentalReference, fitMask, options)
% Fit a single-spin reference with a frequency-offset grid and shared T2.
%
% Inputs
%   params                Previously fitted single-spin MR-SHIFT parameters.
%   expParams             Acquisition and processing parameters.
%   AddLB                 Existing extra broadening setting.
%   ppmGrid               Fixed grid positions, in ppm.
%   experimentalReference Experimental spectrum after subtracting other fits.
%   fitMask               Logical mask for the reference fitting interval.
%   options               Optional fields:
%       T2StartFactor          Initial T2 / original T2 (default 1.2).
%       T2BoundsFactors        T2 bounds / original T2 (default [0.5 3]).
%       MaxIterations          Maximum outer iterations (default 30).
%       RelativeTolerance      Relative change in selected objective RMSE
%                              (default 1e-3, i.e. 0.1%).
%       ConvergencePatience    Consecutive low-improvement iterations required
%                              to stop (default 3).
%       FminbndTolX            Tolerance in log(T2) (default 1e-4).
%       UseNoiseInsensitiveLoss Enable a squared epsilon-insensitive objective
%                               (default false).
%       NoiseSigma             Experimental noise standard deviation.
%       ResidualNoiseSigma     Threshold in noise standard deviations
%                               (default 2.0).
%       AmpMaxIterations       Maximum amplitude MM steps (default 40).
%       AmpRelativeTolerance   Amplitude MM tolerance (default 1e-6).
%
% When enabled, the epsilon-insensitive loss is
%   sum(max(abs(experimental - predicted) - delta, 0).^2),
%   delta = ResidualNoiseSigma * NoiseSigma.
% The same loss is used for amplitude optimisation, T2 optimisation,
% monotonicity checks and convergence. Ordinary RMSE is always reported.
% The stopping rule tracks relative improvement in sqrt(loss / nFitPoints),
% rather than in loss itself, for ConvergencePatience consecutive steps.
%
% GL, frequency-grid positions, and source params remain unchanged.
% Requires generateLineshapeGrid_1p3.m,
% fitLineshapeGridAmplitudes_1p3.m and FastFIDgeneration_1p3.m.

    if nargin < 7 || isempty(options)
        options = struct();
    end

    defaults = struct( ...
        'T2StartFactor', 1.2, ...
        'T2BoundsFactors', [0.5 3], ...
        'MaxIterations', 30, ...
        'RelativeTolerance', 1e-3, ...
        'ConvergencePatience', 3, ...
        'FminbndTolX', 1e-4, ...
        'UseNoiseInsensitiveLoss', false, ...
        'NoiseSigma', NaN, ...
        'ResidualNoiseSigma', 2.0, ...
        'AmpMaxIterations', 40, ...
        'AmpRelativeTolerance', 1e-6);

    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(options, name) || isempty(options.(name))
            options.(name) = defaults.(name);
        end
    end

    validateattributes(options.T2StartFactor, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.T2BoundsFactors, {'numeric'}, ...
        {'vector', 'numel', 2, 'real', 'finite', 'positive'});
    validateattributes(options.MaxIterations, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(options.RelativeTolerance, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(options.ConvergencePatience, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(options.FminbndTolX, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(options.UseNoiseInsensitiveLoss, ...
        {'logical', 'numeric'}, {'scalar'});
    validateattributes(options.ResidualNoiseSigma, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(options.AmpMaxIterations, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(options.AmpRelativeTolerance, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});

    if options.T2BoundsFactors(1) >= options.T2BoundsFactors(2)
        error('T2BoundsFactors must be increasing.');
    end

    useInsensitiveLoss = logical(options.UseNoiseInsensitiveLoss);
    if useInsensitiveLoss
        validateattributes(options.NoiseSigma, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'});
        delta = options.ResidualNoiseSigma * options.NoiseSigma;
    else
        delta = 0;
    end

    if ~isscalar(params.w_CS_J) || ~isscalar(params.T2) || ...
            ~isscalar(params.GL) || ~isscalar(params.Amp)
        error('This function requires a fitted single-spin reference.');
    end
    if ~isfinite(params.T2) || params.T2 <= 0 || ...
            ~isfinite(params.Amp) || params.Amp <= 0
        error('The fitted T2 and amplitude must be positive and finite.');
    end
    if isfield(expParams, 'SimulationShimPhase') && ...
            ~isempty(expParams.SimulationShimPhase)
        error('Disable Shim Effects before estimating the correction.');
    end

    ppmGrid = double(ppmGrid(:));
    experimentalReference = real(double(experimentalReference(:)));
    fitMask = logical(fitMask(:));

    if any(~isfinite(ppmGrid)) || numel(unique(ppmGrid)) ~= numel(ppmGrid)
        error('ppmGrid must contain distinct, finite positions.');
    end
    if numel(experimentalReference) ~= expParams.pts(1) || ...
            numel(fitMask) ~= numel(experimentalReference)
        error('Spectrum and fitting mask must match expParams.pts.');
    end
    if nnz(fitMask) <= numel(ppmGrid)
        error('There are too few fitting points for this grid.');
    end
    if any(~isfinite(experimentalReference(fitMask)))
        error('Nonfinite experimental values were found in the fitting window.');
    end

    if ~isfield(params, 'Skeleton') || isempty(params.Skeleton)
        [~, params.Skeleton] = NSpin_Simulation_1p3( ...
            params, expParams, AddLB);
    end
    if ~isscalar(params.Skeleton) || numel(params.Skeleton.freqAll) ~= 1
        error('The reference must contain a single simulated transition.');
    end

    referenceCS = double(params.w_CS_J);
    offsetsHz = (ppmGrid - referenceCS) * double(expParams.BF1);
    t = (0:expParams.TDeff/2-1) * expParams.dwell;

    originalSpectrum = params.Amp * FastFIDgeneration_1p3( ...
        params, expParams, AddLB);
    originalSpectrum = real(originalSpectrum(:));

    originalT2 = double(params.T2);
    t2Bounds = originalT2 * options.T2BoundsFactors;
    currentT2 = min(max(originalT2 * options.T2StartFactor, ...
        t2Bounds(1)), t2Bounds(2));

    currentParams = params;
    currentParams.T2 = currentT2;
    currentGrid = generateLineshapeGrid_1p3( ...
        currentParams, expParams, AddLB, ppmGrid);

    currentFit = fitAmplitudes( ...
        currentGrid, originalSpectrum, experimentalReference, ...
        fitMask, ppmGrid, referenceCS, delta, options, []);
    currentLoss = currentFit.objective;

    history = makeHistoryRow(0, currentT2, currentFit, offsetsHz, NaN);
    stopReason = 'Maximum iterations reached';
    smallImprovementStreak = 0;
    nFitPoints = nnz(fitMask);
    searchOptions = optimset('Display', 'off', ...
        'TolX', options.FminbndTolX);

    for iteration = 1:options.MaxIterations
        if useInsensitiveLoss && currentLoss <= 1e-12
            stopReason = 'All fitting residuals are within the noise tolerance';
            break;
        end

        weights = currentFit.weights(:);
        modulation = weights.' * exp(1i * 2*pi * (offsetsHz * t));
        expWithKernel = expParams;
        expWithKernel.SimulationShimPhase = modulation;

        objectiveT2 = @(logT2) fixedKernelLoss( ...
            exp(logT2), params, expWithKernel, AddLB, ...
            fitMask, experimentalReference(fitMask), delta);

        previousFixedLoss = objectiveT2(log(currentT2));
        [bestLogT2, bestFixedLoss] = fminbnd( ...
            objectiveT2, log(t2Bounds(1)), log(t2Bounds(2)), ...
            searchOptions);

        if bestFixedLoss < previousFixedLoss - ...
                1e-10 * max(1, previousFixedLoss)
            proposedT2 = exp(bestLogT2);
        else
            proposedT2 = currentT2;
        end

        proposedParams = params;
        proposedParams.T2 = proposedT2;
        proposedGrid = generateLineshapeGrid_1p3( ...
            proposedParams, expParams, AddLB, ppmGrid);

        proposedFit = fitAmplitudes( ...
            proposedGrid, originalSpectrum, experimentalReference, ...
            fitMask, ppmGrid, referenceCS, delta, options, ...
            currentFit.coefficients);
        proposedLoss = proposedFit.objective;

        if proposedLoss > currentLoss + 1e-9 * max(1, currentLoss)
            stopReason = 'No further decrease in fitting objective';
            break;
        end

        previousObjectiveRMSE = sqrt(currentLoss / nFitPoints);
        proposedObjectiveRMSE = sqrt(proposedLoss / nFitPoints);
        improvement = max(0, ...
            (previousObjectiveRMSE - proposedObjectiveRMSE) / ...
            max(previousObjectiveRMSE, eps));

        currentT2 = proposedT2;
        currentFit = proposedFit;
        currentLoss = proposedLoss;
        history(end+1) = makeHistoryRow( ...
            iteration, currentT2, currentFit, offsetsHz, improvement); %#ok<AGROW>

        if useInsensitiveLoss && currentLoss <= 1e-12
            stopReason = 'All fitting residuals are within the noise tolerance';
            break;
        end
        if improvement < options.RelativeTolerance
            smallImprovementStreak = smallImprovementStreak + 1;
        else
            smallImprovementStreak = 0;
        end
        if smallImprovementStreak >= options.ConvergencePatience
            stopReason = sprintf( ...
                'Objective RMSE improvement below %.4g%% for %d consecutive iterations', ...
                100*options.RelativeTolerance, options.ConvergencePatience);
            break;
        end
    end

    result = currentFit;
    result.T2Original = originalT2;
    result.T2Start = history(1).T2;
    result.T2Optimized = currentT2;
    result.T2Bounds = t2Bounds;
    result.GLFixed = params.GL;
    result.offsetsHz = offsetsHz;
    result.meanOffsetHz = sum(currentFit.weights(:) .* offsetsHz);
    result.rmsOffsetHz = sqrt(sum(currentFit.weights(:) .* ...
        (offsetsHz - result.meanOffsetHz).^2));
    result.SimulationShimPhase = currentFit.weights(:).' * ...
        exp(1i * 2*pi * (offsetsHz * t));
    result.history = history;
    result.iterations = numel(history) - 1;
    result.stopReason = stopReason;
    result.relativeTolerance = options.RelativeTolerance;
    result.convergencePatience = options.ConvergencePatience;
    result.smallImprovementStreak = smallImprovementStreak;
    result.rmseAtStart = history(1).rmse;
    result.improvementPercent = 100 * (1 - ...
        result.rmseAfter / max(result.rmseBefore, eps));
    result.nearT2Bound = currentT2 <= 1.01*t2Bounds(1) || ...
        currentT2 >= 0.99*t2Bounds(2);
    result.UseNoiseInsensitiveLoss = useInsensitiveLoss;
    result.noiseSigma = options.NoiseSigma;
    result.noiseThreshold = delta;

    finalParams = params;
    finalParams.T2 = currentT2;
    finalExpParams = expParams;
    finalExpParams.SimulationShimPhase = result.SimulationShimPhase;
    checkSpectrum = params.Amp * result.totalCoefficient * ...
        FastFIDgeneration_1p3(finalParams, finalExpParams, AddLB);
    result.fidConsistencyError = norm( ...
        real(checkSpectrum(:)) - result.fittedSpectrum(:)) / ...
        max(norm(result.fittedSpectrum(:)), eps);
end

function fitResult = fitAmplitudes( ...
    gridSpectra, originalSpectrum, experimentalReference, ...
    fitMask, ppmGrid, referenceCS, delta, options, initialCoefficients)

    if ~logical(options.UseNoiseInsensitiveLoss)
        fitResult = fitLineshapeGridAmplitudes_1p3( ...
            gridSpectra, originalSpectrum, experimentalReference, ...
            fitMask, ppmGrid, referenceCS);
        residual = fitResult.residualAfter(fitMask);
        fitResult.objective = residual(:).' * residual(:);
        fitResult.insensitiveRMSE = fitResult.rmseAfter;
        fitResult.pointsAboveThreshold = nnz(abs(residual) > 0);
        return;
    end

    X = real(double(gridSpectra(fitMask, :)));
    y = real(double(experimentalReference(fitMask)));
    y = y(:);

    if isempty(initialCoefficients)
        coefficients = lsqnonneg(X, y);
    else
        coefficients = max(0, initialCoefficients(:));
    end

    previousResidual = y - X * coefficients;
    previousLoss = insensitiveLoss(previousResidual, delta);

    for j = 1:options.AmpMaxIterations
        clippedResidual = min(max(previousResidual, -delta), delta);
        pseudoTarget = y - clippedResidual;
        proposedCoefficients = lsqnonneg(X, pseudoTarget);
        proposedResidual = y - X * proposedCoefficients;
        proposedLoss = insensitiveLoss(proposedResidual, delta);

        if proposedLoss > previousLoss + 1e-9 * max(1, previousLoss)
            break;
        end

        improvement = max(0, previousLoss - proposedLoss);
        coefficients = proposedCoefficients;
        previousResidual = proposedResidual;
        previousLoss = proposedLoss;

        if previousLoss <= 1e-12 || ...
                improvement <= options.AmpRelativeTolerance * ...
                max(1, previousLoss)
            break;
        end
    end

    amplitudeSum = sum(coefficients);
    if amplitudeSum <= 0
        error('No positive grid amplitudes were recovered.');
    end

    fittedSpectrum = real(double(gridSpectra * coefficients));
    residualBefore = real(double(experimentalReference(:) - ...
        originalSpectrum(:)));
    residualAfter = real(double(experimentalReference(:) - ...
        fittedSpectrum(:)));
    offsetsPPM = ppmGrid(:) - referenceCS;
    weights = coefficients / amplitudeSum;

    fitResult = struct();
    fitResult.coefficients = coefficients;
    fitResult.weights = weights;
    fitResult.ppmGrid = ppmGrid(:);
    fitResult.offsetsPPM = offsetsPPM;
    fitResult.meanOffsetPPM = sum(weights .* offsetsPPM);
    fitResult.fittedSpectrum = fittedSpectrum;
    fitResult.totalCoefficient = amplitudeSum;
    fitResult.residualBefore = residualBefore;
    fitResult.residualAfter = residualAfter;
    fitResult.rmseBefore = sqrt(mean(residualBefore(fitMask).^2));
    fitResult.rmseAfter = sqrt(mean(residualAfter(fitMask).^2));
    fitResult.objective = insensitiveLoss(residualAfter(fitMask), delta);
    fitResult.insensitiveRMSE = sqrt( ...
        fitResult.objective / nnz(fitMask));
    fitResult.pointsAboveThreshold = ...
        nnz(abs(residualAfter(fitMask)) > delta);
end

function value = fixedKernelLoss( ...
    trialT2, params, expWithKernel, AddLB, fitMask, targetWindow, delta)

    trialParams = params;
    trialParams.T2 = trialT2;
    template = params.Amp * FastFIDgeneration_1p3( ...
        trialParams, expWithKernel, AddLB);
    template = real(double(template(fitMask)));
    template = template(:);
    y = targetWindow(:);

    if ~all(isfinite(template)) || dot(template, template) <= eps
        value = inf;
        return;
    end

    scale = bestNonnegativeScale(template, y, delta);
    value = insensitiveLoss(y - scale * template, delta);
end

function scale = bestNonnegativeScale(template, target, delta)

    denominator = dot(template, template);
    ordinaryScale = max(0, dot(template, target) / denominator);
    if delta <= 0
        scale = ordinaryScale;
        return;
    end

    if insensitiveLoss(target - ordinaryScale*template, delta) <= 1e-12
        scale = ordinaryScale;
        return;
    end

    if residualSlope(0, template, target, delta) <= 0
        scale = 0;
        return;
    end

    lower = 0;
    upper = max(1, 2*ordinaryScale);
    for j = 1:40
        if residualSlope(upper, template, target, delta) <= 0
            break;
        end
        upper = upper * 2;
    end

    if residualSlope(upper, template, target, delta) > 0
        scale = upper;
        return;
    end

    for j = 1:50
        middle = (lower + upper)/2;
        if residualSlope(middle, template, target, delta) > 0
            lower = middle;
        else
            upper = middle;
        end
    end
    scale = (lower + upper)/2;
end

function slope = residualSlope(scale, template, target, delta)
    residual = target - scale * template;
    effectiveResidual = sign(residual) .* max(abs(residual) - delta, 0);
    slope = dot(template, effectiveResidual);
end

function value = insensitiveLoss(residual, delta)
    excess = max(abs(residual(:)) - delta, 0);
    value = sum(excess.^2);
end

function row = makeHistoryRow(iteration, T2, fitResult, offsetHz, improvement)
    weights = fitResult.weights(:);
    meanHz = weights.' * offsetHz;
    row = struct( ...
        'iteration', iteration, ...
        'T2', T2, ...
        'rmse', fitResult.rmseAfter, ...
        'insensitiveRMSE', fitResult.insensitiveRMSE, ...
        'objective', fitResult.objective, ...
        'pointsAboveThreshold', fitResult.pointsAboveThreshold, ...
        'relativeImprovement', improvement, ...
        'rmsOffsetHz', sqrt(sum(weights .* (offsetHz - meanHz).^2)), ...
        'amplitudeFactor', fitResult.totalCoefficient);
end
