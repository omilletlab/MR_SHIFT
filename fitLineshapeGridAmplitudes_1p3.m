
function result = fitLineshapeGridAmplitudes_1p3( ...
    gridSpectra, originalSpectrum, experimentalSpectrum, ...
    fitMask, ppmGrid, referenceCS)

    fitMask = logical(fitMask(:));
    experimentalSpectrum = real(experimentalSpectrum(:));
    originalSpectrum = real(originalSpectrum(:));

    if size(gridSpectra, 1) ~= numel(experimentalSpectrum)
        error('Grid and experimental spectrum sizes do not match.');
    end

    if nnz(fitMask) <= numel(ppmGrid)
        error('The fitting window contains too few points.');
    end

    X = real(gridSpectra(fitMask, :));
    y = experimentalSpectrum(fitMask);

    if any(~isfinite(X(:))) || any(~isfinite(y))
        error('The fitting data contain nonfinite values.');
    end

    coefficients = lsqnonneg(X, y);

    totalCoefficient = sum(coefficients);

    if totalCoefficient <= 0
        error('No positive grid amplitudes were recovered.');
    end

    weights = coefficients / totalCoefficient;

    fittedSpectrum = real(gridSpectra * coefficients);

    offsetsPPM = ppmGrid(:) - referenceCS;

    result.coefficients = coefficients;
    result.weights = weights;
    result.ppmGrid = ppmGrid(:);
    result.offsetsPPM = offsetsPPM;
    result.meanOffsetPPM = sum(weights .* offsetsPPM);
    result.fittedSpectrum = fittedSpectrum;
    result.totalCoefficient = totalCoefficient;

    result.residualBefore = experimentalSpectrum - originalSpectrum;
    result.residualAfter = experimentalSpectrum - fittedSpectrum;

    result.rmseBefore = sqrt(mean( ...
        result.residualBefore(fitMask).^2));

    result.rmseAfter = sqrt(mean( ...
        result.residualAfter(fitMask).^2));

end
