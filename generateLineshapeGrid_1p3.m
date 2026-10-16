
function [gridSpectra, originalSpectrum] = generateLineshapeGrid_1p3(params, expParams, AddLB, ppmGrid)

    if ~isfield(params, 'Skeleton') || isempty(params.Skeleton)
        [~, params.Skeleton] = NSpin_Simulation_1p3( ...
            params, expParams, AddLB);
    end

    if ~isscalar(params.w_CS_J) || ...
            ~isscalar(params.T2) || ...
            ~isscalar(params.GL) || ...
            ~isscalar(params.Skeleton) || ...
            numel(params.Skeleton.freqAll) ~= 1

        error('The reference must contain one uncoupled spin.');
    end

    referenceCS = params.w_CS_J;
    originalFreqHz = params.Skeleton.freqAll;

    originalSpectrum = params.Amp * FastFIDgeneration_1p3( ...
        params, expParams, AddLB);

    gridSpectra = zeros(numel(originalSpectrum), numel(ppmGrid));

    for k = 1:numel(ppmGrid)

        offsetHz = (ppmGrid(k) - referenceCS) * expParams.BF1;

        shiftedParams = params;
        shiftedParams.Skeleton.freqAll = originalFreqHz + offsetHz;

        gridSpectra(:, k) = params.Amp * FastFIDgeneration_1p3( ...
            shiftedParams, expParams, AddLB);

    end

end
