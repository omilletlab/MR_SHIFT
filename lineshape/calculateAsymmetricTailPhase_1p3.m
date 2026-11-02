
function modulation = calculateAsymmetricTailPhase_1p3(t,  lambdaHz, direction)

    validateattributes(lambdaHz, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});

    if ~isscalar(direction) || ~ismember(direction, [-1, 1])
        error('Direction must be -1 or +1.');
    end

    % A one-sided exponential frequency-offset distribution has mean
    % direction * lambdaHz. Remove this mean frequency displacement so that
    % the correction changes lineshape/asymmetry without globally shifting
    % the resonance.
    meanOffsetHz = direction * lambdaHz;

    modulation = exp(-1i * 2*pi * meanOffsetHz .* t) ./ ...
        (1 - 1i * direction * 2*pi * lambdaHz .* t);

end
