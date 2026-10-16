
function modulation = calculateAsymmetricTailPhase_1p3(t, eta, lambdaHz, direction)

    validateattributes(eta, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>=', 0, '<=', 1});

    validateattributes(lambdaHz, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});

    if ~isscalar(direction) || ~ismember(direction, [-1, 1])
        error('Direction must be -1 or +1.');
    end

    modulation = (1 - eta) + ...
        eta ./ (1 - 1i * direction * 2*pi * lambdaHz .* t);

end
