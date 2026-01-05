function [prediction, s] = idw_prediction_and_uncertainty(x, known_points, values, p)
    % IDW_PREDICTION_AND_UNCERTAINTY Calculate IDW prediction and uncertainty.
    %   [PREDICTION, S2] = IDW_PREDICTION_AND_UNCERTAINTY(X, KNOWN_POINTS, VALUES, P)
    %   returns the IDW prediction and local variance estimate at point X.
    %
    %   X is a vector representing the coordinates of the prediction point.
    %   KNOWN_POINTS is an NxM matrix where each row represents the coordinates
    %   of a known point.
    %   VALUES is an Nx1 vector of values at the known points.
    %   P is the power parameter for IDW, default is 2.
    
    if nargin < 4
        p = 2; % Default power parameter
    end
    
    % Number of known points
    n = size(known_points, 1);
    n = size(x, 1);
    % Calculate the weights
%     weights = zeros(n, 1);
    for i = 1:n
        weights = zeros(size(known_points, 1), 1);
        A = repmat(x(i,:),size(known_points, 1),1) - known_points;
        distance = sqrt(sum(A.^2, 2));
        weights = 1 ./ (distance.^p);
        weights = weights / sum(weights);
        prediction = sum(weights .* values);
        s(i) = sum(weights .* (values - prediction).^2);
    end
end