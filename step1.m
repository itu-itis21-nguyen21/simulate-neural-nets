function [dN, SPON_SPIKES, search_range] = step1(dN, f0, M, n, SPON_SPIKES)
    % Generate a 4xM matrix of random values between 0 and 1
    randomMatrix = rand(4, M);
    %randomMatrix = rand(2, M);  % 2-neuron network
    %randomMatrix = rand(3, M);

    % Set the specified columns of dN based on the firing probability
    dN(:, n:M+n-1) = randomMatrix <= f0;
    
    % Fix search range to M time bins
    search_range = M;

    % Check if any spikes are generated
    s = sum(dN(:,n:search_range+n-1), "all");
    if s > 0
        SPON_SPIKES = SPON_SPIKES + 1;
    end
end