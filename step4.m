function P = step4(dN, n, M, F, f0, interaction_from_j_to_neurons)
    % Compute H as the minimum of n and M
    H = min(n-1, M);
    
    % Initialize the output matrix for firing probabilities
    P = zeros(4, F);
        
    % Loop over each neuron
    for neuron = 1:4
        % Initialize convolution results
        o = zeros(4, M + F + H - 1);

        % Loop over each neuron j to compute contributions
        for j = 1:4 
            % Extract history and prepare link with future
            history = dN(j, n-H:n-1);
            future = zeros(1, F);
            link = [history, future];
            
            % Compute convolution with the respective interaction matrix
            o(j, :) = conv(link, interaction_from_j_to_neurons{neuron}(j, :));
        end
        
        % Sum the convolutions across neurons
        o = sum(o);
        
        % Compute firing probabilities for the given time steps
        P(neuron, :) = exp(o(H:H+F-1)) * f0;
        P(neuron, :) = 1 - exp(-P(neuron, :));
    end
end