tic
% Define raw spontaneous firing rate (in Hz)
f0_hz = 25;
% Define the time interval Delta t
delta_t = 1e-3;
% Define the time variable
%t = (0:(180*1e2-1)) * delta_t;
t = (0:179) * delta_t;
% Time window (in ms)
%M = 180 * 1e2;
M = 180;
% Spontaneous firing rate
f0 = f0_hz * delta_t;
%f0 = 2.5 * 1e-5;
% Pointer to the current bin number 
n = 1;
% F time bins = inverse of f0
F = 1/f0;
%F = 1 / (f0*1e-3);
% Spiking count
spike_count = zeros(4,1);
% Stopping rule
max_spike_no = 5000;
%max_spike_no = 10000;
% Save after every 1 million time stepsnot
chunk_size = 5e5;
% Generated activity matrix (size 4xK)
dN = zeros(4,chunk_size+M);
% File saving parameters
file_index = 1;
% How far should I search for spikes (in time bins) (M or F)?
search_range = 0;

% COUNTING VARS
% Probabilities (not) exceeding 1 counters
EXCEED1 = zeros(4,1);
NOTEXCEED1 = 0;
% Counter for number of times step 1 generates a spike (spontaneous firing
% rate)
SPON_SPIKES = 0;
% Counter for number of times step 4 generates a spike (spiking based on
% history)
FUNC_SPIKES = 0;

% Interaction functions related variables
au = auto_func(t);
exf = excite_func(t);
exf_delayed = excite_func(t - 3*10^(-2));
ihf = inhibit_func(t);
ihf_delayed = inhibit_func(t - 6*10^(-2));
u = unit_step(t);
u_delayed_3 = unit_step(t - 3*10^(-2));
u_delayed_6 = unit_step(t - 6*10^(-2));

interaction_from_j_to_neuron = cell(1,4);
% Interaction matrix to neuron 1
interaction_from_j_to_neuron{1}(1,:) = au .* u;         % 1->1: auto
%interaction_from_j_to_neuron{1}(1,:) = ihf .* u;
%interaction_from_j_to_neuron{1}(1,:) = 0;
interaction_from_j_to_neuron{1}(3,:) = 0;
interaction_from_j_to_neuron{1}(4,:) = ihf_delayed .* u_delayed_6;  % 4->1: inhibit (delay 60ms)
%interaction_from_j_to_neuron{1}(4,:) = ihf .* u;        % 4->1: inhibit
%interaction_from_j_to_neuron{1}(4,:) = 0;
interaction_from_j_to_neuron{1}(2,:) = 0;
%interaction_from_j_to_neuron{1}(2,:) = exf .* u;
% Interaction matrix to neuron 2
%interaction_from_j_to_neuron{2}(2,:) = ihf .* u;        % 2->2: inhibit
interaction_from_j_to_neuron{2}(2,:) = au .* u;
%interaction_from_j_to_neuron{2}(2,:) = 0;
interaction_from_j_to_neuron{2}(3,:) = 0;
%interaction_from_j_to_neuron{2}(3,:) = exf .* u;
interaction_from_j_to_neuron{2}(4,:) = ihf .* u;        % 4->2: inhibit
%interaction_from_j_to_neuron{2}(4,:) = 0;
interaction_from_j_to_neuron{2}(1,:) = exf .* u;        % 1->2: excite
%interaction_from_j_to_neuron{2}(1,:) = ihf .* u;
%interaction_from_j_to_neuron{2}(1,:) = 0;
% Interaction matrix to neuron 3                   
%interaction_from_j_to_neuron{3}(2,:) = ihf .* u;
%interaction_from_j_to_neuron{3}(2,:) = exf .* u;
%interaction_from_j_to_neuron{3}(3,:) = ihf .* u;        % 3->3: inhibit
%interaction_from_j_to_neuron{3}(3,:) = 0;
interaction_from_j_to_neuron{3}(3,:) = au .* u;
interaction_from_j_to_neuron{3}(4,:) = ihf .* u;        % 4->3: inhibit
%interaction_from_j_to_neuron{3}(4,:) = 0;
interaction_from_j_to_neuron{3}(1,:) = exf_delayed .* u_delayed_3;    % 1->3: excite (delay 30ms)
%interaction_from_j_to_neuron{3}(1,:) = exf .* u;        % 1->3: excite
%interaction_from_j_to_neuron{3}(1,:) = 0;
interaction_from_j_to_neuron{3}(2,:) = 0;
% Interaction matrix to neuron 4
interaction_from_j_to_neuron{4}(2,:) = exf .* u;        % 2->4: excite
interaction_from_j_to_neuron{4}(3,:) = exf .* u;        % 3->4: excite
%interaction_from_j_to_neuron{4}(2,:) = 0;
%interaction_from_j_to_neuron{4}(4,:) = ihf .* u;         % 4->4: auto
interaction_from_j_to_neuron{4}(4,:) = au .* u;
%interaction_from_j_to_neuron{4}(4,:) = 0;
interaction_from_j_to_neuron{4}(1,:) = 0;

%% STEP 1: Generate spontaneous activity for the next M = 180 ms (180 cells)
[dN, SPON_SPIKES, search_range] = step1(dN, f0, M, n, SPON_SPIKES);

%% STEP 2 + 3 + 4:
% If no spikes are observed, back to STEP 1
% If at least 1 spike is observed, point n -> the bin number of the earliest spike.
% then ignore bins n+1:M and compute firing probabilities for n+1:n+F bins
while sum(spike_count >= max_spike_no, "all") < 4
    s = sum(dN(:,n:search_range+n-1), "all");
    if s == 0           % no spikes
        n = search_range + n;
        if n >= chunk_size
            [dN, n, file_index] = periodic_save(dN, n, file_index, M);
        end
        [dN, SPON_SPIKES, search_range] = step1(dN, f0, M, n, SPON_SPIKES);
    elseif s > 0        % spikes
        % find the bin of the earliest spike
        flag = 0;
        for m = n:search_range+n-1
            for i = 1:4
                if dN(i,m) == 1
                    n = m + 1;
                    flag = 1;
                    spike_count = spike_count + dN(:,m);
                    break;
                end
            end
            if flag == 1
                break;
            end
        end
        % ignore bins n+1:M
        %dN(:,n:M) = 0;
        if n >= chunk_size
            [dN, n, file_index] = periodic_save(dN, n, file_index, M);
        end
        P = step4(dN, n, M, F, f0, interaction_from_j_to_neuron);
%% STEP 5: 
% If any probability > 1, set all bins n+1:n+F to 0 and back to STEP 1 with
% n=n+F+1
% If all probabilities <= 1, generate spikes according to the probabilities
% (absolute refractory period) first bins of new window = 0
% Then back to STEP 2
        %s = sum(P > 1, "all");
        s_all = sum(P >= 1, "all");
        %s_all = sum(P > 1, "all");
        if s_all > 0        % probability > 1
            s_rows = sum(P >= 1, 2);
            %s_rows = sum(P > 1, 2);
            s_rows = s_rows >= 1;
            EXCEED1 = EXCEED1 + s_rows;
            %save_file = sprintf("./exceed_1_probabilities/P%d.mat", exceed_1_counter);
            %save(save_file, "P", "-v7.3");
            dN(:,n:n+F-1) = 0;
            n = n + F;
            [dN, SPON_SPIKES, search_range] = step1(dN, f0, M, n, SPON_SPIKES);
            % back to step 2
        else            % all probabilities <=1
            NOTEXCEED1 = NOTEXCEED1 + 1;
            randomMatrix = rand(4, F);
            dN(:,n:n+F-1) = randomMatrix < P;
            for i = 1:4
                if (dN(i,n-1) == 1) && (dN(i,n) == 1)   % refractory period
                    dN(i,n) = 0;
                end
            end   
            ends_rows = sum(P >= 1, 2);
            % update search_range to F (no longer M) since we just
            % generated spikes for F time bins
            search_range = F;
            % Check if any spikes are generated
            s = sum(dN(:,n:search_range+n-1), "all");
            if s > 0
                FUNC_SPIKES = FUNC_SPIKES + 1;      % spikes generated because of functional connectivity
            end
            % back to step 2
        end
    end
end

% Final savehow to find the reflection of a function
[dN, n, file_index] = periodic_save(dN, n, file_index, M);

% Local functions
function u = unit_step(t)
    % UNIT_STEP Generates a unit step function for the input vector t.
    u = double(t >= 0);
end

function f = auto_func(t)
    f = 0 + (log(2500*t + exp(-10))) .* exp(-t/(15*10^(-3)));
    % experiment with new decay auto function
    %f = -5*exp(-500*t);
    %f = log((50/9)*t);
    %f(1) = -10;
    % logistic decay
    %f = -2 ./ (1 + 0.001*exp(70*t));
end

function f = excite_func(t)
    f = 2 * sin((2*pi*t)/(6*10^(-2))) .* exp(-t/(4*10^(-2)));
end

function f = inhibit_func(t)
    f = -3 * sin((2*pi*t)/(12*10^(-2))) .* exp(-t/(4*10^(-2)));
end
toc