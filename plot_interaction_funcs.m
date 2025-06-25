t = (0:179) * 1e-3;                 % Define the time variable
au = auto_func(t);
exf = excite_func(t);
exf_delayed = excite_func(t - 3*10^(-2));
ihf = inhibit_func(t);
ihf_delayed = inhibit_func(t - 6*10^(-2));
u = unit_step(t);
u_delayed_3 = unit_step(t - 3*10^(-2));
u_delayed_6 = unit_step(t - 6*10^(-2));

% Multiply func and unit_step_func element-wise
auto = au .* u;
excite = exf .* u;
delayed_excite = exf_delayed .* u_delayed_3;
inhibit = ihf .* u;
delayed_inhibit = ihf_delayed .* u_delayed_6;

% Plot auto function
%plot(t, auto, "LineWidth",1.5, "Color", [1, 0.5, 0]);     % Plot auto as a function of t
%
plot(t, auto, "LineWidth",1.5);
ylim([-2 2])
xlabel('Time (s)');                 % Label the x-axis
ylabel('auto(t)');                  % Label the y-axis
title('Logistic auto-interaction function'); % Title of the plot
grid on;
%

% Plot excite function
%{
plot(t, excite, "LineWidth",1.5);   % Plot auto as a function of t
ylim([-2 2])
xlabel('Time (s)');                     % Label the x-axis
ylabel('excite(t)');                % Label the y-axis
title('Excitation function');       % Title of the plot
grid on;
%}

% Plot delayed_excite function
%{
plot(t, delayed_excite, "LineWidth",1.5);   % Plot auto as a function of t
ylim([-2 2])
xlabel('Time (s)');                     % Label the x-axis
ylabel('delayed excite(t)');                % Label the y-axis
title('Delayed-excitation function');       % Title of the plot
grid on;
%}

% Plot inhibit function
%{
plot(t, inhibit, "LineWidth",1.5);   % Plot auto as a function of t
ylim([-2 2]);
xlabel('Time (s)');                     % Label the x-axis
ylabel('inhibit(t)');                % Label the y-axis
title('Inhibition function');       % Title of the plot
grid on;
%}

% Plot delayed_inhibit function
%{
plot(t, delayed_inhibit, "LineWidth",1.5);   % Plot auto as a function of t
ylim([-2 2]);
xlabel('Time (s)');                     % Label the x-axis
ylabel('delayed inhibit(t)');                % Label the y-axis
title('Delayed-inhibition function');       % Title of the plot
grid on;
%}


% Local functions
function u = unit_step(t)
    % UNIT_STEP Generates a unit step function for the input vector t.
    u = double(t >= 0);
end

function f = auto_func(t)
    %f = -2 + log(2500*t + exp(-10)) .* exp(-t/(15*10^(-3)));
    f = -2 ./ (1 + 0.001*exp(70*t));
end

function f = excite_func(t)
    f = (2 * sin((2*pi*t)/(6*10^(-2)))) .* exp(-t/(4*10^(-2)));
end

function f = inhibit_func(t)
    f = (-3 * sin((2*pi*t)/(12*10^(-2)))) .* exp(-t/(4*10^(-2)));
end
