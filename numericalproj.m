% Parameters
fs = 10;              % Sampling frequency (Hz)
t = 0:1/fs:2;         % Time vector (0 to 2 seconds, step = 1/fs)
f1 = 1;               % Frequency of the first signal (Hz)
f2 = 4;               % Frequency of the second signal (Hz)

% EEG Signal
x = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t);  % Combined signal

% Step 1: Numerical Differentiation using Forward Difference Method
h = 1/fs;                     % Time step
x_diff = diff(x) / h;         % Numerical differentiation
t_diff = t(1:end-1);          % Time vector for derivative

% Step 2: Discrete Fourier Transform (DFT)
N = length(x);                % Number of samples
X = fft(x);                   % Perform FFT
f = (0:N-1)*(fs/N);           % Frequency vector

% Magnitude Spectrum
X_magnitude = abs(X);

% Step 3: Frequency Band Analysis
delta_band = [0.5 4];          % Delta band range (Hz)
theta_band = [4 8];            % Theta band range (Hz)

% Identify Peaks in Frequency Bands
[~, delta_idx] = max(X_magnitude(f >= delta_band(1) & f <= delta_band(2)));
[~, theta_idx] = max(X_magnitude(f >= theta_band(1) & f <= theta_band(2)));

delta_peak_freq = f(find(f >= delta_band(1) & f <= delta_band(2), 1) + delta_idx - 1);
theta_peak_freq = f(find(f >= theta_band(1) & f <= theta_band(2), 1) + theta_idx - 1);

% Display Results
disp(['Delta peak frequency: ', num2str(delta_peak_freq), ' Hz']);
disp(['Theta peak frequency: ', num2str(theta_peak_freq), ' Hz']);

% Step 4: Plotting
figure;

% Original Signal
subplot(3,1,1);
plot(t, x);
title('EEG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Numerical Derivative
subplot(3,1,2);
plot(t_diff, x_diff);
title('Numerical Derivative of EEG Signal');
xlabel('Time (s)');
ylabel('Derivative');

% Magnitude Spectrum
subplot(3,1,3);
plot(f, X_magnitude);
title('Magnitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 fs/2]);  % Focus on frequencies of interest