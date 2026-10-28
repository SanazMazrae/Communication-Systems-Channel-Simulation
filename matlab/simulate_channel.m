function [ch_output, info] = simulate_channel(ch_input, ch_type, snr, Tm, fd)
%SIMULATE_CHANNEL Discrete-time complex-baseband fading channel simulator.
%   y = simulate_channel(x, channel_type, snr_dB, Ts, fd_Hz)
%   [y, info] = simulate_channel(...) additionally returns channel taps,
%   noiseless receiver samples, noise samples, and model parameters.
%
%   Supported channel_type values:
%     'Ray_TI_FF'  static Rayleigh flat fading
%     'Ray_TI_FS'  static Rayleigh three-tap frequency-selective fading
%     'Ray_TV_FF'  time-varying Rayleigh flat fading
%     'Ray_TV_FS'  time-varying Rayleigh three-tap fading
%     'Ric_TV_FS'  time-varying Rician three-tap fading (K=10 linear)
%     'AWGN'       unfaded additive complex Gaussian noise
%     'Awgn'       accepted alias for the course specification
%
%   Ts is the sampling interval of x (seconds); fd is the maximum
%   diffuse-scattering Doppler frequency (Hz). For sampled symbols,
%   Ts equals the symbol interval only when there is one sample/symbol.
%
%   For all modes, y[n] = sum_l h_l[n] * x[n-l] + w[n].  Delays are
%   integer multiples of Ts. Rayleigh paths have zero-mean unit-variance
%   complex Gaussian gains for static cases. Time-varying diffuse taps
%   use 32 random-angle sinusoids to approximate isotropic scattering
%   and J0(2*pi*fd*lag*Ts) temporal autocorrelation.
%
%   The power-delay profile [0.5,0.3,0.2] is a normalized, illustrative
%   modeling choice, not measured propagation data. Rician LOS power is
%   K/(K+1), and diffuse power is 1/(K+1), shared by the three taps.
%
%   Noise variance = mean(abs(x).^2) / 10^(snr_dB/10). This is an
%   INPUT-REFERENCED expected SNR: the realized received SNR can differ
%   due to fading and finite-length observations. RNG state is controlled
%   by the caller (e.g., rng(2025)).

    supported = {'Ray_TI_FF','Ray_TI_FS','Ray_TV_FF', ...
                 'Ray_TV_FS','Ric_TV_FS','AWGN'};
    if ~(isnumeric(ch_input) && isvector(ch_input) && ~isempty(ch_input) ...
         && all(isfinite(ch_input(:))))
        error('ch_input must be a nonempty finite numeric vector.');
    end
    if ~(ischar(ch_type) || (isstring(ch_type) && isscalar(ch_type)))
        error('ch_type must be a channel-name string.');
    end
    ch_type = char(ch_type);
    if strcmp(ch_type,'Awgn'), ch_type='AWGN'; end
    if ~any(strcmp(ch_type, supported))
        error('Unsupported channel type: %s', ch_type);
    end
    if ~(isscalar(snr) && isnumeric(snr) && isreal(snr) && isfinite(snr))
        error('snr must be a finite scalar in dB.');
    end
    if ~(isscalar(Tm) && isnumeric(Tm) && isreal(Tm) && isfinite(Tm) && Tm>0)
        error('Tm must be a positive sampling interval in seconds.');
    end
    if ~(isscalar(fd) && isnumeric(fd) && isreal(fd) && isfinite(fd) && fd>=0)
        error('fd must be a nonnegative maximum Doppler frequency in Hz.');
    end
    if fd * Tm >= 0.5 && ~strcmp(ch_type, 'AWGN')
        warning('fd*Tm >= 0.5: sampled Doppler process can alias.');
    end

    was_row = isrow(ch_input);
    x = ch_input(:);
    N = numel(x);
    profile = [0.5 0.3 0.2];
    K = NaN;

    if strcmp(ch_type, 'AWGN')
        h = ones(N,1);
    elseif strcmp(ch_type, 'Ray_TI_FF')
        h = repmat(complex_normal() ,N,1);
    elseif strcmp(ch_type, 'Ray_TI_FS')
        h = zeros(N,3);
        for ell=1:3
            h(:,ell) = sqrt(profile(ell)) * complex_normal();
        end
    elseif strcmp(ch_type, 'Ray_TV_FF')
        h = diffuse_fading(N,Tm,fd);
    elseif strcmp(ch_type, 'Ray_TV_FS')
        h = zeros(N,3);
        for ell=1:3
            h(:,ell) = sqrt(profile(ell)) * diffuse_fading(N,Tm,fd);
        end
    else % Ric_TV_FS
        K = 10;  % linear, equivalent to 10 dB
        h = zeros(N,3);
        for ell=1:3
            h(:,ell) = sqrt(profile(ell)/(K+1)) * diffuse_fading(N,Tm,fd);
        end
        t = (0:N-1).' * Tm;
        % LOS occupies the undelayed path; its Doppler is a model choice.
        h(:,1) = h(:,1) + sqrt(K/(K+1)) * exp(1i*2*pi*(fd/2)*t);
    end

    L = size(h,2);
    noiseless = complex(zeros(N,1));
    for ell=1:L
        delay = ell-1;
        if N>delay
            noiseless((delay+1):N) = noiseless((delay+1):N) + ...
                h((delay+1):N,ell) .* x(1:(N-delay));
        end
    end

    input_power = mean(abs(x).^2);
    noise_var = input_power / 10^(snr/10);
    w = sqrt(noise_var/2) * (randn(N,1) + 1i*randn(N,1));
    y = noiseless + w;
    if was_row
        ch_output = y.';
    else
        ch_output = y;
    end
    if nargout>1
        info = struct('channel_type',ch_type,'taps',h, ...
            'delay_samples',0:(L-1),'noiseless',noiseless, ...
            'noise',w,'noise_variance',noise_var, ...
            'sampling_interval_s',Tm,'maximum_doppler_Hz',fd, ...
            'rician_K_linear',K, ...
            'input_reference_snr_dB',snr, ...
            'realized_output_snr_dB',10*log10(mean(abs(noiseless).^2) ...
                                         /mean(abs(w).^2)));
    end
end

function z = complex_normal()
    z = (randn + 1i*randn)/sqrt(2);
end

function h = diffuse_fading(N,Ts,fd)
    % Sum of sinusoids with random arrival angles and independent phases.
    % Ensemble autocorrelation approximates J0(2*pi*fd*tau).
    if fd==0
        h = repmat(complex_normal(),N,1);
        return;
    end
    M=32;
    angles=2*pi*rand(1,M);
    phases=2*pi*rand(1,M);
    time=(0:N-1).' * Ts;
    h=sum(exp(1i*(2*pi*fd*(time*cos(angles)) + phases)),2)/sqrt(M);
end
