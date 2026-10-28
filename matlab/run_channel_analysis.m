%% Communication Systems | Channel simulation and evaluation
% Run this script from the matlab/ folder or from any current directory.
% Requires base MATLAB and the accompanying simulate_channel.m.
% All figures and CSV tables are written to ../figures and ../results.
clear; close all; clc;
rng(2025,'twister');
here = fileparts(mfilename('fullpath'));
addpath(here);
out_fig = fullfile(here,'..','figures');
out_res = fullfile(here,'..','results');
if ~exist(out_fig,'dir'), mkdir(out_fig); end
if ~exist(out_res,'dir'), mkdir(out_res); end

N = 4096; Ts = 1e-3; fd = 50; snr_dB = 20;
ch_types = {'Ray_TI_FF','Ray_TI_FS','Ray_TV_FF', ...
            'Ray_TV_FS','Ric_TV_FS','AWGN'};
I = 2*randi([0,1],N,1)-1;
Q = 2*randi([0,1],N,1)-1;
x = complex(I,Q)/sqrt(2);  % unit-power QPSK input

%% Channel realizations, CIR snapshots, and signal/noise power
hinfo = cell(numel(ch_types),1);
Y = cell(numel(ch_types),1);
summary = zeros(numel(ch_types),5);
for k=1:numel(ch_types)
    rng(2025+k,'twister');
    [Y{k},hinfo{k}] = simulate_channel(x,ch_types{k},snr_dB,Ts,fd);
    m=hinfo{k};
    summary(k,:)=[size(m.taps,2),mean(abs(m.noiseless).^2), ...
        mean(abs(m.noise).^2),m.realized_output_snr_dB, ...
        sqrt(mean(abs(m.noise).^2)/mean(abs(m.noiseless).^2))];
end
fid=fopen(fullfile(out_res,'channel_summary.csv'),'w');
fprintf(fid,'channel,n_taps,received_clean_power,noise_power,measured_SNR_dB,noise_EVM\n');
for k=1:numel(ch_types)
    fprintf(fid,'%s,%d,%.9g,%.9g,%.6f,%.7f\n', ...
        ch_types{k},summary(k,:));
end
fclose(fid);

figure('Color','w','Position',[100 80 1160 610]);
for k=1:numel(ch_types)
    subplot(2,3,k);
    m=hinfo{k};
    stem(m.delay_samples,abs(m.taps(1,:)),'filled','LineWidth',1.2);
    xlim([-0.3 2.4]); ylim([0 max(1.3,1.1*max(abs(m.taps(1,:))))]);
    title(strrep(ch_types{k},'_','\_')); xlabel('Path delay (samples)');
    ylabel('|h_l[0]|'); grid on;
end
saveas(gcf,fullfile(out_fig,'01_channel_impulse_snapshots.png'));

figure('Color','w','Position',[100 100 1100 520]);
keep=[1 2 3 4 5];
for j=1:numel(keep)
    k=keep(j); subplot(numel(keep),1,j);
    plot((0:499)*Ts,abs(hinfo{k}.taps(1:500,1)),'k','LineWidth',1.0);
    ylabel('|h_0[n]|'); title(strrep(ch_types{k},'_','\_'));
    grid on; xlim([0 0.5]);
end
xlabel('Time (s)');
saveas(gcf,fullfile(out_fig,'02_fading_vs_time.png'));

figure('Color','w','Position',[100 100 1050 520]);
selected=[1 2 4 5]; ff=linspace(-0.5,0.5,512);
for j=1:numel(selected)
    k=selected(j); subplot(2,2,j); m=hinfo{k};
    ht=m.taps(1,:).';
    Hf=fftshift(fft(ht,512));
    plot(ff,20*log10(abs(Hf)+1e-10),'k','LineWidth',1.1);
    grid on; xlabel('Normalized frequency (cycles/sample)');
    ylabel('|H(f)| (dB)'); title(strrep(ch_types{k},'_','\_'));
end
saveas(gcf,fullfile(out_fig,'03_frequency_response_snapshots.png'));

%% Input-reference noise: SNR sweeps for Rayleigh, Rician and AWGN
snr_grid=[0 5 10 15 20 25 30]; reps=20;
noise_evm=zeros(numel(ch_types),numel(snr_grid));
for k=1:numel(ch_types)
    for s=1:numel(snr_grid)
        evm_runs=zeros(reps,1);
        for rep=1:reps
            rng(10000+1000*k+100*s+rep,'twister');
            [~,m]=simulate_channel(x,ch_types{k},snr_grid(s),Ts,fd);
            evm_runs(rep)=sqrt(mean(abs(m.noise).^2) / ...
                               mean(abs(m.noiseless).^2));
        end
        noise_evm(k,s)=mean(evm_runs);
    end
end
fid=fopen(fullfile(out_res,'noise_sweep.csv'),'w');
fprintf(fid,'channel,snr_dB,mean_noise_EVM,replicates\n');
for k=1:numel(ch_types)
    for s=1:numel(snr_grid)
        fprintf(fid,'%s,%g,%.8f,%d\n',ch_types{k},snr_grid(s), ...
            noise_evm(k,s),reps);
    end
end
fclose(fid);
figure('Color','w','Position',[100 100 850 500]); hold on;
for k=1:numel(ch_types)
    plot(snr_grid,noise_evm(k,:),'-o','LineWidth',1.25,'MarkerSize',4);
end
set(gca,'YScale','log'); grid on;
xlabel('Input-reference SNR (dB)'); ylabel('Noise-only EVM (RMS)');
legend(strrep(ch_types,'_','\_'),'Location','northeast');
saveas(gcf,fullfile(out_fig,'04_noise_sweep.png'));

%% Doppler effect: temporal correlation of the first Rayleigh path
fd_grid=[0 10 25 50 75 100 150];
correlation=zeros(size(fd_grid)); theo=zeros(size(fd_grid));
reps=15;
for f=1:numel(fd_grid)
    vals=zeros(reps,1);
    for rep=1:reps
        rng(20000+100*f+rep,'twister');
        [~,m]=simulate_channel(x,'Ray_TV_FF',20,Ts,fd_grid(f));
        h=m.taps(:,1);
        vals(rep)=real(sum(h(2:end).*conj(h(1:end-1))) ...
            /sum(abs(h(1:end-1)).^2));
    end
    correlation(f)=mean(vals);
    theo(f)=besselj(0,2*pi*fd_grid(f)*Ts);
end
fid=fopen(fullfile(out_res,'doppler_correlation.csv'),'w');
fprintf(fid,'fd_Hz,empirical_lag1_correlation,J0_prediction\n');
for f=1:numel(fd_grid)
    fprintf(fid,'%g,%.7f,%.7f\n',fd_grid(f),correlation(f),theo(f));
end
fclose(fid);
figure('Color','w','Position',[100 100 820 460]);
plot(fd_grid,correlation,'ko-','LineWidth',1.5); hold on;
plot(fd_grid,theo,'k--','LineWidth',1.4); grid on;
xlabel('Maximum Doppler (Hz)'); ylabel('First-lag correlation');
legend('Simulation','J_0(2\pi f_d T_s)','Location','southwest');
saveas(gcf,fullfile(out_fig,'05_doppler_correlation.png'));

%% Noise and fade visualization for a QPSK signal
figure('Color','w','Position',[100 100 1070 350]);
indices=[6 1 4];
for j=1:3
    k=indices(j);
    subplot(1,3,j); scatter(real(Y{k}(1:500)),imag(Y{k}(1:500)),7,'k','filled');
    axis equal; grid on; xlabel('In-phase'); ylabel('Quadrature');
    title(strrep(ch_types{k},'_','\_'));
end
saveas(gcf,fullfile(out_fig,'06_received_constellations.png'));

fprintf('Finished. Generated six figures and three result tables.\n');
for k=1:numel(ch_types)
    fprintf('%s: taps=%d, measured SNR=%.2f dB\n', ...
            ch_types{k},summary(k,1),summary(k,4));
end
