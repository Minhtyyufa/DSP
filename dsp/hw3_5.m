clear all; clc; close all;
%% Part A
r_s2sig_s = @(r_s) 10^(-r_s/20);
r_p2sig_p = @(r_p) (10^(r_p/20)-1)/(10^(r_p/20)+1);
db2lin = @(db) 10^(db/20);

sample_freq = 6*10^6;
f_cut = [1*10^6 1.2*10^6 1.5*10^6 1.6*10^6];
mags = [0 1 0];
r_p = 2;
r_s = 30;
dev = [r_s2sig_s(r_s) r_p2sig_p(r_p) r_s2sig_s(r_s)];

[n, Wn, beta, ftype] = kaiserord(f_cut, mags, dev, sample_freq);
n = n +rem(n,2);
b_k = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
[H_k, f_k] = freqz(b_k, 1, 512, sample_freq);

[n_pm, f_pm, a_pm, w_pm] = firpmord(f_cut, mags, dev, sample_freq);
b_pm = firpm(n_pm, f_pm, a_pm,w_pm);
[H_pm, f_pm]  = freqz(b_pm, 1,1024,sample_freq);
 
%% b
[h_k t_k] = impz(b_k,1);
[h_pm t_pm] = impz(b_pm, 1);
plot_impulse(h_k, t_k, 'Kaiser');
plot_impulse(h_pm, t_pm, 'Parks-McClellan');


plotter(f_k, H_k, 'Kaiser');
plotter(f_pm, H_pm, 'Parks-McClellan');

%% c
% check_specs returns 0 if specs are not met and 1 if specs are met;

% For Kaiser
k_spec = check_specs(f_k, H_k)

% For Parks-McClellan
pm_spec = check_specs(f_pm, H_pm)

%% d
weight_check = w_pm(1) == max([dev(1)/dev(2) dev(2)/dev(1)])

%% e

b_k = fir1(2+n,Wn,ftype,kaiser(2+n+1,beta),'noscale');
[H_k, f_k] = freqz(b_k, 1, 512, sample_freq);
[n_pm, f_pm, a_pm, w_pm] = firpmord(f_cut, mags, dev, sample_freq);
b_pm = firpm(2+n_pm, f_pm, a_pm,w_pm);
[H_pm, f_pm]  = freqz(b_pm, 1,1024,sample_freq);

k_spec = check_specs(f_k, H_k)
pm_spec = check_specs(f_pm, H_pm)
plotter(f_k, H_k, 'Kaiser 2n');
plotter(f_pm, H_pm, 'Parks-McClellan 2n');
%% functions
function plotter(f, H, g_title)
    figure
    plot(f/(10^6), 20*log10(abs(H)))
    ylim([-50,2]);
    hold on;

    xl = xlim();
    line([1.2, 1.5], [2 2], 'Color','red', ...
        'LineStyle', '--');
    line([1.2, 1.5], [-2 -2], 'Color','red', ...
        'LineStyle', '--');
    line([xl(1), 1], [-30 -30], 'Color','red', ...
        'LineStyle', '--');
    line([1.6, xl(2)], [-30 -30], 'Color','red', ...
        'LineStyle', '--');

    title(g_title);
end
function plot_impulse(h, t, my_title)
    figure
    stem(t, h)
    title(my_title)
end

function check = check_specs(f, H)
    h = 20*log10(abs(H));
    pass = find((f > 1.2*10^6) & (f < 1.5*10^6));
    stop = find((f < 1*10^6) | (f > 1.6*10^6));
    
    pass_band_max_ripple = max(abs(h(pass)))
    stop_band_min_atten = min(abs(h(stop)))
    check = ~(any(h(pass) > 2) | any(h(pass) < -2) | any(h(stop) > -30));
end