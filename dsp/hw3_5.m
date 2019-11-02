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

[n_k, Wn, beta, ftype] = kaiserord(f_cut, mags, dev, sample_freq);
n_k = n_k +rem(n_k,2);
b_k = fir1(n_k,Wn,ftype,kaiser(n_k+1,beta),'noscale');
[H_k, f_k] = freqz(b_k, 1, 512, sample_freq);

[n_pm, f_pm, a_pm, w_pm] = firpmord(f_cut, mags, dev, sample_freq);
b_pm = firpm(n_pm, f_pm, a_pm,w_pm);
[H_pm, f_pm]  = freqz(b_pm, 1,512,sample_freq);
 
%% b
impz(b_k,1);
impz(b_pm, 1);
plotter(f_k, H_k, 'Kaiser');
plotter(f_pm, H_pm, 'Parks-McClellan');

%% c
% check_specs returns 1 if the specs are met and 0 if they are not
k_spec = check_specs(f_k, H_k)
pm_spec = check_specs(f_pm, H_pm)

%% d
weight_check = w_pm(1) == max([dev(1)/dev(2) dev(2)/dev(1)])

%% e

b_k = fir1(2+n_k,Wn,ftype,kaiser(2+n_k+1,beta),'noscale');
[H_k, f_k] = freqz(b_k, 1, 512, sample_freq);
[n_pm, f_pm, a_pm, w_pm] = firpmord(f_cut, mags, dev, sample_freq);
b_pm = firpm(2+n_pm, f_pm, a_pm,w_pm);
[H_pm, f_pm]  = freqz(b_pm, 1,512,sample_freq);

k_spec = check_specs(f_k, H_k)
% The Parks-McClellan doesn't meet specs at n = 58. It has a peak at
% -29.9748 dB in the passband.
pm_spec = check_specs(f_pm, H_pm)

% The Parks-McClellan meets specs at n = 60. 
[n_pm, f_pm, a_pm, w_pm] = firpmord(f_cut, mags, dev, sample_freq);
b_pm = firpm(4+n_pm, f_pm, a_pm,w_pm);
[H_pm, f_pm]  = freqz(b_pm, 1,512,sample_freq);
pm_spec = check_specs(f_pm, H_pm)

plotter(f_k, H_k, 'Kaiser 2n');
plotter(f_pm, H_pm, 'Parks-McClellan 2n');
%% Functions
function plotter(f, H, g_title)
    figure
    plot(f/(10^6), 20*log10(abs(H)))
    ylim([-50,2]);
    hold on;
    xlabel('Frequency (MHz)');

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

% returns 1 if specs are met
function check = check_specs(f, H)

    h = 20*log10(abs(H));
    pass = find((f > 1.2*10^6) & (f < 1.5*10^6));
    stop = find((f < 1*10^6) | (f > 1.6*10^6));
    check = ~(any(h(pass) > 2) | any(h(pass) < -2) | any(h(stop) > -30));
end
