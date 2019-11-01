% Q4
% (a)
close all; clear all; clc;

L = 31;
r = 30;
w_c = chebwin(L, r);
w_c = w_c/sum(w_c);
%%
% (b)
figure;
zplane(roots(w_c));
% (c)
ang_c = abs(angle(roots(w_c)));
n2n_w_c = 2*min(ang_c);
%%
% (d)
n = 1000;
w_c_f= abs(fftshift(fft(w_c, n)));
mag_c = 20*log10(w_c_f);
figure;
plot(mag_c);
title("Magnitude Response of Chebyshev Window");
xlabel("Radian (rad)")
ylabel("Magnitude (dB)")

%%
% (e)
beta = 3.13;
w_k = kaiser(L, beta);
w_k = w_k/sum(w_k);
ang_k = abs(angle(roots(w_k)));
n2n_w_k = 2*min(ang_k);
%%
% (f)
w_k_f = abs(fftshift(fft(w_k, 1000)));
mag_k = 20*log10(w_k_f);
figure;
plot(mag_c);
hold on;
plot(mag_k);
title("Magnitude Response of Chebyshev and Kaiser Window");
legend("Chebyshev Window", "Kaiser Window");
xlabel("Radian (rad)")
ylabel("Magnitude (dB)")

figure;
zplane(roots(w_c), roots(w_k));
title("Zeros of Chebyshev and Kaiser Window");
legend("Chebyshev Window", "Kaiser Window");

%%  
% (h)
[pks_k, loc_k] = findpeaks(mag_c);
[pk, I] = max(pks_k);
sl_level = pk - pks_k(I+1)
%%
% (i)
total_energy = sum(w_c_f.^2);
ind_max_pk = find(w_c_f(loc_k) == pk);
loc_min= islocalmin(w_c_f);
right_null = find(loc_min(loc_k(I):end) == 1 , 1) + loc_k(I);
left_null =  find(loc_min(1:loc_k(I)) == 1);
left_null = left_null(end);
main_lobe_energy = sum(w_c_f(left_null:right_null).^2);
side_lobe_energy_fraction = (total_energy - main_lobe_energy) / total_energy;