%%
clear all; close all; clc;
r_p = 2;
r_s = 30;
fs = 6*1e6;

f2wn = @(f, fs) f/fs;
f2w = @(f) 2*pi*f;
w2wproto_BP = @(w, w_o, B) abs((w^2 - w_o^2)/(B*w));

f_s1 = 10^6;
f_s2 = 1.6*10^6;
f_p1 = 1.2*10^6;
f_p2 = 1.5*10^6;
w_p1 = f2wn(f_p1, fs);
w_p2 = f2wn(f_p2, fs);
w_s1 = f2wn(f_s1, fs);
w_s2 = f2wn(f_s2, fs);

B = w_p2 - w_p1;
w_o = sqrt(w_p2 * w_p1);

w_p = w2wproto_BP(w_p1, w_o, B);
w_s_lp1 = w2wproto_BP(w_s1, w_o, B);
w_s_lp2 = w2wproto_BP(w_s2, w_o, B);
%%
% (b)
[n_b, Wn_b] = buttord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s);
[z_b, p_b, k_b] = butter(n_b, Wn_b);
[n_c1, Wn_c1] = cheb1ord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s);
[z_c1, p_c1, k_c1] = cheby1(n_c1, r_p, [w_p1 w_p2]);
[n_c2, Wn_c2] = cheb2ord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s);
[z_c2, p_c2, k_c2] = cheby2(n_c2, r_s, [w_s1 w_s2]);
[n_e, Wn_e] = ellipord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s);
[z_e, p_e, k_e] = ellip(n_e, r_p, r_s, [w_p1 w_p2]);

%%
% (c)
clc; close all;
[H_b, w_b] = plot_digital(z_b, p_b, k_b, fs, 'butter');
[H_c1, w_c1] = plot_digital(z_c1, p_c1, k_c1, fs, 'cheb1');
[H_c2, w_c2] = plot_digital(z_c2, p_c2, k_c2, fs, 'cheb2');
[H_e, w_e] = plot_digital(z_e, p_e, k_e, fs, 'ellip');
%%
% (d)
close all;
figure;
subplot(2, 2, 1);
zplane(z_b, p_b);
title('butterworth');

subplot(2, 2, 2);
zplane(z_c1, p_c1);
title('cheby1');

subplot(2, 2, 3);
zplane(z_c2, p_c2);
title('cheby2');

subplot(2, 2, 4);
zplane(z_e, p_e);
title('ellip');