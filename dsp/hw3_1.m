%%
% Q1
% (a)
clear all; close all; clc;
r_p = 2;
r_s = 30;

f2w = @(f) 2*pi*f;
w2wproto_BP = @(w, w_o, B) abs((w^2 - w_o^2)/(B*w));

f_s1 = 10^6;
f_s2 = 1.6*10^6;
f_p1 = 1.2*10^6;
f_p2 = 1.5*10^6;
w_p1 = f2w(f_p1);
w_p2 = f2w(f_p2);
w_s1 = f2w(f_s1);
w_s2 = f2w(f_s2);

B = w_p2 - w_p1;
w_o = sqrt(w_p2 * w_p1);

w_p = w2wproto_BP(w_p1, w_o, B);
w_s_lp1 = w2wproto_BP(w_s1, w_o, B);
w_s_lp2 = w2wproto_BP(w_s2, w_o, B);

%%
% (b)
butter_order = @(r_s, r_p, w_s_lp1, w_s_lp2, w_p) (log10((10^(r_s/10)-1)/(10^(r_p/10)-1))/log10(min(w_s_lp1, w_s_lp2)/w_p))/2;
chev_order   = @(r_s, r_p, w_s_lp1, w_s_lp2, w_p) acosh(sqrt((10^(r_s/10)-1)/ (10^(r_p/10)-1)))/ acosh(min(w_s_lp1, w_s_lp2)/w_p);
n_butter = butter_order(r_s, r_p, w_s_lp1, w_s_lp2, w_p);
n_chev   = chev_order(r_s, r_p, w_s_lp1, w_s_lp2, w_p);

n_lp_butter = ceil(n_butter);
n_bp_butter = 2*n_lp_butter;
n_lp_chev = ceil(n_chev);
n_bp_chev = 2*n_lp_chev;

%%
% (c)
clc;
[n_b, Wn_b] = buttord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s, 's');
[z_b, p_b, k_b] = butter(n_b, Wn_b, 's');
[n_c1, Wn_c1] = cheb1ord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s, 's');
[z_c1, p_c1, k_c1] = cheby1(n_c1, r_p, [w_p1 w_p2], 's');
[n_c2, Wn_c2] = cheb2ord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s, 's');
[z_c2, p_c2, k_c2] = cheby2(n_c2, r_s, [w_s1 w_s2], 's');
[n_e, Wn_e] = ellipord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s, 's');
[z_e, p_e, k_e] = ellip(n_e, r_p, r_s, [w_p1 w_p2], 's');
%%
% (d)
% do rest
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

%[z_c1, p_c1, k_c1] = butter(n_c1, Wn_c1, 's');
%%
% (e)
close all;
[n_e_lp, Wn_e_lp] = ellipord(w_p, min([w_s_lp1 w_s_lp2]), r_p, r_s, 's');
[z_e_lp, p_e_lp, k_e_lp] = ellip(n_e_lp, r_p, r_s, min([w_s_lp1 w_s_lp2]), 's');
figure;
zplane(z_e_lp, p_e_lp);
%%
% (f)
clc; close all;
[H_b, w_b] = plot_analog(z_b, p_b, k_b, 'butter');
[H_c1, w_c1] = plot_analog(z_c1, p_c1, k_c1, 'cheb1');
[H_c2, w_c2] = plot_analog(z_c2, p_c2, k_c2, 'cheb2');
[H_e, w_e] = plot_analog(z_e, p_e, k_e, 'ellip');

%%
% (g)
figure;
plot(w_b, H_b, w_c1, H_c1, w_c2, H_c2, w_e, H_e);
legend("butterworth", "cheby1", "cheby2", "ellip");
xlim([w_p1, w_p2]);
title("zoomed in to magnitude response");

%%
% (h)
isb = find(w_b> w_s2, 1);
r_s_b = max(H_b(isb:end));
isc1 = find(w_c1> w_s2, 1);
r_s_c1 = max(H_c1(isc1:end));
isc2 = find(w_c2> w_s2, 1);
r_s_c2 = max(H_c2(isc2:end));
ise = find(w_e> w_s2, 1);
r_s_e = max(H_e(ise:end));

