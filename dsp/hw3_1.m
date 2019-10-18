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
% do rest
clc;
[n_b, Wn_b] = buttord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s, 's');
[z_b, p_b, k_b] = butter(n_b, Wn_b, 's');
[n_c1, Wn_c1] = cheb1ord([w_p1 w_p2], [w_s1 w_s2], r_p, r_s, 's');
[z_c1, p_c1, k_c1] = cheby1(n_c1, r_p, [w_p1 w_p2], 's');
%%
% (d)
% do rest
zplane(z_b, p_b);
[z_c1, p_c1, k_c1] = butter(n_c1, Wn_c1, 's');
%%
% (e)

%%
% (f)
w = 0:100:f2w(3*10^6);
[b_b, a_b] = zp2tf(z_b, p_b, k_b);
figure;
freqs(b_b, a_b, w)

% [h_b, w_out_ b] = freqs(b_b, a_b, w);