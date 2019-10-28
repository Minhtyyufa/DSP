%%% DSP HW3- Do Hyung, Minh, Yuval
clear all; close all; clc;

%% Problem 1
%% (a)
r_p = 2;
r_s = 30;
fs = 6*1e6;
f_end = 3*1e6;

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
w_end = f2w(f_end);

B = w_p2 - w_p1
w_o = sqrt(w_p2 * w_p1)

w_p = w2wproto_BP(w_p1, w_o, B);
w_s_lp1 = w2wproto_BP(w_s1, w_o, B);
w_s_lp2 = w2wproto_BP(w_s2, w_o, B);

figure
subplot(2,1,1)

line([w_p1, w_p2], [0, 0],'Color', 'red', 'LineStyle', '--');
line([w_p1, w_p2], [-r_p, -r_p],'Color', 'red','LineStyle', '--');
line([0,w_s1], [-r_s,-r_s],'Color', 'red','LineStyle', '--');
line([w_s2, w_end], [-r_s,-r_s],'Color', 'red','LineStyle', '--');
xlim([0, w_end])
xticks([w_s1, w_p1, w_p2, w_s2])
xticklabels({['Ws_1 = ', num2str(round(w_s1/1e6, 1))], ['Wp_1 = ', num2str(round(w_p1/1e6, 1))], ['Wp_2 = ', num2str(round(w_p2/1e6, 1))], ['Ws_2 = ', num2str(round(w_s2/1e6, 1))]})
title('Prewarped')
xlabel('frequency (Mrad/s)')
ylim([-(abs(r_s)+1), 1])

subplot(2,1,2)
line([0, w_p],  [0,0], 'LineStyle', '--');
line([0, w_p],  [-r_p,-r_p], 'LineStyle', '--');
line([min(w_s_lp1, w_s_lp2), pi], [-r_s, -r_s], 'LineStyle', '--')
title('LP Prototype')
ylim([-(abs(r_s)+1), 1])
xticks([1, min(w_s_lp1, w_s_lp2), max(w_s_lp1, w_s_lp2)])
xticklabels({'Wp = 1', ['Ws_1 = ',num2str(min(w_s_lp1, w_s_lp2))] , ...
    ['Ws_2 = ',num2str(max(w_s_lp1, w_s_lp2))]})
xlim([0,pi])
sgtitle('Design Specs')

%% (b)
butter_order = @(r_s, r_p, w_s_lp1, w_s_lp2, w_p) ...
    (log10((10^(r_s/10)-1)/(10^(r_p/10)-1))/log10(min(w_s_lp1, w_s_lp2)/w_p))/2;
chev_order   = @(r_s, r_p, w_s_lp1, w_s_lp2, w_p) ...
    acosh(sqrt((10^(r_s/10)-1)/ (10^(r_p/10)-1)))/ acosh(min(w_s_lp1, w_s_lp2)/w_p);
n_butter = butter_order(r_s, r_p, w_s_lp1, w_s_lp2, w_p);
n_chev   = chev_order(r_s, r_p, w_s_lp1, w_s_lp2, w_p);

n_lp_butter = ceil(n_butter)
n_bp_butter = 2*n_lp_butter
n_lp_chev = ceil(n_chev)
n_bp_chev = 2*n_lp_chev

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
title('Ellip for LP')
%%
% (f)
close all;

[H_b, w_b] = plot_analog(z_b, p_b, k_b, 'butter');
plotspecs([w_s1, w_s2], [w_p1, w_p2], w_b(end), r_s, r_p)


[H_c1, w_c1] = plot_analog(z_c1, p_c1, k_c1, 'cheb1');
plotspecs([w_s1, w_s2], [w_p1, w_p2], w_c1(end), r_s, r_p)

[H_c2, w_c2] = plot_analog(z_c2, p_c2, k_c2, 'cheb2');
plotspecs([w_s1, w_s2], [w_p1, w_p2], w_c2(end), r_s, r_p)

[H_e, w_e] = plot_analog(z_e, p_e, k_e, 'ellip');
plotspecs([w_s1, w_s2], [w_p1, w_p2], w_e(end), r_s, r_p)


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
table(r_s_b,r_s_c1, r_s_c2, r_s_e)

%% Problem 2
f2wn = @(f, fs) 2*f/fs;     %nyquist
f_end = 3*1e6;
wn_end =  f2wn(f_end, fs);

wn_p1 = f2wn(f_p1, fs);
wn_p2 = f2wn(f_p2, fs);
wn_s1 = f2wn(f_s1, fs);
wn_s2 = f2wn(f_s2, fs);

Bn = wn_p2 - wn_p1;
wn_o = sqrt(wn_p2 * wn_p1);

wn_p = w2wproto_BP(wn_p1, wn_o, Bn);
wn_s_lp1 = w2wproto_BP(wn_s1, wn_o, Bn);
wn_s_lp2 = w2wproto_BP(wn_s2, wn_o, Bn);
%%  (a)
figure
subplot(2,1,1)
plotspecs([wn_s1, wn_s2]*1e6, [wn_p1, wn_p2]*1e6, wn_end*1e6, r_s, r_p);
title('Prewarped')
xlabel('frequency (rad/s)')
xticks([wn_s1, wn_p1, wn_p2, wn_s2])
xticklabels({['Ws_1 = ', num2str(round(wn_s1, 2))], ['Wp_1 = ', num2str(round(wn_p1, 2))], ['Wp_2 = ', num2str(round(wn_p2, 2))], ['Ws_2 = ', num2str(round(wn_s2, 2))]})
ylim([-(abs(r_s)+1), 1])

subplot(2,1,2)
line([0, wn_p],  [0,0], 'LineStyle', '--');
line([0, wn_p],  [-r_p,-r_p], 'LineStyle', '--');
line([min(wn_s_lp1, wn_s_lp2), pi], [-r_s, -r_s], 'LineStyle', '--')
title('LP Prototype')
ylim([-(abs(r_s)+1), 1])
xticks([1, min(wn_s_lp1, wn_s_lp2), max(wn_s_lp1, wn_s_lp2)])
xticklabels({'Wp = 1', ['Ws_1 = ', num2str(min(wn_s_lp1, wn_s_lp2))], ['Ws_2', num2str(max(wn_s_lp1, wn_s_lp2))]})
xlim([0,pi])
sgtitle('Design Specs')

%%
% (b)
[nD_b, WnD_b] = buttord([wn_p1 wn_p2], [wn_s1 wn_s2], r_p, r_s);
[zD_b, pD_b, kD_b] = butter(nD_b, WnD_b);
[nD_c1, WnD_c1] = cheb1ord([wn_p1 wn_p2], [wn_s1 wn_s2], r_p, r_s);
[zD_c1, pD_c1, kD_c1] = cheby1(nD_c1, r_p, [wn_p1 wn_p2]);
[nD_c2, WnD_c2] = cheb2ord([wn_p1 wn_p2], [wn_s1 wn_s2], r_p, r_s);
[zD_c2, pD_c2, kD_c2] = cheby2(nD_c2, r_s, [wn_s1 wn_s2]);
[nD_e, WnD_e] = ellipord([wn_p1 wn_p2], [wn_s1 wn_s2], r_p, r_s);
[zD_e, pD_e, kD_e] = ellip(nD_e, r_p, r_s, [wn_p1 wn_p2]);
table(nD_b, nD_c1, nD_c2, nD_e)

%%
% (c)
clc; close all;
n = 1e5;
[H_b, w_b] = plot_digital(zD_b, pD_b, kD_b, n ,fs, 'butter');
plotspecs([f_s1, f_s2], [f_p1, f_p2],f_end , r_s, r_p);
[H_c1, w_c1] = plot_digital(zD_c1, pD_c1, kD_c1, n, fs, 'cheb1');
plotspecs([f_s1, f_s2], [f_p1, f_p2],f_end , r_s, r_p);
[H_c2, w_c2] = plot_digital(zD_c2, pD_c2, kD_c2, n, fs, 'cheb2');
plotspecs([f_s1, f_s2], [f_p1, f_p2],f_end , r_s, r_p);
[H_e, w_e] = plot_digital(zD_e, pD_e, kD_e, n,fs, 'ellip');
plotspecs([f_s1, f_s2], [f_p1, f_p2],f_end , r_s, r_p);

%%
% (d)
%close all;
figure;
subplot(2, 2, 1);
zplane(zD_b, pD_b);
title('butterworth');

subplot(2, 2, 2);
zplane(zD_c1, pD_c1);
title('cheby1');

subplot(2, 2, 3);
zplane(zD_c2, pD_c2);
title('cheby2');

subplot(2, 2, 4);
zplane(zD_e, pD_e);
title('ellip');

%% Problem 3

%(a)
[b_c1, a_c1] = zp2tf(z_c1, p_c1, k_c1);
[b_c2, a_c2] = zp2tf(z_c2, p_c2, k_c2);

[bz_c1, az_c1] = impinvar(b_c1, a_c1, fs);
[bz_c2, az_c2] = impinvar(b_c2, a_c2, fs);

% size(az_c1, 2)-1 == n_bp_chev
% size(az_c2, 2)-1 == n_bp_chev
% size(bz_c1, 2)-1 == n_bp_chev
% size(bz_c2, 2)-1 == n_bp_chev
%%% Indeed, they have the order is the same.

%% (b)
[H_ii_c1, f_ii_c1] = plot_digital2(bz_c1, az_c1, n, fs, "Cheby 1");
xlim([f_p1, f_p2]/1e6)
%We see that the passbands isn't equiripple after zooming in
[H_ii_c2, f_ii_c2] = plot_digital2(bz_c2, az_c2,n,  fs, "Cheby 2");
%No equiripple here


%Checking attenuation 
% Cheby 1
%stopband
is1 = find(f_ii_c1 >= f_s1, 1);
rs_ii_c1_1 = abs(max(H_ii_c1(1:is1)));
is2 = find(f_ii_c1 >= f_s2, 1);
rs_ii_c1_2 = abs(max(H_ii_c1(is2:end)));
%passband
ip1 = find(f_ii_c1 >= f_p1, 1);
rp_ii_c1_1 = abs(H_ii_c1(ip1));
ip2 = find(f_ii_c1 >= f_p2, 1);
rp_ii_c1_2 = abs(H_ii_c1(ip2));

% Cheby 2
%stopband
is3 = find(f_ii_c2>= f_s1, 1);
rs_ii_c2_1 = abs(max(H_ii_c2(1:is3)));
is4 = find(f_ii_c2>= f_s2, 1);
rs_ii_c2_2 = abs(max(H_ii_c2(is4:end)));
%passband
ip3 = find(f_ii_c2>= f_p1, 1);
rp_ii_c2_1 = abs(H_ii_c2(ip3));
ip4 = find(f_ii_c2>= f_p2, 1);
rp_ii_c2_2 = abs(H_ii_c2(ip4));

table(rp_ii_c1_1, rp_ii_c1_2, rp_ii_c2_1, rp_ii_c2_1)   %passband check
table(rs_ii_c1_1, rs_ii_c1_2, rs_ii_c2_1, rs_ii_c2_2)   %stopband check

%% (c) ERROR- not sure what he meant with poles/zeros in and on unit circle
[z_ii_c1, p_ii_c1, k_ii_c1] = tf2zp(bz_c1, az_c1);
[z_ii_c2, p_ii_c2, k_ii_c2] = tf2zp(bz_c2, az_c2);
figure
subplot(2, 1,1)
zplane(z_ii_c1, p_ii_c1);
title('Cheby 1')
subplot(2, 1,2)
zplane(z_ii_c2, p_ii_c2);
title('Cheby 2')

% The poles stay within the unit circle
% The zeros on the unit circle are no longer there, some stayed in the
% surroundings of where they were for the not Impulse Invariant filter,
% while others moved a lot.

%% (d)
%analog
[A_c1, B_c1, C_c1, D_c1] = zp2ss(z_c1, p_c1, k_c1);
[A_c2, B_c2, C_c2, D_c2] = zp2ss(z_c2, p_c2, k_c2);

[b_c1, a_c1] = zp2tf(z_c1, p_c1, k_c1);
[y_c1, t_c1_a] = impulse(tf(b_c1, a_c1));

%[y_c1, t_c1_a] = impulse(ss(A_c1, B_c1, C_c1, D_c1));
%[y_c2, t_c2_a] = impulse(ss(A_c2, B_c2, C_c2, D_c2));

n1 = round(t_c1_a(end)*fs);
n2 = round(t_c2_a(end)*fs);

[h_c1, t_c1] = impz(bz_c1, az_c1, n1, fs);
[h_c2, t_c2] = impz(bz_c1, az_c1, n2, fs);

figure
subplot(2,1,1)
plot(t_c1_a, y_c1)
title('IR for Analog Filter')
xlabel('t')
ylabel('magnitude')
subplot(2,1,2)
plot(t_c1, h_c1)
title('IR for Digital IIR Filter')
xlabel('t')
ylabel('magnitude')
sgtitle('Cheby 1')

figure
subplot(2,1,1)
plot(t_c2_a, y_c2)
title('IR for Analog Filter')
xlabel('t')
ylabel('magnitude')
subplot(2,1,2)
plot(t_c2, h_c2)
title('IR for Digital IIR Filter')
xlabel('t')
ylabel('magnitude')
sgtitle('Cheby 2')
