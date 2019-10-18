%%
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
Fs = 6*10^6;
%%
