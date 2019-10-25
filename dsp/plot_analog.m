function [magH, w] = plot_analog(z, p, k, t)
    f2w = @(f) 2*pi*f;
    w = 0:100:f2w(3*10^6);
    [b, a] = zp2tf(z, p, k);
    [H, w] = freqs(b, a, w);
    figure;
    sgtitle(t);
    subplot(2, 1, 1);
    plot(w, unwrap(angle(H))*180/pi);
    title('Phase Response');
    subplot(2, 1, 2);
    magH = 20*log10(abs(H));
    plot(w, magH);
    title('Magnitude Response');
end

