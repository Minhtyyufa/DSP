function [magH, w] = plot_digital(z, p, k, fs, t)
    [b, a] = zp2tf(z, p, k);
    [H, w] = freqz(b, a);
    figure;
    sgtitle(t);
    subplot(2, 1, 1);
    plot(w*fs/1e6, unwrap(angle(H))*180/pi);
    xlabel("Frequency (MHz)")
    ylabel("Phase (degrees)")
    title('Phase Response');
    subplot(2, 1, 2);
    magH = 20*log10(abs(H));
    plot(w*fs/1e6, magH);
    title('Magnitude Response');
    xlabel("Frequency (MHz)")
    ylabel("Magnitude (dB)")
end

