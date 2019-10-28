function [magH, f] = plot_digital2(b, a, n, fs, t)
    [H, f] = freqz(b, a, n, fs);
    figure;
    sgtitle(t);
    subplot(2, 1, 1);
    plot(f/1e6, unwrap(angle(H))*180/pi);
    xlabel("Frequency (MHz)")
    ylabel("Phase (degrees)")
    title('Phase Response');
    subplot(2, 1, 2);
    magH = 20*log10(abs(H));
    plot(f/1e6, magH);
    title('Magnitude Response');
    xlabel("Frequency (MHz)")
    ylabel("Magnitude (dB)")
end
