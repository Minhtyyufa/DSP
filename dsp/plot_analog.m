function [magH, w] = plot_analog(z, p, k, t)
    f2w = @(f) 2*pi*f;
    w = 0:100:f2w(3*10^6);
    [b, a] = zp2tf(z, p, k);
    H = freqs(b, a, w);
    figure;
    
    sgtitle(t);
    subplot(2, 1, 1);
    plot(w/10^6, unwrap(angle(H))*180/pi);
    xlim([0, w(end)/10^6]);
    title('Phase Response');
    xlabel('Frequency (Mrad/s)')
    ylabel('Phase (degrees)')
    
    subplot(2, 1, 2);
    magH = 20*log10(abs(H));
    plot(w/10^6, magH);
    xlim([0, w(end)/10^6]);
    title('Magnitude Response');
    xlabel('Frequency (Mrad/s)')
    ylabel('Magnitude (dB)')
    ylim([-100, 2])

end

