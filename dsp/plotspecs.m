function plotspecs(ws, wp, wb_end, r_s, r_p)
    w_p1 = wp(1);
    w_p2 = wp(2);
    w_s1 = ws(1);
    
    w_s2 = ws(2);
    line([w_p1, w_p2]/1e6, [0, 0],'Color', 'red', 'LineStyle', '--');
    line([w_p1, w_p2]/1e6, [-r_p, -r_p],'Color', 'red','LineStyle', '--');
    line([0,w_s1]/1e6, [-r_s,-r_s],'Color', 'red','LineStyle', '--');
    line([w_s2, wb_end]/1e6, [-r_s,-r_s],'Color', 'red','LineStyle', '--');
end

