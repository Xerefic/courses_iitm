function [wr, dr, phijr_phikr] = peak_picking(H, w)
    H = abs(H);
    [peak, wr] = findpeaks(H, w);
    half_power = 1/sqrt(2) * peak;
    [r,c] = size(wr);
    if (r > 1) | (c > 1)
        wr = wr(end);
        half_power = half_power(end);
    end
    ind_wr = find(w==wr);
    
    ind_wa = find(min(abs(H(1:ind_wr)-half_power)) == abs(H(1:ind_wr)-half_power));
    ind_wb = find(min(abs(H(ind_wr:end)-half_power)) == abs(H(ind_wr:end)-half_power)) + ind_wr-1;
    wa = w(ind_wa);
    wb = w(ind_wb);

    dr = (wb-wa)/wr;
    phijr_phikr = peak * dr * wr^2;
end