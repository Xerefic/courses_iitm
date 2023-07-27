function [wr, dr, phijr_phikr] = circle_fit(H, w)
    rh = real(H);
    ih = imag(H);

    slopes = zeros(size(w));
    for i = 1:1:size(w, 2)
        slopes(i) = abs(ih(i)/rh(i));
    end
    
    rate = zeros(size(w));
    for i = 1:1:size(w, 2)-1
        del_theta = (atan(slopes(i+1)) - atan(slopes(i)))*2;
        del_w2 = w(i+1)^2 - w(i)^2;
        rate(i) = del_theta/del_w2;
    end
    [peak, wr] = findpeaks(rate, w);

    [r,c] = size(wr);
    if (r > 1) | (c > 1)
        wr = wr(end);
    end
    ind_wr = find(w==wr);
    
    yrxr = slopes(ind_wr);
    y1x1 = (yrxr-tan(deg2rad(135/2)))/(1+tan(deg2rad(135/2))*yrxr);
    ind_wa = find(min(abs(slopes(1:ind_wr)-y1x1)) == abs(slopes(1:ind_wr)-y1x1));
    ind_wb = find(min(abs(slopes(ind_wr:end)-y1x1)) == abs(slopes(ind_wr:end)-y1x1)) + ind_wr-1;
    
    x = rh(ind_wa:ind_wb); y = ih(ind_wa:ind_wb);
    x = x(:); y = y(:);
    a = [x y ones(size(x))]\[-(x.^2+y.^2)];
    xc = -.5*a(1);
    yc = -.5*a(2);
    R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
    
    alpha = atan((ih(ind_wr)-yc)/(rh(ind_wr)-xc));
    xa = xc + R*sin(alpha); ya = yc - R*cos(alpha);
    xb = xc - R*sin(alpha); yb = yc + R*cos(alpha);

    ind_wa = find(min(abs(slopes(1:ind_wr)-abs(ya/xa))) == abs(slopes(1:ind_wr)-abs(ya/xa)));
    ind_wb = find(min(abs(slopes(ind_wr:end)-abs(yb/xb))) == abs(slopes(ind_wr:end)-abs(yb/xb))) + ind_wr-1;
    wa = w(ind_wa);
    wb = w(ind_wb);

    dr = (wb^2-wa^2)/(2*wr^2);
    phijr_phikr = 2 * R * wr^2 * dr;
end