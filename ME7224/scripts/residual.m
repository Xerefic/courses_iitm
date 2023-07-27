function [mr, kr] = residual(H_actual, H_generated, w)
    H_actual = abs(H_actual);
    H_generated = abs(H_generated);
    length = floor(size(w, 2)/10);
    w_low = w(1:length); 
    H_actual_low = H_actual(1:length);
    H_actual_low = H_actual_low(:);
    H_generated_low = H_generated(1:length);
    H_generated_low = H_generated_low(:);

    w_high = w(end-length:end);
    H_actual_high = H_actual(end-length:end);
    H_actual_high = H_actual_high(:);
    H_generated_high = H_generated(end-length:end);
    H_generated_high = H_generated_high(:);

    kr = 1/ mean((H_actual_high-H_generated_high));
    mr = mean(-1./w_low.^2) / (mean(H_actual_low-H_generated_low)-1/kr);
end
