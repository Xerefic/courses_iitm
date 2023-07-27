clc; clear;

data = readmatrix('data.xlsx');
% data(:, 2:end) = data(:, 2:end)/4.448222;

modes = data(:, 1);
H11a = data(:, 2) + i .* data(:, 3); H11 = H11a ./ modes.^2;
H12a = data(:, 4) + i .* data(:, 5); H12 = H12a ./ modes.^2;
H13a = data(:, 6) + i .* data(:, 7); H13 = H13a ./ modes.^2;
H14a = data(:, 8) + i .* data(:, 9); H14 = H14a ./ modes.^2;
H15a = data(:, 10) + i .* data(:, 11); H15 = H15a ./ modes.^2;
H16a = data(:, 12) + i .* data(:, 13); H16 = H16a ./ modes.^2;
H17a = data(:, 14) + i .* data(:, 15); H17 = H17a ./ modes.^2;
H18a = data(:, 16) + i .* data(:, 17); H18 = H18a ./ modes.^2;


% start_ind = find(modes==5); end_ind = find(modes==450);
% hold on
% set(gcf,'units','points','position',[0,0,1250,400])
% plot(modes(start_ind:end_ind), imag(H11(start_ind:end_ind)))
% % plot(modes(start_ind:end_ind), real(H14(start_ind:end_ind)))
% xlabel('f')
% % ylabel('Re\{H_{11}(2 \pi f)z\}')
% ylabel('Im\{H_{11}(2 \pi f)\}')
% % ylabel('|H_{11}(2\pi f)|')
% title('Plot of Receptance FRF H_{11}')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mode2 = 92:0.5:97;
start_ind = find(modes==92); end_ind = find(modes==97);

H11_mode2 = H11(start_ind:end_ind);
[w2_11p, d2_11p, ph12_ph12p] = peak_picking(H11_mode2, mode2);
[w2_11c, d2_11c, ph12_ph12c] = circle_fit(H11_mode2, mode2);

H12_mode2 = H12(start_ind:end_ind);
[w2_12p, d2_12p, ph12_ph22p] = peak_picking(H12_mode2, mode2);
[w2_12c, d2_12c, ph12_ph22c] = circle_fit(H12_mode2, mode2);

H13_mode2 = H13(start_ind:end_ind);
[w2_13p, d2_13p, ph12_ph32p] = peak_picking(H13_mode2, mode2);
[w2_13c, d2_13c, ph12_ph32c] = circle_fit(H13_mode2, mode2);

H14_mode2 = H14(start_ind:end_ind);
[w2_14p, d2_14p, ph12_ph42p] = peak_picking(H14_mode2, mode2);
[w2_14c, d2_14c, ph12_ph42c] = circle_fit(H14_mode2, mode2);

H15_mode2 = H15(start_ind:end_ind);
[w2_15p, d2_15p, ph12_ph52p] = peak_picking(H15_mode2, mode2);
[w2_15c, d2_15c, ph12_ph52c] = circle_fit(H15_mode2, mode2);

H16_mode2 = H16(start_ind:end_ind);
[w2_16p, d2_16p, ph12_ph62p] = peak_picking(H16_mode2, mode2);
[w2_16c, d2_16c, ph12_ph62c] = circle_fit(H16_mode2, mode2);

H17_mode2 = H17(start_ind:end_ind);
[w2_17p, d2_17p, ph12_ph72p] = peak_picking(H17_mode2, mode2);
[w2_17c, d2_17c, ph12_ph72c] = circle_fit(H17_mode2, mode2);

H18_mode2 = H18(start_ind:end_ind);
[w2_18p, d2_18p, ph12_ph82p] = peak_picking(H18_mode2, mode2);
[w2_18c, d2_18c, ph12_ph82c] = circle_fit(H18_mode2, mode2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode3 = 285:0.5:290;
start_ind = find(modes==285); end_ind = find(modes==290);

H11_mode3 = H11(start_ind:end_ind);
[w3_11p, d3_11p, ph13_ph13p] = peak_picking(H11_mode3, mode3);
[w3_11c, d3_11c, ph13_ph13c] = circle_fit(H11_mode3, mode3);

H12_mode3 = H12(start_ind:end_ind);
[w3_12p, d3_12p, ph13_ph23p] = peak_picking(H12_mode3, mode3);
[w3_12c, d3_12c, ph13_ph23c] = circle_fit(H12_mode3, mode3);

H13_mode3 = H13(start_ind:end_ind);
[w3_13p, d3_13p, ph13_ph33p] = peak_picking(H13_mode3, mode3);
[w3_13c, d3_13c, ph13_ph33c] = circle_fit(H13_mode3, mode3);

H14_mode3 = H14(start_ind:end_ind);
[w3_14p, d3_14p, ph13_ph43p] = peak_picking(H14_mode3, mode3);
[w3_14c, d3_14c, ph13_ph43c] = circle_fit(H14_mode3, mode3);

H15_mode3 = H15(start_ind:end_ind);
[w3_15p, d3_15p, ph13_ph53p] = peak_picking(H15_mode3, mode3);
[w3_15c, d3_15c, ph13_ph53c] = circle_fit(H15_mode3, mode3);

H16_mode3 = H16(start_ind:end_ind);
[w3_16p, d3_16p, ph13_ph63p] = peak_picking(H16_mode3, mode3);
[w3_16c, d3_16c, ph13_ph63c] = circle_fit(H16_mode3, mode3);

H17_mode3 = H17(start_ind:end_ind);
[w3_17p, d3_17p, ph13_ph73p] = peak_picking(H17_mode3, mode3);
[w3_17c, d3_17c, ph13_ph73c] = circle_fit(H17_mode3, mode3);

H18_mode3 = H18(start_ind:end_ind);
[w3_18p, d3_18p, ph13_ph83p] = peak_picking(H18_mode3, mode3);
[w3_18c, d3_18c, ph13_ph83c] = circle_fit(H18_mode3, mode3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape_2p = [ph12_ph12p, ph12_ph22p, ph12_ph32p, ph12_ph42p, ph12_ph52p, ph12_ph62p, ph12_ph72p, ph12_ph82p] ./ sqrt(ph13_ph13p);
shape_2c = [ph12_ph12c, ph12_ph22c, ph12_ph32c, ph12_ph42c, ph12_ph52c, ph12_ph62c, ph12_ph72c, ph12_ph82c] ./ sqrt(ph13_ph13p);

shape_3p = [ph13_ph13p, ph13_ph23p, ph13_ph33p, ph13_ph43p, ph13_ph53p, ph13_ph63p, ph13_ph73p, ph13_ph83p] ./ sqrt(ph13_ph13p);
shape_3c = [ph13_ph13c, ph13_ph23c, ph13_ph33c, ph13_ph43c, ph13_ph53c, ph13_ph63c, ph13_ph73c, ph13_ph83c] ./ sqrt(ph13_ph13p);

hold on
plot(shape_2c, 'DisplayName', 'r=2')
plot(shape_3c, 'DisplayName', 'r=3')
title('Mode Shape for Circle Fit')
legend()
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 50:0.5:450;
start_ind = find(modes==50); end_ind = find(modes==450);
H11a_mode = H11a(start_ind:end_ind);
H15a_mode = H15a(start_ind:end_ind);

H11a_p = @(x) x^2 * (ph12_ph12p/((w2_11p^2-x^2)+1i*d2_11p*w2_11p*x) + ph13_ph13p/((w3_11p^2-x^2)+1i*d3_11p*w3_11p*x));
H11a_c = @(x) x^2 * (ph12_ph12c/((w2_11c^2-x^2)+1i*d2_11c*w2_11c*x) + ph13_ph13c/((w3_11c^2-x^2)+1i*d3_11c*w3_11c*x));

H15a_p = @(x) x^2 * (ph12_ph52p/((w2_15p^2-x^2)+1i*d2_15p*w2_15p*x) + ph13_ph53p/((w3_15p^2-x^2)+1i*d3_15p*w3_15p*x));
H15a_c = @(x) x^2 * (ph12_ph52c/((w2_15c^2-x^2)+1i*d2_15c*w2_15c*x) + ph13_ph53c/((w3_15c^2-x^2)+1i*d3_15c*w3_15c*x));

% hold on
% set(gcf,'units','points','position',[0,0,1250,400])
% plot(mode, abs(H15a_mode), 'DisplayName', 'Analytical')
% fplot(@(x) abs(H15a_p(x)), [50, 450], 'DisplayName', 'Peak Picking');
% fplot(@(x) abs(H15a_c(x)), [50, 450], 'DisplayName', 'Circle Fit');
% 
% xlabel('f')
% ylabel('|H_{15}^a(2 \pi f)|')
% legend()
% title('Plot of Accelerance FRF H_{15}')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 50:0.5:450;
start_ind = find(modes==50); end_ind = find(modes==450);
H12a_mode = H12a(start_ind:end_ind);
H16a_mode = H16a(start_ind:end_ind);

H12a_p = @(x) x^2 * (ph12_ph22p/((w2_12p^2-x^2)+1i*d2_12p*w2_12p*x) + ph13_ph23p/((w3_12p^2-x^2)+1i*d3_12p*w3_12p*x));
H12a_c = @(x) x^2 * (ph12_ph22c/((w2_12c^2-x^2)+1i*d2_12c*w2_12c*x) + ph13_ph23c/((w3_12c^2-x^2)+1i*d3_12c*w3_12c*x));
i = 1;
H12a_pg = zeros(size(mode));
H12a_cg = zeros(size(mode));
for w = mode
    H12a_pg(i) = H12a_p(w);
    H12a_cg(i) = H12a_c(w);
    i = i + 1;
end
[mr12p, kr12p] = residual(H12a_mode, H12a_pg, mode);
[mr12c, kr12c] = residual(H12a_mode, H12a_cg, mode);
H12a_p_corr = @(x) H12a_p(x) -1./(x.^2*mr12p) + 1/kr12p;
H12a_c_corr = @(x) H12a_c(x) -1./(x.^2*mr12c) + 1/kr12c;

H16a_p = @(x) x^2 * (ph12_ph62p/((w2_16p^2-x^2)+1i*d2_16p*w2_16p*x) + ph13_ph63p/((w3_16p^2-x^2)+1i*d3_16p*w3_16p*x));
H16a_c = @(x) x^2 * (ph12_ph62c/((w2_16c^2-x^2)+1i*d2_16c*w2_16c*x) + ph13_ph63c/((w3_16c^2-x^2)+1i*d3_16c*w3_16c*x));
i = 1;
H16a_pg = zeros(size(mode));
H16a_cg = zeros(size(mode));
for w = mode
    H16a_pg(i) = H16a_p(w);
    H16a_cg(i) = H16a_c(w);
    i = i + 1;
end
[mr16p, kr16p] = residual(H16a_mode, H16a_pg, mode);
[mr16c, kr16c] = residual(H16a_mode, H16a_cg, mode);
H16a_p_corr = @(x) H16a_p(x) -1/(x.^2*mr16p) + 1/kr16p;
H16a_c_corr = @(x) H16a_c(x) -1/(x.^2*mr16c) + 1/kr16c;


% hold on
% set(gcf,'units','points','position',[0,0,1250,400])
% plot(mode, abs(H16a_mode), 'DisplayName', 'Analytical')
% % plot(mode, abs(H12a_pg), 'DisplayName', 'Peak Picking');
% % plot(mode, abs(H12a_cg), 'DisplayName', 'Circle Fit');
% fplot(@(x) abs(H16a_p_corr(x)), [50, 450], 'DisplayName', 'Peak Picking');
% fplot(@(x) abs(H16a_c_corr(x)), [50, 450], 'DisplayName', 'Circle Fit');
% xlabel('f')
% ylabel('|H_{16}^a(2 \pi f)|')
% legend()
% title('Plot of Corrected Accelerance FRF H_{16}')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 5:0.5:125;
start_ind = find(modes==5); end_ind = find(modes==125);
H11_mode = H11(start_ind:end_ind);
H12_mode = H12(start_ind:end_ind);
H13_mode = H13(start_ind:end_ind);
H14_mode = H14(start_ind:end_ind);
H15_mode = H15(start_ind:end_ind);
H16_mode = H16(start_ind:end_ind);
H17_mode = H17(start_ind:end_ind);
H18_mode = H18(start_ind:end_ind);

m = 6;

[H11_rfp, modal_par] = rfp(H11_mode, mode, m);
freq11 = modal_par(:, 1);
damp11 = modal_par(:, 2);
Ci11 = modal_par(:, 3);

[H12_rfp, modal_par] = rfp(H12_mode, mode, m);
freq12 = modal_par(:, 1);
damp12 = modal_par(:, 2);
Ci12 = modal_par(:, 3);

[H13_rfp, modal_par] = rfp(H13_mode, mode, m);
freq13 = modal_par(:, 1);
damp13 = modal_par(:, 2);
Ci13 = modal_par(:, 3);

[H14_rfp, modal_par] = rfp(H14_mode, mode, m);
freq14 = modal_par(:, 1);
damp14 = modal_par(:, 2);
Ci14 = modal_par(:, 3);

[H15_rfp, modal_par] = rfp(H15_mode, mode, m);
freq15 = modal_par(:, 1);
damp15 = modal_par(:, 2);
Ci15 = modal_par(:, 3);

[H16_rfp, modal_par] = rfp(H16_mode, mode, m);
freq16 = modal_par(:, 1);
damp16 = modal_par(:, 2);
Ci16 = modal_par(:, 3);

[H17_rfp, modal_par] = rfp(H17_mode, mode, m);
freq17 = modal_par(:, 1);
damp17 = modal_par(:, 2);
Ci17 = modal_par(:, 3);

[H18_rfp, modal_par] = rfp(H18_mode, mode, m);
freq18 = modal_par(:, 1);
damp18 = modal_par(:, 2);
Ci18 = modal_par(:, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 5:0.5:125;
start_ind = find(modes==5); end_ind = find(modes==125);
H11_mode = H11(start_ind:end_ind);
H17_mode = H17(start_ind:end_ind);

[H11_rfp_2, modal_par] = rfp(H11_mode, mode, 2);
[H11_rfp_3, modal_par] = rfp(H11_mode, mode, 3);
[H11_rfp_4, modal_par] = rfp(H11_mode, mode, 4);
[H11_rfp_5, modal_par] = rfp(H11_mode, mode, 5);
[H11_rfp_6, modal_par] = rfp(H11_mode, mode, 6);

[H17_rfp_2, modal_par] = rfp(H17_mode, mode, 2);
[H17_rfp_3, modal_par] = rfp(H17_mode, mode, 3);
[H17_rfp_4, modal_par] = rfp(H17_mode, mode, 4);
[H17_rfp_5, modal_par] = rfp(H17_mode, mode, 5);
[H17_rfp_6, modal_par] = rfp(H17_mode, mode, 6);


% hold on
% set(gcf,'units','points','position',[0,0,1250,400])
% plot(mode, abs(H17_mode), 'DisplayName', 'Analytical')
% plot(mode, abs(H17_rfp_2), 'DisplayName', 'RFP m=2')
% plot(mode, abs(H17_rfp_3), 'DisplayName', 'RFP m=3')
% plot(mode, abs(H17_rfp_4), 'DisplayName', 'RFP m=4')
% plot(mode, abs(H17_rfp_5), 'DisplayName', 'RFP m=5')
% plot(mode, abs(H17_rfp_6), 'DisplayName', 'RFP m=6')
% 
% xlabel('f')
% ylabel('|H_{17}(2 \pi f)|')
% legend()
% title('Plot of Reactance FRF H_{17}')
% hold off





