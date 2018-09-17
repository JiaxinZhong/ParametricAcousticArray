%% import the data
dir_head = 'kzk/data/';

set_num = 1;
switch set_num
case 1
	% set 1 for P0 = 118 dB
	data = load(strcat(dir_head, 'kzk_cal_gen_fb1k_f160k_a0p1_P0118.mat'));
case 2
	% set 2 for P0 = 135 dB
	data = load(strcat(dir_head, 'kzk_cal_gen_fb4k_f160k_a0p1_P0135.mat'));
end


%% plot the SPL along the axis of the transducer
% set the marker
num_crv = 5;
num_mkr = 8;
stt_mkr = zeros(num_crv,1);
for i = 1:num_crv
	stp_mkr(i) = fix(length(data.z)/num_mkr);
end
stt_mkr = fix((1:num_crv)'/num_crv.*stp_mkr);
i_mkr = cell(num_crv,1);
for i = 1:num_crv
	i_mkr{i} = stt_mkr(i):stp_mkr(i):length(data.z);
end

figure
plot(data.z(i_mkr{1}(1)), data.Lp{1}(1,1), '-o', 'color', ca1);
hold on
plot(data.z(i_mkr{2}(1)), data.Lp{2}(1,1), '-*', 'color', ca2);
plot(data.z(i_mkr{3}(1)), data.Lp{3}(1,1), '-^', 'color', ca3);
plot(data.z(i_mkr{4}(1)), data.Lp{data.N1}(1,1), '-d', 'color', ca4);
plot(data.z(i_mkr{5}(1)), data.Lp{data.N2}(1,1), '-o', ...
	'MarkerFaceColor', ca5, 'color', ca5);

plot(data.z(i_mkr{1}), data.Lp{1}(1,i_mkr{1}), 'o', 'color', ca1);
plot(data.z(i_mkr{2}), data.Lp{2}(1,i_mkr{2}), '*', 'color', ca2);
plot(data.z(i_mkr{3}), data.Lp{3}(1,i_mkr{3}), '^', 'color', ca3);
plot(data.z(i_mkr{4}), data.Lp{data.N1}(1,i_mkr{4}), 'd', 'color', ca4);
plot(data.z(i_mkr{5}), data.Lp{data.N2}(1,i_mkr{5}), 'o', ...
	'MarkerFaceColor', ca5, 'color', ca5);

plot(data.z, data.Lp{1}(1,:), '-', 'color', ca1);
plot(data.z, data.Lp{2}(1,:), '-', 'color', ca2);
plot(data.z, data.Lp{3}(1,:), '-', 'color', ca3);
plot(data.z, data.Lp{data.N1}(1,:), '-', 'color', ca4);
plot(data.z, data.Lp{data.N2}(1,:), '-', 'color', ca5);

xlim([0,1]);
switch set_num
case 1
	ylim([-100,120]);
case 2
	ylim([0,140]);
end
xlabel('$z$ (m)');
ylabel('SPL (dB)');

legend({'1 kHz', '2 kHz', '3 kHz', '60 kHz', '61 kHz'});

 %% save the figure
 print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', ...
 	'-r200');
