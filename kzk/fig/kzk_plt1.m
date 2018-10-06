%% import the data
data = cell(4,1);
dir_head = 'kzk/data/';

set_num = 1;
switch set_num
case 1
	% set 1 for P0 = 118 dB
	data{1} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd500_f160k_a0p1_P0118.mat'));
	data{2} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd1k_f160k_a0p1_P0118.mat'));
	data{3} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd2k_f160k_a0p1_P0118.mat'));
	data{4} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd4k_f160k_a0p1_P0118.mat'));
case 2
	% set 2 for P0 = 135 dB
	data{1} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd500_f160k_a0p1_P0135.mat'));
	data{2} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd1k_f160k_a0p1_P0135.mat'));
	data{3} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd2k_f160k_a0p1_P0135.mat'));
	data{4} = load(strcat(dir_head, ...
		'kzk_cal_gen_fd4k_f160k_a0p1_P0135.mat'));
end

%% plot the SPL along the axis of the transducer
% set the marker
num_crv = 4;
num_mkr = 8;
stt_mkr = zeros(num_crv,1);
for i = 1:num_crv
	stp_mkr(i) = fix(length(data{i}.z)/num_mkr);
end
stt_mkr = fix((1:num_crv)'/num_crv.*stp_mkr);
i_mkr = cell(num_crv,1);
for i = 1:num_crv
	i_mkr{i} = stt_mkr(i):stp_mkr(i):length(data{i}.z);
end

figure
plot(data{1}.z(i_mkr{1}(1)), data{1}.Lp{1}(1,1), '-o', 'color', ca1);
hold on
plot(data{2}.z(i_mkr{2}(1)), data{2}.Lp{1}(1,1), '-*', 'color', ca2);
plot(data{3}.z(i_mkr{3}(1)), data{3}.Lp{1}(1,1), '-^', 'color', ca3);
plot(data{4}.z(i_mkr{4}(1)), data{4}.Lp{1}(1,1), '-d', 'color', ca4);

plot(data{1}.z(i_mkr{1}), data{1}.Lp{1}(1,i_mkr{1}), 'o', 'color', ca1);
plot(data{2}.z(i_mkr{2}), data{2}.Lp{1}(1,i_mkr{2}), '*', 'color', ca2);
plot(data{3}.z(i_mkr{3}), data{3}.Lp{1}(1,i_mkr{3}), '^', 'color', ca3);
plot(data{4}.z(i_mkr{4}), data{4}.Lp{1}(1,i_mkr{4}), 'd', 'color', ca4);

plot(data{1}.z, data{1}.Lp{1}(1,:), '-', 'color', ca1);
plot(data{2}.z, data{2}.Lp{1}(1,:), '-', 'color', ca2);
plot(data{3}.z, data{3}.Lp{1}(1,:), '-', 'color', ca3);
plot(data{4}.z, data{4}.Lp{1}(1,:), '-', 'color', ca4);

xlim([0,1]);
switch set_num
case 1
	ylim([10,70]);
case 2
	ylim([40,100]);
end
xlabel('$z$ (m)');
ylabel('SPL (dB)');

legend({'$f_\mathrm{b} = 500$ Hz', '$f_\mathrm{b} = 1$ kHz'...
 '$f_\mathrm{b} = 2$ kHz', '$f_\mathrm{b} = 4$ kHz'});

 %% save the figure
 print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', ...
 	'-r200');
