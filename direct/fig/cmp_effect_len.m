hr = [30; 50; 80; 100];
f_U = linspace(30000,100000,200);
Tc = 20;
pa = 101.325;

l = zeros(length(f_U),length(hr));
for i = 1:length(hr)
	l(:,i) = cal_effect_len(f_U, hr(i), Tc, pa);
end

figure
num_crv = length(hr);
stp_mkr = ceil(length(f_U)/7);
stt_mkr = ceil((1:num_crv)/num_crv*stp_mkr);
idx_mkr = cell(num_crv,1);
for i = 1:num_crv
    idx_mkr{i} = stt_mkr(i):stp_mkr:length(f_U);
end

type_mkr0 = {'o','*','^','d', 'o'};
type_mkr = strcat('-',type_mkr0);
facecolor_mkr = {'none', 'none', 'none', 'none', cacell{5}};

f_U = f_U/1000;
for i = 1:num_crv
	plot(f_U(stt_mkr(i)), l(stt_mkr(i),i), type_mkr{i}, ...
		'color', cacell{i}, 'MarkerFaceColor', facecolor_mkr{i});
	hold on
end
for i = 1:num_crv
	plot(f_U, l(:,i), '-', 'color', cacell{i});
	plot(f_U(idx_mkr{i}), l(idx_mkr{i},i), type_mkr0{i}, ...
		'color', cacell{i}, 'MarkerFaceColor', facecolor_mkr{i});
end

hr_str = sprintfc('$$h_\\mathrm{r} = $$ %g', hr);
legend(hr_str);
ylabel('Effective absorption length (m)')
xlabel('$f_\mathrm{U}$ (kHz)')

print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', '-r200');

