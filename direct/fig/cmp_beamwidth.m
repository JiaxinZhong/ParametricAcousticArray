f_A = [500;1000;2000;5000];
f_U = linspace(30000,100000,200);

Tc = 20;
pa = 101.325;
hr = 50;
theta_paa = zeros(length(f_U),length(f_A));
for i = 1:length(f_A)
    theta_paa(:,i) = cal_beamwidth(f_U, f_A(i), hr, Tc, pa);
end
theta_paa = theta_paa*180/pi; % in degree


figure
num_crv = length(f_A);
stp_mkr = ceil(length(f_U)/7);
stt_mkr = ceil((1:num_crv)/num_crv*stp_mkr);
idx_mkr = cell(num_crv,1);
for i = 1:num_crv
    idx_mkr{i} = stt_mkr(i):stp_mkr:length(f_U);
end

type_mkr0 = {'o','*','^','d'};
type_mkr = strcat('-',type_mkr0);
facecolor_mkr = {'none', 'none', 'none', 'none'};

f_U = f_U/1000;
for i = 1:num_crv
	plot(f_U(stt_mkr(i)), theta_paa(stt_mkr(i),i), type_mkr{i}, ...
		'color', cacell{i}, 'MarkerFaceColor', facecolor_mkr{i});
	hold on
end
for i = 1:num_crv
	plot(f_U, theta_paa(:,i), '-', 'color', cacell{i});
	plot(f_U(idx_mkr{i}), theta_paa(idx_mkr{i},i), type_mkr0{i}, ...
		'color', cacell{i}, 'MarkerFaceColor', facecolor_mkr{i});
end

f_A_str = sprintfc('$$f_\\mathrm{A} = $$ %g Hz', f_A);
legend(f_A_str);
ylabel('Half power beamwidth ($^\circ$)')
xlabel('$f_\mathrm{U}$ (kHz)')

print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', '-r200');
