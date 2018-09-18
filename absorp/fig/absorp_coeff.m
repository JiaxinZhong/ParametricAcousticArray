%% Initialization
Tc = 20; % Ambient atmospheric temperature, in Celcius
pa = 101.325; % the atmospheric pressure, in kilopscals
hr = [0;10;20;40;60;80;100]; % relative humidity as a percentage
f = logspace(log10(29900),log10(300000),100); % frequency

% pure-tone sound attenuation coeff. in dB per meter, 
%	for atmospheric absorption
alpha = zeros(length(hr),length(f)); 
for i = 1:length(hr)
    hr_now = hr(i);
    [~, alpha(i,:)] = cal_absorp_coeff(f, hr_now, Tc, pa);
end

%% plot the figure
figure;
num_crv = length(hr);
stp_mkr = ceil(length(f)/3);
stt_mkr = ceil((1:num_crv)/num_crv*stp_mkr);
idx_mkr = cell(num_crv,1);
for i = 1:num_crv
    idx_mkr{i} = stt_mkr(i):stp_mkr:length(f);
end

type_mkr0 = {'o','*','^','d','o','^','d'};
type_mkr = strcat('-',type_mkr0);
facecolor_mkr = {'none', 'none', 'none', 'none', cacell{5:num_crv}};
f = f/1000; 
for i = 1:num_crv
    loglog(f(stt_mkr(i)), alpha(i,stt_mkr(i)), type_mkr{i}, ...
        'color', cacell{i}, 'MarkerFaceColor', facecolor_mkr{i});
    hold on
end
for i = 1:num_crv
    loglog(f, alpha(i,:), '-', 'color', cacell{i});
    loglog(f(idx_mkr{i}), alpha(i,idx_mkr{i}), type_mkr0{i}, ...
        'color', cacell{i}, 'MarkerFaceColor', facecolor_mkr{i});
end
xlim([min(f) max(f)]);
xlabel('Frequency (kHz)')
ylabel('Absorption coefficient $\alpha_\mathrm{dB}$ (dB/m)')
ftick = [30;50;80;100;150;200;250;300];
set(gca, 'xtick', ftick);
hr_str = sprintfc('$$h_\\mathrm{r}$$ = %g', hr);
legend(hr_str)
set(gca,'gridlinestyle','-');
set(gca,'MinorGridLineStyle','-')

print(sprintf('%s_cache.jpg',mfilename('fullpath')), '-djpeg', '-r200');
