%% Import the data
fdList = [(400:200:5000)';(500:1000:4500)'];
fdList = sort(fdList);
zList = [0.25;0.5;0.75;1];
splMat = zeros(length(fdList), length(zList));

for i = 1:length(fdList)
	iFn = sprintf('kzk/data/kzk_cal_gen_fd%s_f160k_a0p1_P0128.mat', ...
		num2str(fdList(i)));
	iData = load(iFn);
    iSpl = iData.Lp{iData.Nd}(1,:);
    iZ = iData.z;
    % Interpolate
    splMat(i,:) = interp1(iZ(2:end), iSpl(2:end).', zList, 'spline');
end

%% plot the figure
figure
for i = 1:length(zList)
    plot(fdList, splMat(:,i));
    hold on
end
xlim([400,5000])
xlabel('频率 (Hz)')
ylabel('SPL (dB)')
legend(strcat('$z = ',sprintfc('%g',zList), '$ m'));
figAddMarker(7)


%% Export the figure
print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', '-r300');

