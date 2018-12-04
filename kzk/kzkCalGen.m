%% Physical parameters
tempCelcius = 20;                 
relHumidityPer = 50;             % relative humidity as a percentage
c0 = 331.6 + 0.6*tempCelcius;  
rho0 = 1.2;
radPrj = 0.1;   % radius of the projector
beta = 1.2;
LP0 = 128;
PRS_REF = 20e-6;
c1 = 1;
c2 = 1;
freq1 = 38000; 
freqDiff = 2000;           % difference frequency
freq2 = freq1+freqDiff;
freqBase = gcd(freq1,freq2);
M = fix(freq2/freqBase*2+10);
prsUltr = 10^(LP0/20) * PRS_REF;
velUltr = prsUltr/rho0/c0;
eps0 = velUltr/c0;
freqMid = (freq1+freq2)/2;

NDiff = freqDiff/freqBase;
NMid = freqMid/freqBase;
NU1 = freq1/freqBase;
NU2 = freq2/freqBase;
wavnumMid = 2*pi*freqMid/c0;
rayDistMid = wavnumMid*radPrj^2/2;
lD0 = 1/(beta*wavnumMid*eps0);

freqLst = (1:M).'*freqBase;
sacLst = cal_absorp_coeff(freqLst, relHumidityPer, tempCelcius, 101.325);
absorpLenUltr = 1/2/sacLst(NU1);

%% Numerical parameters
sigDel = 5e-4;
zetaDel = 2e-2;

zetaMax = 20;
% zMax = absorpLenUltr*5;               % the maximal z [m]
zMax = 0.05;
sigMax = zMax/rayDistMid;
% sig_max = 5;
I = fix(sigMax/sigDel);
J = fix(zetaMax/zetaDel);
sigLst = linspace(0, sigMax, I+1)';
zetaLst = linspace(0,zetaMax, J+1)';
xiLst = zetaLst * (sigLst.' +1);
zLst = sigLst * rayDistMid;
n = (1:M)';

[Lp,prsMag,prsPha] = kzkCal(sigDel,sigMax,zetaDel,zetaMax,NU1,NU2,M,c1,c2,sacLst,rayDistMid,...
	lD0,LP0, wavnumMid);

% Plot
figure
plot(sigLst*rayDistMid, Lp{NDiff}(1,:),'-');
xlabel('$z$ (m)')
ylabel('SPL (dB)')

figure
thePrsLvl = Lp{NDiff};
thePrsLvl(:,1) = 0;
zArr = repmat(zLst.', size(xiLst,1), 1);
pcolor(zArr, xiLst, thePrsLvl);
colorbar
MAX_LVL = max(max(thePrsLvl));
caxis([MAX_LVL-30, MAX_LVL]);
shading interp
% ylim([0,15])

print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', '-r300');
save('kzk/data/kzkCalGen_cache.mat','Lp','prsMag','prsPha',...
    'rayDistMid','sigLst',...
    'freq1','freq2','freqDiff','freqMid','freqLst',...
    'NU1','NU2','NMid','NDiff',...
	'xiLst', 'tempCelcius', 'relHumidityPer', 'rho0', 'radPrj',...
    'freqBase',...
	'beta', 'LP0','zetaLst','sigMax','zMax',...
	'sacLst','I','J','sigDel','zetaDel','c0','M','PRS_REF','zLst',...
    'absorpLenUltr');
