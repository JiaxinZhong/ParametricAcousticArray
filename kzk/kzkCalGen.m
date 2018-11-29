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
freq1 = 60000; 
freqDiff = 1000;           % difference frequency
freq2 = freq1+freqDiff;
freqBase = gcd(freq1,freq2);
M = fix(freq2/freqBase*2+10);
prsUltr = 10^(LP0/20) * PRS_REF;
velUltr = prsUltr/rho0/c0;
eps0 = velUltr/c0;
freqMid = (freq1+freq2)/2;

Nd = freqDiff/freqBase;
N0 = freqMid/freqBase;
N1 = freq1/freqBase;
N2 = freq2/freqBase;
wavnumMid = 2*pi*freqMid/c0;
rayDistMid = wavnumMid*radPrj^2/2;
lD0 = 1/(beta*wavnumMid*eps0);

freqLst = (1:M).'*freqBase;
sacLst = cal_absorp_coeff(freqLst, relHumidityPer, tempCelcius, 101.325);
absorpLenUltr = 1/2/sacLst(N1);

%% Numerical parameters
sigDel = 5e-4;
zetaDel = 2e-2;

zetaMax = 15;
zMax = absorpLenUltr*5;               % the maximal z [m]
sigMax = zMax/rayDistMid;
% sig_max = 5;
I = fix(sigMax/sigDel);
J = fix(zetaMax/zetaDel);
sigLst = linspace(0, sigMax, I+1)';
zetaLst = linspace(0,zetaMax, J+1)';
xiLst = zetaLst * (sigLst.' +1);
zLst = sigLst * rayDistMid;
n = (1:M)';

[Lp,prsMag,prsPha] = kzkCal(sigDel,sigMax,zetaDel,zetaMax,N0,N1,N2,M,c1,c2,sacLst,rayDistMid,...
	lD0,LP0, wavnumMid);


% reserveOrderLst = 1:4;
for i = 1:M
%     reserve the ultrasound and their sum frequency
    if (N1~=i && N2~=i && (N1*2)~=i && (N2*2)~=i && (N1+N2)~=i ...
            && i~=Nd && i~=2*Nd && i~=3*Nd && i~=4*Nd && i~=5*Nd)
        Lp{i} = {};
        prsMag{i} = {};
        prsPha{i} = {};
    end
end

% Plot
figure
plot(sigLst*rayDistMid, Lp{Nd}(1,:),'-');
xlabel('$z$ (m)')
ylabel('SPL (dB)')

print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', '-r300');
save('kzk/data/kzkCalGen_cache.mat','Lp','prsMag','prsPha',...
    'rayDistMid','sigLst',...
    'freq1','freq2','freqDiff','freqMid','freqLst',...
    'N1','N2',...
	'N0','xiLst', 'tempCelcius', 'relHumidityPer', 'rho0', 'radPrj',...
    'freqBase','Nd',...
	'beta', 'LP0','zetaLst','sigMax','zMax',...
	'sacLst','I','J','sigDel','zetaDel','c0','M','PRS_REF','zLst',...
    'absorpLenUltr');
