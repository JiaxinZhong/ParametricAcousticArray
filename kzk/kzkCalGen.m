%% Physical parameters
Tc = 20;                    % 
hr = 50;                    % relative humidity as a percentage
c0 = 331.6+0.6*Tc;  
rho0 = 1.2;
a_rad = 0.1;
beta = 1.2;
LP0 = 128;
PREF = 20e-6;
c1 = 1;
c2 = 1;
f1 = 60000; 
fd = 1000;                  % difference frequency
f2 = f1+fd;
fb = gcd(f1,f2);
M = fix(f2/fb*2+5);
P0 = 10^(LP0/20) * PREF;
v0 = P0/rho0/c0;
eps0 = v0/c0;
f0 = (f1+f2)/2;

Nd = fd/fb;
N0 = f0/fb;
N1 = f1/fb;
N2 = f2/fb;
k0 = 2*pi*f0/c0;
RD0 = k0*a_rad^2/2;
lD0 = 1/(beta*k0*eps0);

%% Numerical parameters
Dsig = 3e-4;
Dzeta = 2e-2;

zeta_max = 15;
z_max = 1;                  % the maximal z [m]
sig_max = z_max/RD0;
% sig_max = 5;
I = fix(sig_max/Dsig);
J = fix(zeta_max/Dzeta);
sig = linspace(0, sig_max, I+1)';
zeta = linspace(0,zeta_max, J+1)';
xi = zeta * (sig.' +1);
z = sig * RD0;
n = (1:M)';
f = (1:M).'*fb;
alpha = cal_absorp_coeff(f, hr, Tc, 101.325);
[Lp,p] = kzk_cal(Dsig,sig_max,Dzeta,zeta_max,N0,N1,N2,M,c1,c2,alpha,RD0,...
	lD0,LP0);


for i = 1:M
%     reserve the ultrasound and their sum frequency
    if (N1~=i && N2~=i && (N1*2)~=i && (N2*2)~=i && (N1+N2)~=i ...
            && i~=Nd && i~=2*Nd && i~=3*Nd)
        Lp{i} = {};
        p{i} = {};
    end
end

% Plot
figure
plot(sig*RD0, Lp{Nd}(1,:),'-');
xlabel('$z$ (m)')
ylabel('SPL (dB)')

print(sprintf('%s_cache.jpg', mfilename('fullpath')), '-djpeg', '-r300');
save('kzk/data/kzk_cal_gen_cache.mat','Lp','p','RD0','sig','f0','N1','N2',...
	'N0','f1','f2','xi', 'Tc', 'hr', 'rho0', 'a_rad','fb','fd','Nd',...
	'beta', 'LP0','f','zeta','sig_max','z_max',...
	'alpha','I','J','Dsig','Dzeta','c0','M','PREF','z');
