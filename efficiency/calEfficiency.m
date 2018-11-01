E = 130e9; 
rho = 2300;
nu = 0.278;
fm = 60e3;
Lambda = 10.22;
Qm = 500;
t_a2 = 2*pi*fm/Lambda*(E/12/rho/(1-nu^2))^(-1/2);
N = 0.00474;

eps0 = 1e-9/36/pi;
tand = 0.02;
eps33 = 1400*eps0;
Cf = eps33*pi*5/9/t_a2;
R0 = 1/(2*pi*fm*Cf*tand);

a = (0.001:0.001:0.01);
etaEm = zeros(length(a),1);
etaMa = etaEm;
etaEa = etaEm;
t = etaEm;
electriPower = t;
for iA = 1:length(a)
    aNow = a(iA);
    t(iA) = t_a2*aNow^2;
    Meff = 192/Lambda^2*rho*pi*aNow^2*t(iA);
    Keff = 16*pi^2/aNow^2*E*t(iA)^3/(1-nu^2);
    
    Rm = 1/Qm*32*pi/Lambda*(3*rho*E/(1-nu^2))^0.5*t(iA)^2;
    Rr = 415*pi*aNow^2;
    ka = 2*pi*fm/343*aNow;
    Rr = Rr * (1-2*besselj(1,2*ka)/2/ka);
    R = (Rm+Rr)/N^2;
    
    etaEm(iA) = 1/R / (1/R0+1/R);
    etaMa(iA) = 1/(Rm/Rr+1);
    etaEa = etaEm.*etaMa;
    electriPower(iA) = 10^((130-10*log10(etaEa(iA))+20*log10(0.1)-118.2)/10);
end

etaEaPerHundred = round(etaEa * 100,2); % rounds to 2 digits
etaEmPerHundred = round(etaEm * 100,2);
etaMaPerHundred = round(etaMa * 100,2);

electriPower = round(electriPower,2);