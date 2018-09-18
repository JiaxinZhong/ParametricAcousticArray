f_A = 4000;
f_U = 60000;
a = 0.1; % the radius of the transducer
Tc = 20;
pa = 101.325;
hr = 50;
theta = (-180:0.001:190)*pi/180; 
[L_D, ~, L_D_conv, ~] = cal_direct(f_U, f_A, a, theta, hr, Tc, pa);
MIN_L = -30;
L_D(L_D<MIN_L) = MIN_L;
L_D_conv(L_D_conv<MIN_L) = MIN_L;

figure
polarplot(theta, L_D);
rlim([MIN_L,0])
hold on
polarplot(theta, L_D_conv, '--')
print(sprintf('%s_chache.jpg', mfilename('fullpath')), '-djpeg', '-r200');
