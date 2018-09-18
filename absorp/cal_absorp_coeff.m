function [alpha_Np, alpha_dB] = cal_absorp_coeff(f, hr, Tc, pa)
% =========================================================================
% FUNCTION  Calculate the absorption coefficients by the atmosphere based 
%			on the standard ISO 9613-1.
% REFERENCE
%	[1] ISO Technical Committees: Noise. Acoustics — Attenuation of sound 
%		during propagation outdoors — Part 1: Calculation of the absorption
%		of sound by the atmosphere: ISO 9613-1:1993[S]. Geneva: 
%		International Organization for Standardization, 1993.
%	[2]	National Physical Laboratory. NPL Acoustics: Calculation of 
%		absorption of sound by the atmosphere[EB/OL]. [2018-08-08].
%		http://resource.npl.co.uk/acoustics/techguides/absorption/.
% -------------------------------------------------------------------------
% INPUT
%	f		- Frequency, in Hertz. 
%			- Column vector or scalar
%	hr 		- relative humidity as a percentage. 
%	Tc		- Ambient atmospheric temperature, in Celcius
%	pa 		- the atmospheric pressure, in kilopscals
% -------------------------------------------------------------------------
% OUTPUT
%	alpha_dB	- pure-tone sound attenuation coeff. in dB per meter, for 
%					atmospheric absorption
%				- the same as f in dimension
%	alpha_Np	- pure-tone sound attenutaion oeffi. in Np per meter, for 
%					atmospheric absorption
%				- the same as f in dimension
% -------------------------------------------------------------------------
% VERSION INFO
%	Author					- Jiaxin Zhong
%	Last modified on		- August 7, 2018
%	Version number			- 1.0.20180807
% =========================================================================

	% reference air temperature in Kelvins, i.e. 20 Celcius degree
	T0 = 293.15; 
	% the triple-point isotherm temperature (i.e. +0.01 Celsius degree)
	T01 = 273.16; 
	% Ambient atmospheric temperature, in Kelvins
	T = Tc + 273.15; 

	pr = 101.325; % reference ambient atmospheric pressure, in kilopascals

    C = -6.8346*(T01/T).^1.261 + 4.6151;
    psat = pr * 10.^C; % saturation vapour pressure
	% the molar concentration of water vapoupr as a percentage
    h = hr.*(psat/pr).*(pa/pr); 
    
    % the oxygen relaxation frequency
    f_rO = pa./pr .* (24 + 4.04*10^4 * h .* (0.02+h) ./ (0.391+h) );
    % the nitrogen relaxation frequency
    f_rN = pa./pr .* (T/T0).^(-1/2) .* (9 + 280 * h ...
        .* exp(-4.17.*((T/T0).^(-1/3)-1)));
    
    alpha_Np = f.^2 .* (1.84*10^(-11) .* pr/pa .* (T/T0)^(1/2) ...
        + (T/T0)^(-5/2) * (0.01275*exp(-2239.1/T)./(f_rO+f.^2/f_rO) ...
        + 0.1068*exp(-3352.0/T)./(f_rN+f.^2/f_rN)) );
	alpha_dB = 20/log(10) * alpha_Np;
end
