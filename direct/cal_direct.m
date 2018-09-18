function [L_D, D, L_D_ap, D_ap] = cal_direct_paa(f_U, f_A, a, theta, hr,...
    Tc, pa)
% =========================================================================
% FUNCTION	Calculate the directivity of the PAA
% REFERENCE
%	[1] BERKTAY H O. Possible exploitation of non-linear acoustics in 
%		underwater transmitting applications[J]. Journal of Sound and 
%		Vibration, 1965, 2(4): 435–461.
% -------------------------------------------------------------------------
% INPUT
%	f_U			- the frequency of the ultrasound
%	f_A			- the frequency of the audio sound
%	a			- the radius of the transducer
%	theta		- the angle in radians
%				- can be a column vector
%	hr			- relative humidity as a percentage
% 	pa			- the atmospheric pressure, in kilopascals
% -------------------------------------------------------------------------
% OUTPUT
%	D			- the directivity function
%	L_D			- the directivity function in dB (20*log(D))
%	D_ap		- the beam pattern resulting from the finite cross section 
%                   of the primary beam
%	L_D_ap		- D_ap in dB (20*log(D_ap))
% -------------------------------------------------------------------------
% VERSION INFO
% Author			- Jiaxin Zhong
% Last modified		- 2018-08-08
% Version			- 1.0.20180808
% =========================================================================

	k_A = 2*pi*f_A/343;
	k_U = 2*pi*f_U/343;

	alpha_U = cal_absorp_coeff(f_U, hr, Tc, pa);
	% the composite (total) attenuation coeff.
	alpha_T = 2*alpha_U;

	ka_A = k_A * a;
	sint = (sin(theta));

	D_ap = 2 * abs(besselj(1,ka_A * sint) ./ (ka_A * sint));
	L_D_ap = 20*log10(D_ap);
	D_A = (1 + (2*k_A/alpha_T * (sin(theta/2)).^2).^2).^(-1/2);
	D = D_ap .* D_A;
	L_D = 20*log10(D);
end
