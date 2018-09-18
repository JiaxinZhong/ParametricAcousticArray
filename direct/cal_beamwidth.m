function [theta] = cal_beamwidth(f_U, f_A, hr, Tc, pa)
% =========================================================================
% FUNCTION	Calculate the beamwidth at the far field
% REFERENCE
%	[1] BERKTAY H O. Possible exploitation of non-linear acoustics in 
%		underwater transmitting applications[J]. Journal of Sound and 
%		Vibration, 1965, 2(4): 435–461.
% -------------------------------------------------------------------------
% INPUT
%	f_U			- the frequency of the ultrasound
%	f_A			- the frequency of the audio sound
%	hr			- relative humidity as a percentage
% 	pa			- the atmospheric pressure, in kilopascals
% -------------------------------------------------------------------------
% OUTPUT
%	theta		- the half power beamwidth of the parametric acoustic array
% -------------------------------------------------------------------------
% VERSION INFO
% Author		- Jiaxin Zhong
% Last modified	- 2018-08-08
% Version		- 1.0.20180808
% =========================================================================

	% the wavenumber of audible sound
	k_A = 2*pi*f_A/343;
	% the absorption coefficient in Np
	alpha_U = cal_absorp_coeff(f_U, hr, Tc, pa);
	% the composite (total) attenuation coeff.
	alpha_T = 2*alpha_U;
	% theta_paa = 4 * atan(sqrt(alpha_T ./ 2 / k_A));
	theta = 4 * (sqrt(alpha_T ./ 2 / k_A));
end
