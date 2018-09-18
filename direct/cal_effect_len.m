function [l] = cal_effect_len(f, hr, Tc, pa)
% =========================================================================
% FUNCTION 	calculate the effective length of the parametrci acoustic array
% REFERENCE
%	[1] BERKTAY H O. Possible exploitation of non-linear acoustics in 
%		underwater transmitting applications[J]. Journal of Sound and 
%		Vibration, 1965, 2(4): 435–461.
% -------------------------------------------------------------------------
% Input:
%	f		- Frequency, in Hertz. 
%			- Column vector or scalar
%	hr 		- relative humidity as a percentage. 
%			- scalar
%	Tc		- Ambient atmospheric temperature, in Celcius
%			- scalar
%	pa 		- the atmospheric pressure, in kilopscals
%			- scalar
% -------------------------------------------------------------------------
% Output:
%	l			- the effective absoption length of the PAA
%				- the same as f in dimension
% -------------------------------------------------------------------------
% VERSION INFO
%	Author				- Jiaxin Zhong
%	VERSION				- 1.0.20180918
% =========================================================================

	alpha_Np = cal_absorp_coeff(f, hr, Tc, pa);
	l = 1./(2*alpha_Np);

end
