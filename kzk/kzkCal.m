function [prsLvl,prsMag,prsPha] = kzkCal(sigDel,sigMax,zetaDel,zetaMax,...
    N1,N2,M,c1,c2,sac,rayDistMid,lD0,LP0,wavnumMid)
% =========================================================================
% FUNCTION 	calculate the nonlienar pressure field based on the kzk 
%			equation using the Implicit Backward Finite Difference
% REFERENCE 
%	[1] AANONSEN S I, BARKVE T, NAZE J, et, al. Distortion and harmonic 
%		generation in the nearfield of a finite amplitude sound beam[J]. 
%		The Journal of the Acoustical Society of America, 1989, 75(3): 
%		749–768.
%	[2] TJOTTA N J, Sigva T, VEFRING E H. Propagation and interaction of 
%		two collinear finite amplitude sound beams[J]. The Journal of the 
%		Acoustical Society of America, 1990, 88(6): 2859–2870.
%	[3]	SMITH G D. Numerical solution of partial differential equations – 
%		Finite Difference Methods[M]. 3rd ed. Oxford: Oxford University 
%		Press, 1985.
%	[4] KAMAKURA T, HAMADA N, AOKI K, et, al. Nonlinearly generated 
%		spectral components in the nearfield of a directive sound 
%		source[J]. The Journal of the Acoustical Society of America, 1989,
%		85(6): 2331–2337.
% -------------------------------------------------------------------------
% INPUT
%	Dsig	- step size in the direction of sigma (axis of the transducer)
%	sig_max	- maximal size in the direction of sigma (axis of the trans.)
%	Dzeta	- step size in the direction of zeta (transverse plane)
%	zeta_max- maximal size in the direction of zeta
%	N0		- = f0/fb, where fb is the basic frequency
%	N1		- = f1/fb
%	N2 		- = f2/fb
%	M		- the maxial number of truncated harmonics
%	c1, c2  - distribution functions along in the transverse plane
%	alpha	- the absorption coefficients in Np/m
%			- a column vector
%	RD0		- Rayleigh distance
%	lD0		- shock wave formation distance
%	LP0		- the dB value of the peak amplitude of the source
% -------------------------------------------------------------------------
% OUTPUT
%	Lp		- sound power pressure of the harmonics in dB
%   prsMag  - the magnitude of the pressure
%   prsPha  - the phase of the pressure
% -------------------------------------------------------------------------
% VERSION INFO
%	Author				- Jiaxin Zhong
%	Last modified 		- Nov 9, 2018
%	Version 			- 2.1
% =========================================================================

NMid = (N1+N2)/2;
NDiff = N2-N1;

I = fix(sigMax/sigDel);
n = (1:M).'; 
sigLst = linspace(0, sigMax, I+1)';
J = fix(zetaMax/zetaDel);
zetaLst = linspace(0,zetaMax, J+1)';
j_a = fix(1/zetaDel);

am = (1-0.5./(1:J-1).') * (1./ n.') * sigDel/4/(zetaDel^2)*NMid; % (J-1) x M
ap = (1+0.5./(1:J-1).') * (1./ n.') * sigDel/4/(zetaDel^2)*NMid; % (J-1) x M
b = -(1./n) *sigDel/2/(zetaDel^2)*NMid; % M x 1
A = cell(M,1);
for nn = 1:M
	bb = b(nn);
	amm = am(1:end,nn);
	app = ap(1:end-1,nn);
	A{nn} = sparse(diag([2*bb;bb*ones(J-1,1)]) + diag([-2*bb;app],1) + ...
        diag(amm,-1));
end

G = zeros(J+1, I+1, M);
H = G;

% boundary conditions
G(0+1:j_a+1,0+1,N1) = c1*cos(N1/NMid*zetaLst(0+1:j_a+1).^2);
G(0+1:j_a+1,0+1,N2) = c2*cos(N2/NMid*zetaLst(0+1:j_a+1).^2);
H(0+1:j_a+1,0+1,N1) = c1*sin(N1/NMid*zetaLst(0+1:j_a+1).^2);
H(0+1:j_a+1,0+1,N2) = c2*sin(N2/NMid*zetaLst(0+1:j_a+1).^2);
B = cell(M,1);
for nn = 1:M
	B{nn} = sparse(-sigDel * 1^2 * sac(nn) * rayDistMid* eye(J));
end

U = cell(M,1);
for nn = 1:M
	U{nn} = zeros(J+1,I+1); 
end
V = U;
for iSig = 1:I
	tic % start the stopwatch
	fprintf('i = %s, I = %s\n', num2str(iSig), num2str(I));
	for nn = 1:M
		for jj = 0:J-1
			% use matrix multiplication to accelerate computations
			GG1 = G(jj+1,iSig-1+1,1:nn-1);
			GG1 = GG1(:);
			GG2 = flipud(GG1);
			HH1 = H(jj+1,iSig-1+1,1:nn-1);
			HH1 = HH1(:);
			HH2 = flipud(HH1);
			Un = 1/2 * (GG1.' * GG2 - HH1.' * HH2);
			Vn = HH1.' * GG2;
				
			GG1 = (G(jj+1,iSig-1+1,1:M-nn)); % m-n
			GG1 = GG1(:);
			GG2 = (G(jj+1,iSig-1+1,nn+1:M)); % m
			GG2 = GG2(:);
			HH1 = (H(jj+1,iSig-1+1,1:M-nn)); % m-n
			HH1 = HH1(:);
			HH2 = (H(jj+1,iSig-1+1,nn+1:M)); % m
			HH2 = HH2(:);
			Un = Un - (GG1.' * GG2 + HH1.' * HH2);
			Vn = Vn + (HH1.' * GG2 - GG1.' * HH2);

			Un = 1/NMid*Un*sigDel*nn*rayDistMid/2/lD0;
			Vn = 1/NMid*Vn*sigDel*nn*rayDistMid/2/lD0;

			U{nn}(jj+1,iSig-1+1) = Un/(sigLst(iSig-1+1)+1);
			V{nn}(jj+1,iSig-1+1) = Vn/(sigLst(iSig-1+1)+1);
		end

		AA = A{nn} / (sigLst(iSig+1)+1)^2;
		II = eye(2*(J)) - [B{nn}, AA; -AA, B{nn}];
		% GH_ = II\(eye(2*(J)) * [G{nn}(0+1:J-1+1,ii-1+1);...
		GH_ = II\(eye(2*(J)) * [G(0+1:J-1+1,iSig-1+1,nn);...
			   % H{nn}(0+1:J-1+1,ii-1+1)] + ...
	   		H(0+1:J-1+1,iSig-1+1,nn)] + ...
			[U{nn}(0+1:J-1+1,iSig-1+1);V{nn}(0+1:J-1+1,iSig-1+1)]);
		% G{nn}(0+1:J-1+1,ii+1) = GH_(1:J);
		% H{nn}(0+1:J-1+1,ii+1) = GH_(J+1:end);
		G(0+1:J-1+1,iSig+1,nn) = GH_(1:J);
		H(0+1:J-1+1,iSig+1,nn) = GH_(J+1:end);
	end
	% end of the stopwatch
	fprintf('Elapsed time is: %s\n', datestr(datenum(0,0,0,0,0,toc),...
		'HH:MM:SS'));
end

PRS_REF = 20e-6;
prsMag = cell(M,1);
prsPha = cell(M,1);
prsLvl = cell(M,1);
p0 = 10^(LP0/20) * PRS_REF;

N_RSV_LST = [(1:3)*NDiff, N1, N2, N1*2, N2*2, N1+N2, N1+NDiff, N2+NDiff];
for nn = N_RSV_LST
    % reserve the ultrasound and their sum frequency
%     if (N1~=nn && N2~=nn && (N1*2)~=nn && (N2*2)~=nn && (N1+N2)~=nn ...
%             && nn~=NDiff && nn~=2*NDiff && nn~=3*NDiff && ...
%             nn~=4*NDiff && nn~=5*NDiff)
    prsMag{nn} = sqrt(abs(G(:,:,nn)).^2 + abs(H(:,:,nn)).^2);
    for iSig = 0:I
        prsMag{nn}(:,iSig+1) = p0*prsMag{nn}(:,iSig+1)/(sigLst(iSig+1)+1);
        for jZeta = 0:J
            prsPha{nn}(jZeta+1,iSig+1) = -(wavnumMid * rayDistMid * ...
                sigLst(iSig+1) + (zetaLst(jZeta+1))^2 * ...
                (sigLst(iSig+1)+1)/NMid + atan2(G(jZeta+1,iSig+1,nn),...
                H(jZeta+1,iSig+1,nn)));
        end
    end
    prsLvl{nn} = 20*log10(prsMag{nn}/sqrt(2)/PRS_REF);
end

