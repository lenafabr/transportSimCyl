function [xvals,kavals] = mtlens2xvals(mtlens,mtrad,Dc,bulkcapture,domrad,domlen,options)
% Converts MT lengths to region boundaries
% Inputs:
% mtlens - microtubule lengths (nmt x 1)
% mtrad - mt capture radius
% ka_tip - absorption rate at the tip
% ka_bulk - absorption rate in bulk
% domlen - domain length
%%
% mtlens = [1,2,3]';
% mtrad = 0.4;
% domlen = 5;
opt = struct();
% simple addition for overlapping tips?
opt.addtips = true;
% use one sided convention instead of two sided
% one sided: MT tip region extends from mtlens to mtlens-2*mtrad
% two sided (default): MT tip region extends from mtlens+mtrad to
% mtlens-mtrad
opt.onesided = true;
% add absorbing region near soma?
opt.addAbsRegNearSoma = false;
% use rates from 3D sims?
opt.use3Drate = false;
% max number of mts allowed
maxnmt = 14;
% use anisotropic comet
opt.comlen = -1;

if(exist('options','var')); opt = copyStruct(options,opt); end

if(opt.comlen<0)
	comlen = mtrad;
else
	comlen = opt.comlen;
end

mtlens = sort(mtlens);
if(size(mtlens,2)>size(mtlens,1)); mtlens = mtlens'; end
if(opt.addAbsRegNearSoma)
	mtlens = [mtlens;0];
end
if(length(mtlens)>maxnmt)
	error('not set up for more than %d MTs',maxnmt);
end

captrate1D = @(a) (-8).*Dc.*domrad.^2.*(a.^4+(-4).*a.^2.*domrad.^2+3.*domrad.^4....
																												+4.*domrad.^4.*log(a.*domrad.^(-1))).^(-1);

ka_bulk = 0;

% set up capture rates
if(opt.use3Drate)
	% use capture rates obtained from sims with randomly dispersed MTs
	ka_tip = [0;0.00879656563517817;0.0193577234976811;
						0.0311796674230274;0.0468546970780642;0.0632344474181689;
						0.0830403595891463;0.105341444771911;0.130798742053226;
						0.156539372060372;0.182802162478362;0.214866227206204;
						0.244938712223344;0.281650590636203;0.317736845202610];
elseif(opt.addtips)
	% linearly add capture rates
	ka0 = captrate1D(mtrad);
	ka_tip = [0;(1:14)'*ka0];
else
	% calculate capture rates for effective area
	avals = sqrt(1:14)'*mtrad;
	ka_tip = [0;captrate1D(avals)];
end
if(bulkcapture>0)
		ka_bulk = ka_tip;
end


%%
if(opt.onesided)
	% left end of mts
	xmlist = max(mtlens-2*comlen,0);
	% right end of mts
	xplist = min(mtlens,domlen);
else
	% left end of mts
	xmlist = max(mtlens-comlen,0);
	% right end of mts
	xplist = min(mtlens+comlen,domlen);
end

allchangepts = unique(sort([xmlist;xplist]));
allchangepts(allchangepts==0|allchangepts==domlen) = [];
xvals = [allchangepts;domlen]';
nbulk = zeros(length(xvals),1);
ntip = zeros(length(xvals),1);

for cc = 1:length(xvals)
	nbulk(cc) = nnz(xmlist>=xvals(cc));
	ntip(cc) = nnz(xmlist<xvals(cc))-nnz(xplist<xvals(cc));
end

if(any(ka_bulk>0))
	kavals = (flipud(ka_tip(ntip+1)+ka_bulk(nbulk+1)))';
else
	kavals = (flipud(ka_tip(ntip+1)))';
end

xvals = [domlen-fliplr(xvals(1:end-1)),domlen];

end

