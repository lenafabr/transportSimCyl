%% calculate cargo capture time for random microtubule configurations
%% set parameters for the analytical model
domlenreal = 10; %domain length in microns
domradreal = 1; %domain radius in microns
diffconstreal = 0.01; % diffusion coefficient in um^2/s
mtsizereal = 0.4; % capture region size in microns
vreal = 1; % velocity in um/s

domlen = 1; % dimensionless length unit
tau = domlenreal^2./diffconstreal; %time units in seconds

domrad = domradreal/domlenreal; % dimensionless radius
MTrad = mtsizereal/domlenreal/2; % dimensionless capture radius
diffconst = diffconstreal*tau/domlenreal^2; % dimensionless diffusion coefficient
vel = vreal*(tau/domlenreal); % dimensionless velocity
%% read microtubule lengths from file
% other microtubule configurations can also be specified here
% format: nmt x nconfig
fname = "mtlengths.txt";
mtlengths = readmatrix(fname)'/domlenreal;
%% calculate MFPT
tic

% initial particle position
% default: start at cell tip
% set to true for particles initiated uniformly
startunif = false;

% options for particle capture

% reflecting boundary at cell body?
domrefbound = false; % set to true for reflecting boundary at the soma

% include time to reach cell body?
getTimetoSoma = false;

% where are particles captured?
% default: capture at microtubule tip
% set to true for particles captured anywhere on the microtubule
bulkcapture = false;

% options to convert 3D microtubule configurations to 1D intervals
options_mt = struct();
options_mt.comlen = MTrad; % length of comet (can be changed for anisotropic capture regions)

mfpt = zeros(size(mtlengths,2),1);
for trc = 1:size(mtlengths,2)
	% calculate the positions and absorbance rates in 1D intervals
	[xvals,kavals] = mtlens2xvals(mtlengths(:,trc),MTrad,diffconst,bulkcapture,...
																	domrad,domlen,options_mt);
	% calculate MFPT
	mfpt(trc) = multiAbsReg_full(xvals,kavals,diffconst,startunif,domrefbound,getTimetoSoma,vel);
end
mfpt = mfpt*tau; % convert MFPT to real units
toc
