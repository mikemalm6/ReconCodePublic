%   Modified by Michael Malmberg, Oct 2022, from Brian Hargeaves' gridlut
%	function [dat] = gridfast(ksamps,dcf,gridsize,ixmin,ixmax,iymin,
%                           iymax,wts)
%
%	Function does gridding with a pre-calculated d ensity function,
%	pre-computed weights, and indices. Modified from Brian Hargreaves
%	gridkb function to split it up between precomputation and gridding
%
%	ksamps		K-space data (complex).
%	dcf		Density compensation factors at each sample.
%	gridsize 	Size of grid onto which to grid data.
%   ixmin       low x index of grid affected by each data point
%   ixmax       high x index of grid affected by each data point
%   iymin       low y index of grid affected by each data point
%   iymax       high y index of grid affected by each data point
%	wts         kernel weights for each index for each data point, as given
%               by gridkbprep.m
%
%	dat		Gridded data (grid goes from -.5 to 0.5 in Kx and Ky).
%

function [dat] = gridfast(ksamps,dcf,gridsize,ixmin,ixmax,iymin,iymax,wtsall)


% --------- Check Enough Inputs ------
if (nargin < 8)
	disp('Usage:  dat=gridfast(ksamps,dcf,gridsize,ixmin,ixmax,iymin,iymax,wtsall)'); 
	disp(' ');
	error('Not enough input arguments');
end


% --------- Check k-space samples are complex-valued -------
s = size(ksamps);
if (s(2)==2)
	ksamps = reshape(ksamps(:,1)+1i.*ksamps(:,2),s([1 3:length(s)]));
elseif (s(1)==2)
    ksamps = reshape(ksamps(1,:)+1i.*ksamps(2,:),s(2:end));
end
if isreal(ksamps)
	ksamps = ksamps + 2*eps*1i;
	disp('Warning k-space samples should be complex-valued, or Nx2');
end


% --------- Check DCFs do not exceed 1.0 -------
if ( max(abs(dcf)) > 1.0)
	disp('Warning:  DCF values should be between 0 and 1.0');
	disp('		(Some Recon Programs may fail)');
end

% --------- Call Mex Function for this ---------

[dat] = gridfast_mex(ksamps,dcf,gridsize,ixmin,ixmax,iymin,iymax,wtsall);


