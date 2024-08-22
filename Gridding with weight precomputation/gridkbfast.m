%   Modified by Michael Malmberg from Brian Hargreaves' gridkb function
%   [dat]  = gridkbfast(ksamps,dcf,ixmin,ixmax,iymin,iymax,wtsall)
%
%	Function does gridding with a pre-calculated density function,
%		and Kaiser-Bessel gridding kernel as given through precalculated 
%       weights and corresponding grid indices. Based off of gridkb
%       function by Brian Hargreaves, which utilizes gridding based on
%		Jackson et al, 1991. 
%   This function makes the second part of a pair with gridkbprep.m which 
%   performs the precomputation of the weights and indices.
%
%	ksamps	K-space data (complex).
%	dcf		Density compensation factors at each sample.
%   ixmin   Low x index of grid for each data point
%   ixmax   High x index of grid for each data point
%   iymin   Low y index of grid for each data point
%   iymax   High x index of grid for each data point
%   wtsall     weights for each data point, as determined by gridkbprepMM.m

%	dat		Gridded data (grid goes from -.5 to 0.5 in Kx and Ky).



function [dat] = gridkbfast(varargin)



% ================= First assign input arguments =====================

% 	Get k-space trajectory, DCFs and data, either from parameter 
% 		list or from file.

if (~ischar(varargin{1}) )
	ksamps = varargin{1};
	dcf = varargin{2};
    ixmin = varargin{3};
    ixmax = varargin{4};
    iymin = varargin{5};
    iymax = varargin{6};
    wtsall = varargin{7};
else
    error('Must input all required values');
end

%	---- Now get other arguments, or assign defaults. ----
%
gridsize = max(max(ixmax(:)),max(iymax(:)))+1; % add 1 because it was done in c, and indices start at 0

[dat] = gridfast(ksamps,dcf,gridsize,ixmin,ixmax,iymin,iymax,wtsall);



