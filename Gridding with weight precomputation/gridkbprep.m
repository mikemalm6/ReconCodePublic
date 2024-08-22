%
%   Michael Malmberg  - altered from Brian Hargreaves' gridkb to allow for precomputation of weights
%   University of Utah 2024
%
%   [indxminAll,indxmaxAll, indyminAll, indymaxAll,wts,kbku,kbkval] = 
%       gridkbprep(ktraj, gridsize, kwidth, overgridfactor)
%
%	Function does pre-gridding with a pre-calculated density function,
%		and Kaiser-Bessel gridding kernel, as described in 
%		Jackson et al, 1991.  The k-space locations (ktraj) are 
%		normalized (ie  |real(ktraj)|<0.5 and |imag(k)|<0.5  ).  
%       The rest of the gridding procedure is done in gridkbfast.m
%
%   Inputs:
%	    ktraj		K-space sample locations (complex), normalized so
%	    			that |real(k)|<0.5 and |imag(k)|<0.5.  
%	    gridsize 	Integer size of grid onto which to grid data.
%	    kwidth		Kernel width, w as described in Jackson's paper.
%	    overgridfactor	Factor by which these grid parameters overgrid.
%	        		(This is gridsize/(FOV/resolution) )
%
%   Outputs:
%       indxminAll  Index of minimum x grid location in cartesian space
%                   (nRO x nLin matrix)
%       indxmaxAll  Index of maximum x grid location in cartesian space
%                   (nRO x nLin matrix)
%       indyminAll  Index of minimum y grid location in cartesian space
%                   (nRO x nLin matrix)
%       indymaxAll  Index of maximum y grid location in cartesian space
%                   (nRO x nLin matrix)
%       wts         Weights for multiplication with data points and dcf
%                   based on linear interpolation of kernel
%                   (nWeights x nSamples matrix) (nSamples = nRO x nLin)
%	    kbkval		Kaiser-bessel kernel values.
%	    kbku		Kaiser-bessel radius in grid points.

function [indxminAll,indxmaxAll, indyminAll, indymaxAll,wts,kbku,kbkval] = gridkbprep(varargin)


% ================= First assign input arguments =====================

% 	Get k-space trajectory, DCFs and data, either from parameter 
% 		list or from file.

if (~ischar(varargin{1}) )
	ktraj = varargin{1};
	otherargs = {varargin{2:end}};
end

%	---- Now get other arguments, or assign defaults. ----
%

if (length(otherargs)>=1)		% -- Get gridsize. 
	gridsize = otherargs{1}; 
else
	disp('Warning:  No gridsize passed - using gridsize=256');
	gridsize = 256;
end

if (length(otherargs)>=2)		% -- Get kernel width.
	kwidth = otherargs{2}; 
else
	disp('Warning:  No kernel width passed - using width=1.5');
	kwidth = 1.5;
end

if (length(otherargs)>=3)		% -- Get overgridfactor
	overgridfactor = otherargs{3};
else
	disp('Warning:  No overgridfactor passed.  Using overgridfactor=1.5');
	disp('		This is an important parameter in selecting');
	disp('		the convolution kernel, so if you are actually');
	disp('		using these images, you should specifiy this.');
	overgridfactor = 1.5;
end

% ======= Try to load the parameters from file, if already created =======
try
% create filename that would be saved
tit = recon_makeGridPrepTitle();

% see if filename has already been saved
d = pwd;
cd /v/raid10/users/mmalmberg/GithubCode/ReconCode/mmalmCode/Gridding/GriddingWeights_PreCalc
load(tit);
cd(d);

% load if so
catch
% ========== Calculate Kaiser-Bessel Kernel for given parameters =========

[kerneltable,u] = calckbkernel(kwidth,overgridfactor);



% ========== Do Gridding =========

[indxminAll,indxmaxAll, indyminAll, indymaxAll,wts] = gridprep(ktraj,gridsize,kwidth/2, kerneltable);

end
kbkval = kerneltable;		% ---- Assign output value.
kbku = u;			% ---- Assign output value.


