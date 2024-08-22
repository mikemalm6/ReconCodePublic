%
%   Michael Malmberg modified this function from Brian Hargreaves' gridlut
%   University of Utah 2024
%
%   function [indxminAll,indxmaxAll, indyminAll, indymaxAll,wts] = 
%       gridprep(ktraj, gridsize, convhalfwidth, kerneltable)
%
%	Function does pre-gridding with a pre-calculated density function,
%		and gridding kernel from the given kerneltable (lookup table).
%
%   Inputs:
%   	ktraj		K-space sample locations (complex).
%   	gridsize 	Size of grid onto which to grid data.
%   	convhalfwidth	(optional) width of convolution kernal in grid points.
%   	kerneltable	Convolution kernel lookup table.
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
%   

function [indxminAll,indxmaxAll, indyminAll, indymaxAll,wts] = gridprep(ktraj,gridsize,convwidth, kerneltable)



% --------- Check Enough Inputs ------
if (nargin < 1)
	disp('Usage:  dat=gridprep(ktraj,gridsize,convhalfwidth,kerneltable)'); 
	disp(' ');
	error('Not enough input arguments');
end

% --------- Check k-space trajectory is complex-valued -------
s = size(ktraj);
if (s(2)==2)
	ktraj = ktraj(:,1)+1i.*ktraj(:,2);
elseif (s(1)==2)
    ktraj = ktraj(1,:)+1i.*ktraj(2,:);
end
if isreal(ktraj)
	ktraj = ktraj + 2*eps*1i;
	disp('Warning: k-space trajectory should be complex-valued, or Nx2');
end



% --------- Check k-space file does not exceed 0.5 -------
if ( max(abs(real(ktraj(:))))>=0.5)  || (max(abs(imag(ktraj(:))))>=0.5 )
	disp('Warning:  k-space location radii should not exceed 0.5.');
	disp('		Some samples could be ignored.');
end


% --------- Default Arguments -------------
if (nargin < 2)
	gridsize=256;
end
if (nargin < 3)
	convwidth=2;
end
if (nargin < 4)
	kerneltable=[1 0];	 
end

% --------- Call Mex Function for this ---------

[indxminAll,indxmaxAll, indyminAll, indymaxAll,wts] = gridprep_mex(ktraj,gridsize,convwidth,kerneltable);


