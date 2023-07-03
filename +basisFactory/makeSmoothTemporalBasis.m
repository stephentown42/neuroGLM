function bases = makeSmoothTemporalBasis(shape, duration, nBases, binfun, nlOffset)
% function bases = makeSmoothTemporalBasis(shape, duration, nBases, binfun, nlOffset)
%
% Generic funciton for making smooth temporal basis functions of different
% shapes
%
% Input
%   shape: 'raised cosine', 'boxcar' or 'nonlinearly scaled cosine'
%   duration: two-element array of the time that needs to be covered
%   nBases: number of basis vectors to fill the duration
%   binfun: function that converts time to bins
%   nlOffset: time to add to center of each basis function
%
% Output
%   BBstm: basis object
%
% Example:
%   basisFactory.makeSmoothTemporalBasis(
%       'nonlinearly scaled cosine',...
%       [30 5000],...
%       8,...
%       @(t)(t==0)+ceil(t/10;),...
%       10);

if nargin < 1
    help basisFactory.makeSmoothTemporalBasis
    bases=[];
    return
end

nkbins = binfun(duration); % number of bins for the basis functions
basis_timebin = repmat((1:nkbins)', 1, nBases); % time indices for basis

switch shape
    case 'nonlinearly scaled cosine'
        
        nonlinearity = @(x)(log(x + 1e-20));
        invert_nonlin = @(x)(exp(x) - 1e-20);
        
        if nlOffset < 0
            error('nlOffset must be greater than 0');
        end
        
        % Get non-linear (logarithmic) range of time values
        if numel(duration)==2
            x_range = nonlinearity(duration + nlOffset);
        else
            x_range = nonlinearity([0 duration] + nlOffset);
        end

        binSize = 1e3 / binfun(1e3);
        x_interval = diff(x_range) / (nBases-1);          % spacing between peaks
        x_centers = x_range(1): x_interval : x_range(2);        
        b_centers = invert_nonlin(x_centers);
        max_bin = invert_nonlin(x_range(2)+2*x_interval) - nlOffset;  
        
        iht = (0:binSize:max_bin)';
        basis_timebin = repmat(iht/binSize, 1, nBases)+1;

        % Nonlinear (max, min) Cosine function
        ff = @(x, c, dc) (cos(max(-pi, min(pi, (x-c)*pi/dc/2))) + 1)/2;

        BBstm = ff(...
            repmat(nonlinearity(iht + nlOffset), 1, nBases),...
            repmat(x_centers, numel(iht), 1),...
            x_interval);
        
        
    case 'raised cosine'
        %   ^
        %  / \
        % /   \______
        %      ^
        %     / \
        % ___/   \___
        %         ^
        %        / \
        % ______/   \
        % For raised cosine, the spacing between the centers must be 1/4 of the
        % width of the cosine        
        dbcenter = nkbins / (3 + nBases); % spacing between bumps
        width = 4 * dbcenter; % width of each bump
        b_centers = 2 * dbcenter + dbcenter*(0:nBases-1);
        % location of each bump centers
        % bcenters = round(bcenters); % round to closest ms <-- BAD!!! NEVER DO THIS
        bfun = @(x,period)((abs(x/period)<0.5).*(cos(x*2*pi/period)*.5+.5));
        BBstm = bfun(basis_timebin-repmat(b_centers,nkbins,1), width);

    case 'boxcar'
        width = nkbins / nBases;
        BBstm = zeros(size(basis_timebin));
        b_centers = width * (1:nBases) - width/2;
        for k = 1:nBases
            idx = basis_timebin(:, k) > ceil(width * (k-1)) & basis_timebin(:, k) <= ceil(width * k);
            BBstm(idx, k) = 1 / sum(idx);
        end
    otherwise
        error('Unknown basis shape');
end

bases = Basis();
bases.type = [shape '@' mfilename];
bases.shape = shape;
bases.duration = duration;
bases.nBases = nBases;
bases.binfun = binfun;
bases.B = BBstm;
bases.edim = size(bases.B, 2);
bases.tr = basis_timebin - 1;
bases.centers = b_centers;
