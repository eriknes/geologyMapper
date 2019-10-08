function G = permTopoMapper(varargin)

%PERMTOPOMAPPER Modified Mapper Algorithm For 2D/3D Images
% This method maps the topology of a 3D permeability cube as a weighted,
% undirected graph starting from a point cloud. The Fast Marching 
% Method is used to compute the geodesic distance frome each point.
% The FMM implementation in C is by Gabriel Peyre
%
% INPUT:
% permCube:     Permeability cube/map
% nVertices:    number of vertices/nodes
% plotResults:  create plots
%
% OUTPUT: struct G
% G.Ared       is the reduced, nonweighted adjacency matrix with nV vertices
% G.Wred    is the reduced adjacency matrix with nV vertices
% G.edges   is the edge number associated with Wred
% G.edgeVecs is the edge directional vector
% G.W       is the full adjacency matrix of size nPoints,nPoints
% G.mu       is the matrix with node barycenters
% G.idx     is the cluster number for the points 
% G.points  are the linear indices of the point cloud
% G.cpuTime is the computational time used
%
% [Wred,W, G] = PERMTOPOMAPPER(PermCube) uses 20 vertices by default
%
% [Wred,W, G] = PERMTOPOMAPPER(PermCube, nV) uses nV vertices
%
% [Wred,W, G] = PERMTOPOMAPPER(PermCube, nV, 1) also creates figures
%
% (c) Erik Nesvold (2018)

t1 = tic;
FMtime = 0;
getd = @(p)path(p,path);
getd('toolbox_graph/');
getd('toolbox_general/')


%% Check Input Parameters

switch length(varargin)
    case 0
        error("No input arguments")
    case 1
        permCube = varargin{1};
        plots = 0;
        nVertices = 20;
    case 2
        permCube = varargin{1};
        plots = 0;
        nVertices = varargin{2};
    case 3
        permCube = varargin{1};
        plots = varargin{3};
        nVertices = varargin{2};
end

%% Check permeability cube

if isnan(permCube)
    error('Perm Cube has NaN entries');
elseif isinf(permCube)
    error('Perm Cube has Inf entries');
elseif permCube < 0
    error('Perm Cube has negative entries');
end


%% Parameters

verbose             = 1;            % print info
nSpectralDim        = 5;            % number of spectral dims used for clustering
nPoints             = nVertices*40; % size of point cloud
nNeighbors          = 10;           % number of neighbors for each point
eps                 = 1e-10;        % threshold for Laplacian eigenvalues

%% 1 Find indices of non-zero pixels

[ni,nj,nk]              = size(permCube);

if verbose
    if nk == 1
        fprintf('2D image of size %d x %d\n', ni, nj);
    elseif nk > 1
        fprintf('3D image of size %d x %d x %d\n', ni, nj, nk);
    end
end

II                      = find(permCube); % linear indices of non-zero values
[ii,jj,kk]              = ind2sub(size(permCube),II);  % sub indices

%% 2 Sample a Point Cloud

if nPoints < numel(permCube)
    indSample         = randperm(numel(permCube), nPoints); % sample nPoints
    [~,indSample,~]   = intersect(II,indSample);
else
    indSample         = 1:numel(II); % sample all points
end

if verbose
if numel(II) == numel(indSample)
    fprintf('%d points in cloud: all nonzero pixels\n', numel(indSample));
else
    fprintf('%d points in cloud\n', numel(indSample));
end
end
nPoints   = numel(indSample);

%% 3  Find local adjacency structure

% Start with a small radius around each point and expand

A               = zeros(nPoints); % Adjacency matrix
W               = zeros(nPoints); % Weighted adjacency matrix (connection strength)

PCtemp      = Inf(ni,nj,nk);
scaleVal    = max([ni nj nk]);

meanRad = 0;

for i = 1:nPoints

    % information
    div = nPoints / 10;
    if mod(i, div) < 1
        fprintf('Finished %d %%\n', round(100*i/nPoints));
    end
    
    pointInd = indSample(i);
    
    n_neighbors     = 0;
    if nk > 1
        dr              = round(min([ni nj nk]) / 20); % radius increment
    else
        dr              = round(min([ni nj]) / 20); % radius increment
    end
    
    radius          = 0; % initial search radius
    disconnectedComponent = 0;
    nFinitePixels   = 0;
    
    while n_neighbors < nNeighbors && radius < min([ni/2 nj/2]) && ~disconnectedComponent
        
        radius  = radius + dr; % expand radius
        
        %[ibig,jbig,kbig] = ind2sub(size(permCube), II(pointInd));
        
        [PCsmall, ctrSmall] = extractSubsectionPerm(permCube,II(pointInd),radius);
        %binaryDist = bwdistgeodesic(logical(PCsmall), ctrSmall, 'quasi-euclidean');
        [ic,jc,kc] = ind2sub(size(PCsmall),ctrSmall);
        clear options
        options.constraint_map = -Inf(size(PCsmall));
        options.constraint_map(PCsmall~=0) = 1;
        
        %keepInd = find(binaryDist > radius);
        %PCsmall(binaryDist > radius) = 0;
        PCsmall(PCsmall~=0) = PCsmall(PCsmall~=0);
        %PCsmall(PCsmall==0) = 0;
        startVec = [ic,jc,kc]';
        if nk == 1
            startVec(3) = [];
            %options.end_points(3,:) = [];
        end
        t2 = tic;
        [distCube,~,~] = perform_fast_marching(PCsmall, startVec, options); 
        FMtime = FMtime + toc(t2);
        distCube = distCube*max(size(PCsmall)) / scaleVal;
        distCube = replaceSubsection(distCube, PCtemp, II(pointInd), radius);
        distFinite = distCube; % keep only pixels with real numbers
        distFinite(isinf(distCube)) = 0;

        % Find neighboring points
        [c,ia,~] = intersect(II(indSample), find(distFinite));
        
        % Do not include the point itself
        n_neighbors = length(ia) - 1;
        
        if sum(isfinite(distCube(:))) > nFinitePixels
            disconnectedComponent = 0;
            nFinitePixels = sum(isfinite(distCube(:)));
        else
            disconnectedComponent = 1;
        end
        
    end
    
    meanRad = meanRad + radius;
    
    distVec = distCube(c);
    [~, order] = sort(distVec); % sort in ascending order
    ia = ia(order);
    c = c(order);
    
    for j = 2:min(length(c), nNeighbors+1)
        if i ~= ia(j)
            A(i, ia(j)) = 1; % non-weighted adjacency
            W(i, ia(j)) = 1 /(distCube(c(j))); % inverse of distance
            %Wdist(i,ia(j)) = distCube(c(j));
        end
    end
    
end

fprintf('Mean radius is %f\n', meanRad/nPoints);

% make adjacency matrices symmetric
A           = (A + A');
A(A~=0)     = 1;
W           = (W + W') / 2;

%% 4 Graph Laplacian Filter in Mapper method


W(isinf(W))             = 0;                % Set Inf to 0
L                       = diag(sum(W)) - W; % Graph Laplacian
[S,V]                   = eig(L);           % Find spectrum of L
eigVals                 = diag(V);          % Diagonal of V
nullInd                 = find(abs(eigVals) < eps, 1, 'last'); % last nullspace eigval
Ssub                    = S(:,nullInd+1:nullInd+nSpectralDim); % keep nSpectralDim eigenvectors

idx                     = kmedoids(Ssub, nVertices); % spectral clustering with k-medoids

if plots && nk == 1
    figure
    imshow(1-.3*permCube/max(permCube(:)), 'InitialMagnification', 'fit')
    %colormap gray
    hold on
    classes = 1:nVertices;

    for i = 1:nVertices
        c = classes(i);
        ind = (idx  == c);
        plot(jj(indSample(ind)),ii(indSample(ind)), '.', 'markersize', 12);
        hold on
    end
    %l = legend;
    title('Point Cloud Colored by Cluster')
    set(gca, 'FontSize', 20)
end

%% 5 Compute barycenters and node weights

mu   = zeros(nVertices, 3);     % barycenters
nodeVol  = zeros(nVertices, 1); % proportion of points
    
for i = 1:nVertices
    coords = [ii(indSample(idx == i)) ...
        jj(indSample(idx == i)) kk(indSample(idx == i))];
    if size(coords,1) == 1
        mu(i,:)   = coords;
    else
        mu(i,:)   = round(median(coords));
    end
    nodeVol(i)    = sum(idx == i) / length(idx);
end


%% 6 Compute reduced graph matrices

% initialize adjacency matrices for reduced graph
Ared = zeros(nVertices);
Wred = zeros(nVertices);
Edges = zeros(nVertices);
edgeVecs = {};

for i = 1:nPoints

    %[~, neighbors] = sort(bwDist(i,:)); % sort in ascending order

    neighbors   = find(W(i,:));
    [~, t]      = sort(W(i,neighbors),'desc'); % sort in ascending order
    neighbors   = neighbors(t);
    % keep only points in other vertices
    neighbors   = neighbors (idx(neighbors) ~= idx(i));
    
    if ~isempty(neighbors)
        %j = neighbors(1); % closest point

        % update adjacency matrices
        %Ared(idx(i), idx(j)) = 1;
        %Ared(idx(j), idx(i)) = 1;
        % add connection strength
        %Wred(idx(i), idx(j)) = Wred(idx(i), idx(j)) + radius/bwDist(i,j);

        neighboring_nodes = unique(idx(neighbors));
        for j = neighboring_nodes'
            Ared(idx(i), j) = 1;
            Ared(j, idx(i)) = 1;
            k = find(idx(neighbors) == j);
            Wred(idx(i), j) = Wred(idx(i), j) + sum(W(i,neighbors(k)));
            Wred(j, idx(i)) = Wred(j, idx(i)) + sum(W(i,neighbors(k)));
        end

    end
end

% Make Wred symmetric
Wred = (Wred +  Wred') / 2;

[iw, jw] = find(Wred);

for i = 1:length(iw)
    Edges(iw(i), jw(i)) = i;                % Edge number
    edgeVecs{i} = mu(jw(i),:) - mu(iw(i),:);        % Edge direction
end


%% 7 Compute Laplacian and spectrum of reduced graph

Lred                    = diag(sum(Wred)) - Wred;
[~, Vred]               = eig(Lred);
eigValsRed                      = diag(Vred);
nullInd                 = find(abs(eigValsRed) < eps);


t2 = toc(t1);
cpuTime = t2 ;
fprintf('Fast Marching time is %f s\n', FMtime);
fprintf('Total time is %f s\n',cpuTime);


%% Plot graph

G = graph(Wred);

nodeWeight = 50;
edgeWeight = 30;

if plots && nk == 1
    f = figure;
    hold on;
    for j = 1:nVertices
        h = plot(nan,nan);
        cc{j} = get(h,'color'); % get standard Matlab colors
    end
    close(f)
    
    figure
    imshow(1 - permCube / max(permCube(:)), 'InitialMagnification', 'fit');
    hold on
    
    g = plot(G, 'Ydata', mu(:,1), 'Xdata', mu(:,2), 'NodeLabel', []);
    for j = 1:nVertices
       highlight(g, j, 'markersize', sqrt(nodeVol(j))*nodeWeight); 
       highlight(g, j, 'NodeColor', cc{1}); 
       %text(Y(j,2), Y(j,1), num2str(j));
    end
    [I,J] = find(Wred);
    scale = max(Wred(:));
    for j = 1:length(I)
        highlight(g, I(j), J(j), 'LineWidth', edgeWeight/scale*Wred(I(j), J(j)))
        highlight(g, I(j), J(j), 'EdgeColor', [1 .2 0])
    end
    xlabel('X')
    ylabel('Y')
    tit = sprintf('Reduced Graph with %d Nodes', length(Wred));
    title(tit);
    set(gca, 'fontsize', 20)
    
elseif plots && nk > 1
    f = figure;
    hold on;
    for j = 1:nVertices
        h = plot(nan,nan);
        cc{j} = get(h,'color'); % get standard Matlab colors
    end
    close(f)
    
    figure
    %imshow(1-.3*BW, 'InitialMagnification', 'fit')
    hold on
    
    g = plot(G, 'Ydata', mu(:,1), 'Xdata', mu(:,2), 'Zdata', mu(:,3), 'NodeLabel', []);
    for j = 1:nVertices
       highlight(g, j, 'markersize', nodeVol(j)^(1/3)*nodeWeight); 
       highlight(g, j, 'NodeColor', cc{1}); 
       %text(Y(j,2), Y(j,1), num2str(j));
    end
    [I,J] = find(Wred);
    scale = max(Wred(:));
    for j = 1:length(I)
        highlight(g, I(j), J(j), 'LineWidth', edgeWeight/scale*Wred(I(j), J(j)))
        highlight(g, I(j), J(j), 'EdgeColor', [1 .2 0])
    end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    tit = sprintf('Reduced Graph with %d Nodes', length(Wred));
    title(tit);
    set(gca, 'fontsize', 20)
end

points = II(indSample);

%% Return struct G

clear G;

G.Ared         = Ared;      % Reduced graph adjacency matrix (no weights) 
G.Wred      = Wred;         % Reduced graph adjacency matrix
G.Edges     = Edges;        % Edge number matrix
G.edgeVecs  = edgeVecs;     % Vectors corresponding to edge matrix entries
G.mu         = mu;          % Barycenter of reduced graph nodes
G.W         = W;            % Full point cloud adjacency matrix
G.points    = points;       % Point cloud (linear indices)
G.idx       = idx;          % Cluster number
G.cpuTime   = cpuTime;

end

