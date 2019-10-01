function G = topoMapper(varargin)


%TOPOMAPPER Modified Mapper Algorithm For 2D/3D Images
% This method maps the topology of a 2D/higher-dim cube as a weighted,
% undirected graph starting from a point cloud. The Fast Marching 
% Method is used to compute the geodesic distance frome each point.
% 
% INPUT: 
% BW: Binary 2D/3D input image
% nVertices: number of vertices/nodes
% plotResults: create plots
% 
%
% OUTPUT: struct G
% G.Wred is the reduced matrix with nV vertices
% G.W is the full adjacency matrix of size nPoints,nPoints
% G.mu is the (nv,2) or (nV,3) matrix with vertex barycenters
% G.idx is the class 
% G.points are the linear indices of the points in the point cloud
% G.cpuTime is the computational time used
%
% [Wred,W, G] = TOPOMAPPER(BW) binary input image, uses 20 vertices by default
%
% [Wred,W, G] = TOPOMAPPER(PermCube, nV) uses nV vertices
%
% [Wred,W, G] = TOPOMAPPER(PermCube, nV, 1) also creates figures
%
% (c) Erik Nesvold (2018)

t1 = tic;

%% Check Input Parameters

switch nargin
    case 0
        error("No input arguments")
    case 1
        BW = logical(varargin{1});
        plotResults = 0;   % do not plot results
        nVertices = 20; 
    case 2
        BW = logical(varargin{1});
        plotResults = 0;   % do not plot results
        nVertices = varargin{2};   
    case 3
        BW = logical(varargin{1});
        plotResults = varargin{3}; % do not plot
        nVertices = varargin{2};
end

%% Other Fixed Parameters

nSpectralDim        = 5; % used for spectral clustering
nPoints             = nVertices*30; % number of points in point cloud
nNeighbors          = 10; % number of neighbors for each point in point cloud
weightImage         = 1; % throw point cloud only over 1s 
verbose             = 1;

nodePlotSize        = 50; % node/vertex size for plots
edgeWeight          = 20; % edge weight for plots

eps                 = 1e-10;            % threshold for Laplacian eigenvalues

%% 1 Remove points outside polygon

[ni,nj,nk]              = size(BW);

if verbose
    if nk == 1
        fprintf('2D image of size %d x %d\n', ni, nj);
    elseif nk > 1
        fprintf('3D image of size %d x %d x %d\n', ni, nj, nk);
    end
end

II                      = find(BW); % find indices of pixels = 1, linear indices
[ii,jj,kk]              = ind2sub(size(BW),II);  % full indices


%% 2 Sample a Point Cloud

% Either sample over entire cube or only over pixels = 1
% This amounts to either weighting the graph or not by the fraction of
% water/pore space in the image


% iSample is the point cloud vector
if weightImage == 0
    if nPoints < numel(BW)
        indSample         = randperm(numel(BW), nPoints); % sample nPoints
        [~,indSample,~]   = intersect(II,indSample);
        %iSample         = iSample(II(iSample));
    else
        indSample         = 1:numel(II); % sample all points
    end
else
    if nPoints < numel(II)
        indSample         = randperm(numel(II), nPoints); % sample nPoints
        %iSample         = II(iSample);
    else
        indSample         = 1:length(II); % sample all points
    end
    
end
nPoints   = numel(indSample);

%% 3  Graph representation

% Start with a small radius around each point and expand until
% at least 5 neighboring points are found

A               = zeros(nPoints); % Adjacency matrix
W               = zeros(nPoints); % Weighted adjacency matrix (connection strength)
%Wdist           = zeros(nPoints); % Weighted adjacency matrix (distance)

BWtemp = Inf(ni,nj,nk);

for i = 1:nPoints

    if verbose == 1
    div = nPoints / 10;
    if mod(i, div) < 1
        fprintf('Finished %d %%\n', round(100*i/nPoints));
    end
    end
    
    pointInd = indSample(i);
    
    n_neighbors     = 0;
    dr              = round(min(ni,nj) / 20); % radius increment
    radius          = 0; % initial search radius
    disconnectedComponent = 0;
    nFinitePixels   = 0;
    
    while n_neighbors < nNeighbors && radius < min([ni/2 nj/2]) && ~disconnectedComponent
        
        radius  = radius + dr; % expand radius
        %BWtemp = false(ni,nj,nk); % initialize cube to false
        
        % equal distance metric in all dimensions
        %bwDistVec = [(ii - ii(pointInd)).^2  (jj - jj(pointInd)).^2 ...
        %    (kk - kk(pointInd)).^2];
        %bwDist = sqrt(sum(bwDistVec,2));

        % only keep the part of BW inside the search radius
        %BWtemp(sub2ind([ni,nj,nk], ii(bwDist < radius),jj(bwDist < radius), ...
        %    kk(bwDist < radius))) = 1;

        % Geodesic distance over BWtemp from point
        [BWsmall, ctrSmall] = extractSubsection(BW,II(pointInd),radius);
        BWdist = bwdistgeodesic(BWsmall, ctrSmall, 'quasi-euclidean');
        BWdist = replaceSubsection(BWdist, BWtemp, II(pointInd), radius);
        %BWdist = bwdistgeodesic(BWtemp, II(pointInd), 'quasi-euclidean');
        BWinside = BWdist <= radius; % keep only pixels with real numbers

        % Find neighboring points
        [c,ia,~] = intersect(II(indSample), find(BWinside));
        
        
        % Do not include the point itself
        n_neighbors = length(ia) - 1;
        
        if sum(isfinite(BWdist(:))) > nFinitePixels
            disconnectedComponent = 0;
            nFinitePixels = sum(isfinite(BWdist(:)));
        else
            disconnectedComponent = 1;
        end
        
    end
    
    distVec = BWdist(c);
    [~, order] = sort(distVec); % sort in ascending order
    ia = ia(order);
    c = c(order);
    
    for j = 1:min(length(c), nNeighbors)
        if i ~= ia(j)
            A(i, ia(j)) = 1;
            W(i, ia(j)) = 1/(BWdist(c(j)));
            %Wdist(i,ia(j)) = BWdist(c(j));
        end
    end
    
end

% make adjacency matrices symmetric
A           = (A + A');
A(A~=0)     = 1;
W           = W + W';
%Wdist       = Wdist + Wdist';

% make Matlab graph structure
%Gd          = graph(Wdist);
%bwDist      = Gd.distances;

%% 4 Graph Laplacian Filter

% TODO: compute normalized Laplacian ?

W(isinf(W))             = 0;                % Set Inf to 0
L                       = diag(sum(W)) - W; % Graph Laplacian
[S,V]                   = eig(L);           % Find spectrum of L
eigVals                 = diag(V);          % Diagonal of V
nullInd                 = find(abs(eigVals) < eps, 1, 'last'); % last nullspace eigval
Ssub                    = S(:,nullInd+1:nullInd+nSpectralDim); % keep nSpectralDim eigenvectors
%Ssub    = [ii(iSample) jj(iSample)];

%{
Wt = W;
for i = 1:nPoints 
    Wt(i,:) = W(i,:) / sum(W(i,:)); 
end
Ssub = tsne_p(Wt, [], nSpectralDim);
%}
%Ssub = tsne(Ssub,[], 2, [], 30);

idx                     = kmedoids(Ssub, nVertices); % spectral clustering with k means

if plotResults && nk == 1
    figure
    imshow(1-.3*BW, 'InitialMagnification', 'fit')
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
end

%% 5 Compute barycenters and node weights

mu   = zeros(nVertices, 3);
nodeVol  = zeros(nVertices, 1);
    
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
 
    neighbors = find(W(i,:));
    [~, t] = sort(W(i,neighbors),'desc'); % sort in ascending order
    neighbors = neighbors(t);
    neighbors = neighbors (idx(neighbors) ~= idx(i));
  
    %neighbors = neighbors(idx(neighbors) ~= idx(i)); % keep only points in other vertices
    
    if ~isempty(neighbors)
        j = neighbors(1); % closest point

        % update adjacency matrices
        Ared(idx(i), idx(j)) = 1;
        Ared(idx(j), idx(i)) = 1;
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
Wred = Wred +  Wred';

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

%{
if plotResults
    figure
    plot(sort(diag(Vred)), '--')
    title('Fiedler Eigenvector')
end
%}

%% Plot graph

G = graph(Wred);

nodeWeight = 50;
edgeWeight = 30;

if plotResults && nk == 1
    f = figure;
    hold on;
    for j = 1:nVertices
        h = plot(nan,nan);
        cc{j} = get(h,'color'); % get standard Matlab colors
    end
    close(f)
    
    figure
    imshow(1-.3*BW, 'InitialMagnification', 'fit')
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

    
elseif plotResults && nk > 1
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
       highlight(g, j, 'markersize', (nodeVol(j))^(1/3)*nodeWeight); 
       highlight(g, j, 'NodeColor', cc{1}); 
       %highlight(g, j, 'NodeColor', cc{j}); 
       %text(Y(j,2), Y(j,1), num2str(j));
    end
    [I,J] = find(Wred);
    scale = max(Wred(:));
    for j = 1:length(I)
        highlight(g, I(j), J(j), 'LineWidth', edgeWeight/scale*Wred(I(j), J(j)))
        highlight(g, I(j), J(j), 'EdgeColor', [1 .2 0])
    end
end

points = II(indSample);

clear G;
G.Wred = Wred;
G.Edges     = Edges;        % edge number
G.edgeVecs  = edgeVecs;     % direction of edge
G.W =W; 
G.cpuTime = cpuTime; 
G.mu = mu;
G.points = points;
G.idx = idx;

end

