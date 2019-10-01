function showDistributions(EV, nGroups, eps, cols)

% params
ftsz            = 20;
numDim          = 2;
lw              = 3;
mrkSz           = 10;

m = size(EV,1);
mGroup = m / nGroups;

D               = zeros(1,size(EV,1));


nn = 0;
for i = 1:size(EV,1)
    for j = i+1:size(EV,1)
        nn          = nn+1;
        f1          = EV(i,:) + eps;
        f2          = EV(j,:) + eps;
        %[~,i1]  = find(f1 < eps);
        %[~,i2]  = find(f2 < eps);
        %imin = max([i1 i2]) + 1;
        f1          = log(f1);
        f2          = log(f2);
        D(nn)       = mean(abs(f1 - f2));
    end
end

DM              = squareform(D);
[md, st]        = mdscale(DM, numDim);

fprintf('MDS Stress: %d\n', st);

%{
figure
hold on
for i = 1:nGroups
    sel = (i-1)*mGroup+1 : i*mGroup;
    plot(md(sel,1),md(sel,2), '*', 'col', cols{i}, 'markersize', mrkSz);
end
%}


pdf = {};
xi = {};
for i = 1:nGroups
    sel = (i-1)*mGroup+1 : i*mGroup;
    [pdf{i}, xi{i}] = ksdensity(md(sel,:));
end

figure
hold on

for i = 1:nGroups
    X = reshape(xi{i}(:,1), 30, 30);
    Y = reshape(xi{i}(:,2), 30, 30);
    F = reshape(pdf{i}, 30, 30);
    
    contour(X,Y,F,3, 'col', cols{i}, 'LineStyle', '-', 'LineWidth', 2)
end

xlabel('1st MDS coordinate', 'FontSize', ftsz)
ylabel('2nd MDS coordinate', 'FontSize', ftsz)
%leg = legend('All Data','Tide Dominated','Coarse Sand');
%leg.set('FontSize', ftsz);
set(gca, 'FontSize', ftsz);

end