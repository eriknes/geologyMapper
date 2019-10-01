function [pt1, pt2, ptCtr] = findDimensionIndices(center, n, radius)

if center - radius < 1
    pt1 = 1;
else
    pt1 = center - radius;
end

if center + radius > n
    pt2 = n;
else
    pt2 = center+radius;
end

ptCtr = 1 + (center - pt1);

end