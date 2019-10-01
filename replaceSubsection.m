function BWlarge = replaceSubsection(BWsmall,BWlarge,ctrPixLarge, radius)
    [ni,nj,nk] = size(BWlarge);
    [ic,jc,kc] = ind2sub(size(BWlarge), ctrPixLarge);
    [i1,i2] = findDimensionIndices(ic, ni, radius);
    [j1,j2] = findDimensionIndices(jc, nj, radius);
    if nk > 1
        [k1,k2] = findDimensionIndices(kc, nk, radius);
    else
        k1 = 1;
        k2 = 1;
    end
    BWlarge(i1:i2,j1:j2,k1:k2) = BWsmall;
end