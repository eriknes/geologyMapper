function [BWsmall,ctrPixSmall] = extractSubsection(BW,centerPixInd,radius)
    [ni,nj,nk] = size(BW);
    [ic,jc,kc] = ind2sub(size(BW), centerPixInd);
    [i1,i2,icSmall] = findDimensionIndices(ic, ni, radius);
    [j1,j2,jcSmall] = findDimensionIndices(jc, nj, radius);
    if nk > 1
        [k1,k2,kcSmall] = findDimensionIndices(kc, nk, radius);
    else
        k1 = 1;
        k2 = 1;
        kcSmall = 1;
    end
    BWsmall = BW(i1:i2,j1:j2,k1:k2);
    ctrPixSmall = sub2ind(size(BWsmall), icSmall, jcSmall, kcSmall);
end