function [evL, SL] = computeDistance(W1, W2)



L1                       = diag(sum(W1)) - W1;
[S1, V1]              =  eig(L1);
evL1                     = diag(V1); % eigenvalues

L2                       = diag(sum(W2)) - W2;
[S2, V2]              =  eig(L2);
evL2                     = diag(V2); % eigenvalues

end