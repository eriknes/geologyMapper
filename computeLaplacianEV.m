function [evL, S] = computeLaplacianEV(W1)


L1                       = diag(sum(W1)) - W1;
[S, V]              =  eig(L1);
evL                     = diag(V); % eigenvalues


end