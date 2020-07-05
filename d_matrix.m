% K-segment M-degree pseudospectral d_matrix, nodes and weights
% FIXME: haven't tested the weights
function [ D, nodes, weights, N ] = d_matrix(M, K)
  N = (M - 1) * K + 1;
  [ P, Pnodes, Pweights ] = ps_matrix(M);
  P = P*K;
  Pnodes = Pnodes/K;
  Pweights = Pweights/K;
  D = sparse(N, N);
  % equidistant segments for now
  for i = 0:K-1
    x = ( M - 1) * i + 1;
    if i == 0
      D(x:x+(M-1),x:x+(M-1)) = P;
    else
      D(x+1:x+(M-1),x:x+(M-1)) = P(2:end,:);
    end
    weights(x:x+(M-1)) = Pweights;
    off = 2 / K * i - (K - 1) / K;
    nodes(x:x+(M-1)) = Pnodes + off;
  end
end

% single segment psuedospectral matrix
function [ P, nodes, weights ] = ps_matrix(M)
  [ nodes, weights ] = lglnodes(M-1);
  weights = flip(weights)';
  nodes = flip(nodes)';

  L = legendreP(M-1, nodes);
  for i = 1:M
    for j = 1:M
      if i == 1 && j == 1
        P(i,j) = - M * ( M - 1) / 4;
      elseif i == M && j == M
        P(i,j) = M * (M - 1 ) / 4;
      elseif i ~= j
        P(i,j) = L(i) / L(j) / ( nodes(i) - nodes(j) );
      else
        P(i,j) = 0;
      end
    end
  end
end
