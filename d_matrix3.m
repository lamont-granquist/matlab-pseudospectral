function [ DN, tau, wt, N ] = d_matrix3(M, K)
  N = M;
  [ lgt, lgw ] = legsrd(N-1);
  wt = flip(lgw);
  tau = [ -1 flip(-lgt)' ]';

  for i=1:N
    for k=1:N-1
      DN(k,i) = lagrangeD(N, i, k+1, tau);
    end
  end
  tau = tau';
  wt = wt';
end

% first derivative of lagrangeP
function dy = lagrangeD(N, j, k, tau)
  dy = 0;
  for l=1:N
    if not(l==j)
      c = 1/(tau(j)-tau(l));
      for m=1:N
        if not(m==j) && not(m==l)
          c = c*(tau(k)-tau(m))/(tau(j)-tau(m));
        end
      end
      dy = dy + c;
    end
  end
end

