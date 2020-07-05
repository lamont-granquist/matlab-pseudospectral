
function [ D, nodes, weights, N ] = d_matrix2(M, K, meth)
  N = M;
  [x,w,A] = OCnonsymGLReig(M-2, meth);
  D = A / 2;
%  D(abs(D)<1e-9) = 0;
  nodes = x' * 2 - 1;
  weights = w' * 2;
end

% https://scicomp.stackexchange.com/questions/32918/need-an-example-legendre-gauss-radau-pseudospectral-differentiation-matrix-or-ma/33794#33794
%
function [x,w,A] = OCnonsymGLReig(n,meth)
   % code for nonsymmetric orthogonal collocation applications on 0 < x < 1
   % n - interior points
   % meth = 1,2,3,4 for Gauss, Lobatto, Radau (right), Radau (left)
   % x - collocation points
   % w - quadrature weights
   % A - 1st derivative
   na = [1 0 0 1];  nb = [1 0 1 0];  nt = n + 2;
   a = 1.0 - na(meth);  b = 1.0 - nb(meth);
   ab = r_jacobi(n,a,b);   ab(2:n,2) = sqrt(ab(2:n,2));
   T = diag(ab(2:n,2),-1) + diag(ab(:,1)) + diag(ab(2:n,2),+1);
   x = eig(T);  x=sort(x);
   x=0.5*(x+1.0);  x = [0.0;x;1.0];
   xdif = x-x'+eye(nt);     dpx = prod(xdif,2);
   w = (x.^nb(meth)).*((1.0 - x).^na(meth))./(dpx.*dpx); w = w/sum(w); % quadrature
   A = dpx./(dpx'.*xdif);  A(1:nt+1:nt*nt) = 1.0 - sum(A,2);            % derivative
end

