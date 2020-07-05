 function [varargout]=legsrd(n)

 % x=legsrd(n) returns n Legendre-Gauss-Radau points with x(1)=-1.
%  [x,w]= legsrd(n) returns n Legendre-Gauss-Radau points and weights
%  Eigenmethod is used for computing nodes.
%  Last modified on August 30, 2011

 j=[0:n-2];                  % indices    
 A=diag(1./((2*j+1).*(2*j+3)));   % Main diagonal
 j=[1:n-2];
 A=A+diag(sqrt(j.*(j+1))./(2*j+1),1) ...
  +diag(sqrt(j.*(j+1))./(2*j+1),-1);     %  Create Jacobi matrix
 x= sort(eig(sparse(A)));              %  Compute eigenvalues
 x=[-1;x];
     varargout{1}=x;
  if nargout==1, return; end;
y=lepoly(n-1,x);
varargout{2}= (1-x)./(n^2*y.^2);  % return the weights by (3.179)


 






