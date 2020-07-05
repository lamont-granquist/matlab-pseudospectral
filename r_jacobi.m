% ------- recurrence coefficients code from OPQ suite -------------------
function ab=r_jacobi(N,a,b) % recurrence coefficients for monic Jacobi polynomials.
   %    Supplied by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
   if nargin<2, a=0; end;  if nargin<3, b=a; end
   if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
   nu=(b-a)/(a+b+2);
   mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
   if N==1, ab=[nu mu]; return, end 
   N=N-1; n=1:N; nab=2*n+a+b;
   A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
   n=2:N; nab=nab(n);
   B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
   B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
   ab=[A' [mu; B1; B']];
end


