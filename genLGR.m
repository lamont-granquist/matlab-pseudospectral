function [LGR]=genLGR(Narray)
  for Nindex=1:length(Narray)

    N=Narray(Nindex);

    % calculate the diagonal
    indices=0:N-2;
    maindiagonal=diag(1./((2*indices+1).*(2*indices+3)));

    % calculate the jacobi matrix
    indices=1:N-2;
    Jacobimat=maindiagonal+diag(sqrt(indices.*(indices+1))./(2*indices+1),1) +diag(sqrt(indices.*(indices+1))./(2*indices+1),-1);

    nodes = [-1; sort(eig(sparse(Jacobimat)))];
    weights = (1-nodes)./(N^2*legendreP(N-1,nodes).^2);

    N=N+1;
    tau=[nodes; 1];
    Dsquare = zeros(N,N);
    for id=1:N
      z=tau(id);
      for j=1:N
        temp=0;
        for i=1:N
          if not(i==j)
            k = 1/(tau(j)-tau(i));
            for m=1:N
              if not(m==j) && not(m==i)
                k = k*(z-tau(m))/(tau(j)-tau(m));
              end
            end
            temp = temp + k;
          end
        end
        Dsquare(id,j)=temp;
      end
    end
    D = Dsquare(1:end-1,:);
    Ddiag=zeros(size(D));
    Ddiag(:,1:end-1) = diag(diag(D));

    %Save variables in a structure
    LGR.points{Nindex} = nodes;
    LGR.weights{Nindex} = weights;
    LGR.diff_matrix{Nindex}=D;
    LGR.diff_matrix_square{Nindex}=Dsquare;
    LGR.diag_diff_matrix{Nindex}=Ddiag;

  end
