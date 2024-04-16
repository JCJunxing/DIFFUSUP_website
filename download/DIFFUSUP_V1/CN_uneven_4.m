function G = CN_uneven_4(D,R,t,S,m_s,minerals,Boundary)
            %D: diffusion coefficient or coefficients matrix;
            %R: spot radius to center
            %t: time step;
            %S: starting profile;
            %m_s: multi or single
            %minerals: mineral names
            
            %creat dR1,dR2
            N = length(R);
            
            dR1(1,1) = R(2,1)-R(1,1);
            j(1,1) = t/(dR1(1,1))^2;
            
            for i=2:(N-1)
            dR1(i,1) = R(i+1,1)-R(i,1);
            dR2(i,1) = dR1(i-1,1);
            l(i,1) = t/(R(i,1)*dR1(i,1)*dR2(i,1)*(dR2(i,1)+dR1(i,1)));
            end
            
            if strcmp(minerals,'garnet')
                if strcmp(m_s,'single')
                    %creat matrix
                    A = zeros(N);
                    B = zeros(N);
                    if isempty(Boundary)
                        Boundary_end = S(end);
                    else
                        Boundary_end = Boundary(2,1);
                    end
                    S(end) = Boundary_end;
                    %edge matrix
                    %set r
                    A(1,1) = 1+D(1,1)*j(1,1);
                    A(1,2) = -D(1,1)*j(1,1);
                    B(1,1) = 1-D(1,1)*j(1,1);
                    B(1,2) = D(1,1)*j(1,1);
                    A(N,N) = 1;
                    B(N,N) = 1;
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1) = l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                        B(i,i-1) = -l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1) = -l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                    end
                    G = A\B*S;
                    
                elseif strcmp(m_s,'multi')
                    %DX,DY,DZ: diffusion coefficient;
                    %N: spot numers;
                    %L: step length;
                    %t: time;
                    %S: sample starting material with convolution(import array);
                    C = zeros(3*N,1);
                    E = zeros(N,1);
                    if isempty(Boundary)
                        Boundary_1_end = S(end,1);
                        Boundary_2_end = S(end,2);
                        Boundary_3_end = S(end,3);
                    else
                        Boundary_1_end = Boundary(2,1);
                        Boundary_2_end = Boundary(2,2);
                        Boundary_3_end = Boundary(2,3);
                    end
                        S(end,1) = Boundary_1_end;
                        S(end,2) = Boundary_2_end;
                        S(end,3) = Boundary_3_end;
                        
                    for a=1:N
                        C(a,1) = S(a,1);
                    end
                    for a=N+1:2*N
                        C(a,1) = S(a-N,2);
                    end
                    for a=2*N+1:3*N
                        C(a,1) = S(a-2*N,3);
                    end
                    for a=1:N
                        E(a,1) = S(a,4);
                    end 
                    
                    %edge condition1
                    A = zeros(3*N);
                    A(1,1) = 1+D(1,1)*j(1,1);
                    A(1,2) = -D(1,1)*j(1,1);
                    A(1,N+1) = D(1,2)*j(1,1);
                    A(1,N+2) = -D(1,2)*j(1,1);
                    A(1,2*N+1) = D(1,3)*j(1,1);
                    A(1,2*N+2) = -D(1,3)*j(1,1);
                    
                    A(N+1,1) = D(1,4)*j(1,1);
                    A(N+1,2) = -D(1,4)*j(1,1);
                    A(N+1,N+1) = 1+D(1,5)*j(1,1);
                    A(N+1,N+2) = -D(1,5)*j(1,1);
                    A(N+1,2*N+1) = D(1,6)*j(1,1) ;
                    A(N+1,2*N+2) = -D(1,6)*j(1,1);
                    
                    A(2*N+1,1) = D(1,7)*j(1,1);
                    A(2*N+1,2) = -D(1,7)*j(1,1);
                    A(2*N+1,N+1) = D(1,8)*j(1,1);
                    A(2*N+1,N+2) = -D(1,8)*j(1,1);
                    A(2*N+1,2*N+1) = 1+D(1,9)*j(1,1);
                    A(2*N+1,2*N+2) = -D(1,9)*j(1,1);
                    
                    B = zeros(3*N);
                    B(1,1) = 1-D(1,1)*j(1,1);
                    B(1,2) = D(1,1)*j(1,1);
                    B(1,N+1) = -D(1,2)*j(1,1);
                    B(1,N+2) = D(1,2)*j(1,1);
                    B(1,2*N+1) = -D(1,3)*j(1,1);
                    B(1,2*N+2) = D(1,3)*j(1,1);
                    
                    B(N+1,1) = -D(1,4)*j(1,1);
                    B(N+1,2) = D(1,4)*j(1,1);
                    B(N+1,N+1) = 1-D(1,5)*j(1,1);
                    B(N+1,N+2) = D(1,5)*j(1,1);
                    B(N+1,2*N+1) = -D(1,6)*j(1,1);
                    B(N+1,2*N+2) = D(1,6)*j(1,1);
                    
                    B(2*N+1,1) = -D(1,7)*j(1,1);
                    B(2*N+1,2) = D(1,7)*j(1,1);
                    B(2*N+1,N+1) = -D(1,8)*j(1,1);
                    B(2*N+1,N+2) = D(1,8)*j(1,1);
                    B(2*N+1,2*N+1) = 1-D(1,9)*j(1,1);
                    B(2*N+1,2*N+2) = D(1,9)*j(1,1);
                    %edge matrix N
                    A(N,N) = 1;
                    A(2*N,2*N) = 1;
                    A(3*N,3*N) = 1;
                    B(N,N) = 1;
                    B(2*N,2*N) = 1;
                    B(3*N,3*N) = 1;
                    
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1) = l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                        
                        B(i,i-1) = -l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1) = -l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                        
                        A(i,i-1+N) = l(i,1)*D(i,2)*(-R(i,1)*dR2(i,1));
                        A(i,i+N) = l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1+N) = l(i,1)*D(i,2)*(-R(i,1)*dR1(i,1));
                        
                        B(i,i-1+N) = -l(i,1)*D(i,2)*(-R(i,1)*dR2(i,1));
                        B(i,i+N) = -l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1+N) = -l(i,1)*D(i,2)*(-R(i,1)*dR1(i,1));
                        
                        A(i,i-1+2*N) = l(i,1)*D(i,3)*(-R(i,1)*dR2(i,1));
                        A(i,i+2*N) = l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1+2*N) = l(i,1)*D(i,3)*(-R(i,1)*dR1(i,1));
                        
                        B(i,i-1+2*N) = -l(i,1)*D(i,3)*(-R(i,1)*dR2(i,1));
                        B(i,i+2*N) = -l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1+2*N) = -l(i,1)*D(i,3)*(-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1) = l(i,1)*D(i,4)*(-R(i,1)*dR2(i,1));
                        A(i+N,i) = l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+N,i+1) = l(i,1)*D(i,4)*(-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1) = -l(i,1)*D(i,4)*(-R(i,1)*dR2(i,1));
                        B(i+N,i) = -l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+N,i+1) = -l(i,1)*D(i,4)*(-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1+N) = l(i,1)*D(i,5)*(-R(i,1)*dR2(i,1));
                        A(i+N,i+N) = 1+l(i,1)*D(i,5)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+N,i+1+N) = l(i,1)*D(i,5)*(-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1+N) = -l(i,1)*D(i,5)*(-R(i,1)*dR2(i,1));
                        B(i+N,i+N) =  1-l(i,1)*D(i,5)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+N,i+1+N) = -l(i,1)*D(i,5)*(-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1+2*N) = l(i,1)*D(i,6)*(-R(i,1)*dR2(i,1));
                        A(i+N,i+2*N) = l(i,1)*D(i,6)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+N,i+1+2*N) = l(i,1)*D(i,6)*(-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1+2*N) = -l(i,1)*D(i,6)*(-R(i,1)*dR2(i,1));
                        B(i+N,i+2*N) = -l(i,1)*D(i,6)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+N,i+1+2*N) = -l(i,1)*D(i,6)*(-R(i,1)*dR1(i,1));
                        
                        A(i+2*N,i-1) = l(i,1)*D(i,7)*(-R(i,1)*dR2(i,1));
                        A(i+2*N,i) = l(i,1)*D(i,7)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+2*N,i+1) = l(i,1)*D(i,7)*(-R(i,1)*dR1(i,1));
                        
                        B(i+2*N,i-1) = -l(i,1)*D(i,7)*(-R(i,1)*dR2(i,1));
                        B(i+2*N,i) = -l(i,1)*D(i,7)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+2*N,i+1) = -l(i,1)*D(i,7)*(-R(i,1)*dR1(i,1));
                        
                        A(i+2*N,i-1+N) = l(i,1)*D(i,8)*(-R(i,1)*dR2(i,1));
                        A(i+2*N,i+N) = l(i,1)*D(i,8)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+2*N,i+1+N) = l(i,1)*D(i,8)*(-R(i,1)*dR1(i,1));
                        
                        B(i+2*N,i-1+N) = -l(i,1)*D(i,8)*(-R(i,1)*dR2(i,1));
                        B(i+2*N,i+N) = -l(i,1)*D(i,8)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+2*N,i+1+N) = -l(i,1)*D(i,8)*(-R(i,1)*dR1(i,1));
                        
                        A(i+2*N,i-1+2*N) = l(i,1)*D(i,9)*(-R(i,1)*dR2(i,1));
                        A(i+2*N,i+2*N) = 1+l(i,1)*D(i,9)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+2*N,i+1+2*N) = l(i,1)*D(i,9)*(-R(i,1)*dR1(i,1));
                        
                        B(i+2*N,i-1+2*N) = -l(i,1)*D(i,9)*(-R(i,1)*dR2(i,1));
                        B(i+2*N,i+2*N) = 1-l(i,1)*D(i,9)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+2*N,i+1+2*N) = -l(i,1)*D(i,9)*(-R(i,1)*dR1(i,1));
                        
                    end
                    %matrix
                    C = A\B*C;
                    for q = 1:N
                        G(q,1)=C(q,1);
                        G(q,2)=C(q+N,1);
                        G(q,3)=C(q+2*N,1);
                        G(q,4)=1-C(q,1)-C(q+N,1)-C(q+2*N,1);
                    end
                end
            elseif strcmp(minerals,'carbonate')
                if strcmp(m_s,'single')
                    %creat matrix
                    A = zeros(N);
                    B = zeros(N);
                    %edge condition
                    if isempty(Boundary)
                        Boundary_end = S(end);
                    else
                        Boundary_end = Boundary(2,1);
                    end
                    S(end) = Boundary_end;
                    
                    %edge condition
                    %set r
                    A(1,1) = 1+D(1,1)*j(1,1);
                    A(1,2) = -D(1,1)*j(1,1);
                    B(1,1) = 1-D(1,1)*j(1,1);
                    B(1,2) = D(1,1)*j(1,1);
                    A(N,N) = 1;
                    B(N,N) = 1;
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1) = l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                        B(i,i-1) = -l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1) = -l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                    end
                    G = A\B*S;
                    
                elseif strcmp(m_s,'multi')
                    %DX,DY,DZ: diffusion coefficient;
                    %N: spot numers;
                    %L: step length;
                    %t: time;
                    %S: sample starting material with convolution(import array);
                    C = zeros(2*N,1); 
                    E = zeros(N,1);
                    if isempty(Boundary)
                        Boundary_1_end = S(end,1);
                        Boundary_2_end = S(end,2);
                    else
                        Boundary_1_end = Boundary(2,1);
                        Boundary_2_end = Boundary(2,2);
                    end
                    S(end,1) = Boundary_1_end;
                    S(end,2) = Boundary_2_end;
                    
                    for a=1:N
                        C(a,1) = S(a,1);
                    end
                    for a=N+1:2*N
                        C(a,1) = S(a-N,2);
                    end
                    for a=1:N
                        E(a,1) = S(a,3);
                    end
                        
                    %edge condition1
                    A = zeros(2*N);
                    A(1,1) = 1+D(1,1)*j(1,1);
                    A(1,2) = -D(1,1)*j(1,1);
                    A(1,N+1) = D(1,2)*j(1,1);
                    A(1,N+2) = -D(1,2)*j(1,1);

                    A(N+1,1) = D(1,3)*j(1,1);
                    A(N+1,2) = -D(1,3)*j(1,1);
                    A(N+1,N+1) = 1+D(1,4)*j(1,1);
                    A(N+1,N+2) = -D(1,4)*j(1,1);

                    B = zeros(2*N);
                    B(1,1) = 1-D(1,1)*j(1,1);
                    B(1,2) = D(1,1)*j(1,1);
                    B(1,N+1) = -D(1,2)*j(1,1);
                    B(1,N+2) = D(1,2)*j(1,1);
                    
                    B(N+1,1) = -D(1,3)*j(1,1);
                    B(N+1,2) = D(1,3)*j(1,1);
                    B(N+1,N+1) = 1-D(1,4)*j(1,1);
                    B(N+1,N+2) = D(1,4)*j(1,1);
                    %edge condition N
                    A(N,N) = 1;
                    A(2*N,2*N) = 1;
                    B(N,N) = 1;
                    B(2*N,2*N) = 1;
                    
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1) = l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                        
                        B(i,i-1) = -l(i,1)*D(i,1)*(-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1) = -l(i,1)*D(i,1)*(-R(i,1)*dR1(i,1));
                        
                        A(i,i-1+N) = l(i,1)*D(i,2)*(-R(i,1)*dR2(i,1));
                        A(i,i+N) = l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i,i+1+N) = l(i,1)*D(i,2)*(-R(i,1)*dR1(i,1));
                        
                        B(i,i-1+N) = -l(i,1)*D(i,2)*(-R(i,1)*dR2(i,1));
                        B(i,i+N) = -l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i,i+1+N) = -l(i,1)*D(i,2)*(-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1) = l(i,1)*D(i,3)*(-R(i,1)*dR2(i,1));
                        A(i+N,i) = l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+N,i+1) = l(i,1)*D(i,3)*(-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1) = -l(i,1)*D(i,3)*(-R(i,1)*dR2(i,1));
                        B(i+N,i) = -l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+N,i+1) = -l(i,1)*D(i,3)*(-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1+N) = l(i,1)*D(i,4)*(-R(i,1)*dR2(i,1));
                        A(i+N,i+N) = 1+l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        A(i+N,i+1+N) = l(i,1)*D(i,4)*(-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1+N) = -l(i,1)*D(i,4)*(-R(i,1)*dR2(i,1));
                        B(i+N,i+N) = 1-l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1));
                        B(i+N,i+1+N) = -l(i,1)*D(i,4)*(-R(i,1)*dR1(i,1));
                        
                    end
                    %matrix
                    C = A\B*C;
                    for q = 1:N
                        G(q,1)=C(q,1);
                        G(q,2)=C(q+N,1);
                        G(q,3)=1-C(q,1)-C(q+N,1);
                    end
                end
            end
        end