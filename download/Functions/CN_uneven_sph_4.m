function G = CN_uneven_sph_4(D,R,t,S,m_s,minerals)
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
            
            dR2(N,1) = R(N,1)-R(N-1,1);
            p(N,1) = t/(dR2(N,1))^2;
            
            if strcmp(minerals,'garnet')
                if strcmp(m_s,'single')
                    %creat matrix
                    A = zeros(N);
                    B = zeros(N);
                    %edge condition
                    %set r
                    A(1,1) = 1+3*D(1,1)*j(1,1);
                    A(1,2) = -3*D(1,1)*j(1,1);
                    B(1,1) = 1-3*D(1,1)*j(1,1);
                    B(1,2) = 3*D(1,1)*j(1,1);
                    %set r
                    A(N,N-1) = -D(N,1)*p(N,1);
                    A(N,N) = 1+D(N,1)*p(N,1);
                    B(N,N-1) = D(N,1)*p(N,1);
                    B(N,N) = 1-D(N,1)*p(N,1);
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1) = l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        B(i,i-1) = -l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1) = -l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
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
                    end %starting condition
                    
                    %edge condition1
                    A = zeros(3*N);
                    A(1,1) = 1+3*D(1,1)*j(1,1);
                    A(1,2) = -3*D(1,1)*j(1,1);
                    A(1,N+1) = 3*D(1,2)*j(1,1);
                    A(1,N+2) = -3*D(1,2)*j(1,1);
                    A(1,2*N+1) = 3*D(1,3)*j(1,1);
                    A(1,2*N+2) = -3*D(1,3)*j(1,1);
                    
                    A(N+1,1) = 3*D(1,4)*j(1,1);
                    A(N+1,2) = -3*D(1,4)*j(1,1);
                    A(N+1,N+1) = 1+3*D(1,5)*j(1,1);
                    A(N+1,N+2) = -3*D(1,5)*j(1,1);
                    A(N+1,2*N+1) = 3*D(1,6)*j(1,1) ;
                    A(N+1,2*N+2) = -3*D(1,6)*j(1,1);
                    
                    A(2*N+1,1) = 3*D(1,7)*j(1,1);
                    A(2*N+1,2) = -3*D(1,7)*j(1,1);
                    A(2*N+1,N+1) = 3*D(1,8)*j(1,1);
                    A(2*N+1,N+2) = -3*D(1,8)*j(1,1);
                    A(2*N+1,2*N+1) = 1+3*D(1,9)*j(1,1);
                    A(2*N+1,2*N+2) = -3*D(1,9)*j(1,1);
                    
                    B = zeros(3*N);
                    B(1,1) = 1-3*D(1,1)*j(1,1);
                    B(1,2) = 3*D(1,1)*j(1,1);
                    B(1,N+1) = -3*D(1,2)*j(1,1);
                    B(1,N+2) = 3*D(1,2)*j(1,1);
                    B(1,2*N+1) = -3*D(1,3)*j(1,1);
                    B(1,2*N+2) = 3*D(1,3)*j(1,1);
                    
                    B(N+1,1) = -3*D(1,4)*j(1,1);
                    B(N+1,2) = 3*D(1,4)*j(1,1);
                    B(N+1,N+1) = 1-3*D(1,5)*j(1,1);
                    B(N+1,N+2) = 3*D(1,5)*j(1,1);
                    B(N+1,2*N+1) = -3*D(1,6)*j(1,1);
                    B(N+1,2*N+2) = 3*D(1,6)*j(1,1);
                    
                    B(2*N+1,1) = -3*D(1,7)*j(1,1);
                    B(2*N+1,2) = 3*D(1,7)*j(1,1);
                    B(2*N+1,N+1) = -3*D(1,8)*j(1,1);
                    B(2*N+1,N+2) = 3*D(1,8)*j(1,1);
                    B(2*N+1,2*N+1) = 1-3*D(1,9)*j(1,1);
                    B(2*N+1,2*N+2) = 3*D(1,9)*j(1,1);
                    
                    %edge condition N
                    A(N,N) = 1+D(N,1)*p(N,1);
                    A(N,N-1) = -D(N,1)*p(N,1);
                    A(N,2*N) = D(N,2)*p(N,1);
                    A(N,2*N-1) = -D(N,2)*p(N,1);
                    A(N,3*N) = D(N,3)*p(N,1);
                    A(N,3*N-1) = -D(N,3)*p(N,1);
                    
                    A(2*N,N) = D(N,4)*p(N,1);
                    A(2*N,N-1) = -D(N,4)*p(N,1);
                    A(2*N,2*N) = 1+D(N,5)*p(N,1);
                    A(2*N,2*N-1) = -D(N,5)*p(N,1);
                    A(2*N,3*N) = D(N,6)*p(N,1);
                    A(2*N,3*N-1) = -D(N,6)*p(N,1);
                    
                    A(3*N,N) = D(N,7)*p(N,1);
                    A(3*N,N-1) = -D(N,7)*p(N,1);
                    A(3*N,2*N) = D(N,8)*p(N,1);
                    A(3*N,2*N-1) = -D(N,8)*p(N,1);
                    A(3*N,3*N) = 1+D(N,9)*p(N,1);
                    A(3*N,3*N-1) = -D(N,9)*p(N,1);
                    
                    B(N,N) = 1-D(N,1)*p(N,1);
                    B(N,N-1) = D(N,1)*p(N,1);
                    B(N,2*N) = -D(N,2)*p(N,1);
                    B(N,2*N-1) = D(N,2)*p(N,1);
                    B(N,3*N) = -D(N,3)*p(N,1);
                    B(N,3*N-1) = D(N,3)*p(N,1);
                    
                    B(2*N,N) = -D(N,4)*p(N,1);
                    B(2*N,N-1) = D(N,4)*p(N,1);
                    B(2*N,2*N) = 1-D(N,5)*p(N,1);
                    B(2*N,2*N-1) = D(N,5)*p(N,1);
                    B(2*N,3*N) = -D(N,6)*p(N,1);
                    B(2*N,3*N-1) = D(N,6)*p(N,1);
                    
                    B(3*N,N) = -D(N,7)*p(N,1);
                    B(3*N,N-1) = D(N,7)*p(N,1);
                    B(3*N,2*N) = -D(N,8)*p(N,1);
                    B(3*N,2*N-1) = D(N,8)*p(N,1);
                    B(3*N,3*N) = 1-D(N,9)*p(N,1);
                    B(3*N,3*N-1) = D(N,9)*p(N,1);
                    
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1) = l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i,i-1) = -l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1) = -l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i,i-1+N) = l(i,1)*D(i,2)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i+N) = l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1+N) = l(i,1)*D(i,2)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i,i-1+N) = -l(i,1)*D(i,2)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i+N) = -l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1+N) = -l(i,1)*D(i,2)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i,i-1+2*N) = l(i,1)*D(i,3)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i+2*N) = l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1+2*N) = l(i,1)*D(i,3)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i,i-1+2*N) = -l(i,1)*D(i,3)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i+2*N) = -l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1+2*N) = -l(i,1)*D(i,3)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1) = l(i,1)*D(i,4)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+N,i) = l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+N,i+1) = l(i,1)*D(i,4)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1) = -l(i,1)*D(i,4)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+N,i) = -l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+N,i+1) = -l(i,1)*D(i,4)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1+N) = l(i,1)*D(i,5)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+N,i+N) = 1+l(i,1)*D(i,5)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+N,i+1+N) = l(i,1)*D(i,5)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1+N) = -l(i,1)*D(i,5)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+N,i+N) =  1-l(i,1)*D(i,5)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+N,i+1+N) = -l(i,1)*D(i,5)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1+2*N) = l(i,1)*D(i,6)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+N,i+2*N) = l(i,1)*D(i,6)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+N,i+1+2*N) = l(i,1)*D(i,6)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1+2*N) = -l(i,1)*D(i,6)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+N,i+2*N) = -l(i,1)*D(i,6)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+N,i+1+2*N) = -l(i,1)*D(i,6)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+2*N,i-1) = l(i,1)*D(i,7)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+2*N,i) = l(i,1)*D(i,7)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+2*N,i+1) = l(i,1)*D(i,7)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+2*N,i-1) = -l(i,1)*D(i,7)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+2*N,i) = -l(i,1)*D(i,7)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+2*N,i+1) = -l(i,1)*D(i,7)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+2*N,i-1+N) = l(i,1)*D(i,8)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+2*N,i+N) = l(i,1)*D(i,8)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+2*N,i+1+N) = l(i,1)*D(i,8)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+2*N,i-1+N) = -l(i,1)*D(i,8)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+2*N,i+N) = -l(i,1)*D(i,8)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+2*N,i+1+N) = -l(i,1)*D(i,8)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+2*N,i-1+2*N) = l(i,1)*D(i,9)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+2*N,i+2*N) = 1+l(i,1)*D(i,9)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+2*N,i+1+2*N) = l(i,1)*D(i,9)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+2*N,i-1+2*N) = -l(i,1)*D(i,9)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+2*N,i+2*N) = 1-l(i,1)*D(i,9)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+2*N,i+1+2*N) = -l(i,1)*D(i,9)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
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
                    %set r
                    A(1,1) = 1+3*D(1,1)*j(1,1);
                    A(1,2) = -3*D(1,1)*j(1,1);
                    B(1,1) = 1-3*D(1,1)*j(1,1);
                    B(1,2) = 3*D(1,1)*j(1,1);
                    %set r
                    A(N,N-1) = -D(N,1)*p(N,1);
                    A(N,N) = 1+D(N,1)*p(N,1);
                    B(N,N-1) = D(N,1)*p(N,1);
                    B(N,N) = 1-D(N,1)*p(N,1);
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1) = l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        B(i,i-1) = -l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1) = -l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
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
                    for a=1:N
                        C(a,1) = S(a,1);
                    end
                    for a=N+1:2*N
                        C(a,1) = S(a-N,2);
                    end
                    for a=1:N
                        E(a,1) = S(a,3);
                    end
                    %starting condition
                    
                    %edge condition1
                    A = zeros(2*N);
                    A(1,1) = 1+3*D(1,1)*j(1,1);
                    A(1,2) = -3*D(1,1)*j(1,1);
                    A(1,N+1) = 3*D(1,2)*j(1,1);
                    A(1,N+2) = -3*D(1,2)*j(1,1);

                    A(N+1,1) = 3*D(1,3)*j(1,1);
                    A(N+1,2) = -3*D(1,3)*j(1,1);
                    A(N+1,N+1) = 1+3*D(1,4)*j(1,1);
                    A(N+1,N+2) = -3*D(1,4)*j(1,1);

                    B = zeros(2*N);
                    B(1,1) = 1-3*D(1,1)*j(1,1);
                    B(1,2) = 3*D(1,1)*j(1,1);
                    B(1,N+1) = -3*D(1,2)*j(1,1);
                    B(1,N+2) = 3*D(1,2)*j(1,1);
                    
                    B(N+1,1) = -3*D(1,3)*j(1,1);
                    B(N+1,2) = 3*D(1,3)*j(1,1);
                    B(N+1,N+1) = 1-3*D(1,4)*j(1,1);
                    B(N+1,N+2) = 3*D(1,4)*j(1,1);

                    %edge condition N
                    A(N,N) = 1+D(N,1)*p(N,1);
                    A(N,N-1) = -D(N,1)*p(N,1);
                    A(N,2*N) = D(N,2)*p(N,1);
                    A(N,2*N-1) = -D(N,2)*p(N,1);

                    A(2*N,N) = D(N,3)*p(N,1);
                    A(2*N,N-1) = -D(N,3)*p(N,1);
                    A(2*N,2*N) = 1+D(N,4)*p(N,1);
                    A(2*N,2*N-1) = -D(N,4)*p(N,1);

                    B(N,N) = 1-D(N,1)*p(N,1);
                    B(N,N-1) = D(N,1)*p(N,1);
                    B(N,2*N) = -D(N,2)*p(N,1);
                    B(N,2*N-1) = D(N,2)*p(N,1);
                    
                    B(2*N,N) = -D(N,3)*p(N,1);
                    B(2*N,N-1) = D(N,3)*p(N,1);
                    B(2*N,2*N) = 1-D(N,4)*p(N,1);
                    B(2*N,2*N-1) = D(N,4)*p(N,1);
                  
                    %mid condition
                    for i = 2:N-1
                        A(i,i-1) = l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i) = 1+l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1) = l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i,i-1) = -l(i,1)*D(i,1)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i) = 1-l(i,1)*D(i,1)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1) = -l(i,1)*D(i,1)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i,i-1+N) = l(i,1)*D(i,2)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i,i+N) = l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i,i+1+N) = l(i,1)*D(i,2)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i,i-1+N) = -l(i,1)*D(i,2)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i,i+N) = -l(i,1)*D(i,2)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i,i+1+N) = -l(i,1)*D(i,2)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1) = l(i,1)*D(i,3)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+N,i) = l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+N,i+1) = l(i,1)*D(i,3)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1) = -l(i,1)*D(i,3)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+N,i) = -l(i,1)*D(i,3)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+N,i+1) = -l(i,1)*D(i,3)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        A(i+N,i-1+N) = l(i,1)*D(i,4)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        A(i+N,i+N) = 1+l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        A(i+N,i+1+N) = l(i,1)*D(i,4)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
                        B(i+N,i-1+N) = -l(i,1)*D(i,4)*((dR2(i,1))^2-R(i,1)*dR2(i,1));
                        B(i+N,i+N) = 1-l(i,1)*D(i,4)*(R(i,1)*dR2(i,1)+R(i,1)*dR1(i,1)-(dR2(i,1))^2+(dR1(i,1))^2);
                        B(i+N,i+1+N) = -l(i,1)*D(i,4)*(-(dR1(i,1))^2-R(i,1)*dR1(i,1));
                        
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