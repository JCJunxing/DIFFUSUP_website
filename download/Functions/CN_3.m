 function G = CN_3(D,N,L,t,S,m_s,minerals,Boundary)
            %D: diffusion coefficient or coefficients matrix;
            %N: spot numers;
            %L: step length;
            %t_interval: time step;
            %S: starting profile;
            %m_s: multi or single
            %c_g: carbonate or garnet
            if strcmp(minerals,'garnet')
                if strcmp(m_s,'single')
                    %creat matrix
                    A = zeros(N);
                    B = zeros(N);
                    %edge condition
                    if isempty(Boundary)
                        Boundary_0 = S(1);
                    else
                        Boundary_0 = Boundary(1,1);
                    end
                    S(1) = Boundary_0;
                    
                    %edge matrix
                    A(1,1) = 1;
                    B(1,1) = 1;
                    %set r
                    r = D(N,1)*t/(2*(L^2));
                    A(N,N) = 1+2*r;
                    A(N,N-1) = -2*r;
                    B(N,N) = 1-2*r;
                    B(N,N-1) = 2*r;
                    %mid condition
                    for i = 2:N-1
                        r = D(i,1)*t/(2*(L^2));
                        A(i,i-1) = -r;
                        A(i,i) = 1+2*r;
                        A(i,i+1) = -r;
                        B(i,i-1) = r;
                        B(i,i) = 1-2*r;
                        B(i,i+1) = r;
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
                    %boundary condition
                    if isempty(Boundary)
                        Boundary_1_0 = S(1,1);
                        Boundary_2_0 = S(1,2);
                        Boundary_3_0 = S(1,3);
                    else
                        Boundary_1_0 = Boundary(1,1);
                        Boundary_2_0 = Boundary(1,2);
                        Boundary_3_0 = Boundary(1,3);
                    end
                        S(1,1) = Boundary_1_0;
                        S(1,2) = Boundary_2_0;
                        S(1,3) = Boundary_3_0;
                    
                    
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
                    
                    
                    
                    %edge matrix1 
                    A = zeros(3*N);
                    A(1,1) = 1;
                    A(N+1,N+1) = 1;
                    A(2*N+1,2*N+1) = 1;
                    B = zeros(3*N);
                    B(1,1) = 1;
                    B(N+1,N+1) = 1;
                    B(2*N+1,2*N+1) = 1;
                    %edge condition N
                    r1 = D(N,1)*t/(2*(L^2));
                    r2 = D(N,2)*t/(2*(L^2));
                    r3 = D(N,3)*t/(2*(L^2));
                    r4 = D(N,4)*t/(2*(L^2));
                    r5 = D(N,5)*t/(2*(L^2));
                    r6 = D(N,6)*t/(2*(L^2));
                    r7 = D(N,7)*t/(2*(L^2));
                    r8 = D(N,8)*t/(2*(L^2));
                    r9 = D(N,9)*t/(2*(L^2));
                    A(N,N) = 1+2*r1;
                    A(N,N-1) = -2*r1;
                    A(N,2*N) = 2*r2;
                    A(N,2*N-1) = -2*r2;
                    A(N,3*N) = 2*r3;
                    A(N,3*N-1) = -2*r3;
                    A(2*N,N) = 2*r4;
                    A(2*N,N-1) = -2*r4;
                    A(2*N,2*N) = 1+2*r5;
                    A(2*N,2*N-1) = -2*r5;
                    A(2*N,3*N) = 2*r6;
                    A(2*N,3*N-1) = -2*r6;
                    A(3*N,N) = 2*r7;
                    A(3*N,N-1) = -2*r7;
                    A(3*N,2*N) = 2*r8;
                    A(3*N,2*N-1) = -2*r8;
                    A(3*N,3*N) = 1+2*r9;
                    A(3*N,3*N-1) = -2*r9;
                    B(N,N) = 1-2*r1;
                    B(N,N-1) = 2*r1;
                    B(N,2*N) = -2*r2;
                    B(N,2*N-1) = 2*r2;
                    B(N,3*N) = -2*r3;
                    B(N,3*N-1) = 2*r3;
                    B(2*N,N) = -2*r4;
                    B(2*N,N-1) = 2*r4;
                    B(2*N,2*N) = 1-2*r5;
                    B(2*N,2*N-1) = 2*r5;
                    B(2*N,3*N) = -2*r6;
                    B(2*N,3*N-1) = 2*r6;
                    B(3*N,N) = -2*r7;
                    B(3*N,N-1) = 2*r7;
                    B(3*N,2*N) = -2*r8;
                    B(3*N,2*N-1) = 2*r8;
                    B(3*N,3*N) = 1-2*r9;
                    B(3*N,3*N-1) = 2*r9;
                    %mid condition
                    for i = 2:N-1
                        r1 = D(i,1)*t/(2*(L^2));
                        r2 = D(i,2)*t/(2*(L^2));
                        r3 = D(i,3)*t/(2*(L^2));
                        r4 = D(i,4)*t/(2*(L^2));
                        r5 = D(i,5)*t/(2*(L^2));
                        r6 = D(i,6)*t/(2*(L^2));
                        r7 = D(i,7)*t/(2*(L^2));
                        r8 = D(i,8)*t/(2*(L^2));
                        r9 = D(i,9)*t/(2*(L^2));
                        A(i,i-1) = -r1;
                        A(i,i) = 1+2*r1;
                        A(i,i+1) = -r1;
                        B(i,i-1) = r1;
                        B(i,i) = 1-2*r1;
                        B(i,i+1) = r1;
                        A(i,i-1+N) = -r2;
                        A(i,i+N) = 2*r2;
                        A(i,i+1+N) = -r2;
                        B(i,i-1+N) = r2;
                        B(i,i+N) = -2*r2;
                        B(i,i+1+N) = r2;
                        A(i,i-1+2*N) = -r3;
                        A(i,i+2*N) = 2*r3;
                        A(i,i+1+2*N) = -r3;
                        B(i,i-1+2*N) = r3;
                        B(i,i+2*N) = -2*r3;
                        B(i,i+1+2*N) = r3;
                        
                        A(i+N,i-1) = -r4;
                        A(i+N,i) = 2*r4;
                        A(i+N,i+1) = -r4;
                        B(i+N,i-1) = r4;
                        B(i+N,i) = -2*r4;
                        B(i+N,i+1) = r4;
                        A(i+N,i-1+N) = -r5;
                        A(i+N,i+N) = 1+2*r5;
                        A(i+N,i+1+N) = -r5;
                        B(i+N,i-1+N) = r5;
                        B(i+N,i+N) = 1-2*r5;
                        B(i+N,i+1+N) = r5;
                        A(i+N,i-1+2*N) = -r6;
                        A(i+N,i+2*N) = 2*r6;
                        A(i+N,i+1+2*N) = -r6;
                        B(i+N,i-1+2*N) = r6;
                        B(i+N,i+2*N) = -2*r6;
                        B(i+N,i+1+2*N) = r6;
                        
                        A(i+2*N,i-1) = -r7;
                        A(i+2*N,i) = 2*r7;
                        A(i+2*N,i+1) = -r7;
                        B(i+2*N,i-1) = r7;
                        B(i+2*N,i) = -2*r7;
                        B(i+2*N,i+1) = r7;
                        A(i+2*N,i-1+N) = -r8;
                        A(i+2*N,i+N) = 2*r8;
                        A(i+2*N,i+1+N) = -r8;
                        B(i+2*N,i-1+N) = r8;
                        B(i+2*N,i+N) = -2*r8;
                        B(i+2*N,i+1+N) = r8;
                        A(i+2*N,i-1+2*N) = -r9;
                        A(i+2*N,i+2*N) = 1+2*r9;
                        A(i+2*N,i+1+2*N) = -r9;
                        B(i+2*N,i-1+2*N) = r9;
                        B(i+2*N,i+2*N) = 1-2*r9;
                        B(i+2*N,i+1+2*N) = r9;
                        
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
                        Boundary_0 = S(1);
                    else
                        Boundary_0 = Boundary(1,1);
                    end
                    S(1) = Boundary_0;
                    
                    %edge condition
                    A(1,1) = 1;
                    B(1,1) = 1;
                    %set r
                    r = D(N,1)*t/(2*(L^2));
                    A(N,N) = 1+2*r;
                    A(N,N-1) = -2*r;
                    B(N,N) = 1-2*r;
                    B(N,N-1) = 2*r;
                    %mid condition
                    for i = 2:N-1
                        r = D(i,1)*t/(2*(L^2));
                        A(i,i-1) = -r;
                        A(i,i) = 1+2*r;
                        A(i,i+1) = -r;
                        B(i,i-1) = r;
                        B(i,i) = 1-2*r;
                        B(i,i+1) = r;
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
                        Boundary_1_0 = S(1,1);
                        Boundary_2_0 = S(1,2);
                    else
                        Boundary_1_0 = Boundary(1,1);
                        Boundary_2_0 = Boundary(1,2);
                    end
                    S(1,1) = Boundary_1_0;
                    S(1,2) = Boundary_2_0;
                    
                    for a=1:N
                        C(a,1) = S(a,1);
                    end
                    for a=N+1:2*N
                        C(a,1) = S(a-N,2);
                    end
                    for a=1:N
                        E(a,1) = S(a,3);
                    end 
                    
                    A = zeros(2*N);
                    A(1,1) = 1;
                    A(N+1,N+1) = 1;
                    B = zeros(2*N);
                    B(1,1) = 1;
                    B(N+1,N+1) = 1;
                    %edge condition N
                    r1 = D(N,1)*t/(2*(L^2));
                    r2 = D(N,2)*t/(2*(L^2));
                    r3 = D(N,3)*t/(2*(L^2));
                    r4 = D(N,4)*t/(2*(L^2));
                    A(N,N) = 1+2*r1;
                    A(N,N-1) = -2*r1;
                    A(N,2*N) = 2*r2;
                    A(N,2*N-1) = -2*r2;
                    A(2*N,N) = 2*r3;
                    A(2*N,N-1) = -2*r3;
                    A(2*N,2*N) = 1+2*r4;
                    A(2*N,2*N-1) = -2*r4;
                    B(N,N) = 1-2*r1;
                    B(N,N-1) = 2*r1;
                    B(N,2*N) = -2*r2;
                    B(N,2*N-1) = 2*r2;
                    B(2*N,N) = -2*r3;
                    B(2*N,N-1) = 2*r3;
                    B(2*N,2*N) = 1-2*r4;
                    B(2*N,2*N-1) = 2*r4;
                    %mid condition
                    for i = 2:N-1
                        r1 = D(i,1)*t/(2*(L^2));
                        r2 = D(i,2)*t/(2*(L^2));
                        r3 = D(i,3)*t/(2*(L^2));
                        r4 = D(i,4)*t/(2*(L^2));
                        A(i,i-1) = -r1;
                        A(i,i) = 1+2*r1;
                        A(i,i+1) = -r1;
                        B(i,i-1) = r1;
                        B(i,i) = 1-2*r1;
                        B(i,i+1) = r1;
                        A(i,i-1+N) = -r2;
                        A(i,i+N) = 2*r2;
                        A(i,i+1+N) = -r2;
                        B(i,i-1+N) = r2;
                        B(i,i+N) = -2*r2;
                        B(i,i+1+N) = r2;
                        A(i+N,i-1) = -r3;
                        A(i+N,i) = 2*r3;
                        A(i+N,i+1) = -r3;
                        B(i+N,i-1) = r3;
                        B(i+N,i) = -2*r3;
                        B(i+N,i+1) = r3;
                        A(i+N,i-1+N) = -r4;
                        A(i+N,i+N) = 1+2*r4;
                        A(i+N,i+1+N) = -r4;
                        B(i+N,i-1+N) = r4;
                        B(i+N,i+N) = 1-2*r4;
                        B(i+N,i+1+N) = r4;
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