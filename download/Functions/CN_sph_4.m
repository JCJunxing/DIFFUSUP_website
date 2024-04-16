function G = CN_sph_4(D,R1,N,L,t,S,m_s,c_g)
            %D: diffusion coefficient or coefficients matrix;
            %R1: first spot radius to mid
            %N: spot numers;
            %L: step length;
            %t_interval: time step;
            %S: starting profile;
            %m_s: multi or single
            %c_g: carbonate or garnet
            if strcmp(c_g,'garnet')
                if strcmp(m_s,'single')
                    %creat matrix
                    A = zeros(N);
                    B = zeros(N);
                    %edge condition
                    %set r
                    r = D(1,1)*t/(2*(L^2));
                    %edge condition
                    A(1,1) = 1+6*r;
                    A(1,2) = -6*r;
                    B(1,1) = 1-6*r;
                    B(1,2) = 6*r;
                    %set r
                    r = D(N,1)*t/(2*(L^2));
                    A(N,N) = 1+2*r;
                    A(N,N-1) = -2*r;
                    B(N,N) = 1-2*r;
                    B(N,N-1) = 2*r;
                    %mid condition
                    for i = 2:N-1
                        r = D(i,1)*t/(2*(R1+(i-1)*L)*(L^2));
                        A(i,i-1) = -r*((R1+(i-1)*L)-L);
                        A(i,i) = 1+2*r*(R1+(i-1)*L);
                        A(i,i+1) = -r*((R1+(i-1)*L)+L);
                        B(i,i-1) = r*((R1+(i-1)*L)-L);
                        B(i,i) = 1-2*r*(R1+(i-1)*L);
                        B(i,i+1) = r*((R1+(i-1)*L)+L);
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
                    r1 = D(1,1)*t/(2*(L^2));
                    r2 = D(1,2)*t/(2*(L^2));
                    r3 = D(1,3)*t/(2*(L^2));
                    r4 = D(1,4)*t/(2*(L^2));
                    r5 = D(1,5)*t/(2*(L^2));
                    r6 = D(1,6)*t/(2*(L^2));
                    r7 = D(1,7)*t/(2*(L^2));
                    r8 = D(1,8)*t/(2*(L^2));
                    r9 = D(1,9)*t/(2*(L^2));
                    A = zeros(3*N);
                    A(1,1) = 1+6*r1;
                    A(1,2) = -6*r1;
                    A(1,N+1) = 6*r2;
                    A(1,N+2) = -6*r2;
                    A(1,2*N+1) = 6*r3;
                    A(1,2*N+2) = -6*r3;
                    A(N+1,1) = 6*r4;
                    A(N+1,2) = -6*r4;
                    A(N+1,N+1) = 1+6*r5;
                    A(N+1,N+2) = -6*r5;
                    A(N+1,2*N+1) = 6*r6;
                    A(N+1,2*N+2) = -6*r6;
                    A(2*N+1,1) = 6*r7;
                    A(2*N+1,2) = -6*r7;
                    A(2*N+1,N+1) = 6*r8;
                    A(2*N+1,N+2) = -6*r8;
                    A(2*N+1,2*N+1) = 1+6*r9;
                    A(2*N+1,2*N+2) = -6*r9;
                    B = zeros(3*N);
                    B(1,1) = 1-6*r1;
                    B(1,2) = 6*r1;
                    B(1,N+1) = -6*r2;
                    B(1,N+2) = 6*r2;
                    B(1,2*N+1) = -6*r3;
                    B(1,2*N+2) = 6*r3;
                    B(N+1,1) = -6*r4;
                    B(N+1,2) = 6*r4;
                    B(N+1,N+1) = 1-6*r5;
                    B(N+1,N+2) = 6*r5;
                    B(N+1,2*N+1) = -6*r6;
                    B(N+1,2*N+2) = 6*r6;
                    B(2*N+1,1) = -6*r7;
                    B(2*N+1,2) = 6*r7;
                    B(2*N+1,N+1) = -6*r8;
                    B(2*N+1,N+2) = 6*r8;
                    B(2*N+1,2*N+1) = 1-6*r9;
                    B(2*N+1,2*N+2) = 6*r9;
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
                        R(i,1) = R1+(i-1)*L;
                        r1 = D(i,1)*t/(2*R(i,1)*(L^2));
                        r2 = D(i,2)*t/(2*R(i,1)*(L^2));
                        r3 = D(i,3)*t/(2*R(i,1)*(L^2));
                        r4 = D(i,4)*t/(2*R(i,1)*(L^2));
                        r5 = D(i,5)*t/(2*R(i,1)*(L^2));
                        r6 = D(i,6)*t/(2*R(i,1)*(L^2));
                        r7 = D(i,7)*t/(2*R(i,1)*(L^2));
                        r8 = D(i,8)*t/(2*R(i,1)*(L^2));
                        r9 = D(i,9)*t/(2*R(i,1)*(L^2));
                        
                        A(i,i-1) = -r1*(R(i,1)-L);
                        A(i,i) = 1+2*r1*R(i,1);
                        A(i,i+1) = -r1*(R(i,1)+L);
                        B(i,i-1) = r1*(R(i,1)-L);
                        B(i,i) = 1-2*r1*R(i,1);
                        B(i,i+1) = r1*(R(i,1)+L);
                        A(i,i-1+N) = -r2*(R(i,1)-L);
                        A(i,i+N) = 2*r2*R(i,1);
                        A(i,i+1+N) = -r2*(R(i,1)+L);
                        B(i,i-1+N) = r2*(R(i,1)-L);
                        B(i,i+N) = -2*r2*R(i,1);
                        B(i,i+1+N) = r2*(R(i,1)+L);
                        A(i,i-1+2*N) = -r3*(R(i,1)-L);
                        A(i,i+2*N) = 2*r3*R(i,1);
                        A(i,i+1+2*N) = -r3*(R(i,1)+L);
                        B(i,i-1+2*N) = r3*(R(i,1)-L);
                        B(i,i+2*N) = -2*r3*R(i,1);
                        B(i,i+1+2*N) = r3*(R(i,1)+L);
                        
                        A(i+N,i-1) = -r4*(R(i,1)-L);
                        A(i+N,i) = 2*r4*R(i,1);
                        A(i+N,i+1) = -r4*(R(i,1)+L);
                        B(i+N,i-1) = r4*(R(i,1)-L);
                        B(i+N,i) = -2*r4*R(i,1);
                        B(i+N,i+1) = r4*(R(i,1)+L);
                        A(i+N,i-1+N) = -r5*(R(i,1)-L);
                        A(i+N,i+N) = 1+2*r5*R(i,1);
                        A(i+N,i+1+N) = -r5*(R(i,1)+L);
                        B(i+N,i-1+N) = r5*(R(i,1)-L);
                        B(i+N,i+N) = 1-2*r5*R(i,1);
                        B(i+N,i+1+N) = r5*(R(i,1)+L);
                        A(i+N,i-1+2*N) = -r6*(R(i,1)-L);
                        A(i+N,i+2*N) = 2*r6*R(i,1);
                        A(i+N,i+1+2*N) = -r6*(R(i,1)+L);
                        B(i+N,i-1+2*N) = r6*(R(i,1)-L);
                        B(i+N,i+2*N) = -2*r6*R(i,1);
                        B(i+N,i+1+2*N) = r6*(R(i,1)+L);
                        
                        A(i+2*N,i-1) = -r7*(R(i,1)-L);
                        A(i+2*N,i) = 2*r7*R(i,1);
                        A(i+2*N,i+1) = -r7*(R(i,1)+L);
                        B(i+2*N,i-1) = r7*(R(i,1)-L);
                        B(i+2*N,i) = -2*r7*R(i,1);
                        B(i+2*N,i+1) = r7*(R(i,1)+L);
                        A(i+2*N,i-1+N) = -r8*(R(i,1)-L);
                        A(i+2*N,i+N) = 2*r8*R(i,1);
                        A(i+2*N,i+1+N) = -r8*(R(i,1)+L);
                        B(i+2*N,i-1+N) = r8*(R(i,1)-L);
                        B(i+2*N,i+N) = -2*r8*R(i,1);
                        B(i+2*N,i+1+N) = r8*(R(i,1)+L);
                        A(i+2*N,i-1+2*N) = -r9*(R(i,1)-L);
                        A(i+2*N,i+2*N) = 1+2*r9*R(i,1);
                        A(i+2*N,i+1+2*N) = -r9*(R(i,1)+L);
                        B(i+2*N,i-1+2*N) = r9*(R(i,1)-L);
                        B(i+2*N,i+2*N) = 1-2*r9*R(i,1);
                        B(i+2*N,i+1+2*N) = r9*(R(i,1)+L);
                        
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
            elseif strcmp(c_g,'carbonate')
                if strcmp(m_s,'single')
                    %creat matrix
                    A = zeros(N);
                    B = zeros(N);
                    %edge condition
                    %set r
                    r = D(1,1)*t/(2*(L^2));
                    %edge condition
                    A(1,1) = 1+6*r;
                    A(1,2) = -6*r;
                    B(1,1) = 1-6*r;
                    B(1,2) = 6*r;
                    %set r
                    r = D(N,1)*t/(2*(L^2));
                    A(N,N) = 1+2*r;
                    A(N,N-1) = -2*r;
                    B(N,N) = 1-2*r;
                    B(N,N-1) = 2*r;
                    %mid condition
                    for i = 2:N-1
                        r = D(i,1)*t/(2*(R1+(i-1)*L)*(L^2));
                        A(i,i-1) = -r*((R1+(i-1)*L)-L);
                        A(i,i) = 1+2*r*(R1+(i-1)*L);
                        A(i,i+1) = -r*((R1+(i-1)*L)+L);
                        B(i,i-1) = r*((R1+(i-1)*L)-L);
                        B(i,i) = 1-2*r*(R1+(i-1)*L);
                        B(i,i+1) = r*((R1+(i-1)*L)+L);
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
                    end %starting condition
                    %edge condition1
                    r1 = D(1,1)*t/(2*(L^2));
                    r2 = D(1,2)*t/(2*(L^2));
                    r3 = D(1,3)*t/(2*(L^2));
                    r4 = D(1,4)*t/(2*(L^2));
                    A = zeros(2*N);
                    A(1,1) = 1+6*r1;
                    A(1,2) = -6*r1;
                    A(1,N+1) = 6*r2;
                    A(1,N+2) = -6*r2;
                    A(N+1,1) = 6*r3;
                    A(N+1,2) = -6*r3;
                    A(N+1,N+1) = 1+6*r4;
                    A(N+1,N+2) = -6*r4;
                    B = zeros(2*N);
                    B(1,1) = 1-6*r1;
                    B(1,2) = 6*r1;
                    B(1,N+1) = -6*r2;
                    B(1,N+2) = 6*r2;
                    B(N+1,1) = -6*r3;
                    B(N+1,2) = 6*r3;
                    B(N+1,N+1) = 1-6*r4;
                    B(N+1,N+2) = 6*r4;
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
                        R(i,1) = R1+(i-1)*L;
                        r1 = D(i,1)*t/(2*R(i,1)*(L^2));
                        r2 = D(i,2)*t/(2*R(i,1)*(L^2));
                        r3 = D(i,3)*t/(2*R(i,1)*(L^2));
                        r4 = D(i,4)*t/(2*R(i,1)*(L^2));
                        A(i,i-1) = -r1*(R(i,1)-L);
                        A(i,i) = 1+2*r1*R(i,1);
                        A(i,i+1) = -r1*(R(i,1)+L);
                        B(i,i-1) = r1*(R(i,1)-L);
                        B(i,i) = 1-2*r1*R(i,1);
                        B(i,i+1) = r1*(R(i,1)+L);
                        A(i,i-1+N) = -r2*(R(i,1)-L);
                        A(i,i+N) = 2*r2*R(i,1);
                        A(i,i+1+N) = -r2*(R(i,1)+L);
                        B(i,i-1+N) = r2*(R(i,1)-L);
                        B(i,i+N) = -2*r2*R(i,1);
                        B(i,i+1+N) = r2*(R(i,1)+L);
                        A(i+N,i-1) = -r3*(R(i,1)-L);
                        A(i+N,i) = 2*r3*R(i,1);
                        A(i+N,i+1) = -r3*(R(i,1)+L);
                        B(i+N,i-1) = r3*(R(i,1)-L);
                        B(i+N,i) = -2*r3*R(i,1);
                        B(i+N,i+1) = r3*(R(i,1)+L);
                        A(i+N,i-1+N) = -r4*(R(i,1)-L);
                        A(i+N,i+N) = 1+2*r4*R(i,1);
                        A(i+N,i+1+N) = -r4*(R(i,1)+L);
                        B(i+N,i-1+N) = r4*(R(i,1)-L);
                        B(i+N,i+N) = 1-2*r4*R(i,1);
                        B(i+N,i+1+N) = r4*(R(i,1)+L);
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