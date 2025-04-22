% clc;  clear;      %  Case2_THM1_static.m       ( Case 2 )
%  \tilde{D}--> Sigma=[0.01  0; 0  0.31];   % Threshold matrix  

A1=[ -1.46           0        2.428; 
     0.1643-0.5*5  0.4-5    -0.3788;
     0.3107          0        -2.23];    
    
A2=[ -1.46           0        2.428; 
     0.1643-0.5*6  0.4-6    -0.3788;
     0.3107          0        -2.23];   

A3=[ -1.46           0        2.428; 
     0.1643-0.5*7  0.4-7    -0.3788;
     0.3107          0        -2.23]; 
 
A4=[ -1.46           0        2.428; 
     0.1643-0.5*8  0.4-8    -0.3788;
     0.3107          0        -2.23];

B1=[0; 0; 1];    B2=B1;    B3=B1;   B4=B1; 
   
rank(ctrb(A1,B1)) ; 
   
C1=[1  1  0;  0  0  1];    C2=C1;    C3=C1;    C4=C1; 

D1=[0.2;  0.2; 0.0];   D2=D1;    D3=D1;    D4=D1;

alpha=0.9;   % The parameter of finite-time control
d_max=0.5;    % The upper bound of disturbance
lambda=0.9;    theta=1.2;   % theta >1,   lambda > 1/theta, 
a1=1.2;    % The parameter of Finsler's lemma
rho=0.6;   % g(t)'*g(t) <= rho x(t)'*x(t),  
c1=0.3;    % The initial energy, x(0)Rx(0) <= c1,   % R=0.1*eye(n_x); 
E_kd=0.2;  E_kc=0.3;     % E_kd=kappa_d,    E_kc=kappa_c,    
T=6;  %<-- settling time, parameter of FTC, T=30 is simulation time,
Sigma=[0.01  0; 0  0.31];  % \tilde{D}--> Sigma: threshold matrix,   
phi=0.1;     % eta(0)= phi*x(0)'*x(0),
x(:,1) = [ -1;  0.5;  0.8];   % initial state x(0), 
y_l(1)= phi*x(:,1)'*x(:,1);   % y_l(1)= \eta(0),  \eta(t)=y_l(n),   
E0 = [ 1 0 0 ; 0 1 0 ; 0 0 1 ];

PI=[-2.4    2.2    0.1    0.1;     % transition rate matrix (TRM),
     0.1   -1.9    1.7    0.1; 
     1.2    0.3   -1.8    0.3;
     2.6    0.1    0.1   -2.8];    

n_x = size(A1,1);   n_u = size(B1,2);   n_y = size(C1,1);  n_d = size(D1,2);


% LMI initialization
setlmis([]);

P1=lmivar(1,[n_x,1]);    % (n_x)*(n_x) symmetric matrix
P2=lmivar(1,[n_x,1]);
P3=lmivar(1,[n_x,1]);
P4=lmivar(1,[n_x,1]);

Omega=lmivar(1,[n_y,1]);

V1=lmivar(1,[n_u,1]);
V2=lmivar(1,[n_u,1]);
V3=lmivar(1,[n_u,1]);
V4=lmivar(1,[n_u,1]);

W1=lmivar(2,[n_u,n_y]);
W2=lmivar(2,[n_u,n_y]);
W3=lmivar(2,[n_u,n_y]);
W4=lmivar(2,[n_u,n_y]);

beta=lmivar(1,[1,0]);
c2=lmivar(1,[1,0]);
gamma=lmivar(1,[1,0]);

% The first mode
lmiterm([3 1 1 P1],1,A1,'s');
lmiterm([3 1 1 W1],a1*(1-E_kd)*(1-E_kc)*B1,C1,'s');
lmiterm([3 1 1 P1],PI(1,1),1);
lmiterm([3 1 1 P2],PI(1,2),1);
lmiterm([3 1 1 P3],PI(1,3),1);
lmiterm([3 1 1 P4],PI(1,4),1);
lmiterm([3 1 1 P1],-alpha,1);
lmiterm([3 1 1 0],rho*eye(n_x));
lmiterm([3 1 1 Omega],(1+1/theta)*C1'*Sigma,C1);    % \tilde{D}--> Sigma,  
lmiterm([3 1 2 W1],a1*(1-E_kd)*(1-E_kc)*B1,1);
lmiterm([3 1 2 Omega],(1+1/theta)*C1'*Sigma,1);   
lmiterm([3 1 3 P1],1,D1);
lmiterm([3 1 4 W1],a1*(1-E_kd)*(E_kc)*B1,1); 
lmiterm([3 1 5 P1],1,B1);
lmiterm([3 1 5 V1],-a1*B1,1);
lmiterm([3 1 5 -W1],a1*(1-E_kd)*(1-E_kc)*C1',1);
lmiterm([3 2 2 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);
lmiterm([3 2 5 -W1],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([3 3 3 0],-eye(n_d));
lmiterm([3 4 4 0],-eye(n_y));
lmiterm([3 4 5 -W1],a1*(1-E_kd)*(E_kc),1);
lmiterm([3 5 5 V1],-a1,1,'s');

lmiterm([4 1 1 beta],-1,1);   % beta > 1,  i.e.  -beta + 1 < 0,
lmiterm([4 1 1 0],1);

lmiterm([5 1 1 gamma],-1,1);   % gamma > 0,

lmiterm([-6 1 1 P1],1,1);    % (7): R < P_{i} < beta R,   % -P1 + R < 0,
lmiterm([6 1 1 0],0.1*eye(n_x));    % R = 0.1*eye(n_x),
lmiterm([7 1 1 P1],1,1);     % P1-beta R < 0, 
lmiterm([-7 1 1 beta],0.1*eye(n_x),1);  


% The second mode
lmiterm([8 1 1 P2],1,A2,'s');
lmiterm([8 1 1 W2],a1*(1-E_kd)*(1-E_kc)*B2,C2,'s');
lmiterm([8 1 1 P1],PI(2,1),1);
lmiterm([8 1 1 P2],PI(2,2),1);
lmiterm([8 1 1 P3],PI(2,3),1);
lmiterm([8 1 1 P4],PI(2,4),1);
lmiterm([8 1 1 P2],-alpha,1);
lmiterm([8 1 1 0],rho*eye(n_x));
lmiterm([8 1 1 Omega],(1+1/theta)*C2'*Sigma,C2);    % \tilde{D}--> Sigma,
lmiterm([8 1 2 W2],a1*(1-E_kd)*(1-E_kc)*B2,1);
lmiterm([8 1 2 Omega],(1+1/theta)*C2'*Sigma,1);     
lmiterm([8 1 3 P2],1,D2);
lmiterm([8 1 4 W2],a1*(1-E_kd)*(E_kc)*B2,1);
lmiterm([8 1 5 P2],1,B2);
lmiterm([8 1 5 V2],-a1*B2,1);
lmiterm([8 1 5 -W2],a1*(1-E_kd)*(1-E_kc)*C2',1);
lmiterm([8 2 2 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);
lmiterm([8 2 5 -W2],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([8 3 3 0],-eye(n_d));
lmiterm([8 4 4 0],-eye(n_y));
lmiterm([8 4 5 -W2],a1*(1-E_kd)*(E_kc),1);
lmiterm([8 5 5 V2],-a1,1,'s');

lmiterm([-10 1 1 P2],1,1);  % (7): R < P_{i} < beta R,   % -P2 + R < 0,  
lmiterm([10 1 1 0],0.1*eye(n_x));    % R = 0.1*eye(n_x),
lmiterm([11 1 1 P2],1,1);     % P2-beta R < 0, 
lmiterm([-11 1 1 beta],0.1*eye(n_x),1);  


% The third mode
lmiterm([12 1 1 P3],1,A3,'s');
lmiterm([12 1 1 W3],a1*(1-E_kd)*(1-E_kc)*B3,C3,'s');
lmiterm([12 1 1 P1],PI(3,1),1);
lmiterm([12 1 1 P2],PI(3,2),1);
lmiterm([12 1 1 P3],PI(3,3),1);
lmiterm([12 1 1 P4],PI(3,4),1);
lmiterm([12 1 1 P3],-alpha,1);
lmiterm([12 1 1 0],rho*eye(n_x));
lmiterm([12 1 1 Omega],(1+1/theta)*C3'*Sigma,C3);    % \tilde{D}--> Sigma, 
lmiterm([12 1 2 W3],a1*(1-E_kd)*(1-E_kc)*B3,1);
lmiterm([12 1 2 Omega],(1+1/theta)*C3'*Sigma,1);     
lmiterm([12 1 3 P3],1,D3);
lmiterm([12 1 4 W3],a1*(1-E_kd)*(E_kc)*B3,1);
lmiterm([12 1 5 P3],1,B3);
lmiterm([12 1 5 V3],-a1*B3,1);
lmiterm([12 1 5 -W3],a1*(1-E_kd)*(1-E_kc)*C3',1);
lmiterm([12 2 2 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);  
lmiterm([12 2 5 -W3],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([12 3 3 0],-eye(n_d));
lmiterm([12 4 4 0],-eye(n_y));
lmiterm([12 4 5 -W3],a1*(1-E_kd)*(E_kc),1);
lmiterm([12 5 5 V3],-a1,1,'s');

lmiterm([-14 1 1 P3],1,1);  % (7): R < P_{i} < beta R,   % -P3 + R < 0, 
lmiterm([14 1 1 0],0.1*eye(n_x));     % R = 0.1*eye(n_x),
lmiterm([15 1 1 P3],1,1);     % P3-beta R < 0, 
lmiterm([-15 1 1 beta],0.1*eye(n_x),1); 


% The fourth mode
lmiterm([17 1 1 P4],1,A4,'s');
lmiterm([17 1 1 W4],a1*(1-E_kd)*(1-E_kc)*B4,C4,'s');
lmiterm([17 1 1 P1],PI(4,1),1);
lmiterm([17 1 1 P2],PI(4,2),1);
lmiterm([17 1 1 P3],PI(4,3),1);
lmiterm([17 1 1 P4],PI(4,4),1);
lmiterm([17 1 1 P4],-alpha,1);
lmiterm([17 1 1 0],rho*eye(n_x));
lmiterm([17 1 1 Omega],(1+1/theta)*C4'*Sigma,C4);    % \tilde{D}--> Sigma, 
lmiterm([17 1 2 W4],a1*(1-E_kd)*(1-E_kc)*B4,1);
lmiterm([17 1 2 Omega],(1+1/theta)*C4'*Sigma,1);     
lmiterm([17 1 3 P4],1,D4);
lmiterm([17 1 4 W4],a1*(1-E_kd)*(E_kc)*B4,1);
lmiterm([17 1 5 P4],1,B4);
lmiterm([17 1 5 V4],-a1*B4,1);
lmiterm([17 1 5 -W4],a1*(1-E_kd)*(1-E_kc)*C4',1);
lmiterm([17 2 2 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);  
lmiterm([17 2 5 -W4],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([17 3 3 0],-eye(n_d));
lmiterm([17 4 4 0],-eye(n_y));
lmiterm([17 4 5 -W4],a1*(1-E_kd)*(E_kc),1);
lmiterm([17 5 5 V4],-a1,1,'s');

lmiterm([-19 1 1 P4],1,1); % (7): R < P_{i} < beta R,   % -P4 + R < 0,
lmiterm([19 1 1 0],0.1*eye(n_x));     % R = 0.1*eye(n_x),
lmiterm([20 1 1 P4],1,1);     % P4-beta R < 0, 
lmiterm([-20 1 1 beta],0.1*eye(n_x),1); 

lmiterm([25 1 1 beta],c1,exp(T*alpha));     % (8): c1 * beta * exp(T*alpha) + 
lmiterm([25 1 1 0],d_max*exp(T*alpha));     %  d * exp(T*alpha) + 
lmiterm([25 1 1 0],y_l*exp(T*alpha));      %  \eta(0) * exp(T*alpha) -c2 < 0,
lmiterm([25 1 1 c2],-1,1);    


lmiterm([-27 1 1 P1],1,1);   % P1 >0, 
lmiterm([-28 1 1 P2],1,1);   % P2 >0,
lmiterm([-29 1 1 P3],1,1);   % P3 >0,
lmiterm([-30 1 1 P4],1,1);   % P4 >0,
lmiterm([-31 1 1 c2],1,1);   % c2 >0, c2 is a positive variable,  
lmiterm([-32 1 1 c2],1,1);    % c1 < c2, 
lmiterm([32 1 1 0],c1);

lmiterm([35 1 1 P1],-1,1);       % (9)  C1'*C1 <= gamma^2 * P1,  
lmiterm([35 1 2 0],C1');
lmiterm([35 2 2 gamma],-1,1);
lmiterm([36 1 1 P1],1,1);       % (9)  P1 <= E0 - phi*I, 
lmiterm([-36 1 1 0],E0);
lmiterm([36 1 1 0],phi*eye(n_x));

lmiterm([37 1 1 P2],-1,1);        % (9)  C2'*C2 <= gamma^2 * P2,
lmiterm([37 1 2 0],C2');
lmiterm([37 2 2 gamma],-1,1);
lmiterm([38 1 1 P2],1,1);        % (9)  P2 <= E0 - phi*I, 
lmiterm([-38 1 1 0],E0);
lmiterm([38 1 1 0],phi*eye(n_x));
 
lmiterm([39 1 1 P3],-1,1);       % (9)  C3'*C3 <= gamma^2 * P3, 
lmiterm([39 1 2 0],C3');
lmiterm([39 2 2 gamma],-1,1);
lmiterm([40 1 1 P3],1,1);       % (9)  P3 <= E0 - phi*I, 
lmiterm([-40 1 1 0],E0);
lmiterm([40 1 1 0],phi*eye(n_x));

lmiterm([41 1 1 P4],-1,1);       % (9)  C4'*C4 <= gamma^2 * P4, 
lmiterm([41 1 2 0],C4');
lmiterm([41 2 2 gamma],-1,1);
lmiterm([42 1 1 P4],1,1);       % (9)  P4 <= E0 - phi*I, 
lmiterm([-42 1 1 0],E0);
lmiterm([42 1 1 0],phi*eye(n_x));


LMIs=getlmis;
n=decnbr(LMIs);
c=zeros(n,1);
for k=1:n
[ci,gammai]=defcx(LMIs,k,c2,gamma);
c(k)=ci+gammai;     % optimization objective: c2^2 + gamma^2,
end

options=[0,0,0,0,0];
[copt,xopt]=mincx(LMIs,c,options);


c2=sqrt(dec2mat(LMIs,xopt,c2))
gamma=sqrt(dec2mat(LMIs,xopt,gamma))   % sqrt(A)=root(A),

P1=dec2mat(LMIs,xopt,P1);
P2=dec2mat(LMIs,xopt,P2);
P3=dec2mat(LMIs,xopt,P3);
P4=dec2mat(LMIs,xopt,P4);

beta=sqrt(dec2mat(LMIs,xopt,beta));

V1=dec2mat(LMIs,xopt,V1);
V2=dec2mat(LMIs,xopt,V2);
V3=dec2mat(LMIs,xopt,V3);
V4=dec2mat(LMIs,xopt,V4);

W1=dec2mat(LMIs,xopt,W1);
W2=dec2mat(LMIs,xopt,W2);
W3=dec2mat(LMIs,xopt,W3);
W4=dec2mat(LMIs,xopt,W4);

K1=inv(V1)*W1;
K2=inv(V2)*W2;
K3=inv(V3)*W3;
K4=inv(V4)*W4;