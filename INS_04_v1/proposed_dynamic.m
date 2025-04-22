% clc;  clear;      %  Case2_THM2_dynamic.m      ( Case 2 )
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
d_max=0.5;   % The upper bound of disturbance
lambda=0.9;    theta=1.2;   % theta >1,   lambda > 1/theta, 
a1=1.2;    % Parameter of Finsler's lemma
rho=0.6;   % g(t)'*g(t) <= rho x(t)'*x(t),  
c1=0.3;   % The initial energy, x(0)Rx(0) <= c1,   % R1=0.1*eye(n_x); 
E_kd=0.2;  E_kc=0.3;     % E_kd=kappa_d,    E_kc=kappa_c,    
T=6;   %<-- settling time, parameter of FTC,  T=30 is simulation time,
Sigma=[0.01  0;  0  0.31];   % Threshold parameter,  % \tilde{D}--> Sigma,
phi=0.1;     % eta(0)= phi*bar_x(0)'*bar_x(0),
x(:,1)= [ -1;  0.5;  0.8];   % initial state x(0), 
xc(:,1)= [ 0; 0; 0];   % initial value of controller state xc(0), 
% bar_x(:,1)= [ -1;  0.5;  0.8; 0; 0;  0];   % initial state bar_x(0), 
% y_l(1)= phi* bar_x(:,1)'*bar_x(:,1);  % y_l(1)= \eta(0),  y_l(n)= \eta(t),
y_l(1)= phi*x(:,1)'*x(:,1);   % y_l(1)= \eta(0),   \eta(t)=y_l(n),
E01 = [ 1 0 0 ; 0 1 0 ; 0 0 1 ];   E02 = [ 1 0 0 ; 0 1 0 ; 0 0 1 ];

PI=[-2.4    2.2    0.1    0.1;     % transition rate matrix (TRM),
     0.1   -1.9    1.7    0.1; 
     1.2    0.3   -1.8    0.3;      %  R1 = 0.1*eye(n_x); 
     2.6    0.1    0.1   -2.8];     %  R2 = 0.1*eye(n_x); 

n_x = size(A1,1);   n_u = size(B1,2);   n_y = size(C1,1);  n_d = size(D1,2);


% LMI initialization
setlmis([]);

P11=lmivar(1,[n_x,1]);    % (n_x)*(n_x) symmetric matrix
P12=lmivar(1,[n_x,1]);    P13=lmivar(1,[n_x,1]);    P14=lmivar(1,[n_x,1]);

P21=lmivar(1,[n_x,1]);    P22=lmivar(1,[n_x,1]);    
P23=lmivar(1,[n_x,1]);    P24=lmivar(1,[n_x,1]);

Omega=lmivar(1,[n_y,1]);

M1=lmivar(1,[n_x,1]);     M2=lmivar(1,[n_x,1]);
M3=lmivar(1,[n_x,1]);     M4=lmivar(1,[n_x,1]);

[N11,n,sN11]=lmivar(1,[1,1]);  
[N12,n,sN12]=lmivar(1,[1,1]);  
[N13,n,sN13]=lmivar(1,[1,1]);  
[N14,n,sN14]=lmivar(1,[1,1]);  
[N15,n,sN15]=lmivar(1,[1,1]);  
[N16,n,sN16]=lmivar(1,[1,1]);  

N1=lmivar(3,[sN11,sN12;sN13,sN14;sN15,sN16]);

[N21,n,sN21]=lmivar(1,[1,1]);  
[N22,n,sN22]=lmivar(1,[1,1]);  
[N23,n,sN23]=lmivar(1,[1,1]);  
[N24,n,sN24]=lmivar(1,[1,1]);  
[N25,n,sN25]=lmivar(1,[1,1]);  
[N26,n,sN26]=lmivar(1,[1,1]);  

N2=lmivar(3,[sN21,sN22;sN23,sN24;sN25,sN26]);

[N31,n,sN31]=lmivar(1,[1,1]);  
[N32,n,sN32]=lmivar(1,[1,1]);  
[N33,n,sN33]=lmivar(1,[1,1]);  
[N34,n,sN34]=lmivar(1,[1,1]);  
[N35,n,sN35]=lmivar(1,[1,1]);  
[N36,n,sN36]=lmivar(1,[1,1]);  

N3=lmivar(3,[sN31,sN32;sN33,sN34;sN35,sN36]);

[N41,n,sN41]=lmivar(1,[1,1]);  
[N42,n,sN42]=lmivar(1,[1,1]);  
[N43,n,sN43]=lmivar(1,[1,1]);  
[N44,n,sN44]=lmivar(1,[1,1]);  
[N45,n,sN45]=lmivar(1,[1,1]);  
[N46,n,sN46]=lmivar(1,[1,1]);  

N4=lmivar(3,[sN41,sN42;sN43,sN44;sN45,sN46]);

V1=lmivar(1,[n_u,1]);      V2=lmivar(1,[n_u,1]);
V3=lmivar(1,[n_u,1]);      V4=lmivar(1,[n_u,1]);

W11=lmivar(2,[n_u,n_y]);      W12=lmivar(2,[n_u,n_y]);
W13=lmivar(2,[n_u,n_y]);      W14=lmivar(2,[n_u,n_y]);

W21=lmivar(2,[n_u,n_x]);      W22=lmivar(2,[n_u,n_x]);
W23=lmivar(2,[n_u,n_x]);      W24=lmivar(2,[n_u,n_x]);

beta=lmivar(1,[1,0]);
c2=lmivar(1,[1,0]);
gamma=lmivar(1,[1,0]);


% The first mode
lmiterm([3 1 1 P11],1,A1,'s');
lmiterm([3 1 1 W11],a1*(1-E_kd)*(1-E_kc)*B1,C1,'s');
lmiterm([3 1 1 P11],PI(1,1),1);
lmiterm([3 1 1 P12],PI(1,2),1);
lmiterm([3 1 1 P13],PI(1,3),1);
lmiterm([3 1 1 P14],PI(1,4),1);
lmiterm([3 1 1 P11],-alpha,1);
lmiterm([3 1 1 0],rho*eye(n_x));
lmiterm([3 1 1 Omega],(1+1/theta)*C1'*Sigma,C1);   % \tilde{D}--> Sigma,  
lmiterm([3 1 2 W21],a1*B1,1);
lmiterm([3 1 2 -N1],(1-E_kd)*(1-E_kc)*C1',1);
lmiterm([3 1 3 W11],a1*(1-E_kd)*(1-E_kc)*B1,1);
lmiterm([3 1 3 Omega],(1+1/theta)*C1'*Sigma,1);   
lmiterm([3 1 4 P11],1,D1);
lmiterm([3 1 5 W11],a1*(1-E_kd)*E_kc*B1,1); 
lmiterm([3 1 6 P11],1,B1);
lmiterm([3 1 6 V1],-a1*B1,1);
lmiterm([3 1 6 -W11],a1*(1-E_kd)*(1-E_kc)*C1',1);
lmiterm([3 2 2 M1],1,1,'s');
lmiterm([3 2 2 P21],PI(1,1),1);
lmiterm([3 2 2 P22],PI(1,2),1);
lmiterm([3 2 2 P23],PI(1,3),1);
lmiterm([3 2 2 P24],PI(1,4),1);
lmiterm([3 2 2 P21],-alpha,1);
lmiterm([3 2 3 N1],(1-E_kd)*(1-E_kc),1);
lmiterm([3 2 5 N1],(1-E_kd)*E_kc,1);
lmiterm([3 2 6 -W21],a1,1);
lmiterm([3 3 3 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);
lmiterm([3 3 6 -W11],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([3 4 4 0],-eye(n_d));
lmiterm([3 5 5 0],-eye(n_y));
lmiterm([3 5 6 -W11],a1*(1-E_kd)*E_kc,1); 
lmiterm([3 6 6 V1],-a1,1,'s');

lmiterm([4 1 1 gamma],-1,1);   % gamma > 0
lmiterm([5 1 1 beta],-1,1);    % beta > 1,  i.e.  -beta + 1 < 0
lmiterm([5 1 1 0],1);
lmiterm([-6 1 1 c2],1,1);   % c2 > 0, (or)  lmiterm([6 1 1 c2],-1,1);

lmiterm([-8 1 1 P11],1,1);    % P11 >0, 
lmiterm([-9 1 1 P12],1,1);    % P12 >0,
lmiterm([-10 1 1 P13],1,1);   % P13 >0,
lmiterm([-11 1 1 P14],1,1);   % P14 >0,

lmiterm([-12 1 1 P21],1,1);   % P21 >0, 
lmiterm([-13 1 1 P22],1,1);   % P22 >0,
lmiterm([-14 1 1 P23],1,1);   % P23 >0,
lmiterm([-15 1 1 P24],1,1);   % P24 >0,

lmiterm([-17 1 1 c2],1,1);    % c1 < c2, 
lmiterm([17 1 1 0],c1);


% The second mode
lmiterm([20 1 1 P12],1,A2,'s');
lmiterm([20 1 1 W12],a1*(1-E_kd)*(1-E_kc)*B2,C2,'s');
lmiterm([20 1 1 P11],PI(2,1),1);
lmiterm([20 1 1 P12],PI(2,2),1);
lmiterm([20 1 1 P13],PI(2,3),1);
lmiterm([20 1 1 P14],PI(2,4),1);
lmiterm([20 1 1 P12],-alpha,1);
lmiterm([20 1 1 0],rho*eye(n_x));
lmiterm([20 1 1 Omega],(1+1/theta)*C2'*Sigma,C2);  % \tilde{D}--> Sigma,  
lmiterm([20 1 2 W22],a1*B2,1);
lmiterm([20 1 2 -N2],(1-E_kd)*(1-E_kc)*C2',1);
lmiterm([20 1 3 W12],a1*(1-E_kd)*(1-E_kc)*B2,1);
lmiterm([20 1 3 Omega],(1+1/theta)*C2'*Sigma,1);   
lmiterm([20 1 4 P12],1,D2);
lmiterm([20 1 5 W12],a1*(1-E_kd)*E_kc*B2,1); 
lmiterm([20 1 6 P12],1,B2);
lmiterm([20 1 6 V2],-a1*B2,1);
lmiterm([20 1 6 -W12],a1*(1-E_kd)*(1-E_kc)*C2',1);
lmiterm([20 2 2 M2],1,1,'s');
lmiterm([20 2 2 P21],PI(2,1),1);
lmiterm([20 2 2 P22],PI(2,2),1);
lmiterm([20 2 2 P23],PI(2,3),1);
lmiterm([20 2 2 P24],PI(2,4),1);
lmiterm([20 2 2 P22],-alpha,1);
lmiterm([20 2 3 N2],(1-E_kd)*(1-E_kc),1);
lmiterm([20 2 5 N2],(1-E_kd)*E_kc,1);
lmiterm([20 2 6 -W22],a1,1);
lmiterm([20 3 3 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);
lmiterm([20 3 6 -W12],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([20 4 4 0],-eye(n_d));
lmiterm([20 5 5 0],-eye(n_y));
lmiterm([20 5 6 -W12],a1*(1-E_kd)*E_kc,1); 
lmiterm([20 6 6 V2],-a1,1,'s');


% The third mode
lmiterm([21 1 1 P13],1,A3,'s');
lmiterm([21 1 1 W13],a1*(1-E_kd)*(1-E_kc)*B3,C3,'s');
lmiterm([21 1 1 P11],PI(3,1),1);
lmiterm([21 1 1 P12],PI(3,2),1);
lmiterm([21 1 1 P13],PI(3,3),1);
lmiterm([21 1 1 P14],PI(3,4),1);
lmiterm([21 1 1 P13],-alpha,1);
lmiterm([21 1 1 0],rho*eye(n_x));
lmiterm([21 1 1 Omega],(1+1/theta)*C3'*Sigma,C3);  % \tilde{D}--> Sigma,  
lmiterm([21 1 2 W23],a1*B3,1);
lmiterm([21 1 2 -N3],(1-E_kd)*(1-E_kc)*C3',1);
lmiterm([21 1 3 W13],a1*(1-E_kd)*(1-E_kc)*B3,1);
lmiterm([21 1 3 Omega],(1+1/theta)*C3'*Sigma,1);   
lmiterm([21 1 4 P13],1,D3);
lmiterm([21 1 5 W13],a1*(1-E_kd)*E_kc*B3,1); 
lmiterm([21 1 6 P13],1,B3);
lmiterm([21 1 6 V3],-a1*B3,1);
lmiterm([21 1 6 -W13],a1*(1-E_kd)*(1-E_kc)*C3',1);
lmiterm([21 2 2 M3],1,1,'s');
lmiterm([21 2 2 P21],PI(3,1),1);
lmiterm([21 2 2 P22],PI(3,2),1);
lmiterm([21 2 2 P23],PI(3,3),1);
lmiterm([21 2 2 P24],PI(3,4),1);
lmiterm([21 2 2 P23],-alpha,1);
lmiterm([21 2 3 N3],(1-E_kd)*(1-E_kc),1);
lmiterm([21 2 5 N3],(1-E_kd)*E_kc,1);
lmiterm([21 2 6 -W23],a1,1);
lmiterm([21 3 3 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);
lmiterm([21 3 6 -W13],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([21 4 4 0],-eye(n_d));
lmiterm([21 5 5 0],-eye(n_y));
lmiterm([21 5 6 -W13],a1*(1-E_kd)*E_kc,1); 
lmiterm([21 6 6 V3],-a1,1,'s');


% The fourth mode
lmiterm([22 1 1 P14],1,A4,'s');
lmiterm([22 1 1 W14],a1*(1-E_kd)*(1-E_kc)*B4,C4,'s');
lmiterm([22 1 1 P11],PI(4,1),1);
lmiterm([22 1 1 P12],PI(4,2),1);
lmiterm([22 1 1 P13],PI(4,3),1);
lmiterm([22 1 1 P14],PI(4,4),1);
lmiterm([22 1 1 P14],-alpha,1);
lmiterm([22 1 1 0],rho*eye(n_x));
lmiterm([22 1 1 Omega],(1+1/theta)*C4'*Sigma,C4);  % \tilde{D}--> Sigma,  
lmiterm([22 1 2 W24],a1*B4,1);
lmiterm([22 1 2 -N4],(1-E_kd)*(1-E_kc)*C4',1);
lmiterm([22 1 3 W14],a1*(1-E_kd)*(1-E_kc)*B4,1);
lmiterm([22 1 3 Omega],(1+1/theta)*C4'*Sigma,1);   
lmiterm([22 1 4 P14],1,D4);
lmiterm([22 1 5 W14],a1*(1-E_kd)*E_kc*B4,1); 
lmiterm([22 1 6 P14],1,B4);
lmiterm([22 1 6 V4],-a1*B4,1);
lmiterm([22 1 6 -W14],a1*(1-E_kd)*(1-E_kc)*C4',1);
lmiterm([22 2 2 M4],1,1,'s');
lmiterm([22 2 2 P21],PI(4,1),1);
lmiterm([22 2 2 P22],PI(4,2),1);
lmiterm([22 2 2 P23],PI(4,3),1);
lmiterm([22 2 2 P24],PI(4,4),1);
lmiterm([22 2 2 P24],-alpha,1);
lmiterm([22 2 3 N4],(1-E_kd)*(1-E_kc),1);
lmiterm([22 2 5 N4],(1-E_kd)*E_kc,1);
lmiterm([22 2 6 -W24],a1,1);
lmiterm([22 3 3 Omega],-(1+1/theta)*(eye(n_y)-Sigma),1);
lmiterm([22 3 6 -W14],a1*(1-E_kd)*(1-E_kc),1);
lmiterm([22 4 4 0],-eye(n_d));
lmiterm([22 5 5 0],-eye(n_y));
lmiterm([22 5 6 -W14],a1*(1-E_kd)*E_kc,1); 
lmiterm([22 6 6 V4],-a1,1,'s');

lmiterm([-23 1 1 P11],1,1);    % (17):  {R1} < P_{1i} < beta {R1},
lmiterm([23 1 1 0],0.1*eye(n_x));   % -P_{1i}+ R1 < 0,  R1 = 0.1*eye(n_x),
lmiterm([24 1 1 P11],1,1);     % P_{1i} - beta R1 < 0, 
lmiterm([-24 1 1 beta],0.1*eye(n_x),1); 

lmiterm([-26 1 1 P12],1,1);    % (17):  {R1} < P_{1i} < beta {R1},
lmiterm([26 1 1 0],0.1*eye(n_x));   % -P_{1i}+ R1 < 0,  R1 = 0.1*eye(n_x),
lmiterm([28 1 1 P12],1,1);     % P_{1i} - beta R1 < 0, 
lmiterm([-28 1 1 beta],0.1*eye(n_x),1); 

lmiterm([-30 1 1 P13],1,1);    % (17):  {R1} < P_{1i} < beta {R1},
lmiterm([30 1 1 0],0.1*eye(n_x));   % -P_{1i}+ R1 < 0,  R1 = 0.1*eye(n_x),
lmiterm([32 1 1 P13],1,1);     % P_{1i} - beta R1 < 0, 
lmiterm([-32 1 1 beta],0.1*eye(n_x),1);

lmiterm([-33 1 1 P14],1,1);    % (17):  {R1} < P_{1i} < beta {R1},
lmiterm([33 1 1 0],0.1*eye(n_x));   % -P_{1i}+ R1 < 0,  R1 = 0.1*eye(n_x),
lmiterm([34 1 1 P14],1,1);     % P_{1i} - beta R1 < 0, 
lmiterm([-34 1 1 beta],0.1*eye(n_x),1);

lmiterm([-36 1 1 P21],1,1);    % (17):  {R2} < P_{2i} < beta {R2},
lmiterm([36 1 1 0],0.1*eye(n_x));   % -P_{2i}+ R2 < 0,  R2 = 0.1*eye(n_x),
lmiterm([37 1 1 P21],1,1);     % P_{2i} - beta R2 < 0, 
lmiterm([-37 1 1 beta],0.1*eye(n_x),1); 

lmiterm([-38 1 1 P22],1,1);    % (17):  {R2} < P_{2i} < beta {R2},
lmiterm([38 1 1 0],0.1*eye(n_x));   % -P_{2i}+ R2 < 0,  R2 = 0.1*eye(n_x),
lmiterm([39 1 1 P22],1,1);     % P_{2i} - beta R2 < 0, 
lmiterm([-39 1 1 beta],0.1*eye(n_x),1);

lmiterm([-40 1 1 P23],1,1);    % (17):  {R2} < P_{2i} < beta {R2},
lmiterm([40 1 1 0],0.1*eye(n_x));   % -P_{2i} +R2 < 0,  R2 = 0.1*eye(n_x),
lmiterm([42 1 1 P23],1,1);     % P_{2i} - beta R2 < 0, 
lmiterm([-42 1 1 beta],0.1*eye(n_x),1);

lmiterm([-43 1 1 P24],1,1);    % (17):  {R2} < P_{2i} < beta {R2},
lmiterm([43 1 1 0],0.1*eye(n_x));   % -P_{2i}+ R2 < 0,  R2 = 0.1*eye(n_x),
lmiterm([45 1 1 P24],1,1);     % P_{2i} - beta R2 < 0, 
lmiterm([-45 1 1 beta],0.1*eye(n_x),1);

lmiterm([46 1 1 beta],c1,exp(T*alpha));    % (18): c1*beta* exp(T*alpha) + 
lmiterm([46 1 1 0],d_max*exp(T*alpha));    % d* exp(T*alpha) + 
lmiterm([46 1 1 0],y_l*exp(T*alpha));     % \eta(0)* exp(T*alpha)-c2 < 0
lmiterm([46 1 1 c2],-1,1);  

lmiterm([47 1 1 P11],-1,1);     % (19) 1st:  C1'*C1 <= gamma^2 * P11,  
lmiterm([47 1 2 0],C1');
lmiterm([47 2 2 gamma],-1,1);
lmiterm([48 1 1 P11],1,1);       % (19) 2nd:  P11 - E01 - phi*I <=0, 
lmiterm([-48 1 1 0],E01);
lmiterm([48 1 1 0],phi*eye(n_x));
lmiterm([49 1 1 P21],1,1);       % (19) 2nd:  P21 - E02 - phi*I <=0, 
lmiterm([-49 1 1 0],E02);
lmiterm([49 1 1 0],phi*eye(n_x));

lmiterm([50 1 1 P12],-1,1);     % (19) 1st:  C2'*C2 <= gamma^2 * P12,  
lmiterm([50 1 2 0],C2');
lmiterm([50 2 2 gamma],-1,1);
lmiterm([52 1 1 P12],1,1);       % (19) 2nd:  P12 - E01 - phi*I <=0, 
lmiterm([-52 1 1 0],E01);
lmiterm([52 1 1 0],phi*eye(n_x));
lmiterm([53 1 1 P22],1,1);       % (19) 2nd:  P22 - E02 - phi*I <=0, 
lmiterm([-53 1 1 0],E02);
lmiterm([53 1 1 0],phi*eye(n_x));

lmiterm([55 1 1 P13],-1,1);     % (19) 1st:  C3'*C3 <= gamma^2 * P13,  
lmiterm([55 1 2 0],C3');
lmiterm([55 2 2 gamma],-1,1);
lmiterm([56 1 1 P13],1,1);       % (19) 2nd:  P13 - E01 - phi*I <=0, 
lmiterm([-56 1 1 0],E01);
lmiterm([56 1 1 0],phi*eye(n_x));
lmiterm([58 1 1 P23],1,1);       % (19) 2nd:  P23 - E02 - phi*I <=0, 
lmiterm([-58 1 1 0],E02);
lmiterm([58 1 1 0],phi*eye(n_x));

lmiterm([59 1 1 P14],-1,1);     % (19) 1st:  C4'*C4 <= gamma^2 * P13,  
lmiterm([59 1 2 0],C4');
lmiterm([59 2 2 gamma],-1,1);
lmiterm([60 1 1 P14],1,1);       % (19) 2nd:  P14 - E01 - phi*I <=0, 
lmiterm([-60 1 1 0],E01);
lmiterm([60 1 1 0],phi*eye(n_x));
lmiterm([61 1 1 P24],1,1);       % (19) 2nd:  P24 - E02 - phi*I <=0, 
lmiterm([-61 1 1 0],E02);
lmiterm([61 1 1 0],phi*eye(n_x));

lmiterm([62 1 1 N11],1,1);      
lmiterm([-63 1 1 N11],1,1);    
lmiterm([63 1 1 0],-1);

lmiterm([64 1 1 N13],1,1);      
lmiterm([-65 1 1 N13],1,1);      
lmiterm([65 1 1 0],-1.5);

lmiterm([66 1 1 N15],1,1);
lmiterm([-67 1 1 N15],1,1);
lmiterm([67 1 1 0],-2);

lmiterm([68 1 1 N21],1,1);
lmiterm([-69 1 1 N21],1,1);
lmiterm([69 1 1 0],-1);

lmiterm([70 1 1 N23],1,1);
lmiterm([-71 1 1 N23],1,1);
lmiterm([71 1 1 0],-1.5);

lmiterm([72 1 1 N25],1,1);
lmiterm([-73 1 1 N25],1,1);
lmiterm([73 1 1 0],-2);

lmiterm([74 1 1 N31],1,1);
lmiterm([-75 1 1 N31],1,1);
lmiterm([75 1 1 0],-1);

lmiterm([76 1 1 N33],1,1);
lmiterm([-77 1 1 N33],1,1);
lmiterm([77 1 1 0],-1.5);

lmiterm([78 1 1 N35],1,1);
lmiterm([-79 1 1 N35],1,1);
lmiterm([79 1 1 0],-2);

lmiterm([80 1 1 N41],1,1);
lmiterm([-81 1 1 N41],1,1);
lmiterm([81 1 1 0],-1);

lmiterm([82 1 1 N43],1,1);
lmiterm([-83 1 1 N43],1,1);
lmiterm([83 1 1 0],-1.5);

lmiterm([84 1 1 N45],1,1);
lmiterm([-85 1 1 N45],1,1);
lmiterm([85 1 1 0],-2);

LMIs=getlmis;
n=decnbr(LMIs);
c=zeros(n,1);
for k=1:n
[ci,gammai]=defcx(LMIs,k,c2,gamma);
c(k)=ci+gammai;     % optimization objective: c2^2 + gamma^2
end

options=[0,10,0,0,0];
[copt,xopt]=mincx(LMIs,c,options);


c2=sqrt(dec2mat(LMIs,xopt,c2))
gamma=sqrt(dec2mat(LMIs,xopt,gamma))

P11=dec2mat(LMIs,xopt,P11);
P12=dec2mat(LMIs,xopt,P12);
P13=dec2mat(LMIs,xopt,P13);
P14=dec2mat(LMIs,xopt,P14);

P21=dec2mat(LMIs,xopt,P21);
P22=dec2mat(LMIs,xopt,P22);
P23=dec2mat(LMIs,xopt,P23);
P24=dec2mat(LMIs,xopt,P24);

beta=sqrt(dec2mat(LMIs,xopt,beta));

V1=dec2mat(LMIs,xopt,V1);
V2=dec2mat(LMIs,xopt,V2);
V3=dec2mat(LMIs,xopt,V3);
V4=dec2mat(LMIs,xopt,V4);

W11=dec2mat(LMIs,xopt,W11);
W12=dec2mat(LMIs,xopt,W12);
W13=dec2mat(LMIs,xopt,W13);
W14=dec2mat(LMIs,xopt,W14);

W21=dec2mat(LMIs,xopt,W21);
W22=dec2mat(LMIs,xopt,W22);
W23=dec2mat(LMIs,xopt,W23);
W24=dec2mat(LMIs,xopt,W24);

M1=dec2mat(LMIs,xopt,M1);
M2=dec2mat(LMIs,xopt,M2);
M3=dec2mat(LMIs,xopt,M3);
M4=dec2mat(LMIs,xopt,M4);

N1=dec2mat(LMIs,xopt,N1);
N2=dec2mat(LMIs,xopt,N2);
N3=dec2mat(LMIs,xopt,N3);
N4=dec2mat(LMIs,xopt,N4);

K11=inv(P21)*M1;   
K12=inv(P22)*M2;    
K13=inv(P23)*M3;   
K14=inv(P24)*M4;

K21=inv(P21)*N1;     
K22=inv(P22)*N2;  
K23=inv(P23)*N3;  
K24=inv(P24)*N4;

K31=inv(V1)*W21;      
K32=inv(V1)*W22;  
K33=inv(V1)*W23;  
K34=inv(V1)*W24;

K41=inv(V1)*W11;      
K42=inv(V1)*W12;   
K43=inv(V1)*W13;   
K44=inv(V1)*W14;    