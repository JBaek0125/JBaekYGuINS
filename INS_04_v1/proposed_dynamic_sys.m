%  Case2_Test_dynamic.m        ( Case 2 )
%  \tilde{D}-->Sigma=[0.01  0;  0  0.31];  % Threshold matrix

tm=30;  h=0.01;  t=0:h:(tm-h);   % running time, sampling interval,
n=1;   lambda=0.9;

PI=[-2.4    2.2    0.1    0.1;    % transition rate matrix (TRM) 
     0.1   -1.9    1.7    0.1; 
     1.2    0.3   -1.8    0.3;
     2.6    0.1    0.1   -2.8];  

nsim=1;   % See the function simCTMC.m
instt=1;
[JumpTime,JumpList] = simCTMC(PI,tm,nsim,instt); % Generate jump times and corresponding modes
JumpTime=[JumpTime,tm];
JumpList=[JumpList,JumpList(end)];
k=1;
% The flags of triggering times and triggering interval
Count2=0;r1=1;

while n<=(tm/h)
    if n==1
       %% Initialization
        switch JumpList(1)
            case 1
                A=A1;B=B1;C=C1;D=D1;K1=K11;K2=K21;K3=K31;K4=K41;
            case 2
                A=A2;B=B2;C=C2;D=D2;K1=K12;K2=K22;K3=K32;K4=K42;
            case 3
                A=A3;B=B3;C=C3;D=D3;K1=K13;K2=K23;K3=K33;K4=K43;
            case 4
                A=A4;B=B4;C=C4;D=D4;K1=K14;K2=K24;K3=K34;K4=K44;
        end 
        x2(:,1)= [ -1; 0.5; 0.8 ];   % initial state x(0)=x(:,1), 
        xc2(:,1)= [ 0; 0; 0 ];   % initial controller state xc(0)=xc(:,1), 
        d_12(:,1)= 0.05 * sin(1*h);   % disturbance, 
        y2(:,1)= C * x2(:,1);     % y(0)=C*x(0),
        y_e2(:,1)= 0 * y2(:,1);    % y_e(:,n) = y(t_k),  
        y_a2(:,1) = 0 * y_e2(:,1); 
        y_l2(1)= phi*x2(:,1)'*x2(:,1);     % eta(0) = phi*x(0)'*x(0), 
        u2(:,1)= K3 * xc2(:,1) + K4 * y_a2(:,1);
%         kd=rand(1,tm/h+1);
%         kc=rand(1,tm/h+1);
        intervals2(1) = -1;   % triggering interval,         
    else
        t1=h*(n-1);
        t2=h*n;
       %% Mode
        if n*h>JumpTime(k)   % Switching System Modes
            flag=JumpList(k);
        switch flag
            case 1
                A=A1;B=B1;C=C1;D=D1;K1=K11;K2=K21;K3=K31;K4=K41;
            case 2
                A=A2;B=B2;C=C2;D=D2;K1=K12;K2=K22;K3=K32;K4=K42;
            case 3
%               A=A1;B=B1;C=C1;D=D1;K1=K11;K2=K21;K3=K31;K4=K41;
                 A=A3;B=B3;C=C3;D=D3;K1=K13;K2=K23;K3=K33;K4=K43;
            case 4
%               A=A2;B=B2;C=C2;D=D2;K1=K12;K2=K22;K3=K32;K4=K42;
                 A=A4;B=B4;C=C4;D=D4;K1=K14;K2=K24;K3=K34;K4=K44;
        end
            k=k+1;
        end 
        
       %% Disturbance,  d_1(:,1) = w(t),
        if t1>=12 && t1<=18
            d_12(:,n)=0.05*sin(t1);
        else
            d_12(:,n)=0.00*sin(t1);
        end
        
       %% Plant
        x2(:,n) = h*( A * x2(:,n-1) + B * u2(:,n-1) + D * d_12(:,n-1)) + x2(:,n-1);
        y2(:,n) = C * x2(:,n);
        
       %% Controller system
       xc2(:,n) = h*( K1 * xc2(:,n-1) + K2 * y_a2(:,n-1) ) + xc2(:,n-1);
       
       %% ET  
        ey2 = y_e2(:,n-1) - y2(:,n);   % the triggering error
        y_l2(n)= h * ( -lambda* y_l2(n-1) + y2(:,n)'*Sigma * y2(:,n) - ey2' * ey2 ) + y_l2(n-1);
        
       if ey2' * ey2 >= y2(:,n)'*Sigma * y2(:,n) + y_l2(n)
        y_e2(:,n)    = y2(:,n);
        Count2       = Count2 + 1;
        intervals2(n)= (n - r1) * h;
        r1          = n;
       else
        y_e2(:,n)    = y_e2(:,n-1);
        intervals2(n)= -1;
       end 
       
       g(:,n)= 0.7 * y_e2(:,n);
       
       y_a2(:,n)=((kd(n)<=E_kd)*0*y_e2(:,n)) + ( (1-(kd(n)<=E_kd))*( (kc(n)<=E_kc)*g(:,n)+(1-(kc(n)<=E_kc))*y_e2(:,n) ) );
       
       
       %% Controller
        u2(:,n) = K3*xc2(:,n) + K4*y_a2(:,n);
    end
    n=n+1;
end

alpha1=(kd<=E_kd);
alpha2=(kc<=E_kc);


% figure(1)
% stem(t,intervals,'b-','linewidth',1)
% axis([0 tm 0 1.5]) 
% xlabel('Time(s)');
% ylabel('Release intervals');
% mlstr = {'The proposed general dynamic ETS'};
% legend(mlstr,'interpreter','latex')
% legend(mlstr,'interpreter','latex')
% 
% 
% figure(2)
% plot(t,x2(1,:),'b-','linewidth',1.5)
% hold on
% plot(t,x2(2,:),'r--','linewidth',1.5)
% hold on
% plot(t,x2(3,:),'g-','linewidth',1.5)
% hold off
% % axis([0 tm -1 1])
% xlabel('Time(s)');
% ylabel('System state');
% legend('x_1(t)','x_2(t)','x_3(t)')
