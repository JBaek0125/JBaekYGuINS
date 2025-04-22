%  Test_static.m     ( Solving  THM1_static.m )    (Case 1)

tm=30;  h=0.01;  t=0:h:(tm-h);  % running time, sampling interval,
n=1;  

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
Count3=0;r1=1;

lambda1=0.9;  

while n<=(tm/h)
    if n==1
       %% Initialization
        switch JumpList(1)
            case 1
                A=A1;B=B1;C=C1;D=D1;K=K1;
            case 2
                A=A2;B=B2;C=C2;D=D2;K=K2;
            case 3
                A=A3;B=B3;C=C3;D=D3;K=K3; 
            case 4
                A=A4;B=B4;C=C4;D=D4;K=K4;
        end 
        x3(:,1)= [ -1;  0.5;  0.8];   % initial state x(0)=x(:,1), 
        d_13(:,1)= 0.05 * sin(1*h);   % disturbance, 
        y3(:,1) = C * x3(:,1);     % y(0)=C*x(0),
        y_e3(:,1)= 0 * y3(:,1);    % y_e(:,n) = y(t_k),  
        y_l3(1) = phi1*x3(:,1)'*x3(:,1);     % eta(0) = phi1*x(0)'*x(0), 
        u3(:,1) = K * y_e3(:,1);
%          kd=rand(1,tm/h+1);
%          kc=rand(1,tm/h+1);
        intervals3(1) = -1;   % triggering interval,         
    else
        t1=h*(n-1);
        t2=h*n;
       %% Mode
        if n*h>JumpTime(k)   % Switching System Modes
            flag=JumpList(k);
        switch flag
            case 1
                A=A1;B=B1;C=C1;D=D1;K=K1;
            case 2
                A=A2;B=B2;C=C2;D=D2;K=K2;
            case 3
                A=A3;B=B3;C=C3;D=D3;K=K3; 
            case 4
                A=A4;B=B4;C=C4;D=D4;K=K4;
        end
            k=k+1;
        end 
        
       %% Disturbance,  d_1(:,1) = w(t),
        if t1>=12 && t1<=18
            d_13(:,n)=0.05*sin(t1);
        else
            d_13(:,n)=0.00*sin(t1);
        end
        
       %% Plant
        x3(:,n) = h*( A * x3(:,n-1) + B * u3(:,n-1) + D * d_13(:,n-1)) + x3(:,n-1);
        y3(:,n) = C * x3(:,n);
        
       %% ET  
        ey3 = y_e3(:,n-1) - y3(:,n);   % the triggering error
        y_l3(n)= h * ( -lambda1* y_l3(n-1) +  sigma1 * y3(:,n)' * y3(:,n) - ey3' * ey3 ) + y_l3(n-1);
        
       if ey3' * ey3 >= 0*sigma1 * y3(:,n)' * y3(:,n) + y_l3(n)
        y_e3(:,n)    = y3(:,n);
        Count3       = Count3 + 1;
        intervals3(n)= (n - r1) * h;
        r1          = n;
       else
        y_e3(:,n)    = y_e3(:,n-1);
        intervals3(n)= -1;
       end 
       
       g(:,n)= 0.7 * y_e3(:,n);
       
       y_a(:,n)=((kd(n)<=E_kd)*0*y_e3(:,n)) + (  (1-(kd(n)<=E_kd))*( (kc(n)<=E_kc)*g(:,n) + (1-(kc(n)<=E_kc) )*y_e3(:,n)) );
       
       
       %% Controller
        u3(:,n) = K * y_a(:,n);
    end
    n=n+1;
end

% 
% figure(1)
% stairs(JumpTime,JumpList,'b-','linewidth',1.0)
% hold off
% axis([0 tm 0.5 5])
% xlabel('Time(s)');
% ylabel('Modes');
% legend('Sysem modes')
% 
% 
% 
% 
% figure(2)
% subplot(2,1,1);
% stairs(t,(kd(1:3000)<=E_kd),'b-','linewidth',0.5)
% axis([0 tm -0.5 1.5])
% % xlabel('Time(s)');
% ylabel('\kappa_{d}(t)');
% mlstr = {'The attack flag $\kappa_{d}(t)$ of DoS attack'};
% legend(mlstr,'interpreter','latex')
% subplot(2,1,2);
% stairs(t,(kc(1:3000)<=E_kc),'b-','linewidth',0.5)
% axis([0 tm -0.5 1.5])
% xlabel('Time(s)');
% ylabel('\kappa_{c}(t)');
% mlstr = {'The attack flag $\kappa_{c}(t)$ of deception attack'};
% legend(mlstr,'interpreter','latex')
% 
% 
% figure(3)
% subplot(2,1,1);
% plot(t,y_ld,'b-','linewidth',1) 
% axis([0 tm 0 0.25]) 
% xlabel('Time(s)');
% ylabel('Internal dynamic function');
% mlstr = {'Internal dynamic function'};
% legend(mlstr,'interpreter','latex')
% subplot(2,1,2);
% stem(t,intervals,'b-','linewidth',1)
% axis([0 tm 0 1.5]) 
% xlabel('Time(s)');
% ylabel('Release intervals');
% mlstr = {'Dynamic ETS'};
% legend(mlstr,'interpreter','latex')
% 
% 
% figure(4)
% plot(t,xd(1,:),'b-','linewidth',1.5)
% hold on
% plot(t,xd(2,:),'r--','linewidth',1.5)
% hold on
% plot(t,xd(3,:),'g-','linewidth',1.5)
% hold off
% axis([0 tm -1 1])
% xlabel('Time(s)');
% ylabel('System state');
% legend('x_1(t)','x_2(t)','x_3(t)')
