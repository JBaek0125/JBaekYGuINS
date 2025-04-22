%  case3_Test_static.m      
%  \tilde{D}--> sigma=0.01;  with  \eta(t)=0;  ( static ETS ) 

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

lambda=0.9;  

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
        xs(:,1)= [ -1;  0.5;  0.8];   % initial state x(0)=x(:,1), 
        d_1s(:,1)= 0.05 * sin(1*h);   % disturbance, 
        ys(:,1) = C * xs(:,1);     % y(0)=C*x(0),
        y_es(:,1)= 0 * ys(:,1);    % y_e(:,n) = y(t_k),  
        us(:,1) = K * y_es(:,1);
%         kd=rand(1,tm/h+1);
%         kc=rand(1,tm/h+1);
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
            d_1s(:,n)=0.05*sin(t1);
        else
            d_1s(:,n)=0.00*sin(t1);
        end
        
       %% Plant
        xs(:,n) = h*( A * xs(:,n-1) + B * us(:,n-1) + D * d_1s(:,n-1)) + xs(:,n-1);
        ys(:,n) = C * xs(:,n);
        
       %% ET  
        eys = y_es(:,n-1) - ys(:,n);   % Triggering error   % eta(t)= y_l(n)= 0,
        % y_l(n)= h * ( -lambda* y_l(n-1) +  sigma * y(:,n)' * y(:,n) - ey' * ey ) + y_l(n-1);
        
       if eys' * eys >= sigmas*ys(:,n)'* ys(:,n)     % eta(t)= 0,
        y_es(:,n)    = ys(:,n);
        Count3       = Count3 + 1;
        intervals3(n)= (n - r1) * h;
        r1          = n;
       else
        y_es(:,n)    = y_es(:,n-1);
        intervals3(n)= -1;
       end 
       
       g(:,n)= 0.7 * y_es(:,n);
       
       y_a(:,n)=((kd(n)<=E_kd)*0*y_es(:,n)) + ( (1-(kd(n)<=E_kd))*( (kc(n)<=E_kc)*g(:,n) + (1-(kc(n)<=E_kc) )*y_es(:,n) ) );
       
       
       %% Controller
        us(:,n) = K * y_a(:,n);
    end
    n=n+1;
end


% figure(1)
% subplot(2,1,1);
% stem(t,intervals,'b-','linewidth',1)
% axis([0 tm 0 0.9]) 
% xlabel('Time(s)');
% ylabel('Release intervals');
% mlstr = {'Static ETS'};
% legend(mlstr,'interpreter','latex')
% subplot(2,1,2);
% plot(t,xs(1,:),'b-','linewidth',1.5)
% hold on
% plot(t,xs(2,:),'r--','linewidth',1.5)
% hold on
% plot(t,xs(3,:),'g-','linewidth',1.5)
% hold off
% axis([0 tm -1 1])
% xlabel('Time(s)');
% ylabel('System state');
% legend('x_1(t)','x_2(t)','x_3(t)')
