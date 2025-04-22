figure(1)
plot(t,x1(1,:),'b-','linewidth',1.5)
hold on
plot(t,x2(1,:),'r--','linewidth',1.5)
hold on
plot(t,x3(1,:),'g-.','linewidth',1.5)
hold off
% axis([0 tm -1 1])
xlabel('Time(s)');
ylabel('System state x_1(t)');
legend('Case 1','Case 2','Case 3')

figure(2)
plot(t,x1(2,:),'b-','linewidth',1.5)
hold on
plot(t,x2(2,:),'r--','linewidth',1.5)
hold on
plot(t,x3(2,:),'g-.','linewidth',1.5)
hold off
% axis([0 tm -1 1])
xlabel('Time(s)');
ylabel('System state x_2(t)');
legend('Case 1','Case 2','Case 3')

figure(3)
plot(t,x1(3,:),'b-','linewidth',1.5)
hold on
plot(t,x3(3,:),'r--','linewidth',1.5)
hold on
plot(t,x2(3,:),'g-.','linewidth',1.5)
hold off
% axis([0 tm -1 1])
xlabel('Time(s)');
ylabel('System state x_3(t)');
legend('Case 1','Case 2','Case 3')

figure(4)
stem(t,intervals1,'b-','linewidth',1)
hold on
stem(t,intervals2,'r-','linewidth',1)
hold on
stem(t,intervals3,'g-','linewidth',1)
hold off
axis([0 tm 0 1.5]) 
xlabel('Time(s)');
ylabel('Release intervals');
legend('Case 1','Case 2','Case 3')


