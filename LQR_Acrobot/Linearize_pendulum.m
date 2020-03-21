init_point = [pi,0,0,0];

tspan = 0:0.001:10;
x0 = [0,0,0.0,0];
% x0 = init_point;
%% Linearizing the Model
syms q1 q2 q1dot q2dot 

% Parameters
m1 = 8;
m2 = 8;
l1 = 0.5;
l2 = 1;
g  = 9.81;
I1 = m1*l1^2;
I2 = m2*l2^2;


%%
M = [I1+I2+m2*l1^2 + 2*m2*l1*l2*cos(q2),  I2+ m2*l1*l2*cos(q2);
           I2 + m2*l1*l2*cos(q2),                      I2];
       

C = [-2*m2*l1*l2*sin(q2)*q2dot,  -m2*l1*l2*sin(q2)*q2dot;
     m2*l1*l2*sin(q2)*q1dot,              0];
 
T = [-m1*g*l1*sin(q1)-m2*g*(l1*sin(q1) + l2*sin(q1+q2));
            -m2*g*l2*sin(q1+q2)];
        
f2 = inv(M)*{T - C*[q1dot;q2dot]};

f1 = [q1dot;q2dot];

f = [f1;f2];
f0 = double(subs(f,[q1,q2,q1dot,q2dot], init_point));

%%
diff_f = jacobian(f,[q1,q2,q1dot,q2dot]);
A = double(subs(diff_f,[q1,q2,q1dot,q2dot], init_point));


%%
del_f_del_u = [0;
     0;
     inv(M)*[0;1]];
B = double(subs(del_f_del_u, [q1,q2,q1dot,q2dot], init_point));


% x


%% ODE
[t,state] = ode45(@(t,x) dynamics(t,x,A,B,f0,init_point),tspan,x0);
% state(:,1:2) = state(:,1:2)/pi;
x1 = l1*sin(state(:,1));
y1 = -l1*cos(state(:,1));

x2 = x1 + l2*sin(state(:,1) + state(:,2));
y2 = y1 - l2*cos(state(:,1) + state(:,2));

plot(t,state);
legend('q1','q2','q1dot','q2dot');

figure(2)
hold on;
plot(x1,y1,'-*','MarkerIndices',1:100:length(y1));
plot(x1(1),y1(1),'-p','MarkerFaceColor','green',...
    'MarkerSize',15);
plot(x1(end),y1(end),'-p','MarkerFaceColor','red',...
    'MarkerSize',15);

figure(3)
hold on;
plot(x2,y2);
plot(x2(1),y2(1),'-p','MarkerFaceColor','green',...
    'MarkerSize',15);
plot(x2(end),y2(end),'-p','MarkerFaceColor','red',...
    'MarkerSize',15);

%% Animation
% Y = state;
% T = t;
% L = l1+l2;
% x = [ l1*sin(Y(:,1)),  l1*sin(Y(:,1))+l2*sin(Y(:,2)+Y(:,1))];
% y = [-l1*cos(Y(:,1)), -l1*cos(Y(:,1))-l2*cos(Y(:,2)+Y(:,1))];
% 
% % Convert radians to degrees
% ang = Y(:,1:2)*180/pi;
% 
% figure;
% subplot(2,1,1);
% xlabel('time (sec)'); ylabel('angle (\circ)');
% 
% tic;    % start timing
% for id = 1:length(T)
%    % The top plot shows a time series of link angles
%    subplot(2,1,1);
%    plot(T,ang, 'LineWidth', 2);
%    line(T(id), ang(id,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
%    line(T(id), ang(id,2), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
%    xlabel('time (sec)'); ylabel('angle (deg)');
% 
%    % The bottom plot shows the animation of the double pendulum
%    subplot(2,1,2);
%    plot([0, x(id,1);x(id,1), x(id,2)], [0, y(id,1);y(id,1), y(id,2)], ...
%       '.-', 'MarkerSize', 20, 'LineWidth', 2);
%    axis equal; axis([-2*L 2*L -2*L 2*L]);
%    title(sprintf('Time: %0.2f sec', T(id)));
% 
%    drawnow;
% end
% fprintf('Animation (Regular): %0.2f sec\n', toc);



function dxdt = dynamics(t,x,A,B,f0,x0)
Q = 100*eye(4);
R = 100;

[X,K,L] = icare(A,B,Q,R);
u = -K*(x-x0');
% u = -K*x;
% u = 0;
dxdt = A*(x-x0') + B*u +f0;

end

        
