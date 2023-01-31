addpath(genpath('test'));
EPS = 1e-6;

% point
P = [5,5,100];

% surface function
func = @(x,y) x.^2 + y.^2;
syms x y;
func_x = matlabFunction(diff(func,x),'Vars',{x,y});
func_y = matlabFunction(diff(func,y),'Vars',{x,y});

% solve the non-linear eqs
Q = fsolve(@(Q) eqs2solve(Q,P,func),P(1:2));
Q(3) = func(Q(1),Q(2));
dist = norm(P - Q);

% plot the result
[ptsMeshX,ptsMeshY] = meshgrid(linspace(-10,10),linspace(-10,10));
ptsMeshZ = func(ptsMeshX,ptsMeshY);
figure;
surf(ptsMeshX,ptsMeshY,ptsMeshZ,'FaceAlpha',0.3,'EdgeColor','none');
hold on;
scatter3(P(1),P(2),P(3));
scatter3(Q(1),Q(2),Q(3));
quiver3(Q(1),Q(2),Q(3),func_x(Q(1),Q(2)),func_y(Q(1),Q(2)),-1);
line([P(1),Q(1)],[P(2),Q(2)],[P(3),Q(3)]);
xlabel('x'); 
ylabel('y'); 
zlabel('z');



function F = eqs2solve(Q,P,func)
syms x y;
func_x = matlabFunction(diff(func,x),'Vars',{x,y});
func_y = matlabFunction(diff(func,y),'Vars',{x,y});
F(1) = (P(1) - Q(1)) + func_x(Q(1),Q(2)) .* (P(3) - func(Q(1),Q(2)));
F(2) = (P(2) - Q(2)) + func_y(Q(1),Q(2)) .* (P(3) - func(Q(1),Q(2)));
end