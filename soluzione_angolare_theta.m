function [x y] = soluzione_angolare_theta(l,m,Xmin,Xmax)
    
    syms y(t)
    [V] = odeToVectorField(diff(y, 2) == ( m^2 / ((sin(t))^2) - l*(l+1))*y - cos(t)/sin(t) * diff(y) );
    M = matlabFunction(V,'vars', {'t','Y'});
    sol = ode45(M,[Xmin Xmax],[0; 1]);

    x = linspace(Xmin,Xmax,200);
    y = deval(sol,x);

    figure();
    plot(x,(y(1,:)).^2);
    hold on
    plot(x,zeros(200,1));
    %plot(t,y(:,1));

    %sol = ode45(M,[0 1e-10],[0 1e-4]);
    %fplot(@(x)deval(sol,x,1), [0, 1e-10]);

end