function [x y] = soluzione_radiale(E,l,Xmin,Xmax)
    
    e = -1.602e-19;
    eps0 = 8.854e-12;
    m = 9.109e-31;
    ht = 1.054e-34;%6.582e-16; %eV*s % 1.054e-34; %J*s
    
    syms y(t)
    [V] = odeToVectorField(diff(y, 2) == -((2*m /(ht)^2) * E + 2*m*e^2 /(4*pi*eps0*ht^2*t) - l*(l+1)/(t^2) )*y);
    M = matlabFunction(V,'vars', {'t','Y'});
    sol = ode45(M,[Xmin Xmax],[0; 1]);

    x = linspace(Xmin,Xmax,200);
    y = deval(sol,x);

    %figure();
    %plot(x,(y(1,:)).^2);
    %hold on
    %plot(x,zeros(200,1));
%plot(t,y(:,1));

%sol = ode45(M,[0 1e-10],[0 1e-4]);
%fplot(@(x)deval(sol,x,1), [0, 1e-10]);

end