clear all
close all
%% Plot the radial solution of the hydrogen atom as a function of n (Principal quantum number) and l (Azimuthal quantum number)

%En = [-2.17274657657552e-18 -5.42707536104199e-19 -2.41210424395473e-19 -1.35687544168024e-19 -8.68363988024612e-20 -6.03061006404460e-20 -4.43065166473389e-20 -3.39221954345703e-20 -2.66075134277344e-20 -3.39221954345703e-20 -1.14440917968750e-20];

% constants
m = 9.109e-31; %Kg
e = -1.602e-19; %C
h= 6.626e-34; % Js
eps0 = 8.854e-12; %F/m

% set n and l
n=1;
l=0;

E = - (m * e^4)/(8 * eps0^2 * n^2 * h^2);

% set the radius max and min of the graph
Xmin = 1e-20;
Xmax = 3e-10;

%4e-10 8e-10 2e-9 3e-9 4e-9 5e-9 6.5e-9 7.5e-9 8.5e-9 9.5e-9 10e-9
%3.5e-10 8.4e-10 2e-9 3e-9

[r, R] = soluzione_radiale(E,l,Xmin,Xmax);

Y = R.^2;

% normalization 

% integration
A = integ(Y(1,:),r);

PSI = (R.*1/(sqrt(A))).^2; %normalized

RAD = (R.*1/(sqrt(A)))./r;

figure();
p = plot(r./5.3e-11,PSI(1,:),'b','linewidth',2);
st = strcat('$R(r)^2$ normalized to $\int |R(r)|^2 =1$  con n=',num2str(n),' l=',num2str(l));
set(title(st,'FontSize',20),'Interpreter','latex');
set(xlabel('r [m]','FontSize',20),'Interpreter','latex');
set(ylabel('$|R(r)|^2$','FontSize',20),'Interpreter','latex');

figure();
p = plot(r./5.3e-11,RAD(1,:).*(5.3e-11)^(3/2),'b','linewidth',2);
set(title('R(r)/r normalized','FontSize',20),'Interpreter','latex');
set(xlabel('r [m]','FontSize',20),'Interpreter','latex');
set(ylabel('$R(r)^2$','FontSize',20),'Interpreter','latex');

%% Plot the radial solution of the hydrogen atom for orbitals s,p,d,f l=0,1,2,3

%En = [-2.17274657657552e-18 -5.42707536104199e-19 -2.41210424395473e-19 -1.35687544168024e-19 -8.68363988024612e-20 -6.03061006404460e-20 -4.43065166473389e-20 -3.39221954345703e-20 -2.66075134277344e-20 -3.39221954345703e-20 -1.14440917968750e-20];

% constants
m = 9.109e-31; %Kg
e = -1.602e-19; %C
h= 6.626e-34; % Js
eps0 = 8.854e-12; %F/m

% energy spectrum 
for n=1:10
   En(n) = - (m * e^4)/(8 * eps0^2 * n^2 * h^2);
end

%%%%%%%%%%%%%%%%%%%%---- s ----%%%%%%%%%%%%%%%%%%%%

% orbital s for 4 energy levels
l=0; 

% set the radius max and min of the graph
Xmmm = [1e-20 1e-20 1e-20 1e-20];
XMMM = [3e-10 9e-10 1.5e-9 2.5e-9];


for n = 1:4
    
    Xmin = Xmmm(n);
    Xmax = XMMM(n);
    
    E = En(n);

    [r, R] = soluzione_radiale(E,l,Xmin,Xmax);

    Y = R.^2;
    
    % normalization
    
    % integration
    A = integ(Y(1,:),r);

    PSI = (R.*1/(sqrt(A))).^2; %normalized

    RAD = (R.*1/(sqrt(A)))./r;

    yyyy_prob(n,:) = PSI(1,:);
    yyyy_rad(n,:) = RAD(1,:);

    xxxx(n,:) = r(:);
end

figure();
p1 = plot(xxxx(1,:),yyyy_prob(1,:),'b-','linewidth',2);
hold on
p2 = plot(xxxx(2,:),yyyy_prob(2,:),'r--','linewidth',2);
p3 = plot(xxxx(3,:),yyyy_prob(3,:),'g:','linewidth',2);
p4 = plot(xxxx(4,:),yyyy_prob(4,:),'c-.','linewidth',2);


set(xlabel('r [m]','FontSize',20),'Interpreter','latex');
set(ylabel('$|R(r)|^2$','FontSize',20),'Interpreter','latex');
t = title('Radial probability function for orbital s','FontSize',20);
set(t,'Interpreter','latex');

l =legend([p1 p2 p3 p4],'1s','2s','3s','4s');
set(l,'Interpreter','latex','FontSize',20);

%%%%%%%%%%%%%%%%%%%%---- p ----%%%%%%%%%%%%%%%%%%%%

% orbital p for 3 energy levels
l=1; 

% set the radius max and min of the graph
Xmmm = [1e-20 1e-20 1e-20];
XMMM = [9e-10 1.5e-9 2.5e-9];


for n = 1:3
    
    Xmin = Xmmm(n);
    Xmax = XMMM(n);
    
    E = En(n+1);

    [r, R] = soluzione_radiale(E,l,Xmin,Xmax);

    Y = R.^2;
    
    % normalization
    
    % integration
    A = integ(Y(1,:),r);

    PSI = (R.*1/(sqrt(A))).^2; %normalized

    RAD = (R.*1/(sqrt(A)))./r;

    yyyy_prob(n,:) = PSI(1,:);
    yyyy_rad(n,:) = RAD(1,:);

    xxxx(n,:) = r(:);
end

figure();
p1 = plot(xxxx(1,:),yyyy_prob(1,:),'b-','linewidth',2);
hold on
p2 = plot(xxxx(2,:),yyyy_prob(2,:),'r--','linewidth',2);
p3 = plot(xxxx(3,:),yyyy_prob(3,:),'g:','linewidth',2);


set(xlabel('r [m]','FontSize',20),'Interpreter','latex');
set(ylabel('$|R(r)|^2$','FontSize',20),'Interpreter','latex');
t = title('Radial probability function for orbital p','FontSize',20);
set(t,'Interpreter','latex');

l =legend([p1 p2 p3],'2p','3p','4p');
set(l,'Interpreter','latex','FontSize',20);



%%%%%%%%%%%%%%%%%%%%---- d,f ----%%%%%%%%%%%%%%%%%%%%

% orbital d,f for 2 energy levels
l=2; 

% set the radius max and min of the graph
Xmmm = [1e-20 1e-20];
XMMM = [1.5e-9 2.5e-9];


for n = 1:2
    
    Xmin = Xmmm(n);
    Xmax = XMMM(n);
    
    E = En(n+2);

    [r, R] = soluzione_radiale(E,l,Xmin,Xmax);

    Y = R.^2;
    
    % normalization
    
    % integration
    A = integ(Y(1,:),r);

    PSI = (R.*1/(sqrt(A))).^2; %normalized

    RAD = (R.*1/(sqrt(A)))./r;

    yyyy_prob(n,:) = PSI(1,:);
    yyyy_rad(n,:) = RAD(1,:);

    xxxx(n,:) = r(:);
end

n=4;
l=3; % orbital 4f

E = En(n);

% set the radius max and min of the graph
Xmin = 1e-20;
Xmax = 2.4e-9;

[r, R] = soluzione_radiale(E,l,Xmin,Xmax);

Y = R.^2;

% normalization 

% integration
A = integ(Y(1,:),r);

PSI = (R.*1/(sqrt(A))).^2; %normalized

RAD = (R.*1/(sqrt(A)))./r;

yyyy_prob(3,:) = PSI(1,:);
yyyy_rad(3,:) = RAD(1,:);

xxxx(3,:) = r(:);

figure();
p1 = plot(xxxx(1,:),yyyy_prob(1,:),'b-','linewidth',2);
hold on
p2 = plot(xxxx(2,:),yyyy_prob(2,:),'r--','linewidth',2);
p3 = plot(xxxx(3,:),yyyy_prob(3,:),'g:','linewidth',2);


set(xlabel('r [m]','FontSize',20),'Interpreter','latex');
set(ylabel('$|R(r)|^2$','FontSize',20),'Interpreter','latex');
t = title('Radial probability function for orbital d-f','FontSize',20);
set(t,'Interpreter','latex');

l =legend([p1 p2 p3],'3d','4d','4f');
set(l,'Interpreter','latex','FontSize',20);