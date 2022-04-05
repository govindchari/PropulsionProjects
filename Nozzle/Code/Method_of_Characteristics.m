%Govind Chari (gmc93)
%December 31, 2019


%This script creates a minimum length supersonic nozzle, given a target
%exit Mach number. The nozzle directs flow to be perfectly axial and
%cancels all expansion shockwaves

%ISSUE: THE INVERSE PM ONLY WORKS WITH GAMMA=1.4. CHANGE IT.
clc; close all; clear; tic;

%Input Block
Me = input('Exit Mach number: ');
gamma = input('Ratio of specific heats: ');
Rt = input('Throat radius: ');
n = input('Number of characteristics lines: ');

%Calculates the maximum flow deflection angle in degrees
theta_max=rad2deg(0.5*PM(Me,gamma));
nodes=0.5*n*(n+1)+n;
%Initializes a matrix of (x,y) points for the wall, including the throat point
%and sets the first point to (0,Rt)
points=zeros(n+1,2); 
points(1,1)=0; 
points(1,2)=Rt;

%Creates x and y arrays for the x and y coordinates of the nodes
x=zeros(1,nodes);
y=zeros(1,nodes);

%Initializes the arrays that store the node numbers for the centerline nodes and
%wall nodes
centerline_nodes=zeros(1,n);
centerline_nodes(1)=1;
wall_nodes=zeros(1,n);
wall_nodes(1)=n+1;

[v,CL,CR,theta,M,mu]=charSolver(theta_max,n);

%Calculates which nodes are on the centerline
for i=2:n
    centerline_nodes(i)=centerline_nodes(i-1)+n+3-i;
end
%Calculates which nodes are on the wall
for i=2:n
    wall_nodes(i)=wall_nodes(i-1)+n+2-i;
end
%Calculates the coordinates of every node
for i=1:nodes
    
    
    %Calculates coordinate of node 1
    if i==1
        x(i)=-Rt/(tand(theta(i)-mu(i)));
        y(i)=0;
        plot([0 x(i)],[Rt 0],'Color',[0,0.7,0.9]);
        hold on
    
    %Calculates coodinates of interior nodes of first characteristic line
    elseif i<=n 
     syms a b
     eq1= b == a*tand(theta(i)-mu(i))+Rt;
     eq2 = b == (a-x(i-1))*tand(theta(i-1)+mu(i-1))+y(i-1);
     [x(i),y(i)] = vpasolve([eq1,eq2],[a,b]);
     plot([x(i-1) x(i)],[y(i-1) y(i)],'Color',[0,0.7,0.9]);
     hold on

    %Calculates coordinate of first wall point
    elseif i==n+1
        syms a b
        eq1= b == (a-x(i-1))*tand(theta(i)+mu(i))+y(i-1);
        eq2= b == a*tand(theta_max)+Rt;
        [x(i),y(i)] = vpasolve([eq1,eq2],[a,b]);
        points(2,1)=x(i);
        points(2,2)=y(i);
        plot([x(i-1) x(i)],[y(i-1) y(i)],'Color',[0,0.7,0.9]);
        plot([0 x(i)],[Rt y(i)],'Color',[0,0,0]);
        hold on
    
    % Calculates coordinates of centerline nodes    
    elseif(ismember(i,centerline_nodes))
        z=find(centerline_nodes==i);
        x(i)=x(i-n+z-2)-(y(i-n+z-2)/tand(theta(i-n+z-2)-mu(i-n+z-2)));
        y(i)=0;
        plot([0 x(i)],[Rt 0],'Color',[0,0.7,0.9]);
        hold on
            
    %Calculates coordinates of wall nodes
    elseif (ismember(i, wall_nodes))
        z=find(wall_nodes==i);
        prev_wall=wall_nodes(z-1);
        syms a b
        eq1= b == (a-x(prev_wall))*tand(theta(prev_wall))+y(prev_wall);
        eq2 = b == (a-x(i-1))*tand(theta(i-1)+mu(i-1))+y(i-1);
        [x(i),y(i)] = vpasolve([eq1,eq2],[a,b]);
        points(z+1,1)=x(i);
        points(z+1,2)=y(i);
        plot([x(i-1) x(i)],[y(i-1) y(i)],'Color',[0,0.7,0.9]);
        plot([x(prev_wall) x(i)],[y(prev_wall) y(i)],'Color',[0,0,0]);
        hold on
        
    %Calculates the coordinates of all other nodes
    else
        for j=1:n
            if(ismember(i-j,centerline_nodes))
                z=find(centerline_nodes==i-j);
                break;
            end
        end
        syms a b
        eq1= b == (a-x(i-n+z-2))*tand(theta(i-n+z-2)-mu(i-n+z-2))+y(i-n+z-2);
        eq2 = b == (a-x(i-1))*tand(theta(i-1)+mu(i-1))+y(i-1);
        [x(i),y(i)] = vpasolve([eq1,eq2],[a,b]);
        plot([x(i-1) x(i)],[y(i-1) y(i)],'Color',[0,0.7,0.9]);
        hold on
    end
end

%Formats the graph
title('MOC Nozzle Contour')
xlabel('Length (in)')
ylabel('Radius (in)')
axis equal
%xlim([0 x(nodes)]+10)
%ylim([0 y(nodes)]+10)

%Export
saveas(gcf, 'MOC_contour.pdf'); % Exports as PDF
filename='moc_data.xlsx'; % Defines Excel file name
xlswrite(filename, points); % Writes points() to the specified file
toc % Prints time elapsed during code execution

%Given a Mach number and ratio of specific heats, calculates the
%Prandtl-Meyer angle in degrees
function v=PM(M,gamma)
    v=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*((M^2-1))))-atan(sqrt(M^2-1));
end

%Given a Prandtl-Meyer angle in degrees, calculates the
%Mach number, this only works for gamma=1.4
function mach = invPM(v)
    A = 1.3604; 
    B = 0.0962;
    C = -0.5127;
    D = -0.6722;
    E = -0.3278;
    v = deg2rad(v);
    v_0 = (pi/2)*(sqrt(6)-1);
    y = (v/v_0)^(2/3);
    mach = (1 + A*y + B*(y^2) + C*(y^3))/(1 + D*y + E*(y^2));
end

%Given a Mach number, calculates the Mach angle in degrees
function mu=MachAngle(M)
    mu=asind(1/M);
end

%Solves for Prandtl-Meyer angle, left and right
%running characteristic constants, flow deflection angle, Mach number, and
%Mach angle
function [v,CL,CR,theta, M, mu]=charSolver(theta_max,n)
    nodes=0.5*n*(n+1)+n;
    v=zeros(1,nodes);
    CL=zeros(1,nodes);
    CR=zeros(1,nodes);
    theta=zeros(1,nodes);
    M=zeros(1,nodes);
    mu=zeros(1,nodes);
    centerline_nodes=zeros(1,n);
    centerline_nodes(1)=1;
    wall_nodes=zeros(1,n);
    wall_nodes(1)=n+1;
    delta_theta=theta_max/n;
    
    %Calculates which nodes are on the centerline
    for i=2:n
        centerline_nodes(i)=centerline_nodes(i-1)+n+3-i;
    end
    %Calculates which nodes are on the wall
    for i=2:n
        wall_nodes(i)=wall_nodes(i-1)+n+2-i;
    end
    %Calculates the characteristics at each node
    for i=1:nodes
        %Calculates the characteristics at each node of the first wave and
        %its reflection
        if i>=1 && i<=n
            theta(i)=i*delta_theta;
            v(i)=theta(i);
            CR(i)=2*theta(i);
            CL(i)=0;
            M(i)=invPM(v(i));
            mu(i)=MachAngle(M(i));
        %Calculates characteristics at all nodes on centerline
        elseif ismember(i,centerline_nodes) && i>2
            theta(i)=0;
            CR(i)=CR(i-(n+2-find(centerline_nodes==i)));
            v(i)=CR(i);
            CL(i)=-v(i);
            M(i)=invPM(v(i));
            mu(i)=MachAngle(M(i));
        
        %Calculates characteristics at all nodes on wall
        elseif ismember(i,wall_nodes)
            v(i)=v(i-1);
            CL(i)=CL(i-1);
            CR(i)=CR(i-1);
            theta(i)=theta(i-1);
            M(i)=M(i-1);
            mu(i)=mu(i-1);
        %Calculates characteristics at the rest of the nodes
        else
            for j=1:n
                if ismember((i-j),centerline_nodes)
                    x=find(centerline_nodes==i-j);
                    break
                end
            end
            CR(i)=CR(i-n+x-2);
            CL(i)=CL(i-1);
            theta(i)=0.5*(CR(i)+CL(i));
            v(i)=0.5*(CR(i)-CL(i));
            M(i)=invPM(v(i));
            mu(i)=MachAngle(M(i));
        end
        
    end
end
