%Govind Chari (gmc93)
%January 13, 2020


%Given input parameters, this script plots the contour of the Thrust
%Optimized Parabolic nozzle and saves the coordinates into an Excel File
clear;clc;close all
%% Input Parameters
Rt = 1.443*0.5; %Throat Radius
epsilon = 4.9712; %Expansion Ratio
bell = 0.90; %Percent Bell
theta_n = 20;
theta_e = 11;

%Combustion Parameters
p0=350;
T0=2955;
g=1.1647;

%Step Size
delta_theta = 1;
delta_t = 0.1;

%Initializing all arrays
theta=[-135:delta_theta:theta_n-90];
t=[0:delta_t:1];
x=zeros(1,length(theta)+length(t));
y=zeros(1,length(theta)+length(t));

%% Throat Entrant Section
for i=1:45/delta_theta
    x(i)=1.5*Rt*cosd(theta(i));
    y(i)=1.5*Rt*sind(theta(i))+1.5*Rt+Rt;
end

%% Throat Exit Section
for i=45/delta_theta+1:length(theta)
    x(i)=0.382*Rt*cosd(theta(i));
    y(i)=0.382*Rt*sind(theta(i))+0.382*Rt+Rt;
end

%% Solves for the control points of the Bezier Curve
Nx=0.382*Rt*cosd(theta_n-90);
Ny=0.382*Rt*sind(theta_n-90)+0.382*Rt+Rt;
Ex= bell*((sqrt(epsilon)-1)*Rt)/tand(15);
Ey= sqrt(epsilon)*Rt;
m1= tand(theta_n);
m2= tand(theta_e);
c1= Ny-m1*Nx;
c2= Ey-m2*Ex;
Qx= (c2-c1)/(m1-m2);
Qy= (m1*c2-m2*c1)/(m1-m2);

%% Calculates Coordinates of Parabolic Bell
for i=1:length(t)
    x(i+1+(theta_n-90+135)/delta_theta)=((1-t(i))^2)*Nx+2*(1-t(i))*t(i)*Qx+t(i)^2*Ex;
    y(i+1+(theta_n-90+135)/delta_theta)=((1-t(i))^2)*Ny+2*(1-t(i))*t(i)*Qy+t(i)^2*Ey;
end

%% Calculates pressure, temperature, and Mach along the nozzle
for i=1:length(x)
    A_ratio(i)=y(i)^2/Rt^2;
    if x(i)<0
        mach(i)=getSubMach(A_ratio(i),g);
    elseif A_ratio(i)==1
        mach(i)=1;
    else
        mach(i)=getSupMach(A_ratio(i),g);
    end
    temp(i)=getTempRatio(mach(i),g)*T0;
    press(i)=getPressureRatio(mach(i),g)*p0;
    
end

%% Plots Data
plot(x,y)
title('THRUST OPTIMIZED PARABOLIC NOZZLE PROFILE')
xlabel('Length (in)')
ylabel('RADIUS (in)')
axis equal

% plot(x,mach)
% title('Mach Number vs Distance From Throat')
% xlabel('Distance (in)')
% ylabel('Mach Number')

M=[x',y'];
writematrix(M,'nozzle_profile.csv') 


%% Gets Supersonic Mach from Area Ratio and gamma
function mach=getSupMach(e,g)
    h1=(g+1)/2;
    h2=(g+1)/(2*g-2);
    h3=(g-1)/2;
    fun=@(x)((h1^-h2)*(1/x)*(1+h3*x^2)^h2-e);
    if e<1.00006
        mach=fzero(fun,1.0011);  
    elseif e<1.1
        mach=fzero(fun,1.2);
    else
        mach=fzero(fun,1.65);
    end
end

%% Gets Subsonic Mach from Area Ratio and gamma
function mach=getSubMach(e,g)
    h1=(g+1)/2;
    h2=(g+1)/(2*g-2);
    h3=(g-1)/2;
    fun=@(x)((h1^-h2)*(1/x)*(1+h3*x^2)^h2-e);
    mach=fzero(fun,0.9);
end

%% Gets Pressure Ratio from Mach and gamma
function pratio=getPressureRatio(M,g)
    h1=(g-1)/2;
    h2=-g/(g-1);
    pratio=(1+h1*M^2)^h2;
end

%% Gets Temperature Ratio from Mach and gamma
function Tratio=getTempRatio(M,g)
    h1=(g-1)/2;
    Tratio=(1+h1*M^2)^-1;
end
