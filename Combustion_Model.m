% Govind Chari (gmc93)
%  July 3, 2020

%This is a combustion model for BATES grains. Given grain geometry, fuel
%chemistry, and nozzle properties, a graph of thrust vs time and chamber
%pressure vs time is generated, and total impulse and burntime are printed
%out



function Govind_Combustion_Model
    hold off;clc;clear;close all;
    global T;
    global g;
    global R;
    global rho;
    global no;
    global At;
    global a;
    global n;
    
    %% Input Parameters
    %Fuel Parameters
    T=2819;     %Combustion Temperature (in K)
    g=1.21;     %Ratio of Specific Heats
    R=349.18;   %Specific Gas Constant (in J/(kg*K))
    rho=1690;   %Propellant Density (in kg/m^3)
    a=.01907;   %From St.Robert's Law
    n=0.369146; %From St.Robert's Law

    %Grain Geometry
    do=4   /39.3701;    %ID of motor casing and OD of grains (the division is to convert in to m)
    dp0=2  /39.3701;    %Initial Port diameter of grains  (the division is to convert in to m)
    L0=6   /39.3701;    %Initial Grain Length  (the division is to convert in to m)
    no=6;               %Number of grains

    %Nozzle Properties
    At=0.00036516056;    %Throat Area of Nozzle in m^2
    ep=23;           %Expansion Ratio
   
    %Simulation Runtime
    run_time=8; % Sim Runtime
    
    %% Solving for State Vectors
    
    %Calculating Initial Conditions
    p0=344738;    %Initial Pressure that gets iterated
    V0=no*L0*pi*0.5*(do^2-dp0^2);    %Initial Free Volume
    Ab0=no*(0.5*pi*(do^2-dp0^2)+pi*dp0*L0);   %Initial Burn Area

    tspan=[0 run_time];
    z0=[p0;V0;Ab0;dp0;L0];

    [t,z] = ode45(@(t,z) diffeq(t,z),tspan,z0);
    
    %Calculates Stopping Point of Simulation
    for i=1:length(z)
        if z(i,4)>= do || z(i,5)<=0
            stop=i;
            break;
        end
    end
    
    %Eliminates all data after grain has burned out
    z(stop+1,1)=0;
    z(stop+1,3)=0;
    z(stop+2:length(t),:)=[];
    t(stop+2:length(t))=[];

    z(:,1)=z(:,1)+101300;  %Converting gauge pressure to absolute pressure
    
    %Creating a matrix of the state vectors in english units
    z_eng(:,1)=z(:,1)/6894.76;
    z_eng(:,2)=z(:,2)*61023.7;
    z_eng(:,3)=z(:,3)*1550;
    z_eng(:,4)=z(:,4)*39.3701;
    z_eng(:,5)=z(:,5)*39.3701;
    
    %% Calculates Thrust and Impulse
    h1=2*g^2/(g-1);
    h2=2/(g+1);
    h3=(g+1)/(g-1);
    h4=(g-1)/g;
    mach=getExitMach(ep,g);             %Calculates Exit Mach
    pratio=getPressureRatio(mach,g);    %Calclates Pressure Ratio
    Cf=sqrt(h1*(h2^h3)*(1-pratio^h4));  %Calculates Thrust Coefficient
    Th=Cf*At*z(:,1);                    %Calculates Thrust
    Th_eng=Th/4.45;                     %Converts Thrust into lbf from N
    Kn=z(:,3)/At;                       %Calculates Klumming Number
    
    %Calculates Impulse
    impulse=0;
    for i=1:length(t)
        if i==1
        else
            impulse=impulse+Th(i-1)*(t(i)-t(i-1));
        end
    end

    %% Plotting and Output
    hold on
    plot(t,z_eng(:,1),'LineWidth',2)
    plot(t,Th_eng,'LineWidth',2)
    plot(t,Kn,'LineWidth',2)
    grid on
    
    legend('Pressure','Thrust','Klumming Number')
    xlabel('Time (s)')
    ylabel('Pressure/Thrust (psi/lbf)')
    title('Pressure and Thrust vs Time')
    
    fprintf("\nImpulse= %0.0f Ns",impulse);
    fprintf("\nBurn Time= %0.2f secs",t(length(t)));
    fprintf("\nMax Thrust= %0.2f lbf",max(Th_eng));
    fprintf("\nMax Pressure= %0.2f psi",max(z_eng(:,1)));
    
    

end

%% Diffeq to be solved
function zdot=diffeq(t,z)
    global T;
    global g;
    global R;
    global rho;
    global no;
    global At;
    global a;
    global n;

    %Extracting Information from state vector
    p=z(1);
    V=z(2);
    Ab=z(3);
    dp=z(4);
    L=z(5);
    %Calculation of Burn Rate
    r=(a*(p/6894.76)^n)/39.3701;

    %Calculation of p_dot
    h1= g/(R*T);
    h2=2/(g+1);
    h3=(g+1)/(g-1);
    x=r*Ab*rho;
    y=V/(R*T);
    z=At*p*sqrt(h1*(h2^h3));
    p_dot=(x-z)/y;
    
    
    %Calculation of the rest of zdot
    V_dot=r*Ab;
    Ab_dot=2*pi*r*no*(L-2*dp);
    dp_dot=2*r;
    L_dot=-2*r;
    
    zdot=[p_dot;V_dot;Ab_dot;dp_dot;L_dot];
end

%% Gets Exit Mach from Expansion Ratio and gamma
function mach=getExitMach(e,g)
    h1=(g+1)/2;
    h2=(g+1)/(2*g-2);
    h3=(g-1)/2;
    fun=@(x)((h1^-h2)*(1/x)*(1+h3*x^2)^h2-e);
    mach=fzero(fun,2);
end

%% Gets Pressure Ratio from Exit Mach and gamma
function pratio=getPressureRatio(M,g)
    h1=(g-1)/2;
    h2=-g/(g-1);
    pratio=(1+h1*M^2)^h2;
end

