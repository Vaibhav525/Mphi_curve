%%Strain softening homogenous

%%Input
%Geometry(Rectangular)
b=200;        %mm
d=500;        %mm
N=500;        %No. of disretized strips along depth
d_phi=1e-8;  %Delta phi, mm^-1
%Material model
E=2e5;        %Young's Modulus N/mm2
sig_y=250;    %Yield stress:   Mpa
e_y=sig_y/E;  %Yield strain
s_ratio=0.1;  %Strain softening ratio E_soft=-E*s_ratio
E_soft=-E*s_ratio;
e_max=0.0035; %Maximum strain after which unloading occurs
%kinematic softening 

%%Initializing M-phi curve array values
%At each step we will increse phi and calculate M and store at
%Moment[counter] and Curvature[counter], then counter++
counter=1;
curr_phi=0;
Moment=zeros(100000,1);
Curvature=zeros(100000,1);

%Dicretizing the section along depth
Fiber_dist=zeros(N,1);
Fiber_dist(1,1)=d/(2*N);
Fiber_depth=(d/N);
for i=2:length(Fiber_dist)
    Fiber_dist(i,1)=Fiber_dist(i-1,1)+Fiber_depth;
end
Fiber_dist=Fiber_dist-(d/2);   %Distances from neutral axis

Strain=zeros(N,1);      %Strain along depth
Stress=zeros(N,1);      %Stress along depth

%Case-1: Elastic region, e<e_y
e=curr_phi*(d/2); %Strain in extreme fiber
while(e<=e_y)
    %Calculate strains
    Strain=curr_phi.*Fiber_dist;
    Stress=E*Strain;
    Moment(counter,1)=sum(Stress.*Fiber_dist*b*Fiber_depth*1e-6);
    Curvature(counter,1)=curr_phi;
    curr_phi=curr_phi+d_phi;
    e=curr_phi*(d/2);
    counter=counter+1;
    
end
figure();
plot(Stress,Fiber_dist);
xlabel("Stress--->")
ylabel("x----->")
axis ( [-2*max(Stress) 2*max(Stress) -d d] )
grid on;
axis square;

figure;
plot(Curvature(1:counter-1),Moment(1:counter-1));
xlabel("Curvature--->")
xlabel("Moment(kN/m)----->")

%Case-2: Strain softening, e_y<e<e_max
e=curr_phi*(d/2); %Strain in extreme fiber
while(e<=e_y)
    %Calculate strains
    Strain=curr_phi.*Fiber_dist;
    Stress=E*Strain;
    Moment(counter,1)=sum(Stress.*Fiber_dist*b*Fiber_depth*1e-6);
    Curvature(counter,1)=curr_phi;
    curr_phi=curr_phi+d_phi;
    e=curr_phi*(d/2);
    counter=counter+1;
    
end
figure();
plot(Stress,Fiber_dist);
xlabel("Stress--->")
ylabel("x----->")
axis ( [-2*max(Stress) 2*max(Stress) -d d] )
grid on;
axis square;

figure;
plot(Curvature(1:counter-1),Moment(1:counter-1));
xlabel("Curvature--->")
xlabel("Moment(kN/m)----->")