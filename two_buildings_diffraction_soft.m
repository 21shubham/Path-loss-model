

clc
clear all

E0=1;alpha=2*pi/180;
w1=10;d1=2;

n=3/2;

global f
f=28*10^9;
lambda=(3*10^8)/f;

global k
k=2*pi/lambda;

%---calculating effective field at the left edge of second building---

%---first calculating direct field component---

Einc_1L=E0; % direct field at the left edge of first building
Einc_1R=E0*exp(-j*k*d1*cos(alpha)); % direct field at the right edge of first building
Einc_2L=E0*exp(-j*k*w1*cos(alpha)); % direct field at the left edge of second building

%---now calculating diffracted field component---

%---calculating D_alpha_L1---

incident_angle=alpha+(pi/2);
diffracted_angle=3*(pi/2);
L=d1;

D_alpha_L1=diff_coef_soft(n,incident_angle,diffracted_angle,L);

A1=1/sqrt(w1);
dummy5=E0*D_alpha_L1*A1*exp(-j*k*w1);

%---calculating D_alpha_R1---

incident_angle=alpha;
diffracted_angle=pi;
L=w1-d1;

D_alpha_R1=diff_coef_soft(n,incident_angle,diffracted_angle,L);

A1=1/sqrt(w1-d1);
dummy7=E0*D_alpha_R1*A1*exp(-j*k*(w1-d1));

E_diff=dummy5+dummy7;

E_total=abs(Einc_2L+E_diff)

hold on



