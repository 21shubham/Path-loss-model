

function [D]=diff_coef_hard(n,incident_angle,diffracted_angle,L)


global k
global f

A=(-exp(-j.*pi./4))./(2.*n.*sqrt(2.*pi.*k));

%========= To compute D1===========

beta1=diffracted_angle-incident_angle;
gamma1=(pi+(beta1))./(2.*n);
x=2.*k.*L.*(n.^2).*(sin(gamma1).^2);

F(x<100)=sqrt(j.*pi.*x(x<100)).*exp(j.*x(x<100)).*(1-erfz(sqrt(j.*x(x<100))));
F(x>100)=1+(j.*(1./2./x(x>100)))-((3./4).*(1./x(x>100).^2))-(j.*(15./8).*(1./x(x>100).^3))+((75./16).*(1./x(x>100).^4));

D1=A.*cot(gamma1).*F;

%========= To compute D2===========

beta2=diffracted_angle-incident_angle;
gamma2=(pi-(beta2))./(2.*n);
x=2.*k.*L.*(n.^2).*(sin(gamma2).^2);

F(x<100)=sqrt(j.*pi.*x(x<100)).*exp(j.*x(x<100)).*(1-erfz(sqrt(j.*x(x<100))));
F(x>100)=1+(j.*(1./2./x(x>100)))-((3./4).*(1./x(x>100).^2))-(j.*(15./8).*(1./x(x>100).^3))+((75./16).*(1./x(x>100).^4));

D2=A.*cot(gamma2).*F;

%========= To compute D3===========

beta3=diffracted_angle+incident_angle;
gamma3=(pi-(beta3))./(2.*n);
x=2.*k.*L.*(n.^2).*(sin(gamma3).^2);

F(x<100)=sqrt(j.*pi.*x(x<100)).*exp(j.*x(x<100)).*(1-erfz(sqrt(j.*x(x<100))));
F(x>100)=1+(j.*(1./2./x(x>100)))-((3./4).*(1./x(x>100).^2))-(j.*(15./8).*(1./x(x>100).^3))+((75./16).*(1./x(x>100).^4));

D3=A.*cot(gamma3).*F;

%========= To compute D4===========

beta4=diffracted_angle+incident_angle;
gamma4=(pi+(beta4))./(2.*n);
x=2.*k.*L.*(n.^2).*(sin(gamma4).^2);

F(x<100)=sqrt(j.*pi.*x(x<100)).*exp(j.*x(x<100)).*(1-erfz(sqrt(j.*x(x<100))));
F(x>100)=1+(j.*(1./2./x(x>100)))-((3./4).*(1./x(x>100).^2))-(j.*(15./8).*(1./x(x>100).^3))+((75./16).*(1./x(x>100).^4));

D4=A.*cot(gamma4).*F;

epsilon_0=8.85419e-12; %in F/m,free space permittivitty
epsilon_r=5; % relative permittivity
conductivity=0.016;   % its the conductivity for lossy dielectric in S/m

psi_0=incident_angle;
psi_n=(n.*pi)-diffracted_angle;

epsilon_com=epsilon_r-(j.*conductivity./(2.*pi.*f.*epsilon_0));
dummy1=(epsilon_com.*sin(psi_0))-sqrt(epsilon_com-cos(psi_0).^2);
dummy2=(epsilon_com.*sin(psi_0))+sqrt(epsilon_com-cos(psi_0).^2);
r_0_hard=(dummy1./dummy2);

dummy3=(epsilon_com.*sin(psi_n))-sqrt(epsilon_com-cos(psi_n).^2);
dummy4=(epsilon_com.*sin(psi_n))+sqrt(epsilon_com-cos(psi_n).^2);
r_n_hard=(dummy3./dummy4);

D=(r_0_hard.*r_n_hard.*D1)+D2+(r_0_hard.*D3)+(r_n_hard.*D4);