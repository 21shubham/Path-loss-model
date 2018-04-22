

function [D]=diff_coef_hard_pec(n,incident_angle,diffracted_angle,L)


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

D=D1+D2+D3+D4;