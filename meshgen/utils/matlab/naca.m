function y=naca(x)

t=0.12;
c=1.00;
x=x/c;
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;

finite_TE=0;

if finite_TE
  a4 = -0.1015; % For finite thick TW
else
  a4 = -0.1036; % For zero thick TE
end

y=t/0.2*c*(a0*sqrt(x)+a1*x+a2*x.^2+ ...
	   a3*x.^3+a4*x.^4);
