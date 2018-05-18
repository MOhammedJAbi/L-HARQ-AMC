
function [solu,xs,n] = Brent(f,a,b,tol_step, tol_abs,nmax)

% find the solution of f(x)=0 or abs(f(x))< tol=10^-alpha

fa=f(a); fb=f(b);

n=0;
if (abs(fa)< tol_abs )
	solu = a;
    xs=a;
	return;
elseif ( abs(fb)< tol_abs )
	solu = b;
    xs=b;
	return;
end

if sign(fa) == sign(fb)
error( 'f(a) and f(b) do not have opposite signs' );
end
%
if abs(fa) < abs(fb)
tmp=fa;
fa=fb;
fb=tmp;
tmp=a;
a=b;
b=tmp;
end
%
c=a;
fc=fa;
mflag=1;

while ((abs(b-a)>tol_step) && (n<nmax) && (abs(fb)>tol_abs))
n=n+1;
lint(n)=abs(b-a);
%fc=f(c);
if ((fa~=fc) && (fb~=fc))
xs(n)=(a*fb*fc)/((fa-fb)*(fa-fc))+(b*fa*fc)/((fb-fa)*(fb-fc))+(c*fa*fb)/((fc-fa)*(fc-fb));
else
xs(n)=b-fb*(b-a)/(fb-fa);
end

if (xs(n)<(3*a+b)/4) || (xs(n)>b) || (mflag && (abs(xs(n)-b)>=(abs(b-c)/2))) || (~mflag && (abs(xs(n)-b)>=(abs(c-d)/2)))
xs(n)=(a+b)/2;
mflag=1;
else
mflag=0;
end
fs=f(xs(n));
d=c;
c=b;
fc=fb;
if fa*fs<0
b=xs(n);
fb=fs;
else
a=xs(n);
fa=fs;%f(xs(n));
end
%fa=f(a);
%fb=f(b);
if abs(fa)<abs(fb)
tmp=fa;
fa=fb;
fb=tmp;
tmp=a;
a=b;
b=tmp;
end

solu=b;
end

if((abs(fb)>tol_abs))
     error( 'the method did not converge' );
end
