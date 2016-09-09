
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this program takes a discrete nozzle geometry with radius over z-axis and
%approximates by a monomial polynomial

% discrete r-z profile of the boundary
fileID = fopen('boundary_coords.dat','r');
A =fscanf(fileID,'%f %f',[2 Inf]);
Nin=size(A,1);
zdata=A(1,:);
rdata=A(2,:);



%number of points for plotting
Np=1000; 

%number of coefficients for the polynomial expansion
ncoefs=7

%zmin, zmax (should be close to input data!)
zmin=0.
zmax=3.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_plot=linspace(zmin-1,zmax+1,Np);

%interpolate r-z-profile
r_coefs=polyfit(zdata,rdata,ncoefs-1);


r_fit=polyval(r_coefs,z_plot);



%plot(zdata,rdata,'kx',z_plot,r_fit,'b');

maxerr=max(abs(rdata-polyval(r_coefs,zdata)));

%flip polynomial coefficients to smallest exponent first
r_coefs_flip=fliplr(r_coefs);

% %improved interpolation with shift (X-MU(1))/MU(2)
% [r_coefs2,s,mu]=polyfit(zdata,rdata,ncoefs-1);
% 
% r_fit2=polyval(r_coefs,z_plot,s,mu);
% 
% maxerr2=max(abs(rdata-polyval(r_coefs2,zdata,s,mu)));
% 
% %flip polynomial coefficients to smallest exponent first
% r_coefs2_flip=fliplr(r_coefs2);

%OUTPUT
fprintf('==== OUTPUT ==== \n');

fprintf('Max. error on input points  %21.15E \n \n',maxerr);

fprintf('==== PROFILE COEFFICIENTS ==== \n');
fprintf('!R_COEFS (MONOMIAL, zero exponent first): \n');

for i=1:ncoefs
  if(mod(i,4)==0)
    fprintf('\n');
  end
  fprintf('  %21.15E',r_coefs_flip(i));
end
fprintf('\n');
fprintf('!Z-pos (Eq.dist. points): \n');
zIP=linspace(zmin,zmax,ncoefs);
rIP=polyval(r_coefs,zIP);
for i=1:ncoefs
  if(mod(i,4)==0)
    fprintf('\n');
  end
  fprintf('  %21.15E',zIP(i));
end
fprintf('\n');
fprintf('\n');
fprintf('!R-pos (Eq.dist. points): \n');
for i=1:ncoefs
  if(mod(i,4)==0)
    fprintf('\n');
  end
  fprintf('  %21.15E',rIP(i));
end
fprintf('\n');

fprintf('!FORTRAN EVALUATION: \n');
fprintf('rz=');
for i=1:ncoefs
   fprintf('(');
end
fprintf(' %21.15E',r_coefs(1));
for i=2:ncoefs
  if(mod(i,4)==0)
    fprintf(' & \n  ');
  end
  fprintf(')*z + %21.15E',r_coefs(i));
end
fprintf(')\n');


 icheck=0*z_plot+r_coefs_flip(1);
 z_exp=z_plot;
 for i=2:ncoefs
     icheck = icheck+r_coefs_flip(i)*z_exp;
     z_exp=z_exp.*z_plot;
 end
 plot(zdata,rdata,'kx',z_plot,icheck,'r',zIP,rIP,'+');
 hold on;
 plot(zdata,-rdata,'kx',z_plot,-icheck,'r',zIP,-rIP,'+');
 axis([zmin-1 zmax+1 -2 2]);
 axis equal;
 hold off;