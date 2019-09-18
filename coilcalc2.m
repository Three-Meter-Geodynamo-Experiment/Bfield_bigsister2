function [Bx,By,Bz,Br,Bphi,Btheta] = coilcalc2(X,Y,Z,str)
%function [Br,Bphi,Btheta] = coilcalc2(X,Y,Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This it a routine for the coil calculations based on the Biot-Savart
% laws written by Sarah Burnett 2019/09. 
% Specifically, this is for two coils with current in opposite directions.
% It uses elements of single coil calculation coilcalc.m written by Axl.
% This is used for the purpose of providing a nondimensionlized external
% magnetic field for simulations in xshells. 
% Input -- X,Y,Z coordinate to compute the magnetic field, str = 'q' or 'd'
% for quadrupole or dipole respectively
% Output -- B_x, B_y, B_z, B_r,B_phi,B_theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Fixed variables
mu0 = 0.0125663706; %permeability of free space in GAUSS*meter/amp
I = 300; %current (amperes)

N = 400;
t = [0:N]*2*pi/N; %parameterization
R = 1.9; %approximate magnet radius (m)
HT = .631825; %approximate height of top magnet (m)
HB = -0.53975; %approximate height of lower magnet (m)

%Top magnet location
x = R*cos(t); 
y = R*sin(t); 
zT = HT*ones(size(t));
 
% %Bottom magnet location
 zB = HB*ones(size(t));

%Calculate the midpoints
mpx = (circshift(x,-1)+x)/2;
mpx = mpx(1:end-1);
mpy = (circshift(y,-1)+y)/2;
mpy = mpy(1:end-1);
mpzT = (circshift(zT,-1)+zT)/2;
mpzT = mpzT(1:end-1);
mpzB = (circshift(zB,-1)+zB)/2;
mpzB = mpzB(1:end-1);

%calculate displacement vector components
rx = X-mpx;
ry = Y-mpy;
rzT = Z-mpzT;
rzB = Z-mpzB;

%calculate wire end displacement vector dL
dlx = (circshift(x,-1)-x);
dlx = dlx(1:end-1);
dly = (circshift(y,-1)-y);
dly = dly(1:end-1);
dlzT = (circshift(zT,-1)-zT);
dlzT = dlzT(1:end-1);
dlzB = (circshift(zB,-1)-zB);
dlzB = dlzB(1:end-1);

%calculate denominator
magr3T = (rx.^2+ry.^2+rzT.^2).^(3/2);
magr3B = (rx.^2+ry.^2+rzB.^2).^(3/2);

%calculate integrand
%three components of integrand: dl cross r over r3
%Top magnet
dBxT = (dly.*rzT-dlzT.*ry)./magr3T;
dByT = (dlzT.*rx-dlx.*rzT)./magr3T; 
dBzT = (dlx.*ry-dly.*rx)./magr3T; 

%Bottom magnet
dBxB = (dly.*rzB-dlzB.*ry)./magr3B;
dByB = (dlzB.*rx-dlx.*rzB)./magr3B; 
dBzB = (dlx.*ry-dly.*rx)./magr3B; 

%Take the difference between the two magnets
if str == 'q'
    Bx = mu0/(4*pi)*I*(sum(dBxT)-sum(dBxB));
    By = mu0/(4*pi)*I*(sum(dByT)-sum(dByB));
    Bz = mu0/(4*pi)*I*(sum(dBzT)-sum(dBzB));
elseif str == 'd'
    Bx = mu0/(4*pi)*I*(sum(dBxT)+sum(dBxB));
    By = mu0/(4*pi)*I*(sum(dByT)+sum(dByB));
    Bz = mu0/(4*pi)*I*(sum(dBzT)+sum(dBzB));
else
    error('Please choose q for quadrupole or d for dipole.')
end

 [phi,theta,~] = cart2sph(X,Y,Z); %doesn't work correctly so Br, Bphi, Bthet might be off...d 
 
 Br = Bx*sin(theta)*cos(phi) +By*sin(theta)*sin(phi)+Bz*cos(theta);
 Bphi = -Bx*sin(theta)+By*cos(theta);
 Btheta = Bx*cos(theta)*cos(phi)+By*cos(theta)*sin(phi)-Bz*sin(theta);

end