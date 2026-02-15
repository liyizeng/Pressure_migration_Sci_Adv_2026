% This code is used in: Ayuba Akinpelu et al.
% Hydrostatic Pressure Reduces the Mechanosensitivity of Cell Migration
% Journal: Science Advances.
% Year: 2026

% If you have any questions on the code or the math model, please contact:
% Yizeng Li, liyizeng52@hotmail.com

% F-actin & G-actin
% Steady-state
% Unknowns: vn, theta_n, theta_c, v0
% Velocities and equations are written in the frame of the moving cell
% Parameters are in units: μm, s, Pa & mM
% Loop over k_st: substrate stiffness maps to eta_st


clc
clear
close all


FontSize = 18;

RankCheck = 0;
Iter = 21;
N  = 101;

%% Parameters are in units: nm, s, Pa & mM

L = 50;                 % (μm) cell length
dx = L/(N-1);
x  = linspace(0,L,N);
R = 8.31451;            % (J/mol K) Ideal gas constant
T = 310;                % (K) absolute temperature
vc = 0;

Thetac = 0.1;           % (mM) reference value of G-actin
Thetan = 0.2;           % (mM) reference value of F-actin

Jactin0 = 4d-4;         % (mM/s) Coefficient of actin polymerization
thetacc = 0.2d-3;       % (mM) G-actin constant in actin polymerization

gamma0 = 2.3d-3;        % (1/s) Baseline actin depolymerization rate


ksigman = 1d4;          % (Pa /mM) Coefficient of passive actin network pressure

etast0 = 2.0d2;         % (Pa s/μm^2/mM) Baseline focal adhesion strength
etastm = 2.0d4;         % (Pa s/μm^2/mM) Conversion factor between kst and etast
kst_0  = 10;            % (kPa) reference stiffness, for normalization
eta    = 1d-2;          % (Pa s/μm^2/mM)
kad    = 2d3;           % (Pa s/μm) adhesive force, Fad^b = kad*v0

Dtc = 1.d1;             % (μm^2/s) diffusion constant for theta_c



%% Equations

% number of variables of the model
n_var = 4; % v_n, theta_n, theta_c, v0
var_legh = [N, N, N, 1]; % length of eavh variable
var_strt = zeros(1,n_var); % start position of each variable
var_strt(1) = 1;
for in = 2:n_var
    var_strt(in) = var_strt(in-1) + var_legh(in-1);
end
N_var = var_strt(n_var) + var_legh(n_var) - 1; % total size of the matrix

s_vn   = var_strt(1);
s_tn   = var_strt(2);
s_tc   = var_strt(3);
s_v0   = var_strt(4);

N2 = 31;
KST = logspace(-1,1,N2)*3;
Axis_m   = KST; % kPa

V0  = zeros(1,N2);

for loop2 = 1:N2
    fprintf('Outer looop #: %d\n',loop2);
    kst = KST(loop2);

    etast   = (etast0 + (kst/kst_0)^2 * etastm ) * linspace(1,1,N)';

    % initial guess
    thetan = linspace(Thetan*1,Thetan*1,N)';
    thetac = linspace(Thetac*1,Thetac*1,N)';
    v0 = (etast(1)*Jactin0*L)/(Thetan*etast(1)+kad/L);
    vn = linspace(v0,v0,N)';


    X = [vn; thetan; thetac; v0];

    iter = 0;
    ITER = true;
    while ITER
        iter = iter + 1;

        DF = zeros(N_var,N_var);
        Fn = zeros(N_var,1);

        sigma_n = ksigman*thetan;     % N by 1 vector
        sigma   = sigma_n;
        dsigmadtn =  ksigman;

        gamma = gamma0*logspace(-0,0,N)';

        Jactin = Jactin0*linspace(0,1,N)'.*thetac./(thetacc+thetac);
        dJactindtc   = Jactin0*thetacc./(thetacc+thetac).^2;
        dJactindRaca = zeros(N,1);


        %% Equations for vn and the Derivatives
        Fn(s_vn) = - (-3*sigma(1) + 4*sigma(2) - sigma(3)) ...
            + 2*dx*eta*thetan(1)*(vc - vn(1)) ...
            - 2*dx*etast(1)*thetan(1)*(vn(1) + v0);
        DF(s_vn,s_vn)   = -2*dx*eta*thetan(1) - 2*dx*etast(1)*thetan(1);
        DF(s_vn,s_tn)   =  3*dsigmadtn + 2*dx*eta*(vc - vn(1)) ...
            - 2*dx*etast(1)*(vn(1) + v0);
        DF(s_vn,s_tn+1) = -4*dsigmadtn;
        DF(s_vn,s_tn+2) =    dsigmadtn;
        DF(s_vn,s_v0) = -2*dx*etast(1)*thetan(1);

        Fn(s_vn+1:s_vn+N-2) = -(sigma(3:N)-sigma(1:N-2)) ...
            + 2*dx*eta*thetan(2:N-1).*(vc-vn(2:N-1)) ...
            - 2*dx*etast(2:N-1).*thetan(2:N-1).*(vn(2:N-1)+v0);
        for i = 2:N-1
            DF(s_vn+i-1,s_vn+i-1) = -2*dx*thetan(i)*(eta+etast(i));
            DF(s_vn+i-1,[s_tn+i-2,s_tn+i]) = [1,-1]*dsigmadtn;
            DF(s_vn+i-1,s_tn+i-1) = 2*dx*eta*(vc-vn(i)) - 2*dx*etast(i)*(vn(i)+v0);
            DF(s_vn+i-1,s_v0) = -2*dx*etast(i)*thetan(i);
        end

        Fn(s_vn+N-1) = thetan(N)*vn(N);
        DF(s_vn+N-1,s_vn+N-1) = thetan(N);
        DF(s_vn+N-1,s_tn+N-1) = vn(N);

        %% Equations for thetan and the Derivatives
        Fn(s_tn) = thetan(1)*vn(1);
        DF(s_tn,s_vn) = thetan(1);
        DF(s_tn,s_tn) = vn(1);

        Fn(s_tn+1:s_tn+N-2) = (thetan(3:N).*vn(3:N)-thetan(1:N-2).*vn(1:N-2))...
            + 2*dx*gamma(2:N-1).*thetan(2:N-1) - 2*dx*Jactin(2:N-1);
        for i = 2:N-1
            DF(s_tn+i-1,[s_vn+i-2,s_vn+i]) = [-thetan(i-1), thetan(i+1)];
            DF(s_tn+i-1,s_tn+i-2) = -vn(i-1);
            DF(s_tn+i-1,s_tn+i-1) = 2*dx*gamma(i);
            DF(s_tn+i-1,s_tn+i)   =  vn(i+1);
            DF(s_tn+i-1,s_tc+i-1) = -2*dx*dJactindtc(i);
        end

        Fn(s_tn+N-1) = (3*thetan(N)*vn(N) - 4*thetan(N-1)*vn(N-1) ...
            + thetan(N-2)*vn(N-2)) + 2*dx*gamma(N)*thetan(N) - 2*dx*Jactin(N);
        DF(s_tn+N-1,s_vn+N-3:s_vn+N-1) = [thetan(N-2), -4*thetan(N-1), 3*thetan(N)];
        DF(s_tn+N-1,s_tn+N-3) =    vn(N-2);
        DF(s_tn+N-1,s_tn+N-2) = -4*vn(N-1);
        DF(s_tn+N-1,s_tn+N-1) =  3*vn(N) + 2*dx*gamma(N);
        DF(s_tn+N-1,s_tc+N-1) = -2*dx*dJactindtc(N);

        %% Equations for thetac and the Derivatives
        Fn(s_tc) = thetac(1)*vc - Dtc/2/dx*(-3*thetac(1) + 4*thetac(2) - thetac(3));
        DF(s_tc,s_tc)   = vc + 3*Dtc/2/dx;
        DF(s_tc,s_tc+1) =    - 4*Dtc/2/dx;
        DF(s_tc,s_tc+2) =        Dtc/2/dx;

        Fn(s_tc+1:s_tc+N-2) = vc * (thetac(3:N)-thetac(1:N-2)) ...
            - 2*Dtc/dx*(thetac(1:N-2) - 2*thetac(2:N-1) + thetac(3:N))...
            - 2*dx*gamma(2:N-1).*thetan(2:N-1) + 2*dx*Jactin(2:N-1);
        for i = 2:N-1
            DF(s_tc+i-1,s_tn+i-1) = -2*dx*gamma(i);
            DF(s_tc+i-1,s_tc+i-2) = -vc - 2*Dtc/dx;
            DF(s_tc+i-1,s_tc+i-1) = 4*Dtc/dx + 2*dx*dJactindtc(i);
            DF(s_tc+i-1,s_tc+i)   =  vc - 2*Dtc/dx;
        end

        Fn(s_tc+N-1) = 1/2*(sum(thetan(1:N-1)+thetac(1:N-1)) + sum(thetan(2:N)+thetac(2:N)))...
            - L*(Thetan + Thetac)/dx;
        DF(s_tc+N-1,[s_tn,s_tn+N-1]) = 1/2;
        DF(s_tc+N-1,s_tn+1:s_tn+N-2) = 1;
        DF(s_tc+N-1,[s_tc,s_tc+N-1]) = 1/2;
        DF(s_tc+N-1,s_tc+1:s_tc+N-2) = 1;



        %% Equation for v0
        Fn(s_v0) = kad*v0 ...
            + dx/2*(sum(etast(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
            + sum(etast(2:N).*thetan(2:N).*vn(2:N)))...
            + dx/2*v0*(sum(etast(1:N-1).*thetan(1:N-1)) + sum(etast(2:N).*thetan(2:N)));
        DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[etast(1)*thetan(1),etast(N)*thetan(N)];
        DF(s_v0,s_vn+1:s_vn+N-2) = dx*etast(2:N-1).*thetan(2:N-1);
        DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[etast(1)*(vn(1)+v0),etast(N)*(vn(N)+v0)];
        DF(s_v0,s_tn+1:s_tn+N-2) = dx*etast(2:N-1).*(vn(2:N-1)+v0);
        DF(s_v0,s_v0) = kad ...
            + dx/2*(sum(etast(1:N-1).*thetan(1:N-1)) + sum(etast(2:N).*thetan(2:N)));

        %% Solve for the matrix
        DF = sparse(DF);
        X = X - DF\Fn;

        vn     = X(s_vn:s_vn+N-1);
        thetan = X(s_tn:s_tn+N-1);
        thetac = X(s_tc:s_tc+N-1);
        v0     = X(s_v0);

        if iter > 1
            error = abs((X-temp_X)./X);
            error = sum(error)/(N_var);
            if error < 1d-5 || iter == Iter
                ITER = false;
            end
        end
        temp_X = X;
    end

    V0(loop2) = v0;
end

%% Plotting
figure
plot(KST,V0*3600,'linewidth',2); hold on
xlabel('Substrate Stiffness (kPa)')
ylabel('v_{cell} (µm/h)');
set(gca, 'xscale','log','FontSize', FontSize);
xlim([1d-1,100]);
xticks([0.1,1,10,100]);
box off

