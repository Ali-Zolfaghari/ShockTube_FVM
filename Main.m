
%================================ Roe =====================================

clear,clc

gamma = 1.4;
CFL = 0.5;
L = 1;
N = 100;
xd = 0.5;
time = 0;
FinalTime = 0.5;
FinalStep = 100;
step = 1;

dx = L/N;
Nd = fix(xd/dx)+1;
X = 0:dx:L;

cond = 1; % sod : cond = 0; Lax : cond = 1

if (cond==0)
    Rho_Left = 1.0;
    Rho_Right = 0.125;
    Velocity_Left = 0;
    Velocity_Right = 0;
    Pressure_Left = 1;
    Pressure_Right = 0.1;
    
else
    Rho_Left = .445;
    Rho_Right = 0.5;
    Velocity_Left = .698;
    Velocity_Right = 0;
    Pressure_Left = 3.528;
    Pressure_Right = 0.571;
end

r(1:Nd) = Rho_Left;
r(Nd+1:N+1) = Rho_Right;
u(1:Nd) = Velocity_Left;
u(Nd+1:N+1) = Velocity_Right;
p(1:Nd) = Pressure_Left;
p(Nd+1:N+1) = Pressure_Right;

q = zeros(N+1,3);
ua = zeros(N+1,1);
alpha = zeros(3,N);
F = zeros(3,N);
lambda = zeros(3,N);

while ( time<FinalTime && step<FinalStep )
    
    for i = 1:N+1
        ua(i) = max(abs(u(i))+sqrt((gamma*p(i))/r(i)));
    end
    dt = (dx*CFL)/max(ua);
    
    h = (gamma/(gamma-1))*(p./r)+0.5*(u.^2);
    
    for i = 1:N
        RHO = sqrt(r(i+1)*r(i));
        VEL = (sqrt(r(i+1))*u(i+1)+sqrt(r(i))*u(i))/(sqrt(r(i+1))+sqrt(r(i)));
        ENG = (sqrt(r(i+1))*h(i+1)+sqrt(r(i))*h(i))/(sqrt(r(i+1))+sqrt(r(i)));
        SPD = sqrt((gamma-1)*(ENG-0.5*(VEL^2)));
        [ F(:,i),R(:,3*(i-1)+1:3*i),alpha(:,i),lambda(:,i) ] = Flux( RHO,VEL,ENG,SPD,r(i+1),r(i),u(i+1),u(i),p(i+1),p(i),h(i+1),h(i) );
    end
    
    for i = 2:N
        Q(1,1) = r(i);
        Q(2,1) = r(i)*u(i);
        Q(3,1) = (p(i)/(gamma-1))+0.5*r(i)*(u(i)^2);
        Fluxp_Roe = zeros(3,1);
        Fluxn_Roe = zeros(3,1);
        for j = 1:3
            SUMp_Roe = 0;
            SUMn_Roe = 0;
            for k = 1:3
                SUMp_Roe = SUMp_Roe+lambda(k,i)*alpha(k,i)*R(j,3*(i-1)+k);
                SUMn_Roe = SUMn_Roe+lambda(k,i-1)*alpha(k,i-1)*R(j,3*(i-2)+k);
            end
            Fluxp_Roe(j,1) = F(j,i)-0.5*SUMp_Roe;
            Fluxn_Roe(j,1) = F(j,i-1)-0.5*SUMn_Roe;
        end
        Fluxp = Fluxp_Roe;
        Fluxn = Fluxn_Roe;
        Q = Q-(dt/dx)*(Fluxp-Fluxn);
        q(i,1:3) = Q;
    end
    
    r = q(:,1);
    u = q(:,2)./q(:,1);
    p = (gamma-1)*(q(:,3)-0.5*r.*u.^2);
    
    r(1,:) = r(2,:);u(1,:) = u(2,:);p(1,:) = p(2,:);
    r(N+1,:) = r(N,:);u(N+1,:) = u(N,:);p(N+1,:) = p(N,:);
    
    time = time+dt;
    step = step+1;
    
end

p_1 = Pressure_Left;
rho_1 = Rho_Left;
u_1 = Velocity_Left;

p_4 = Pressure_Right;
rho_4 = Rho_Right;
u_4 = Velocity_Right;

a_1 = sqrt((gamma*p_1)/rho_1);
a_4 = sqrt((gamma*p_4)/rho_4);

p_14 = p_1/p_4;
p_34 = 0.00001;

beta = (1+(0.5/a_1)*(u_1-u_4)*(gamma-1)-(gamma-1)*(a_4/a_1)*(p_34-1)*((2*gamma)*((gamma-1)+(gamma+1)*p_34))^(-0.5))^((2*gamma)/(1-gamma));
while ((p_14-p_34*beta)>0.001)
    p_34 = p_34+0.00001;
    beta = (1+(0.5/a_1)*(u_1-u_4)*(gamma-1)-(gamma-1)*(a_4/a_1)*(p_34-1)*((2*gamma)*((gamma-1)+(gamma+1)*p_34))^(-0.5))^((2*gamma)/(1-gamma));
end

p_3=p_4*p_34;
w = (a_4*sqrt(1+((1+gamma)/(2*gamma))*(p_34-1)))+u_4;
up = u_4+(a_4/gamma)*(p_34-1)*sqrt(((2*gamma)/(1+gamma))/(p_34-((1-gamma)/(1+gamma))));
u_2 = up;
u_3 = up;
p_2 = p_3;
rho_3 = rho_4/((((gamma+1)/(gamma-1))+p_34)/(1+p_34*((gamma+1)/(gamma-1))));
rho_2 = rho_1*((p_2/p_1)^(1/gamma));
a_3 = sqrt((gamma*p_3)/rho_3);
a_2 = sqrt((gamma*p_2)/rho_2);

conpos = xd+up*time;
spos = xd+w*time;
pos_1 = xd+(u_1-a_1)*time;
pos_2 = xd+(u_2-a_2)*time;

pos = pos_1:(pos_2-pos_1)/20:pos_2;
x = [0,pos_1,pos,pos_2,conpos,conpos,spos,spos,L,L];

ue = (2/(gamma+1))*((0.5*(gamma-1)*u_1)+a_1+((pos-xd)/time));
pe = p_1*(1-(0.5*(gamma-1))*((ue-u_1)./a_1)).^((2*gamma)/(gamma-1));
re = rho_1*(1-(0.5*(gamma-1))*((ue-u_1)./a_1)).^(2/(gamma-1));
P = [p_1,p_1,pe,p_2,p_2,p_3,p_3,p_4,p_4,0];
U = [u_1,u_1,ue,u_2,u_2,u_3,u_3,u_4,u_4,0];
R = [rho_1,rho_1,re,rho_2,rho_2,rho_3,rho_3,rho_4,rho_4,0];
H = (gamma/(gamma-1))*(P./R)+0.5*(U.^2);

fig = figure(1);
plot(x,R,'--',X,r);
title('density');
legend('exact','Roe');
saveas(fig,['Roe_r',num2str(cond),'.jpg'],'jpeg');

fig = figure(2);
plot(x,U,'--',X,u);
title('velocity');
legend('exact','Roe');
saveas(fig,['Roe_v',num2str(cond),'.jpg'],'jpeg');

fig = figure(3);
plot(x,P,'--',X,p);
title('pressure');
legend('exact','Roe');
saveas(fig,['Roe_p',num2str(cond),'.jpg'],'jpeg');

fig = figure(4);
plot(x,H,'--',X,h);
title('enthalpy');
legend('exact','Roe')
saveas(fig,['Roe_h',num2str(cond),'.jpg'],'jpeg');

close all

%================================ RoeLax ==================================

clear,clc

gamma = 1.4;
CFL = 0.5;
L = 1;
N = 100;
xd = 0.5;
time = 0;
FinalTime = 0.5;
FinalStep = 100;
step = 1;
delta = 0.000001;

dx = L/N;
Nd = fix(xd/dx)+1;
X = 0:dx:L;

cond = 1; % sod : cond = 0; Lax : cond = 1

if (cond==0)
    Rho_Left = 1.0;
    Rho_Right = 0.125;
    Velocity_Left = 0;
    Velocity_Right = 0;
    Pressure_Left = 1;
    Pressure_Right = 0.1;
    
else
    Rho_Left = .445;
    Rho_Right = 0.5;
    Velocity_Left = .698;
    Velocity_Right = 0;
    Pressure_Left = 3.528;
    Pressure_Right = 0.571;
end

r(1:Nd) = Rho_Left;
r(Nd+1:N+1) = Rho_Right;
u(1:Nd) = Velocity_Left;
u(Nd+1:N+1) = Velocity_Right;
p(1:Nd) = Pressure_Left;
p(Nd+1:N+1) = Pressure_Right;

q = zeros(N+1,3);
ua = zeros(N+1,1);
alpha = zeros(3,N);
teta = zeros(3,N);
phi = zeros(3,N);
F = zeros(3,N);
lambda = zeros(3,N);

while ( time<FinalTime && step<FinalStep )
    
    for i = 1:N+1
        ua(i) = max(abs(u(i))+sqrt((gamma*p(i))/r(i)));
    end
    dt = (dx*CFL)/max(ua);
    
    h = (gamma/(gamma-1))*(p./r)+0.5*(u.^2);
    
    for i = 1:N
        RHO = sqrt(r(i+1)*r(i));
        VEL = (sqrt(r(i+1))*u(i+1)+sqrt(r(i))*u(i))/(sqrt(r(i+1))+sqrt(r(i)));
        ENG = (sqrt(r(i+1))*h(i+1)+sqrt(r(i))*h(i))/(sqrt(r(i+1))+sqrt(r(i)));
        SPD = sqrt((gamma-1)*(ENG-0.5*(VEL^2)));
        [ F(:,i),R(:,3*(i-1)+1:3*i),alpha(:,i),lambda(:,i) ] = Flux( RHO,VEL,ENG,SPD,r(i+1),r(i),u(i+1),u(i),p(i+1),p(i),h(i+1),h(i) );
    end
    
    for i = 1:N
        for j = 1:3
            k = i-sign(lambda(j,i));
            if (k==0 || k==N+1)
                teta(j,i) = 1;
            else
                if (alpha(j,i)~=0)
                    teta(j,i) = alpha(j,k)/alpha(j,i);
                else
                    teta(j,i) = alpha(j,k)/delta;
                end
            end
            
        end
    end
    
    for i = 1:N
        for j = 1:3
            phi(j,i) = max([0,min([1,2*teta(j,i)]),min([2,teta(j,i)])]); %Superbee
            %phi(j,i) = (abs(teta(j,i))+teta(j,i))/(1+abs(teta(j,i))); %VanLeer
        end
    end
    
    for i = 2:N
        Q(1,1) = r(i);
        Q(2,1) = r(i)*u(i);
        Q(3,1) = (p(i)/(gamma-1))+0.5*r(i)*(u(i)^2);
        Fluxp_Roe = zeros(3,1);
        Fluxn_Roe = zeros(3,1);
        Fluxp_RoeLax = zeros(3,1);
        Fluxn_RoeLax = zeros(3,1);
        for j = 1:3
            SUMp_Roe = 0;
            SUMn_Roe = 0;
            SUMp_RoeLax = 0;
            SUMn_RoeLax = 0;
            for k = 1:3
                SUMp_Roe = SUMp_Roe+lambda(k,i)*alpha(k,i)*R(j,3*(i-1)+k);
                SUMn_Roe = SUMn_Roe+lambda(k,i-1)*alpha(k,i-1)*R(j,3*(i-2)+k);
                SUMp_RoeLax = SUMp_RoeLax+phi(k,i)*(sign(lambda(k,i))-(dt/dx)*lambda(k,i))*lambda(k,i)*alpha(k,i)*R(j,3*(i-1)+k);
                SUMn_RoeLax = SUMn_RoeLax+phi(k,i-1)*(sign(lambda(k,i-1))-(dt/dx)*lambda(k,i-1))*lambda(k,i-1)*alpha(k,i-1)*R(j,3*(i-2)+k);
            end
            Fluxp_Roe(j,1) = F(j,i)-0.5*SUMp_Roe;
            Fluxn_Roe(j,1) = F(j,i-1)-0.5*SUMn_Roe;
            Fluxp_RoeLax(j,1) = 0.5*SUMp_RoeLax;
            Fluxn_RoeLax(j,1) = 0.5*SUMn_RoeLax;
        end
        Fluxp = Fluxp_Roe+Fluxp_RoeLax;
        Fluxn = Fluxn_Roe+Fluxn_RoeLax;
        Q = Q-(dt/dx)*(Fluxp-Fluxn);
        q(i,1:3) = Q;
    end
    
    r = q(:,1);
    u = q(:,2)./q(:,1);
    p = (gamma-1)*(q(:,3)-0.5*r.*u.^2);
    
    r(1,:) = r(2,:);u(1,:) = u(2,:);p(1,:) = p(2,:);
    r(N+1,:) = r(N,:);u(N+1,:) = u(N,:);p(N+1,:) = p(N,:);
    
    time = time+dt;
    step = step+1;
    
end

p_1 = Pressure_Left;
rho_1 = Rho_Left;
u_1 = Velocity_Left;

p_4 = Pressure_Right;
rho_4 = Rho_Right;
u_4 = Velocity_Right;

a_1 = sqrt((gamma*p_1)/rho_1);
a_4 = sqrt((gamma*p_4)/rho_4);

p_14 = p_1/p_4;
p_34 = 0.00001;

beta = (1+(0.5/a_1)*(u_1-u_4)*(gamma-1)-(gamma-1)*(a_4/a_1)*(p_34-1)*((2*gamma)*((gamma-1)+(gamma+1)*p_34))^(-0.5))^((2*gamma)/(1-gamma));
while ((p_14-p_34*beta)>0.001)
    p_34 = p_34+0.00001;
    beta = (1+(0.5/a_1)*(u_1-u_4)*(gamma-1)-(gamma-1)*(a_4/a_1)*(p_34-1)*((2*gamma)*((gamma-1)+(gamma+1)*p_34))^(-0.5))^((2*gamma)/(1-gamma));
end

p_3=p_4*p_34;
w = (a_4*sqrt(1+((1+gamma)/(2*gamma))*(p_34-1)))+u_4;
up = u_4+(a_4/gamma)*(p_34-1)*sqrt(((2*gamma)/(1+gamma))/(p_34-((1-gamma)/(1+gamma))));
u_2 = up;
u_3 = up;
p_2 = p_3;
rho_3 = rho_4/((((gamma+1)/(gamma-1))+p_34)/(1+p_34*((gamma+1)/(gamma-1))));
rho_2 = rho_1*((p_2/p_1)^(1/gamma));
a_3 = sqrt((gamma*p_3)/rho_3);
a_2 = sqrt((gamma*p_2)/rho_2);

conpos = xd+up*time;
spos = xd+w*time;
pos_1 = xd+(u_1-a_1)*time;
pos_2 = xd+(u_2-a_2)*time;

pos = pos_1:(pos_2-pos_1)/20:pos_2;
x = [0,pos_1,pos,pos_2,conpos,conpos,spos,spos,L,L];

ue = (2/(gamma+1))*((0.5*(gamma-1)*u_1)+a_1+((pos-xd)/time));
pe = p_1*(1-(0.5*(gamma-1))*((ue-u_1)./a_1)).^((2*gamma)/(gamma-1));
re = rho_1*(1-(0.5*(gamma-1))*((ue-u_1)./a_1)).^(2/(gamma-1));
P = [p_1,p_1,pe,p_2,p_2,p_3,p_3,p_4,p_4,0];
U = [u_1,u_1,ue,u_2,u_2,u_3,u_3,u_4,u_4,0];
R = [rho_1,rho_1,re,rho_2,rho_2,rho_3,rho_3,rho_4,rho_4,0];
H = (gamma/(gamma-1))*(P./R)+0.5*(U.^2);

fig = figure(1);
plot(x,R,'--',X,r);
title('density');
legend('exact','Roe & Lax');
saveas(fig,['RoeLax_r',num2str(cond),'.jpg'],'jpeg');

fig = figure(2);
plot(x,U,'--',X,u);
title('velocity');
legend('exact','Roe & Lax');
saveas(fig,['RoeLax_v',num2str(cond),'.jpg'],'jpeg');

fig = figure(3);
plot(x,P,'--',X,p);
title('pressure');
legend('exact','Roe & Lax');
saveas(fig,['RoeLax_p',num2str(cond),'.jpg'],'jpeg');

fig = figure(4);
plot(x,H,'--',X,h);
title('enthalpy');
legend('exact','Roe & Lax');
saveas(fig,['RoeLax_h',num2str(cond),'.jpg'],'jpeg');

close all

%================================ Jamson ==================================

clear,clc

gamma = 1.4;
CFL = 0.3;
L = 1;
N = 100;
xd = 0.5;
time = 0;
FinalTime = 0.5;
FinalStep = 100;
step = 1;
k2 = 1/2;
k4 = 1/32;

dx = L/N;
Nd = fix(xd/dx)+1;
X = 0:dx:L;

cond = 1; % sod : cond = 0; Lax : cond = 1

if (cond==0)
    Rho_Left = 1.0;
    Rho_Right = 0.125;
    Velocity_Left = 0;
    Velocity_Right = 0;
    Pressure_Left = 1;
    Pressure_Right = 0.1;
    
else
    Rho_Left = .445;
    Rho_Right = 0.5;
    Velocity_Left = .698;
    Velocity_Right = 0;
    Pressure_Left = 3.528;
    Pressure_Right = 0.571;
end

r(1:Nd) = Rho_Left;
r(Nd+1:N+1) = Rho_Right;
u(1:Nd) = Velocity_Left;
u(Nd+1:N+1) = Velocity_Right;
p(1:Nd) = Pressure_Left;
p(Nd+1:N+1) = Pressure_Right;

q = zeros(N+1,3);
ua = zeros(N+1,1);
v = zeros(3,N);
e2 = zeros(3,N);
e4 = zeros(3,N);
d = zeros(3,N);
d2 = zeros(3,N);
d4 = zeros(3,N);
F = zeros(3,N);
lambda = zeros(3,N);

while ( time<FinalTime && step<FinalStep )
    
    for i = 1:N+1
        ua(i) = max(abs(u(i))+sqrt((gamma*p(i))/r(i)));
    end
    dt = (dx*CFL)/max(ua);
    
    h = (gamma/(gamma-1))*(p./r)+0.5*(u.^2);
    
    for i = 2:N
        v(i) = abs(p(i+1)-2*p(i)-p(i-1))/abs(p(i+1)+2*p(i)+p(i-1));
    end
    v(1) = 0;v(N+1) = 0;
    
    for i = 1:N
        e2(i) = k2*max([v(i),v(i+1)]);
        e4(i) = max([0,k4-e2(i)]);
        lambda(i) = 0.5*(abs(u(i))+sqrt((gamma*p(i))/r(i))+abs(u(i+1))+sqrt((gamma*p(i+1))/r(i+1)));
    end
    
    for i = 1:N
        q0 = [r(i);r(i)*u(i);(p(i)/(gamma-1))+0.5*r(i)*(u(i)^2)];
        q1 = [r(i+1);r(i+1)*u(i+1);(p(i+1)/(gamma-1))+0.5*r(i+1)*(u(i+1)^2)];
        d2(:,i) = e2(i)*lambda(i)*(q1-q0);
    end
    
    for i = 1:N
        if (i==1)
            q0 = [r(i);r(i)*u(i);(p(i)/(gamma-1))+0.5*r(i)*(u(i)^2)];
            q3 = [r(i+2);r(i+2)*u(i+2);(p(i+2)/(gamma-1))+0.5*r(i+2)*(u(i+2)^2)];
        elseif (i==N)
            q0 = [r(i-1);r(i-1)*u(i-1);(p(i-1)/(gamma-1))+0.5*r(i-1)*(u(i-1)^2)];
            q3 = [r(i+1);r(i+1)*u(i+1);(p(i+1)/(gamma-1))+0.5*r(i+1)*(u(i+1)^2)];
        else
            q0 = [r(i-1);r(i-1)*u(i-1);(p(i-1)/(gamma-1))+0.5*r(i-1)*(u(i-1)^2)];
            q3 = [r(i+2);r(i+2)*u(i+2);(p(i+2)/(gamma-1))+0.5*r(i+2)*(u(i+2)^2)];
        end
        q1 = [r(i);r(i)*u(i);(p(i)/(gamma-1))+0.5*r(i)*(u(i)^2)];
        q2 = [r(i+1);r(i+1)*u(i+1);(p(i+1)/(gamma-1))+0.5*r(i+1)*(u(i+1)^2)];
        d4(:,i) = -e4(i)*lambda(i)*(q3-3*q2+3*q1-q0);
    end
    
    for i = 1:N
        d(:,i) = d2(:,i)+d4(:,i);
        FL = [r(i)*u(i);p(i)+r(i)*(u(i)^2);r(i)*u(i)*h(i)];
        FR = [r(i+1)*u(i+1);p(i+1)+r(i+1)*(u(i+1)^2);r(i+1)*u(i+1)*h(i+1)];
        F(:,i) = 0.5*(FL+FR)-d(:,i);
    end
    
    for i = 2:N
        Q(1,1) = r(i);
        Q(2,1) = r(i)*u(i);
        Q(3,1) = (p(i)/(gamma-1))+0.5*r(i)*(u(i)^2);
        Q = Q-(dt/dx)*(F(:,i)-F(:,i-1));
        q(i,1:3) = Q;
    end
    
    r = q(:,1);
    u = q(:,2)./q(:,1);
    p = (gamma-1)*(q(:,3)-0.5*r.*u.^2);
    
    r(1,:) = r(2,:);u(1,:) = u(2,:);p(1,:) = p(2,:);
    r(N+1,:) = r(N,:);u(N+1,:) = u(N,:);p(N+1,:) = p(N,:);
    
    time = time+dt;
    step = step+1;
    
end

p_1 = Pressure_Left;
rho_1 = Rho_Left;
u_1 = Velocity_Left;

p_4 = Pressure_Right;
rho_4 = Rho_Right;
u_4 = Velocity_Right;

a_1 = sqrt((gamma*p_1)/rho_1);
a_4 = sqrt((gamma*p_4)/rho_4);

p_14 = p_1/p_4;
p_34 = 0.00001;

beta = (1+(0.5/a_1)*(u_1-u_4)*(gamma-1)-(gamma-1)*(a_4/a_1)*(p_34-1)*((2*gamma)*((gamma-1)+(gamma+1)*p_34))^(-0.5))^((2*gamma)/(1-gamma));
while ((p_14-p_34*beta)>0.001)
    p_34 = p_34+0.00001;
    beta = (1+(0.5/a_1)*(u_1-u_4)*(gamma-1)-(gamma-1)*(a_4/a_1)*(p_34-1)*((2*gamma)*((gamma-1)+(gamma+1)*p_34))^(-0.5))^((2*gamma)/(1-gamma));
end

p_3=p_4*p_34;
w = (a_4*sqrt(1+((1+gamma)/(2*gamma))*(p_34-1)))+u_4;
up = u_4+(a_4/gamma)*(p_34-1)*sqrt(((2*gamma)/(1+gamma))/(p_34-((1-gamma)/(1+gamma))));
u_2 = up;
u_3 = up;
p_2 = p_3;
rho_3 = rho_4/((((gamma+1)/(gamma-1))+p_34)/(1+p_34*((gamma+1)/(gamma-1))));
rho_2 = rho_1*((p_2/p_1)^(1/gamma));
a_3 = sqrt((gamma*p_3)/rho_3);
a_2 = sqrt((gamma*p_2)/rho_2);

conpos = xd+up*time;
spos = xd+w*time;
pos_1 = xd+(u_1-a_1)*time;
pos_2 = xd+(u_2-a_2)*time;

pos = pos_1:(pos_2-pos_1)/20:pos_2;
x = [0,pos_1,pos,pos_2,conpos,conpos,spos,spos,L,L];

ue = (2/(gamma+1))*((0.5*(gamma-1)*u_1)+a_1+((pos-xd)/time));
pe = p_1*(1-(0.5*(gamma-1))*((ue-u_1)./a_1)).^((2*gamma)/(gamma-1));
re = rho_1*(1-(0.5*(gamma-1))*((ue-u_1)./a_1)).^(2/(gamma-1));
P = [p_1,p_1,pe,p_2,p_2,p_3,p_3,p_4,p_4,0];
U = [u_1,u_1,ue,u_2,u_2,u_3,u_3,u_4,u_4,0];
R = [rho_1,rho_1,re,rho_2,rho_2,rho_3,rho_3,rho_4,rho_4,0];
H = (gamma/(gamma-1))*(P./R)+0.5*(U.^2);

fig = figure(1);
plot(x,R,'--',X,r);
title('density');
legend('exact','Jamson');
saveas(fig,['Jamson_r',num2str(cond),'.jpg'],'jpeg');

fig = figure(2);
plot(x,U,'--',X,u);
title('velocity');
legend('exact','Jamson');
saveas(fig,['Jamson_v',num2str(cond),'.jpg'],'jpeg');

fig = figure(3);
plot(x,P,'--',X,p);
title('pressure');
legend('exact','Jamson');
saveas(fig,['Jamson_p',num2str(cond),'.jpg'],'jpeg');

fig = figure(4);
plot(x,H,'--',X,h);
title('enthalpy');
legend('exact','Jamson');
saveas(fig,['Jamson_h',num2str(cond),'.jpg'],'jpeg');

close all








