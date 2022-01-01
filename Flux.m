function [ F,R,alpha,lambda ] = Flux( r,u,h,a,rr,rl,ur,ul,pr,pl,hr,hl )

dr = rr-rl;
du = ur-ul;
dp = pr-pl;

R(1,1) = 1;
R(2,1) = u-a; 
R(3,1) = h-u*a;
R(1,2) = 1;
R(2,2) = u; 
R(3,2) = 0.5*(u^2);
R(1,3) = 1;
R(2,3) = u+a; 
R(3,3) = h+u*a;

alpha(1,1) = (0.5/(a^2))*(dp-a*r*du);
alpha(2,1) = dr-(dp/(a^2));
alpha(3,1) = (0.5/(a^2))*(dp+a*r*du);

lambda(1,1) = abs(u-a);
lambda(2,1) = abs(u);
lambda(3,1) = abs(u+a);

F(1,1) = (rr*ur)+(rl*ul);
F(2,1) = (rr*(ur^2)+pr)+(rl*(ul^2)+pl);
F(3,1) = (rr*ur*hr)+(rl*ul*hl);
F = 0.5*F;

end

