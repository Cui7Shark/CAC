function dudt=transmem(t,u)
global a a1 K Res Rus Cap Muli;

T = 300;
k_Boltzmann = 1.38e-23;
q = 2.46;
qe = 1.65e-19;
X = qe/(k_Boltzmann*T);%%95.0725
rm=0.80e-9;
sigma= 1.2;
dm=5e-9;
w0=2.65;
yita=0.15;
N0=1.5e9;
Vep=0.258;
Th = 10;
alpha = 1e9;
dudt = zeros(a+a1, 1);
%%100ns   10kv/cm  plus
Vapp=  0.*(t<0)+(2*1e10.*t).*(t<10e-9 & t>=0) + 2e2.*(t>=10e-9 & t<=90e-9) + (-2e10.*t+2e3).*(t>90e-9 & t<100e-9)+ 0.*(t>=100e-9);
disp(t)
%%
for i = 1:(a+a1)
    [row,col]=find(K(:,1:2)==i);                                           %%判断点i是否为膜上节点
    if size(col,1)==0 && i<=a                                              %%非膜上节点的电流
        dudt(i) = Muli*Res(i,:)*u(1:a) + Muli*Rus(i)*Vapp;                           %求电流
    elseif  size(col,1)==1 && i <= a
        Vm=@(u) u(K(row,2))-u(K(row,1));                                   %外电势-内电势
        vm=@(u) Vm(u)*X;                                                   %X=q*e/(k_Boltzmann*T);%%95.0725

        Gp=@(u) (sigma*pi*rm*rm/dm) * (exp(vm(u))-1)/((w0*exp(w0-yita*vm(u))-yita*vm(u))*exp(vm(u))/(w0-yita*vm(u))- ...
            ((w0*exp(w0+yita*vm(u))+yita*vm(u))/(w0+yita*vm(u))));
     
%%计算跨膜电流
        %细胞膜内节点
        if col==1 
            dudt(i) = Muli*Res(i,:)*u(1:a) - Muli*Gp(u)*u(row+a)*Vm(u)*K(row,3)*Th*1e-12;      %1e-6*1e-6*;
        %细胞膜外-节点
        elseif col==2 
            dudt(i) = Muli* Res(i,:)*u(1:a) + Muli*Gp(u)*u(row+a)*Vm(u)*K(row,3)*Th*1e-12;      %1e-6*1e-6;
        end
    end
%%计算孔隙密度
    if  i > a
        Vm1=@(u) (u(K(i-a,2))-u(K(i-a,1)));
        dudt(i) = alpha * (exp((Vm1(u)/Vep)^2) * (1-((u(i)/N0)*exp(-q*(Vm1(u)/Vep)^2))));
    end
end
