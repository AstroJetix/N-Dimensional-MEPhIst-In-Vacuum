function [V,Es,Accept] = IdealGas2(Vel,T,dVmax,NTrial,Np)
% �������, ������������
% ���������� �������� ������� �������� (V), 
% ������� (E), 
% ����� �������� (Accept)
% Vel - ��������� ��������
% T - ����������� ����
% dVmax - ������������ ��������� ��������
% NTrial - ����� ���������
% Np - ����� ������
% Copyright 2003-2009 S. V. Porshnev

beta=1/T;
i=1:Np;
v(i)=Vel;
Eold=dot(v,v)/2;
V(1)=mean(v);
Es=Eold;
Accept=0;
k=2;
for i=1:NTrial
  for j=1:Np
     N=floor(Np*rand(1)+1); % ���������
                            % ����� �������
    dV=(2*rand(1)-1)*dVmax; % ��������� ��������� 
                            % ��������
     Vtemp=v(j);
     v(j)=v(j)+dV;
     Enew=0.5*dot(v,v);
     de=Enew-Eold;% ������� 
                  % ��������� �������
    if de>0
       % ��� ����������� 
       if exp(-beta*de)>rand(1)
         Accept=Accept+1;
         Eold=Enew;
       else
         % ��� �� �����������  
         v(j)=Vtemp;
       end
    else
       Accept=Accept+1;
       Eold=Enew;
    end
    V(k)=mean(v); 
    Es(k)=Eold/Np;
    k=k+1;
  end
end
Accept=Accept/(NTrial*Np);