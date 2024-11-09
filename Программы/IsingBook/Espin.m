function z = Espin(Nspin,s,h,J)
% �������, ������������ ������� 
% �������� ������������ ������ 
% Nspin - ����� ������ �������
% s - �������, ���������� ���������� 
% � ��������� ������
% h - ������������� �������� ���������� ����
% J - ��������� ��������� ��������������
% Copyright 2003-2009 S. V. Porshnev

Ns=Nspin.^0.5;
E=0;
for i=1:Ns
   for j=1:Ns
     % �������� ������������� ��������� �������
     if j==1
        Left=Ns;
     else
        Left==j-1;  
     end
     if j==Ns
        Right=1;
     else
        Right=j+1;
     end
     if i==Ns
       Up=1;
     else
       Up=i+1;
     end
     if i==1
        Down=Ns;
     else
        Down=i-1;
     end
     % ������� �������
     E=E-s(i,j)*(h+J*(s(i,Left)+s(i,Right)+s(Down,j)+s(Up,j)));
   end
end
z=E;
