function Keo = stiffness_element_2d_3node(Ysi,Possi,thicki,xin,yin)
%   Ysi—杨氏模量；
%   Possi—泊松比；
%   thick—单元厚度，仅用于平面应力情况下；
%   xin—存放单元节点的x坐标
%   yin—存放单元节点的y坐标
%   输出刚度矩阵keo(6*6)
%单元面积为A
A=0.05*0.05/2;
%生成弹性矩阵D
D=Ysi/(1-Possi*Possi)*[1 Possi 0; Possi 1 0; 0 0 (1-Possi)/2];
at=[xin(2)*yin(3)-xin(3)*yin(2) xin(3)*yin(1)-xin(1)*yin(3) xin(1)*yin(2)-xin(2)*yin(1)];
bt=[yin(2)-yin(3) yin(3)-yin(1) yin(1)-yin(2)];
ct=[xin(3)-xin(2) xin(1)-xin(3) xin(2)-xin(1)];
%Keij存储单元刚度矩阵中的分块矩阵
Keijt=zeros(2,2);
%Keo存储单元刚度矩阵
Keo=zeros(6,6);
for ii=1:3
   for jj=1:3 
       bii=1/(2*A)*[bt(ii) 0; 0 ct(ii); ct(ii) bt(ii)];
       bjj=1/(2*A)*[bt(jj) 0; 0 ct(jj); ct(jj) bt(jj)];
       Keij=((bii)')*D*bjj*thicki*A;
       Keo(ii*2-1:ii*2,jj*2-1:jj*2)=Keij;
   end
end
end