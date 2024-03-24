function Stresso = stress_element(Ysi,Possi,xin,yin,qin)
%   Ysi—杨氏模量；
%   Possi—泊松比；
%   thick—单元厚度，仅用于平面应力情况下；
%   xin—存放单元节点的x坐标
%   yin—存放单元节点的y坐标
%   qin—存放单元节点的x位移
%单元面积为A
A=0.05*0.05/2;
%生成弹性矩阵D
at=[xin(2)*yin(3)-xin(3)*yin(2) xin(3)*yin(1)-xin(1)*yin(3) xin(1)*yin(2)-xin(2)*yin(1)];
bt=[yin(2)-yin(3) yin(3)-yin(1) yin(1)-yin(2)];
ct=[xin(3)-xin(2) xin(1)-xin(3) xin(2)-xin(1)];
ii=1;jj=2;mm=3;
sii=Ysi/(2*(1-Possi*Possi)*A)*[bt(ii) Possi*ct(ii);Possi*bt(ii) ct(ii);(1-Possi)/2*ct(ii) (1-Possi)/2*bt(ii)];
sjj=Ysi/(2*(1-Possi*Possi)*A)*[bt(jj) Possi*ct(jj);Possi*bt(jj) ct(jj);(1-Possi)/2*ct(jj) (1-Possi)/2*bt(jj)];
smm=Ysi/(2*(1-Possi*Possi)*A)*[bt(mm) Possi*ct(mm);Possi*bt(mm) ct(mm);(1-Possi)/2*ct(mm) (1-Possi)/2*bt(mm)];
St=[sii sjj smm];
Stresso=St*qin;
end

