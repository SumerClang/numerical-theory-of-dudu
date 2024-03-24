function KKo = stiffness_element_2d_3node_assemble(KKi,Kei,elementi)
%该函数进行单元刚度矩阵的组装
%输入单元刚度矩阵Kei 
%输入单元的节点编号elementi
%输出整体刚度矩阵 KKo 
KKo=KKi;
DOF(1)=2*elementi(1)-1; 
DOF(2)=2*elementi(1); 
DOF(3)=2*elementi(2)-1; 
DOF(4)=2*elementi(2); 
DOF(5)=2*elementi(3)-1; 
DOF(6)=2*elementi(3); 
for n1=1:6 
    for n2=1:6 
        KKo(DOF(n1),DOF(n2))= KKo(DOF(n1),DOF(n2))+Kei(n1,n2); 
    end 
end
end