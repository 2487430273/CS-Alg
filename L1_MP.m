%数据清除
clear all;
clc;

%  生成欠采样矩阵XOmmiga；
X_1 = rand(11,3);
X_2 = rand(3,11);
X = X_1 * X_2;                %生成秩为1的4*6的矩阵X
sparsity =0.9;               %Ommig  稀疏度（矩阵中的非零元素占比）
K = rank(X);
[m,n] = size(X);
Ommiga = generateSparseBinaryMatrix(m, n, sparsity);
XOmmiga = X .* Ommiga;       %生成欠采样矩阵XOmmiga
X1 = zeros(size(XOmmiga));
x0=0;
x1=1;
sigma = 1 ;
elta = 1;
u_1 = X1(:,1)';                      %
v_1 = X1(1,:);
result_1 = [];
result_2 = [];

%  初始化
R = XOmmiga;
X2 = rand(size(XOmmiga));
v0 = (X2(1,:))';          %v向量随机化

% 计算
for k = 1:K
    v = v0;
    if ((sigma < 1*10^(-5)) && (elta <= 5*10^(-4)))  %uk和vk迭代完成的条件
        break;
    else
        E1=((norm((R-(u_1)'*v_1),1))/(norm(R,1)));
        for q = 1:10                        %迭代次数
            for i = 1:m
                r_i = R(i,:);                    %定义R的行向量r(i)
                I_i = Ommiga(i,:);
                a_i = r_i;                            %
                b = v .* ((I_i)');                    %
                fun = @(x)norm(((a_i)' - x*(b)),1);
                %                 u(i)=fminsearch(fun,x0);            %文中式10
                u(i)=fminbnd(fun,x0,x1);
                %                 u(i)=fminunc(fun,x0);
                u_1(i) = u(i);
            end
            %根据u  优化v
            for j = 1:n
                r_j = R(:,j);
                J_j = Ommiga(:,j);
                c_j = r_j;
                d = u_1 .* J_j;
                fun = @(x)norm((c_j - x*d),1);
                %                 v(j)=fminsearch(fun,x0);   %文中式14 fminsearch 函数用于无约束的单目标函数最小化问题，它使用了 Nelder-Mead 搜索算法。
                v(j)=fminbnd(fun,x0,x1);%函数用于在给定的区间内求解单变量函数的最小值点。它使用了一种被称为黄金分割搜索的算法。
                %                 v(j)=fminunc(fun,x0);%fminunc
                %                 函数用于无约束的多目标函数最小化问题，它使用了拟牛顿法或共轭梯度法等优化算法  在 MATLAB 中，用于解决有约束的函数最小化问题的函数包括 fmincon 和 fminimax。
                v_1(j) = v(j);
            end
        end
        E2 = ((norm((R-(u_1)'*v_1),1))/(norm(R,1)));
        Elta_1 = ((norm(R,1))/(norm(XOmmiga,1)));
        sigma = E1 -E2;
        R = R - (((u_1)'*v_1) .* Ommiga);      %计算R（k+1）
        Elta_2 = ((norm(R,1))/(norm(XOmmiga,1)));
        elta = Elta_1 - Elta_2;
    end
    result_1 = [result_1;u];
    result_2 = [result_2,v];
end
u_2 = result_1;
v_2 = result_2;
%输出
M = u_2'*v_2'
X