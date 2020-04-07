clc,close all,clear all;
x=xlsread('data1.xlsx', 'A1:B5401');
xx = xlsread('data2.xlsx','A1:B1081')
temperature45list = x(:,2);%读取皮肤外侧数据，作为模型另一边的温度情况
temperature01list = xx(:,2);%读取I层靠近环境数据，作为模型另一边的温度情况
numberOfPDE = 1;
pdem1 = createpde(numberOfPDE);
numberOfPDE = 1;
pdem2 = createpde(numberOfPDE);
numberOfPDE = 1;
pdem3 = createpde(numberOfPDE);
numberOfPDE = 1;
pdem4 = createpde(numberOfPDE);
rho1 = 300; % 密度
Cp1 = 1377; % 比热
k1 = 0.082; % 热导率
h1 = 0.6;

rho2 = 862; % 密度
Cp2 = 2100; % 比热
k2 = 0.37; % 热导率
h2 = 6;

rho3 = 74.2; % 密度
Cp3 = 1726; % 比热
k3 = 0.045; % 热导率
h3 = 3.6;

rho4 = 1.18; % 密度
Cp4 = 1005; % 比热
k4 = 0.028; % 热导率
h4 = 5;

R1 = [3 4 0 h1+h2+h3+h4 h1+h2+h3+h4 0  0 0 100 100]';
R2 = [3 4 h1 h1+h2+h3+h4 h1+h2+h3+h4 h1  0 0 100 100]';
R3 = [3 4 h1+h2 h1+h2+h3+h4 h1+h2+h3+h4 h1+h2  0 0 100 100]';
R4 = [3 4 h1+h2+h3 h1+h2+h3+h4 h1+h2+h3+h4 h1+h2+h3  0 0 100 100]';

g1 = decsg(R1);
geometryFromEdges(pdem1 ,g1);
%figure(1)
%pdegplot(pdem1,'EdgeLabels','on'); 
%axis([-.1 15.3 49 51]);

g2 = decsg(R2);
geometryFromEdges(pdem2 ,g2);
%figure(2)
%pdegplot(pdem2,'EdgeLabels','on'); 
%axis([-.1 15.3 5 6]);

g3 = decsg(R3);
geometryFromEdges(pdem3 ,g3);
%figure(3)
%pdegplot(pdem3,'EdgeLabels','on'); 
%axis([-.1 15.3 5 6]);

g4 = decsg(R4);
geometryFromEdges(pdem4 ,g4);
%figure(4)
%pdegplot(pdem4,'EdgeLabels','on'); 
%axis([-.1 15.3 5 6]);


%虽然希望三角网格的划分越精细越好，但是由于上下边界无穷远的设定，导致这部分无用网格太多，
%所以注意大小
hmax = .2; % element size
msh1=generateMesh(pdem1,'Hmax',hmax);
%figure
%pdeplot(pdem1); 
%axis([-.1 .7 5 6]);

hmax = .5; % element size
msh2=generateMesh(pdem2,'Hmax',hmax);
%figure
%pdeplot(pdem2); 
%axis([.5 6.7 2 8]);

hmax = .5; % element size
msh3=generateMesh(pdem3,'Hmax',hmax);
%figure
%pdeplot(pdem3); 
%axis([6.5 10.3 3 7]);

hmax = .5; % element size
msh4=generateMesh(pdem4,'Hmax',hmax);
%figure
%pdeplot(pdem4); 
%axis([10.1 15.3 2 8]);



temperature = 75%环境温度已经在另一个模型中应用
temperature01 = 37
temperature12 = 37
temperature23 = 37
temperature34 = 37
time = 1;
temperature12list = ones(1,1200)
temperature23list = ones(1,1200)
temperature34list = ones(1,1200)

for i2 = 1:1100

	if (temperature01list(time+1) - temperature01list(time) >0.001)
		%第一层
		temperature01 = temperature01list(time,1);
		temperature45 = temperature45list(time*5,1);%每隔5秒作为一次迭代，所以选取5秒后的皮肤外侧温度
		applyBoundaryCondition(pdem1,'dirichlet','Edge',4,'u',temperature01);
		applyBoundaryCondition(pdem1,'dirichlet','Edge',2,'u',temperature45);
		applyBoundaryCondition(pdem1,'dirichlet','Edge',[1,3],'u',37);
		c = k1;
		a = 0;
		f = 0;
		d = rho1*Cp1; 
		specifyCoefficients(pdem1,'m',0,'d',0,'c',c,'a',a,'f',f);
		tlist = 0:.1:5;
		setInitialConditions(pdem1, 0);
		R = solvepde(pdem1,tlist);
		u = R.NodalSolution;%反解每个节点的温度
		p = pdem1.Mesh.Nodes;%得出每个mesh后的节点坐标
		x = p(1,:);
		y = p(2,:);
		num1 = 1;
		%figure 
		%pdeplot(pdem1,'XYData',u(:,1),'Contour','off','ColorMap','hot'); 
		%axis([-.1 15 49 51]);
		nodechoose =zeros(1,1000)%选取的作为该界面的坐标值
		num2 = length(x)
		for i = 1:num2
			if ((x(i)>0.5) && (x(i)<0.7) && (y(i)>35) && (y(i) <65))
				nodechoose(num1) = i;
				num1 = num1 + 1;
			end
		end
		num3 = 0;
		temperature12 = 0
		nodechoose(1)
		for i = 1:1000
			if((nodechoose(i)) ~=0)
				temperature12 = temperature12 + u(nodechoose(i));
				num3 = num3 +1;
			end
		end
		temperature12 = temperature12/num3;
		temperature12list(time) = temperature12;

		%第二层
		applyBoundaryCondition(pdem2,'dirichlet','Edge',4,'u',temperature12);
		applyBoundaryCondition(pdem2,'dirichlet','Edge',2,'u',temperature45);
		applyBoundaryCondition(pdem2,'dirichlet','Edge',[1,3],'u',37);
		c = k2;
		a = 0;
		f = 0;
		d = rho2*Cp2; 
		specifyCoefficients(pdem2,'m',0,'d',0,'c',c,'a',a,'f',f);
		tlist = 0:.1:5;
		setInitialConditions(pdem2, 0);
		R = solvepde(pdem2,tlist);
		u = R.NodalSolution;%反解每个节点的温度
		p = pdem2.Mesh.Nodes;%得出每个mesh后的节点坐标
		x = p(1,:);
		y = p(2,:);
		num1 = 1;
		%figure 
		%pdeplot(pdem2,'XYData',u(:,1),'Contour','off','ColorMap','hot'); 
		%axis([-.1 15 49 51]);
		nodechoose =zeros(1,1000)%选取的作为该界面的坐标值
		num2 = length(x)
		for i = 1:num2
			if ((x(i)>6) && (x(i)<7.5) && (y(i)>35) && (y(i) <65))
				nodechoose(num1) = i;
				num1 = num1 + 1;
			end
		end
		num3 = 0;
		temperature23 = 0
		nodechoose(1)
		for i = 1:1000
			if((nodechoose(i)) ~=0)
				temperature23 = temperature23 + u(nodechoose(i));
				num3 = num3 +1;
			end
		end
		temperature23 = temperature23/num3;
		temperature23list(time) = temperature23;
	%temperature23

	%第三层
		applyBoundaryCondition(pdem3,'dirichlet','Edge',4,'u',temperature23);
		applyBoundaryCondition(pdem3,'dirichlet','Edge',2,'u',temperature45);
		applyBoundaryCondition(pdem3,'dirichlet','Edge',[1,3],'u',37);
		c = k3;
		a = 0;
		f = 0;
		d = rho3*Cp3; 
		specifyCoefficients(pdem3,'m',0,'d',0,'c',c,'a',a,'f',f);
		tlist = 0:.1:5;
		setInitialConditions(pdem3, 0);
		R = solvepde(pdem3,tlist);
		u = R.NodalSolution;%反解每个节点的温度
		p = pdem3.Mesh.Nodes;%得出每个mesh后的节点坐标
		x = p(1,:);
		y = p(2,:);
		num1 = 1;
	%figure 
	%pdeplot(pdem3,'XYData',u(:,1),'Contour','off','ColorMap','hot'); 
	%axis([-.1 15 49 51]);
		nodechoose =zeros(1,1000)%选取的作为该界面的坐标值
		num2 = length(x)
		for i = 1:num2
			if ((x(i)>9.5) && (x(i)<11) && (y(i)>35) && (y(i) <65))
				nodechoose(num1) = i;
				num1 = num1 + 1;
			end
		end
		num3 = 0;
		temperature34 = 0
		nodechoose(1)
		for i = 1:1000
			if((nodechoose(i)) ~=0)
				temperature34 = temperature34 + u(nodechoose(i));
				num3 = num3 +1;
			end
		end
		temperature34 = temperature34/num3;
		temperature34list(time) = temperature34
	%temperature34

	%第四层
		applyBoundaryCondition(pdem4,'dirichlet','Edge',4,'u',temperature34);
		applyBoundaryCondition(pdem4,'dirichlet','Edge',2,'u',temperature45);
		applyBoundaryCondition(pdem4,'dirichlet','Edge',[1,3],'u',37);
		c = k4;
		a = 0;
		f = 0;
		d = rho4*Cp4; 
		specifyCoefficients(pdem4,'m',0,'d',0,'c',c,'a',a,'f',f);
		tlist = 0:.1:5;
		setInitialConditions(pdem4, 0);
		R = solvepde(pdem4,tlist);
		u = R.NodalSolution;%反解每个节点的温度
		p = pdem4.Mesh.Nodes;%得出每个mesh后的节点坐标
		x = p(1,:);
		y = p(2,:);
		num1 = 1;
	%figure 
	%pdeplot(pdem4,'XYData',u(:,1),'Contour','off','ColorMap','hot'); 
	%axis([-.1 15 49 51]);
		nodechoose =zeros(1,1000)%选取的作为该界面的坐标值
		num2 = length(x)
		for i = 1:num2
			if ((x(i)>14.5) && (x(i)<15.2) && (y(i)>35) && (y(i) <65))
				nodechoose(num1) = i;
				num1 = num1 + 1;
			end
		end
		num3 = 0;
		temperature45 = 0
		nodechoose(1)
		for i = 1:1000
			if((nodechoose(i)) ~=0)
				temperature45 = temperature45 + u(nodechoose(i));
				num3 = num3 +1;
			end
		end
		temperature45 = temperature45/num3;
		%temperature45


		time = time+1

	
	else
		i2 = 1100;
	end
end














































