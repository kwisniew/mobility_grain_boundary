 fileID = fopen('mobility_100K.txt','r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 m_100K = fscanf(fileID,formatSpec,sizeA);
 fclose(fileID);
% B=A'
 loglog(m_100K(1,:),m_100K(2,:))
 
 fileID = fopen('mobility_150K.txt','r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 m_150K = fscanf(fileID,formatSpec,sizeA);
 fclose(fileID);
 loglog(m_150K(1,:),m_150K(2,:))
 
  
 fileID = fopen('mobility_200K.txt','r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 m_200K = fscanf(fileID,formatSpec,sizeA);
 fclose(fileID);
 loglog(m_200K(1,:),m_200K(2,:))
  
 fileID = fopen('mobility_250K.txt','r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 m_250K = fscanf(fileID,formatSpec,sizeA);
 fclose(fileID);
 loglog(m_250K(1,:),m_250K(2,:))
  
 fileID = fopen('mobility_300K.txt','r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 m_300K = fscanf(fileID,formatSpec,sizeA);
 fclose(fileID);
 loglog(m_300K(1,:),m_300K(2,:))
 
 max_Na = zeros(1,5);
 x = size(m_100K);
 max_Na(1) = x(2);
 x = size(m_150K);
 max_Na(2) = x(2);
 x = size(m_200K);
 max_Na(3) = x(2);
  x = size(m_250K);
 max_Na(4) = x(2);
  x = size(m_300K);
 max_Na(5) = x(2)
 
 max_Na_all = min(max_Na)
 all_mobilities = zeros(5,max_Na_all);
 all_mobilities(1,1:max_Na_all) = m_100K(2,1:max_Na_all);
  all_mobilities(2,1:max_Na_all) = m_150K(2,1:max_Na_all);
   all_mobilities(3,1:max_Na_all) = m_200K(2,1:max_Na_all);
    all_mobilities(4,1:max_Na_all) = m_250K(2,1:max_Na_all);
     all_mobilities(5,1:max_Na_all) = m_300K(2,1:max_Na_all);
 
Ea = zeros(2,max_Na_all);
T = [100,150,200,250,300]

plot(1./T, log(all_mobilities(:,30))', 'b--*','MarkerSize',10)
c = polyfit(1./T,log( all_mobilities(:,30) )',1);
% c(1)*1./T(1)+c(2)
y_est = polyval(c,1./T);
hold on
plot(1./T,y_est,'r--','LineWidth',2)
hold off


 for i=1:max_Na_all
     Ea(:,i) = polyfit(1./T,log( all_mobilities(:,i)' ),1);
 end
 plot(1./T, log(all_mobilities(:,30))', 'b--*','MarkerSize',10)
 y_est = polyval(Ea(:,30),1./T);
 hold on
 plot(1./T,y_est,'r--','LineWidth',2)
 hold off
 
 
 q=1.6021766208*10^(-19); %C
 k = 1.38064852*10^(-23); %J/K
 Ea(1,:) = -Ea(1,:)*k/q;
 semilogx(m_100K(1,:),Ea(1,:))
 semilogx(m_100K(1,:),Ea(2,:))
 
 