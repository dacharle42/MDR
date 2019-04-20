clearvars; close all;

% load results from Dizzy output: 'rxy_z' (tab-delimited text)
% where x=p01NF y=p10NF z=p01NF/p01PF=p10NF/p10PF
r=load('r21_20.txt');

t=r(:,1);   %time(h)
ANF=r(:,2); % active NF promoter
APF=r(:,3); % active PF promoter
INF=r(:,4); % inactive NF promoter
IPF=r(:,5); % inactive PF promoter
MNF=r(:,6); % NF mRNA
MPF=r(:,7); % PF mRNA
P1NF=r(:,8); % NF protein1
P1PF=r(:,9); % NF protein2
P2NF=r(:,10); % PF protein1
P2PF=r(:,11); % PF protein2

% plot time courses
figure;hold on;
plot(t,P1NF,'c','LineWidth',2);
plot(t,P2NF,'b','LineWidth',2);
plot(t,P1PF,'m','LineWidth',2);
plot(t,P2PF,'r','LineWidth',2);
xlabel('time (h)','FontSize',20);
ylabel('protein copy numbers','FontSize',20);
set(gca,'FontSize',20,'XLim',[0 t(end)+1]);
legend('P1_{NF}','P2_{NF}','P1_{PF}','P2_{PF}','Orientation','Horizontal','Location','North');

%plot histograms
binEdges=linspace(0,4000,100);
figure;hold on;
h1NF=histogram(P1NF,binEdges);
h2NF=histogram(P2NF,binEdges);
h1PF=histogram(P1PF,binEdges);
h2PF=histogram(P2PF,binEdges);
set(h1NF,'FaceColor','c');
set(h2NF,'FaceColor','b');
set(h1PF,'FaceColor','m');
set(h2PF,'FaceColor','r');
xlabel('protein copy numbers','FontSize',20);
ylabel('cell counts','FontSize',20);
set(gca,'FontSize',20,'XLim',[0 4000]);
legend('P1_{NF}','P2_{NF}','P1_{PF}','P2_{PF}');

% calculate statistics
m1NF=mean(P1NF);
m2NF=mean(P2NF);
m1PF=mean(P1PF);
m2PF=mean(P2PF);

s1NF=std(P1NF);
s2NF=std(P2NF);
s1PF=std(P1PF);
s2PF=std(P2PF);

c1NF=s1NF/m1NF;
c2NF=s2NF/m2NF;
c1PF=s1PF/m1PF;
c2PF=s2PF/m2PF;

% plot CVs versus means
figure;hold on;
plot(m1NF,c1NF,'c>','LineWidth',2,'MarkerSize',10);
plot(m2NF,c2NF,'bs','LineWidth',2,'MarkerSize',10);
plot(m1PF,c1PF,'m>','LineWidth',2,'MarkerSize',10);
plot(m2PF,c2PF,'rs','LineWidth',2,'MarkerSize',10);
xlabel('mean copy nr.','FontSize',20);
ylabel('CV','FontSize',20);
set(gca,'FontSize',20,'XLim',[m1NF-s1NF m2NF+s2NF]);
legend('P1_{NF}','P2_{NF}','P1_{PF}','P2_{PF}','Location','EastOutside');
grid on;box on;
