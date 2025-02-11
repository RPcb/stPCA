% clc;
clear;
close all;

rawdata=importdata('chu_time_200genes.csv');
newdata_union.textdata(:,1)=rawdata.textdata(:,1);
profile=rawdata.data;
[gene_num, cell_num]=size(profile);

c=zeros(6,cell_num);
clu={'00h','12h','24h','36h','72h','96h'};
figure;
for i=1:6
    tp_num=0;
    for j=1:cell_num
        if isempty(cell2mat(strfind(rawdata.textdata(1,1+j),clu(i))))
            continue;
        else
            tp_num=tp_num+1;
            c(i,tp_num)=j;
        end
    end
    c_num(i)=tp_num;
    
    tmp_profile=profile(:,c(i,1:tp_num));
    dist=squareform(pdist(tmp_profile','euclidean'));
    [sorted_dist,idx]=sort(dist,2);
    k=10;
    KNN_net_idx=idx(:,2:k+1);
    KNN_net=zeros(c_num(i),c_num(i));
    for j=1:c_num(i)
        KNN_net(j,KNN_net_idx(j,:))=dist(j,KNN_net_idx(j,:));
    end
    KNN_net_v=mean(KNN_net,2);
    sources=[];
    targets=[];
    weights=[];
    edge_num=0;
    for j1=1:c_num(i)
        for j2=1:k
            edge_num=edge_num+1;
            sources(edge_num)=j1;
            targets(edge_num)=KNN_net_idx(j1,j2);
            weights(edge_num)=dist(i,KNN_net_idx(j1,j2));
        end
    end
    head=1;
    G=graph(sources,targets,weights);
    TR = shortestpathtree(G,head);
%     figure;
    subplot(2,3,i);
    p=plot(TR,'Layout','auto');
    
    colormap=[1.0000    0.4118    0.1608;0.6000    0.6000    0.2353; 0.3922    0.8314    0.0745;...
        0.0588    1.0000    1.0000;0.0745    0.6235    1.0000;1     0     1];
    highlight(p,1:c_num(i),'NodeColor',colormap(i,:),'MarkerSize',6);
    highlight(p,1,'NodeColor',[0 0 0],'MarkerSize',10);
    %     legend(clu{i});
    % legend('00h','12h','24h','36h','72h','96h','FontSize',10);
    highlight(p,TR,'EdgeColor',[0.8020 0.8020 0.8020],'LineWidth',0.2)
    title(['Tree, stage=',num2str(i)]);
    set(gca,'FontSize',20);
    
    
    clear Ddist pTime
    for j=1:c_num(i)
        [tmp,Ddist(j)]=shortestpath(G,head,j);
        pTime(j)=Ddist(j)/(KNN_net_v(head)+KNN_net_v(j));
    end
    
    tmpi=c(i,1:c_num(i));
    [pTime,idx]=sort(pTime);
    newdata_token(i).data=tmp_profile(:,idx);
    newdata_token(i).textdata(:,1)=rawdata.textdata(:,1);
    newdata_token(i).textdata(:,2:c_num(i)+1)=rawdata.textdata(:,tmpi(idx)+1);
    
    newdata_union.data(:,sum(c_num(1:i-1))+1:sum(c_num(1:i)))=tmp_profile(:,idx);
    newdata_union.textdata(1:gene_num+1,sum(c_num(1:i-1))+2:sum(c_num(1:i))+1)=rawdata.textdata(:,tmpi(idx)+1);
    
end
save newdata_token newdata_token
save newdata_union newdata_union

