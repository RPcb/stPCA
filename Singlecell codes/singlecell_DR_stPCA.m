% clc;
clear;
close all;

addpath("..");
% rawdata=importdata('chu_time_200genes.csv');
% hline = textscan(fpi, '%s', 1, 'delimiter', '\n');
% clear format;
% format='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s';
% cellinfo =textscan(fpi, format,100000,'delimiter', ',');
% rawdata.celltype=cellinfo{1,5};
load newdata_union.mat
newdata=newdata_union;
[gene_num, cell_num]=size(newdata.data);

c=zeros(6,cell_num);
c_num=zeros(6,1);
clu={'00h','12h','24h','36h','72h','96h'};
for i=1:6
    tp_num=0;
    for j=1:cell_num
        if isempty(cell2mat(strfind(newdata.textdata(1,1+j),clu(i))))
            continue;
        else
            tp_num=tp_num+1;
            c(i,tp_num)=j;
            c_stage(j)=i;
        end
    end
    c_num(i)=tp_num;
end

start=1;
ii=start-1;
window=60;
tid=0;
all_flat_z=zeros(floor(cell_num/window),window);

while ii+window<=cell_num
    tid=tid+1
    trainlength=window;
    xx=newdata.data(:,ii+1:ii+trainlength);

    [input_dimensions, time_point]=size(xx);
    
    noisestrength=0;          % noise could be added
    xx_noise=xx+noisestrength*rand(size(xx));
    
    %critical zone
    flat_y=zeros(trainlength,1);
    
    embeddings_num=20;    % 20,  dimensions after dimension reduction
    
    traindata=xx_noise(:,1:trainlength);
    traindata=traindata-mean(traindata,2);

    X=traindata';      

    m=trainlength;
    n=input_dimensions;
    L=embeddings_num;
    
    
    P=X(2:end,:);
    Q=X(1:end-1,:);
    a=0.01;  %10^(-11);
    b=1-a;
    
    
    
    %%  solve Z  %%
    
    H=zeros(n*L,n*L);
    H(1:n, 1:n)=a*X'*X-b*P'*P;
    H(1:n, n+1:2*n)=b*P'*Q;
    
    for j=2:L-1
        H(n*(j-1)+1:n*j, n*(j-1)+1-n:n*j-n)=b*Q'*P;
        H(n*(j-1)+1:n*j, n*(j-1)+1:n*j)=a*X'*X-b*P'*P-b*Q'*Q;
        H(n*(j-1)+1:n*j, n*(j-1)+1+n:n*j+n)=b*P'*Q;
    end
    
    H(n*(L-1)+1:n*L, n*(L-2)+1:n*(L-1))=b*Q'*P;
    H(n*(L-1)+1:n*L, n*(L-1)+1:n*L)=a*X'*X-b*Q'*Q;
    
    
    [V,D]=eig(H);
    
    
    ao=(diag(D));
    ao=real(ao);
    [aa, eigvIdx]=sort(ao,'descend');
    V=V(:,eigvIdx);
    %
    cW=V(:,1);

    aW=null(H);
    max_Z=0;
    ti=1;
    for ci=1:size(aW,2)
        cW=aW(:,ci);
        W=reshape(cW, n, L);
        Z=X*W;
        tmp=norm(Z,'fro');
        if tmp>max_Z
            max_Z=tmp;
            ti=ci;
        end
    end
   

    for ci=1:length(aa)
        if aa(ci)>0
            break;
        end
    end

    
    cW=V(:,ci);
    %          cW=V(:,8);
    W=reshape(cW, n, L);
    
    Z=X*W*max(abs(aa));
    
    weight_W=W*max(abs(aa));

    
    if ii==start
        weight_Z=norm(Z,'fro');
    end
    
    ii
    %     Z=Z/weight_Z              % normalization Z
    
    %  flat Z
    clear flat_z flat_z_pred
    sd_flat_z=0;
    sd_flat_z_pred=0;
    for zi=1:size(Z,1)
        num=0;
        for zj=1:size(Z,2)
            if zi-zj+1<1
                break;
            end
            num=num+1;
            tmp(num)=Z(zi-zj+1,zj);
            %                 [num2str(zi-zj+1),', ', num2str(zj)]
        end
        tmp=tmp(1:num);
        flat_z(zi)=mean(tmp);
        %         flat_z(zi)=max(tmp);
        sd_flat_z=sd_flat_z+std(tmp);
    end
    sd_flat_z=sd_flat_z/size(Z,1);
    
    for zj=2:size(Z,2)
        num=size(Z,2)-zj+1;
        for ni=1:num
            zi=size(Z,1)-ni+1;
            tmp(ni)=Z(zi,zj+ni-1);
            %                  [num2str(zi),', ', num2str(zj+ni-1)]
        end
        tmp=tmp(1:num);
        flat_z_pred(zj-1)=mean(tmp);
        %         flat_z_pred(zj-1)=max(tmp);
        sd_flat_z_pred=sd_flat_z_pred+std(tmp);
        %             num
    end
    sd_flat_z_pred=sd_flat_z_pred/(size(Z,2)-1);
    
    all_flat_z(tid,1:window)=flat_z;
    all_flat_z_pred(tid,1:embeddings_num-1)=flat_z_pred;

    
    hz=hankel(flat_z,[flat_z(m) flat_z_pred])
    %% temporal component exaction ,  Y=U1S1V1, part 1
    [U1,S1,V1]=svd(hz');
    pc=V1';
    pc_num=3;
    
    variance=abs(diag(S1));
    %     variance=var(pc');
    variance=variance/(sum(variance)+0.0001);
    
%     figure;
%     for i=1:pc_num
%         subplot(1,pc_num,i);
%         plot(pc(i,:),'r-*','LineWidth',2,'MarkerSize',4);
%         title(['DPCA, PC: #', num2str(i), ' , ( ', num2str(variance(i)*100),'% )'], 'FontSize', 20);
%         set(gca,'FontSize',20);
%     end
%     
%     figure;
%     plot3(pc(1,:),pc(2,:),pc(3,:),'r-*', 'LineWidth',1,'MarkerSize',2);
%     title('DPCA, 3 PCs reconstruction');
%     set(gca,'FontSize',20);
    
    
    colormap=[1.0000    0.4118    0.1608;0.6000    0.6000    0.2353; 0.3922    0.8314    0.0745;...
        0.0588    1.0000    1.0000;0.0745    0.6235    1.0000;1     0     1];
    
    
    
    ii=ii+window;
    
  
end

 colormap=[1.0000 0.4118 0.1608;0.6000 0.6000 0.2353;0.3922 0.8314 0.0745;...
        0.0588    1.0000    1.0000;0.0745    0.6235    1.0000;1     0     1];
 
 clu={'00h','12h','24h','36h','72h','96h'};
for i=1:6
    tp_num=0;
    for j=1:cell_num
        if isempty(cell2mat(strfind(newdata.textdata(1,1+j),clu(i))))
            continue;
        else
            tp_num=tp_num+1;
            c(i,tp_num)=j;
        end
    end
    c_num(i)=tp_num;
end   

% all_flat_z(tid,1:window)=flat_z;
all_flat_z=all_flat_z(1:tid,:)';
union_z=all_flat_z(:);
all_flat_z_pred=all_flat_z_pred(1:tid,:)';
union_z_pred=all_flat_z_pred(:);
% close all
figure;
plot(union_z,'Color',[0.8 0.8 0.8],'LineWidth',0.1);
hold on;
for i=1:length(union_z)
    scatter(i,union_z(i),'MarkerFaceColor',colormap(c_stage(i),:),'MarkerEdgeColor',colormap(c_stage(i),:));
    hold on;
end
title('flat z');
set(gca,'FontSize',20);


win_len=10;
for i=1:size(union_z,1)-win_len+1
    tmp_var_flat_z(i)=std(union_z(i:i+win_len-1));
end
figure;
plot(tmp_var_flat_z,'b-*', 'LineWidth',2);
title('SD of z');
set(gca,'FontSize',20);
 