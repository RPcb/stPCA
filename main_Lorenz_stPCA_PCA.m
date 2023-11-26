% clc;
clear;
close all;
Y=mylorenz(20);     % dimension of X, n=60
noisestrength=0;
XX=Y+noisestrength*rand(size(Y));% noise could be added

start=1;
ii=158;             % case ID=158
    xx=XX(3000+ii:size(XX,1), :)';
    [input_dimensions, time_point]=size(xx);
    
    noisestrength=0;          % noise could be added
    xx_noise=xx+noisestrength*rand(size(xx));
    trainlength=50;             %  the known length, m
    flat_y=zeros(trainlength,1);
    
    embeddings_num=20;    % embedding dimension, L
    
    traindata=xx_noise(:,1:trainlength);
    traindata=traindata-mean(traindata,2);

    X=traindata';      
    
    m=trainlength;
    n=input_dimensions;
    L=embeddings_num;
    
    
    P=X(2:end,:);
    Q=X(1:end-1,:);
    a=0.1;  %10^(-11);
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


    
    for ci=1:length(aa)
        if aa(ci)>0
            break;
        end
    end
    ci

    
    cW=V(:,ci);
%          cW=V(:,8);
    W=reshape(cW, n, L);
    
    Z=X*W*max(abs(aa))
    
    weight_W=W*max(abs(aa));
%     figure;
%     heatmap(weight_W);
%     title('DPCA, weights W');
%     set(gca,'FontSize',20);
    
    
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
        end
        tmp=tmp(1:num);
        flat_z(zi)=mean(tmp);
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
    
    
    figure;
    plot([1:length(flat_z)],flat_z,'b-*', 'LineWidth',3);
    hold on;
    plot([length(flat_z)+1:length(flat_z)+length(flat_z_pred)],flat_z_pred, 'c-*', 'LineWidth',3);
    title('stPCA, flat z');
    set(gca,'FontSize',20);
    
    

    hz=hankel(flat_z,[flat_z(m) flat_z_pred]);
    %% temporal component exaction ,  Y=U1S1V1, part 1
    [U1,S1,V1]=svd(hz');
    pc_dpca=V1';
    pc_num=3;
    
    variance=abs(diag(S1));
    %     variance=var(pc');
    variance_dpca=variance/(sum(variance)+0.0001);
 
    [pc,expvar]=eig(hz'*hz);
    [expvar,idx]=sort(real(diag(expvar)),'descend');
    pc_dpca=real(pc(:,idx));
    
    Zdpca=hz*pc_dpca;
    pc_dpca=Zdpca';
    
    variance=abs(diag(expvar));
    variance_dpca=variance/(sum(variance)+0.0001);
    
    
    % compare with PCA 
    [pc,expvar]=eig(X'*X);
    [expvar,idx]=sort(real(diag(expvar)),'descend');
    pc2=real(pc(:,idx));
    
    Zpca=X*pc2;
    
    variance=abs(diag(expvar));
    variance_pca=variance/(sum(variance)+0.0001);
    
    


    figure(1);
    for i=1:pc_num
        subplot(2,pc_num,i);
        plot(pc_dpca(i,:),'r-*','LineWidth',2,'MarkerSize',4);
        title(['stPCA, PC: #', num2str(i)], 'FontSize', 20);
        set(gca,'FontSize',15);
    end
     

    for i=1:pc_num
        subplot(2,pc_num,pc_num+i);
        plot(Zpca(:,i),'b-*','LineWidth',2,'MarkerSize',4);
        title(['PCA, PC: #', num2str(i)], 'FontSize', 20);
        set(gca,'FontSize',15);
    end


    figure(2);
    subplot(2,2,1);
    plot(pc_dpca(1,:),pc_dpca(2,:),'r-*', 'LineWidth',1,'MarkerSize',2);
    title('stPCA, 2 PCs reconstruction');
    xlabel('PC1');
    ylabel('PC2');
    set(gca,'FontSize',20);
    
    subplot(2,2,2);
    plot3(pc_dpca(1,:),pc_dpca(2,:),pc_dpca(3,:),'r-*', 'LineWidth',1,'MarkerSize',2);
    title('stPCA, 3 PCs reconstruction');
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    set(gca,'FontSize',20);

    subplot(2,2,3);
    plot(Zpca(:,1),Zpca(:,2),'b-*', 'LineWidth',1,'MarkerSize',2);
    title('PCA, 2 PCs reconstruction');
    xlabel('PC1');
    ylabel('PC2');
    set(gca,'FontSize',20);
    
    subplot(2,2,4);
    plot3(Zpca(:,1),Zpca(:,2),Zpca(:,3),'b-*', 'LineWidth',1,'MarkerSize',2);
    title('PCA, 3 PCs reconstruction');
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    set(gca,'FontSize',20);
     
   
%     saveas(gcf,['PCA_DPCA_reconstruct_fig/PC_reconstruct_m=',num2str(m),'_L=',num2str(L),'_noise=',num2str(noisestrength),'_ii=',num2str(ii),'.fig']);
%     saveas(gcf,['PCA_DPCA_reconstruct_fig/PC_reconstruct_m=',num2str(m),'_L=',num2str(L),'_noise=',num2str(noisestrength),'_ii=',num2str(ii),'.jpg']);

    figure(3);
    subplot(1,3,1);
    plot([1:length(flat_z)],flat_z,'b-*', 'LineWidth',3);
    hold on;
    plot([length(flat_z)+1:length(flat_z)+length(flat_z_pred)],flat_z_pred, 'c-*', 'LineWidth',3);
    title('stPCA, flat z');
    xlabel('time points');
    ylabel('state value');
    set(gca,'FontSize',20);
    
    subplot(1,3,2);
    plot3(flat_z(1:end-2),flat_z(2:end-1),flat_z(3:end),'r-*', 'LineWidth',1,'MarkerSize',2);
    title('stPCA, flat z phase');
    xlabel('Z^1');
    ylabel('Z^2');
    zlabel('Z^3');
    set(gca,'FontSize',20);
    
    subplot(1,3,3);
    plot3(Zpca(1:end-2,1),Zpca(2:end-1,1),Zpca(3:end,1),'b-*', 'LineWidth',1,'MarkerSize',2);
    title('PCA, PC1 phase');
    xlabel('Z^1');
    ylabel('Z^2');
    zlabel('Z^3');
    set(gca,'FontSize',20);

