% clc;
warning off;
clear;
close all;

load UU_18nodes.mat;

win_m=10;
for ni=1:180
    myzones(ni,:)=150+ni: 150+ni+win_m-1;
end

flag=1;
clear temp_var_y sd_x  ac
for ii=1:length(myzones)      %1:size(myzones,1)


    xx=UU(1:18,myzones(ii,:));

   sd_x(ii)=mean(std(xx,0,2));

    noisestrength=0;          % noise could be added

    xx_noise=xx+noisestrength*rand(size(xx));
    [input_dimensions, time_points]=size(xx);
    trainlength=time_points;
    %critical zone
    flat_y=zeros(trainlength,1);

    embeddings_num=5;    % dimensions after dimension reduction

    traindata=xx_noise(:,1:trainlength);
    traindata=traindata-mean(traindata,2);
    %         traindata=traindata-ones(size(traindata,1),1)*mean(traindata);
    %         [coeff,score,latent] = pca(traindata');
    %         var_pca(ii)=norm(score(:,1:6)*score(:,1:6)','fro');
    %
    %        testSD(ii)=mean(std(traindata,0,2));
    % %
    X=traindata';       % m=7, L=6, n=16

    m=trainlength;
    n=input_dimensions;
    L=embeddings_num;
    P=X(2:end,:);
    Q=X(1:end-1,:);

    a=0.01;     %10^(-11);    %  alpha
    b=1-a;                            %  beta



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
    eigv(ii)=aa(ci);
    cW=V(:,ci);

    W=reshape(cW, n, L);

    Z=X*W*max(abs(aa));


    if ii==1
        weight_Z=norm(Z,'fro');
    end

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
        sd_flat_z_pred=sd_flat_z_pred+std(tmp);
        %             num
    end
    sd_flat_z_pred=sd_flat_z_pred/(size(Z,2)-1);


    temp_var_y(ii)=norm(Z,'fro');
    
    temp_var_flat_y(ii)=std(flat_z);
    temp_var_flat_y_pred(ii)=std(flat_z_pred);
    tmp_r1(ii)=max(abs(aa));


    all_flat_z(ii,:)=flat_z;
    all_flat_z_pred(ii,:)=flat_z_pred;

    all_ac=autocorr([flat_z flat_z_pred]);
    all_ac=autocorr([flat_z]);
    ac(ii)=all_ac(2);
    continue;


    %% SVD for the Hankel matrix Z

    [U1,S1,V1]=svd(Z');
    pc=V1';
    pc_num=5;
    variance=diag(S1);
    variance=variance/(sum(variance));
    sd_pc1(ii)=std(pc(1,:));

    figure(1);
    for i=1:3
        subplot(1,3,i);
        plot(pc(i,:),'r-*','LineWidth',2,'MarkerSize',4);
        title(['PC: #', num2str(i), ' , ( ', num2str(variance(i)*100),'% )'], 'FontSize', 20);
        set(gca,'FontSize',20);
    end
    suptitle(['ZONE: ', num2str(ii)]);

%     saveas(gcf,['fig_res3/ZONE_ ',num2str(floor(ii)),'.fig']);
%     saveas(gcf,['fig_res3/ZONE_ ',num2str(floor(ii)),'.jpg']);


    continue;
    %% temporal component exaction ,  Y=U1S1V1, part 1
    for i=1:3 %pc_num
        figure;
        plot(pc(i,:),'r-*','LineWidth',2,'MarkerSize',4);
        title(['PC: #', num2str(i), ' , ( ', num2str(variance(i)*100),'% )'], 'FontSize', 20);
        set(gca,'FontSize',20);
    end

    %% spatial component exaction,  Y=U1S1V1,   part 1
    Xpc=X'/V1';
    Xpc=reshape(Xpc', trainlength, 4, 4);
    map=cmocean('balance');
    for i=1:3  %pc_num
        xr=real(squeeze(Xpc(i,:,:)));
        figure;
        contourf(xr);
        colormap(map);
        colorbar;
        title(['PC: #', num2str(i), ' , ( ', num2str(variance(i)*100),'% )'], 'FontSize', 20);
        set(gca,'FontSize',20);
    end

    %% clustering of temporal components, V1, Y
    pc=V1';
    figure;
    % sz = 25;
    c = linspace(1,0.2,length(pc(1,:)));
    % scatter(pc(1,:), pc(2,:), 80, c, 'filled');
    plot(pc(1,:), pc(2,:), 'r-*', 'MarkerSize', 6,'LineWidth',1);
    title('Clustering of temporal components from Y', 'FontSize', 20);
    xlabel('PC1');
    ylabel('PC2');
    set(gca,'FontSize',20);

    figure;
    plot3(pc(1,:), pc(2,:), pc(3,:), 'r-*', 'MarkerSize', 6,'LineWidth',1);
    % scatter3(pc(1,:), pc(2,:), pc(3,:), 80*ones(1,nt) ,  pc(3,:), 'filled');
    title('Clustering of temporal components from Y', 'FontSize', 30);
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    set(gca,'FontSize',20);
    grid on;

end




for ii=1:size(UU,1)
    figure;
    plot(UU(ii,150+1:340),'*-','LineWidth',3);
    title(['X',num2str(ii)], 'FontSize', 20);
    %hold on;
    set(gca,'FontSize',20);
end


figure;
plot(temp_var_y,'*-', 'LineWidth',2);
title('SD of Z');
set(gca,'FontSize',20);



