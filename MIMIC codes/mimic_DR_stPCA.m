%clc;
warning off;
clear;
close all;

folderName = 'subject_z_results_all';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

start=1;
fileDir=dir(['../subject_items_at_dhours_all/' '*_filtered_events_hadm.csv']);
for fi=1:length(fileDir)
    fileDir(fi).name
    textdata=importdata(['../subject_items_at_dhours_all/' fileDir(fi).name]);
    allxx=textdata.data;

    [input_dimensions, time_point]=size(allxx);
    time_point,fi

    trainlength=6;

    clear tmp_svd_proj_sd;
    for ii=1:time_point-trainlength+1
        xx=allxx(:,ii:ii+trainlength-1);

        noisestrength=0;          % noise could be added
        xx_noise=xx+noisestrength*rand(size(xx));



        %critical zone
        flat_y=zeros(trainlength,1);

        embeddings_num=4;    % dimensions after dimension reduction

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


        for ci=1:length(aa)
            if aa(ci)>0
                break;
            end
        end

        cW=V(:,ci);
        W=reshape(cW, n, L);

        Z=X*W*max(abs(aa));

        weight_W=W*max(abs(aa));


        if ii==start
            weight_Z=norm(Z,'fro');
        end

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

            end
            tmp=tmp(1:num);
            flat_z_pred(zj-1)=mean(tmp);
            sd_flat_z_pred=sd_flat_z_pred+std(tmp);
            %             num
        end
        sd_flat_z_pred=sd_flat_z_pred/(size(Z,2)-1);

        all_flat_z(ii,1:trainlength)=flat_z;
        all_flat_z_pred(ii,1:embeddings_num-1)=flat_z_pred;


        hz=hankel(flat_z,[flat_z(m) flat_z_pred]);
        [U1,S1,V1]=svd(hz');
        pc_dpca=V1';
        pc_num=3;
        variance=abs(diag(S1));
        variance_dpca=variance/(sum(variance)+0.0001);

        [pc,expvar]=eig(hz'*hz);
        [expvar,idx]=sort(real(diag(expvar)),'descend');
        pc_dpca=real(pc(:,idx));

        Zdpca=hz*pc_dpca;
        pc_dpca=Zdpca';


        tmp_svd_proj_sd(:,ii)=std(pc_dpca(1:3,:),0,2);

    end
    svd_projection_sd{fi}=tmp_svd_proj_sd;

    clear tmp_z_sd tmp_z_pred_sd
    for i=1:time_point-trainlength+1
        tmp_z_sd(i)=std(all_flat_z(i,:));
        tmp_z_pred_sd(i)=std(all_flat_z_pred(i,:));
    end
    window_z_sd{fi}=tmp_z_sd;
    window_z_pred_sd{fi}=tmp_z_pred_sd;

end
save window_z_sd window_z_sd
save  window_z_pred_sd  window_z_pred_sd
save svd_projection_sd svd_projection_sd

for i=1:length(window_z_sd)
    h1=figure(1);
    plot(window_z_sd{i},'b' ,'LineWidth',2);
    hold on;
    plot(window_z_pred_sd{i},'r', 'LineWidth',2);


    token=strfind(fileDir(i).name,'_');
    title(['fid=',num2str(i),', ICUID=',fileDir(i).name(1:token-1)]);
    xlabel('Hours');
    ylabel('SD of z');
    set(gca,'FontSize',20);
    saveas(h1,['subject_z_results_all/ICUID=',fileDir(i).name(1:token-1),'.fig']);
    saveas(h1,['subject_z_results_all/ICUID=',fileDir(i).name(1:token-1),'.jpg']);
    hold off;

    h2=figure(2);
    tmp_pc=svd_projection_sd{i};
    plot(tmp_pc(1,:),tmp_pc(2,:),'Color',[0.8 0.8 0.8],'LineWidth',0.1);
    hold on;
    s=scatter(tmp_pc(1,:),tmp_pc(2,:),'MarkerFaceColor','r','MarkerEdgeColor',[0 0 0]);
    hold on;
    s=scatter(tmp_pc(1,1),tmp_pc(2,1),'MarkerFaceColor','y','MarkerEdgeColor',[0 0 0]);
    s.SizeData=180;


    title(['fid=',num2str(i),', ICUID=',fileDir(i).name(1:token-1)]);
    xlabel('Hours');
    ylabel('SD of SVD PC projection');
    saveas(h2,['subject_z_results_all/PC_ICUID=',fileDir(i).name(1:token-1),'.fig']);
    saveas(h2,['subject_z_results_all/PC_ICUID=',fileDir(i).name(1:token-1),'.jpg']);
    set(gca,'FontSize',20);
    hold off;
end

clear ICUSTAY_IDs

fileDir=dir(['../subject_items_at_dhours_all/' '*_filtered_events_hadm.csv']);
for i=1:length(window_z_sd)
    token=strfind(fileDir(i).name,'_');
    ICUSTAY_IDs{i}=fileDir(i).name(1:token-1);
end
save ICUSTAY_IDs ICUSTAY_IDs