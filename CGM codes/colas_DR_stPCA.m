clc
clear
close all
warning off;

  
load ColasData_sub_1
my_data=ColasData_sub_1;
line=1:1512;

for k=1:189
clear temp_var_y
wlen=20;
for ii=1:size(my_data,2)-wlen      
    
    k1=k*8-1
    xx=my_data(line(k1-5:k1),ii:ii+wlen-1);
    noisestrength=0;          % noise could be added
    
    xx_noise=xx+noisestrength*rand(size(xx));
    [input_dimensions, time_points]=size(xx);
    trainlength=time_points;
    %critical zone
    flat_y=zeros(trainlength,1);
    
    embeddings_num=4;    % dimensions after dimension reduction
    
    traindata=xx_noise(:,1:trainlength);
    traindata=traindata-mean(traindata,2);
    
    %%  PCA
    [coeff,score,latent] = pca(traindata');
    flat_Z_pca(ii,:)=traindata'*coeff(:,1);
    var_pca(ii)=norm(traindata'*coeff(:,1),'fro');
    
    %% DPCA
    X=traindata';
    m=trainlength;
    n=input_dimensions;
    L=embeddings_num;
    P=X(2:end,:);
    Q=X(1:end-1,:);
    a=0.01;    
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
    ti
 
    for ci=1:length(aa)
        if aa(ci)>0
            break;
        end
    end
    ci
    aa(ci)
    
    cW=V(:,ci);
    W=reshape(cW, n, L);
    
    Z=X*W*max(abs(aa));
    
    weight_W=W*max(abs(aa));
   
    ii
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
     end
    sd_flat_z_pred=sd_flat_z_pred/(size(Z,2)-1);
    
    temp_var_y(ii)=norm(Z,'fro');
    temp_var_flat_y(ii)=std(flat_z);
    temp_var_flat_y_pred(ii)=std(flat_z_pred);
    
    tmp_r1(ii)=max(abs(aa));
       
    x(k)={temp_var_y};
    mean_SD(k)=mean(temp_var_y);
    
    
    hz=hankel(flat_z,[flat_z(m) flat_z_pred]);
    
    
    
    %% temporal component exaction ,  Y=U1S1V1, part 1
    [U1,S1,V1]=svd(hz');
    pc=V1';
    pc_num=3;
    
    variance=abs(diag(S1));
    variance=variance/(sum(variance)+0.0001);
     
    ii
     
end

end
set(0,'defaultfigurecolor','w')
load ColasLabel1
%Plot the SD curve for each patient
figure;
for p1=1:25
    subplot(5,5,p1)
    plot(x{1,p1},'b','LineWidth',2);
    title( ColasLabel1(p1))
end
figure;
for p2=26:50
    subplot(5,5,p2-25)
    plot(x{1,p2},'b','LineWidth',2);
    title( ColasLabel1(p2))
end
figure;
for p3=51:75
    subplot(5,5,p3-50)
    plot(x{1,p3},'b','LineWidth',2);
    title( ColasLabel1(p3))
end
figure;
for p4=76:100
    subplot(5,5,p4-75)
    plot(x{1,p4},'b','LineWidth',2);
    title( ColasLabel1(p4))
end
figure;
for p5=101:125
    subplot(5,5,p5-100)
    plot(x{1,p5},'b','LineWidth',2);
    title( ColasLabel1(p5))
end
figure;
for p6=126:150
    subplot(5,5,p6-125)
    plot(x{1,p6},'b','LineWidth',2);
    title( ColasLabel1(p6))
end
figure;
for p7=151:172
    subplot(5,5,p7-150)
    plot(x{1,p7},'b','LineWidth',2);
    title( ColasLabel1(p7))
end
figure;
for p8=173:189
    subplot(5,5,p8-172)
    plot(x{1,p8},'r','LineWidth',2);
    title( ColasLabel1(p8))
end


%Plot the mean curve of the patient's SD
figure;
plot(mean_SD,'-*','LineWidth',2)
hold on;
plot([0,190],[mean(mean_SD),mean(mean_SD)],'LineWidth',1)
hold off;