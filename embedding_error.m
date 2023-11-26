function [total_embedding_error] = embedding_error(Z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[trainlength, embeddings_num]=size(Z);

norm_Z=reshape(normalize(Z(:)),trainlength,embeddings_num);

%clear flat_z flat_z_pred
sd_flat_z=0;
sd_flat_z_pred=0;
for zi=1:size(norm_Z,1)
    num=0;
    for zj=1:size(norm_Z,2)
        if zi-zj+1<1
            break;
        end
        num=num+1;
        tmp(num)=norm_Z(zi-zj+1,zj);
        %                 [num2str(zi-zj+1),', ', num2str(zj)]
    end
    tmp=tmp(1:num);
    %flat_z(zi)=mean(tmp);
    sd_flat_z=sd_flat_z+std(tmp);
end
%sd_flat_z=sd_flat_z/size(norm_Z,1);

for zj=2:size(norm_Z,2)
    num=size(norm_Z,2)-zj+1;
    for ni=1:num
        zi=size(norm_Z,1)-ni+1;
        tmp(ni)=norm_Z(zi,zj+ni-1);
        %                  [num2str(zi),', ', num2str(zj+ni-1)]
    end
    tmp=tmp(1:num);
    %flat_z_pred(zj-1)=mean(tmp);
    sd_flat_z_pred=sd_flat_z_pred+std(tmp);
    %             num
end
%sd_flat_z_pred=sd_flat_z_pred/(size(norm_Z,2)-1);
total_embedding_error=sd_flat_z+sd_flat_z_pred;
end