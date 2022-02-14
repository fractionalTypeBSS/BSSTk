load example_data

stim=data0(371,:);c=double(stim);e(2:size(stim,2))=c(1:size(stim,2)-1);ee=c-e;
stim_info_1=find(ee(:)==512);    stim_info_2=find(ee(:)==1024);
stim_info_3=find(ee(:)==2048);   stim_info_4=find(ee(:)==4096);
stim_info_5=find(ee(:)==8192);
stim_info_stan=[stim_info_1; stim_info_2; stim_info_3; stim_info_4];
stim_info_sort=sort(stim_info_stan);
stim_info_dev=stim_info_5;
save stim_info_rt stim_info_sort stim_info_dev

data=data0; data(307:end,:)=[]; data(3:3:306,:)=[];

%BSSTk
display('BSSTk');
BSS_func_kishida('example_data_rf2',data,1000)

%wightedBSSTk
display('wighted BSSTk');
data0=[];for k=1:207  data0=[data0 data(:,stim_info_dev(k):stim_info_dev(k)+95)*0.2 data(:,stim_info_dev(k)+96:stim_info_dev(k)+276) data(:,stim_info_dev(k)+277:stim_info_dev(k)+499)*0.2];end
BSS_func_kishida('example_data_wavelet_rf2',data0,1000)
