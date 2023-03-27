clc;
clear;
close all;
%initial values
fall_min=inf;
non_fall_max=0;
resultant_array=[];
droll_arr=[];
dpitch_arr=[];
nf_resultant_array=[];
nf_droll_arr=[];
nf_dpitch_arr=[];
%read data
dataset='D:\TARP\dataset';
folders=dir(dataset);
for i=1:length(folders)
    folder_name=folders(i).name;
    if startsWith(folder_name,'.')
        continue
    end
    folder_path=fullfile(dataset,folder_name);
    sub_folders=dir(folder_path);
    for j=1:length(sub_folders)
        sub_folder_name=sub_folders(j).name;
        if startsWith(sub_folder_name,'.')
            continue
        end
        sub_folder_path = fullfile(folder_path, sub_folder_name);
        files=dir(sub_folder_path);
        for j=1:length(files)
            file_name=files(j).name;
            if startsWith(file_name,'.')
                continue
            end
            file_path = fullfile(sub_folder_path, file_name);
            data = load(file_path);
            if strcmp(sub_folder_name, 'fall')
                resultant=sqrt(data.droll.^2+data.dpitch.^2);
                fall_min=min(fall_min,max(resultant));
                resultant_array=[resultant_array;resultant];
                droll_arr=[droll_arr;data.droll];
                dpitch_arr=[dpitch_arr;data.dpitch];
            elseif strcmp(sub_folder_name, 'non-fall')
                resultant=sqrt(data.droll.^2+data.dpitch.^2);
                non_fall_max=max(non_fall_max,max(resultant));
                nf_resultant_array=[nf_resultant_array;resultant];
                nf_droll_arr=[nf_droll_arr;data.droll];
                nf_dpitch_arr=[nf_dpitch_arr;data.dpitch];
            end
        end
    end
end
if(fall_min<non_fall_max)
    disp(['There is overlap. FT1=',num2str(fall_min),' rad/s'])
    FT1=fall_min;
else
    disp('There is no overlap')
end

%%% FT2 calculation %%%
window_size = 50; %cuz sampling rate = 50Hz
min_alpha_res = Inf;
ft1_exceed_idx=[];
for i = window_size/2+1:length(resultant_array)-window_size/2
    if(resultant_array(i)>FT1)
        idx=i;
        ft1_exceed_idx=[ft1_exceed_idx;idx];
        start_index = idx - window_size/2; % find start index for window
        end_index = idx + window_size/2; % find end index for window
        alpha_p=diff(dpitch_arr(start_index:end_index));
        alpha_r=diff(droll_arr(start_index:end_index));
        alpha_res=sqrt(alpha_p.^2+alpha_r.^2);
        min_alpha_res=min(min_alpha_res,max(alpha_res));
    end
end
disp(['FT2=',num2str(min_alpha_res),' rad/s2'])
FT2=min_alpha_res;

%%% FT3 calculation %%%
min_theta_res=inf;
for i = 1:length(ft1_exceed_idx)
    % Calculate the start and end indices for the 1.7 second window
    start_idx = max(1, ft1_exceed_idx(i) - 1.2 * 50);
    end_idx = min(length(resultant), ft1_exceed_idx(i) + 0.5 * 50);

    % Calculate the time array for the window
    t = ((start_idx:end_idx) - ft1_exceed_idx(i)) / 50;

    % Calculate the trunk angle using cumulative integration
    theta_r(start_idx:end_idx) = cumtrapz(t, droll_arr(start_idx:end_idx));
    theta_p(start_idx:end_idx) = cumtrapz(t, dpitch_arr(start_idx:end_idx));

    %Resultant trunk angle
    theta_res=sqrt(theta_p.^2+theta_r.^2);
    min_theta_res=min(min_theta_res,max(theta_res));
end
disp(['FT3=',num2str(min_theta_res),' rad'])
FT3=min_theta_res;

nf_ft1_exceed_idx = find(nf_resultant_array > FT1);
false_positives=0;
nf_alpha_res=0;
nf_theta_res=0;
for i=1:length(nf_ft1_exceed_idx)
        idx=nf_ft1_exceed_idx(i);
        start_index = max(1,idx - window_size/2); % find start index for window
        end_index = min(length(nf_dpitch_arr),idx + window_size/2); % find end index for window
        alpha_p=diff(nf_dpitch_arr(start_index:end_index));
        alpha_r=diff(nf_droll_arr(start_index:end_index));
        nf_alpha_res=sqrt(alpha_p.^2+alpha_r.^2);

        start_idx = max(1, idx - 1.2 * 50);
        end_idx = min(length(nf_droll_arr), idx + 0.5 * 50);

        % Calculate the time array for the window
        t = ((start_idx:end_idx) - idx) / 50;

        % Calculate the trunk angle using cumulative integration
        theta_r(start_idx:end_idx) = cumtrapz(t, nf_droll_arr(start_idx:end_idx));
        theta_p(start_idx:end_idx) = cumtrapz(t, nf_dpitch_arr(start_idx:end_idx));

        %Resultant trunk angle
        nf_theta_res=sqrt(theta_p.^2+theta_r.^2);

        if((nf_alpha_res>FT2) & (nf_theta_res>FT3))
            false_positives=false_positives+1;
        end
end
disp(['False Positives: ',num2str(false_positives)])
