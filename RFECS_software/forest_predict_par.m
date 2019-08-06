function forest_predict_par(forest_all,nvar,ntrees,len_var,mod_vec,extn,peak_filt_dist,input_path,output_path,chr_set)
% enhancer predictions using p300 vs TSS-bg

%forest_all=forest_p300_tssbg_all; 
%nvar=15;                 
%ntrees=65; 
%len_var=len_var; 
%mod_vec=[1:15]; 
%extn=extn; 
%peak_filt_dist=20;
%input_path='test_set/all_mods_'; 
%output_path='test_set/';
%chr_set=chr_set;





%output_path='hESC/H1_diff_rpkm/hg18_chr_mods/imr90_comb/';




for x=1:length(chr_set)
        %x
        all_mods=load(char(strcat(input_path,chr_set(x)))); %2492506 x 15
        all_mods=all_mods(:,mod_vec);
        index_map=[];
        mean_bg_all=[];

        for k=1:len_var
            limit=floor(length(all_mods)/len_var)*len_var+k-1; %2492500
            if(limit>length(all_mods))
                limit=limit-len_var;
            end
            index_map=[index_map;[k:len_var:limit]']; %1:20:1492481
            cd4_set=[];
            for j=1:nvar
                cd4_set=[cd4_set reshape(all_mods(k:limit,j)',len_var,[])']; %124625 x 300
            end
            %predicting background score
            y_bg=[];

            parfor i=1:ntrees %for i=1:ntrees %%%CHANGED parfor to for
                y_bg(:,i)=str2num(char(eval(forest_all{i},cd4_set)));
            end

            mean_bg=mean(y_bg')';
            mean_bg_all=[mean_bg_all;mean_bg];
        end
    file=strcat(output_path,chr_set(x),extn);
    [s sind]=sort(index_map);
    mean_bg_all=mean_bg_all(sind); %
    %save(char(file),'mean_bg_all');
    
    %save scores of all the bins to a separate file
    %note that we are using (bin location + len_var/2) as in chr_enhancers.txt files
    filename=strcat(output_path,chr_set(x),extn,'_allbins','.txt');
    fid=fopen(char(filename),'w');
    for i=1:length(mean_bg_all)
        fprintf(fid,'%s\t%d\t%f\n',char(chr_set(x)),(i+(len_var/2))*100,mean_bg_all(i));
    end
    fclose(fid);
    
    index_sel=find(mean_bg_all>0.5); %mean_bg_all 2492487 x 1
    ind = intersect(find(index_sel > len_var), find(index_sel<=length(index_map) - len_var));
    index_sel=index_sel(ind); %locations that are predicted as enhancers, not too close to the chromosome end

    mean_bg_sel=mean_bg_all(index_sel); %1098877 x 1

    patches2=[];

    start=index_sel(1); %21
    prev=index_sel(1); %21 %debug from here
    for i=2:length(index_sel) %2:1098877
        if(index_sel(i)-prev==1)
            prev=index_sel(i);
        else
            patches2=[patches2;start prev];
            start=index_sel(i);
            prev=start;
        end
    end
    patches2=[patches2;start prev]; %109897 x 2 tells which indices contain adjacent enhancer predictions

    patches2_bg={};
    for i=1:length(patches2) %109897
        abc=[];
        for j=patches2(i,1):patches2(i,2)
            abc=[abc mean_bg_all(j)];
        end
        patches2_bg{i}=abc;
    end
    
    
    patches2_bg_lengths=zeros(1,length(patches2_bg));
    for i=1:length(patches2_bg_lengths)
       patches2_bg_lengths(i)=length(patches2_bg{i});
    end    
    
    patches2_bg_hist=hist(patches2_bg_lengths,2:max(patches2_bg_lengths));
    % max(patches2_bg_lengths) 210494
    %patches2_bg gives the scores of different enhancer groups (adjacent enhancers form a group)
    mean_y_bg_uniq=[0:1:ntrees]./ntrees; %1x66 66 values between 0 and 1, equally separated
    peaks=[];
    peaks_prob=[];
    %wrong_locs=[];
    %filename=strcat('when_recursion_fails.txt');
    %fid=fopen(char(filename),'w');
    for i=1:length(patches2) %1:109897
        %fprintf(fid,'%f\t%f\n',i,patches2_bg_lengths(i));
        %i=37605
        %vals1=patches2_bg{i}; %first is 5614 length vector
        %[pks locs]=peak_call(vals1,peak_filt_dist); %if peak_dist is 2kb, it should be 20
        %start_temp=patches2(i,1);
        %[pks locs]=peak_call_without_recursion_version3(vals1, peak_filt_dist, start_temp);
        %locs_sorted=sort(locs);
        %locs2_sorted=sort( patches2(i,1)+locs2-1);
        %[pks locs]=peak_call_without_recursion(vals1, peak_filt_dist, start_temp);
        %wrong_locs=[wrong_locs,length(locs)-length(find(sort(locs)==sort(patches2(i,1)+locs2-1)))];
        %find(sort(locs)==sort(patches2(i,1)+locs2-1))
        %if( length(locs)~= length( patches2(i,1)+locs2-1)
        %    print1='false';
        %    i
        %    print1
        %end
        %locs and pks 249 and 249 length vectors
        %peaks=[peaks, locs];
        %peaks=[peaks patches2(i,1)+locs-1]; % change condition to patches2_2kb
        %converts a large region to 249 individual enhancers
        %peaks_prob=[peaks_prob pks];
	vals1=patches2_bg{i};
        [pks locs]=peak_call(vals1,peak_filt_dist); %if peak_dist is 2kb, it should be 20
        peaks=[peaks patches2(i,1)+locs-1]; % change condition to patches2_2kb
        peaks_prob=[peaks_prob pks];



    end
    %fclose(fid);

        [s sind]=sort(peaks);
        peaks=peaks(sind);
        peaks_prob=peaks_prob(sind);
        prev=1;
        peaks_filt=[];
        peaks_prob_filt=[];

        for i=2:length(peaks)
            if(peaks(i)-peaks(prev)<=peak_filt_dist)
                if(peaks_prob(i)>peaks_prob(prev))
                    prev=i;
                end
            else
                peaks_filt=[peaks_filt;peaks(prev)];
                peaks_prob_filt=[peaks_prob_filt;peaks_prob(prev)];
                prev=i;
            end
        end
        mean_bg_filt=zeros(length(peaks_filt),1);
        for i=1:length(peaks_filt)
            ind1=find(index_sel==peaks_filt(i));
            mean_bg_filt(i)=mean_bg_sel(ind1);
        end
        loc=(peaks_filt+(len_var/2))*100;
        filename=strcat(output_path,chr_set(x),extn,'.txt');
        fid=fopen(char(filename),'w');
        for i=1:length(loc)
            fprintf(fid,'%s\t%d\t%f\n',char(chr_set(x)),loc(i),mean_bg_filt(i));
        end
        fclose(fid);
end
