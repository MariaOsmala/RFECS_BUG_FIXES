function [indiv_enhancers_pks indiv_enhancers_locs]=peak_call(prob_dist,peak_dist, start)
%prob_dist are the prediction scores of a single enhancer region, composed
%of adjacent regions having an enhancer score larger than 0.5
%prob_dist is a vector of enhancer scores between 0.5 and 1
%peak_dist=20
%peak_call(vals1,peak_filt_dist);

indiv_enhancers_pks=[];
indiv_enhancers_locs=[];

if(length(prob_dist)<=peak_dist)
     %finds the maximum value -> pks
    %locs is the first occurrence of the maximal value
    [pks locs]=max(prob_dist);
    indiv_enhancers_pks=[indiv_enhancers_pks pks];
    indiv_enhancers_locs=[indiv_enhancers_locs start+locs-1];
    
else
    [pks locs]=max(prob_dist);
    indiv_enhancers_pks=[indiv_enhancers_pks pks];
    indiv_enhancers_locs=[indiv_enhancers_locs start+locs-1]; %21
  
    first_vector_end=locs-peak_dist-1; %-19
    second_vector_start=locs+peak_dist+1; %21
    clear('vectors');
    vectors.vector=prob_dist(1:first_vector_end);
    vectors(2).vector=prob_dist(second_vector_start:end);
    vectors(1).start=start;
    vectors(2).start=start+second_vector_start;
    vectors(1).length=length(vectors(1).vector);
    vectors(2).length=length(vectors(2).vector);
    
    
    %remove zero-length vectors
    
    nonempty=extractfield(vectors, 'length')~=0;
    vectors=vectors(nonempty);
    
    %it can be that the vector start is less than start,
    problem=find((extractfield(vectors, 'start')<start) ==1);
    if(~isempty(problem))
        
        vectors(problem).vector=vectors(problem).vector(start-vectors(problem).start+1:end);
        vectors(problem).start=start;
        vectors(problem).length=length(vectors(problem).vector);
    end    
    
    clear('vectors_temp');
    while(~isempty(vectors)) %1
       %if(length(indiv_enhancers_locs)>1)
       %    indiv_enhancers_locs(end-1:end)
       %end
       %are there vectors with length less or equal to 20 -> add to indiv_enhancers
       peak_dist_length= find((extractfield(vectors, 'length')<=(peak_dist)+1) ==1);
       if(~isempty(peak_dist_length))
           for i=1:length(peak_dist_length)
              % short='short';
               %short
               %vectors(peak_dist_length(i)).length
               %length(vectors)-length(peak_dist_length)
               [pks locs]=max( vectors(peak_dist_length(i)).vector );
               indiv_enhancers_pks=[indiv_enhancers_pks pks];
               indiv_enhancers_locs=[indiv_enhancers_locs vectors(peak_dist_length(i)).start+locs-1];
           
           end
           %remove peak_dist_length from vectors
           vectors(peak_dist_length)=[];
       end
       if(~isempty(vectors))
           peak_dist_length= find((extractfield(vectors, 'length')>peak_dist) ==1); %error here
       else
           peak_dist_length=[];
       end
       
       if(~isempty(peak_dist_length))
           %length(peak_dist_length)
           for i=1:length(peak_dist_length)
                   %i
               [pks locs]=max( vectors(peak_dist_length(i)).vector );
               indiv_enhancers_pks=[indiv_enhancers_pks pks];
               indiv_enhancers_locs=[indiv_enhancers_locs vectors(peak_dist_length(i)).start+locs-1];
  
               first_vector_end=locs-peak_dist-1;
               second_vector_start=locs+peak_dist+1;
               clear('vectors_temp');
               vectors_temp.vector=vectors(peak_dist_length(i)).vector(1:first_vector_end); %error occurs here
               vectors_temp(2).vector=vectors(peak_dist_length(i)).vector(second_vector_start:end);
               
               vectors_temp(1).start=vectors(peak_dist_length(i)).start;
               vectors_temp(2).start=vectors(peak_dist_length(i)).start+second_vector_start;
               
               vectors_temp(1).length=length(vectors_temp(1).vector);
               vectors_temp(2).length=length(vectors_temp(2).vector);
               %remove zero-length vectors
    
               nonempty=extractfield(vectors_temp, 'length')~=0;
               vectors_temp=vectors_temp(nonempty);
               
               %add vectors_temp to vectors
               if(length(vectors_temp)~=0)
                   vectors=[vectors, vectors_temp];
               end
           
           end
               %remove peak_dist_length from vectors
               vectors(peak_dist_length)=[];
       end
       
       
       
          
        
    end
end
    
end






%if(test-1>peak_dist)
%   test    
%end