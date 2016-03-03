function profiles = getradialprofiles(clouds)
    for i=1:size(clouds,3)
        profiles(:,i) = mean(clouds(140:150,:,i),1);
        profiles(:,i) = profiles(:,i)-mean(profiles(1:5,i));
    end
    
end