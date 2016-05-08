function [ integral, error ] = RadialIntegral_temp( dat, center, radius )

if radius <= sqrt(2)
    disp(['WARNING! For radius smaller than ',num2str(sqrt(2)),' this program gives zero for integral. This will be taken care of in version 2']);
end

size1 = size(dat,1);
size2 = size(dat,2);
i0 = center(1); j0 = center(2);

%% Calculate N+1 Corners Matrix
corners = -1*ones(size1+1,size2+1);
for i = 1:size1+1
    for j = 1:size2+1
        % for each element calculate if the point is inside or not
        corners(i,j) = ( (i-0.5-i0)^2 + (j-0.5-j0)^2 ) <= radius^2;
    end
end

%% Determine number of corners inside circle
noCorners = -1*ones(size1,size2);
for i = 1:size1
    for j = 1:size2
        % for each element add up the four corners
        noCorners(i,j) = corners(i,j) + corners(i,j+1) + corners(i+1,j) + corners(i+1,j+1);
    end
end

%% Fill out the NxN matrix weights
weights = zeros(size1,size2);
for i = 1:size1
    for j = 1:size2
        
        % If the noCorners is 0
        if noCorners(i,j) == 0
            weights(i,j) = 0;
        
        % If the noCorners is 4
        elseif noCorners(i,j) == 4
            weights(i,j) = 1;
            
        % If the noCorners is 1
        elseif noCorners(i,j) == 1
            % bottom left
            if corners(i,j)
                x1 = [i-0.5, j0 + sqrt(radius^2 - (i-0.5-i0)^2)];
                x2 = [i0 + sqrt(radius^2 - (j-0.5-j0)^2), j-0.5];
                ht = x1(2) - (j-0.5);
                ln = x2(1) - (i-0.5);
            % bottom right
            elseif corners(i+1,j)
                x1 = [i+0.5, j0 + sqrt(radius^2 - (i+0.5-i0)^2)];
                x2 = [i0 - sqrt(radius^2 - (j-0.5-j0)^2), j-0.5];
                ht = x1(2) - (j-0.5);
                ln = (i+0.5) - x2(1);
            % top right
            elseif corners(i+1,j+1)
                x1 = [i+0.5, j0 - sqrt(radius^2 - (i+0.5-i0)^2)];
                x2 = [i0 - sqrt(radius^2 - (j+0.5-j0)^2), j+0.5];
                ht = (j+0.5) - x1(2);
                ln = (i+0.5) - x2(1);
            % top left
            elseif corners(i,j+1)
                x1 = [i-0.5, j0 - sqrt(radius^2 - (i-0.5-i0)^2)];
                x2 = [i0 + sqrt(radius^2 - (j+0.5-j0)^2), j+0.5];
                ht = (j+0.5) - x1(2);
                ln = x2(1) - (i-0.5);
            end
            % Compute area
            weights(i,j) = 0.5 * ht * ln + segment(radius, x1 - center, x2 - center);
            
        % If the noCorners is 2
        elseif noCorners(i,j) == 2
            % Side left
            if corners(i,j) && corners(i,j+1)
                x1 = [i0 + sqrt(radius^2 - (j+0.5-j0)^2), j+0.5];
                x2 = [i0 + sqrt(radius^2 - (j-0.5-j0)^2), j-0.5];
                ln1 = x1(1) - (i-0.5);
                ln2 = x2(1) - (i-0.5);
            % Side bottom
            elseif corners(i,j) && corners(i+1,j)
                x1 = [i-0.5, j0 + sqrt(radius^2 - (i-0.5-i0)^2)];
                x2 = [i+0.5, j0 + sqrt(radius^2 - (i+0.5-i0)^2)];
                ln1 = x1(2) - (j-0.5);
                ln2 = x2(2) - (j-0.5);
            % Side right
            elseif corners(i+1,j) && corners(i+1,j+1)
                x1 = [i0 - sqrt(radius^2 - (j+0.5-j0)^2), j+0.5];
                x2 = [i0 - sqrt(radius^2 - (j-0.5-j0)^2), j-0.5];
                ln1 = (i+0.5) - x1(1);
                ln2 = (i+0.5) - x2(1);
            % Side top
            elseif corners(i,j+1) && corners(i+1,j+1)
                x1 = [i-0.5, j0 - sqrt(radius^2 - (i-0.5-i0)^2)];
                x2 = [i+0.5, j0 - sqrt(radius^2 - (i+0.5-i0)^2)];
                ln1 = (j+0.5) - x1(2);
                ln2 = (j+0.5) - x2(2);
            end
            % Compute area
            if ln1 == ln2
                weights(i,j) = ln1 + segment(radius, x1 - center, x2 - center);
            elseif ln1 > ln2
                weights(i,j) = ln2 + 0.5 * (ln1 - ln2) + segment(radius, x1 - center, x2 - center);
            elseif ln1 < ln2
                weights(i,j) = ln1 + 0.5 * (ln2 - ln1) + segment(radius, x1 - center, x2 - center);
            end
            
        % If the noCorners is 3
        elseif noCorners(i,j) == 3
            % bottom left
            if ~corners(i,j)
                x1 = [i-0.5, j0 - sqrt(radius^2 - (i-0.5-i0)^2)];
                x2 = [i0 - sqrt(radius^2 - (j-0.5-j0)^2), j-0.5];
                ht = x1(2) - (j-0.5);
                ln = x2(1) - (i-0.5);
            % bottom right
            elseif ~corners(i+1,j)
                x1 = [i+0.5, j0 - sqrt(radius^2 - (i+0.5-i0)^2)];
                x2 = [i0 + sqrt(radius^2 - (j-0.5-j0)^2), j-0.5];
                ht = x1(2) - (j-0.5);
                ln = (i+0.5) - x2(1);
            % top right
            elseif ~corners(i+1,j+1)
                x1 = [i+0.5, j0 + sqrt(radius^2 - (i+0.5-i0)^2)];
                x2 = [i0 + sqrt(radius^2 - (j+0.5-j0)^2), j+0.5];
                ht = (j+0.5) - x1(2);
                ln = (i+0.5) - x2(1);
            % top left
            elseif ~corners(i,j+1)
                x1 = [i-0.5, j0 + sqrt(radius^2 - (i-0.5-i0)^2)];
                x2 = [i0 - sqrt(radius^2 - (j+0.5-j0)^2), j+0.5];
                ht = (j+0.5) - x1(2);
                ln = x2(1) - (i-0.5);
            end
            % Compute area
            weights(i,j) = 1 - 0.5 * ht * ln + segment(radius, x1 - center, x2 - center);
%             disp(['hi ',num2str(i),' ',num2str(j),' ',num2str(ln),' ',num2str(ht),' ',num2str(segment(rad, x1 - cen, x2 - cen)),' ',num2str(weights(i,j)),' '])
        end

    end
end


%% Calculate circular integral
integral = sum(sum(dat .* weights));
error =  pi*radius^2 - sum(sum(weights));

end

function [ area ] = segment( rad, r1, r2 )
    % Find the two angles, given between (-pi,pi]
    theta1 = atan2(r1(2),r1(1)); theta2 = atan2(r2(2),r2(1));
    % Make sure theta2 is greater than theta 1
    if theta1 > theta2
        temp = theta2; theta2 = theta1; theta1 = temp;
    end
    theta = theta2 - theta1;
    % Make sure that the difference is smaller than pi/2. If it isn't, you
    % need to add 2*pi to theta 1 making theta1 larger than pi.
    if theta > 1.6
        theta1 = theta1 + 2*pi;
        theta = theta1 - theta2;
    end
    
    area = (0.5 * theta - 0.5 * sin(theta) ) * rad^2;
end