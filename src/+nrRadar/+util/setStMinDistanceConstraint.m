function index = setStMinDistanceConstraint(stPosition, nPerGroup, minDist)
% Make sequential groups of nPerGroup points with pairwise spacing >= minDist.
% Points can reappear in later groups.

N = size(stPosition,1);
groups = {};
i = 1;

while i <= N
    this = zeros(1, nPerGroup);
    this(1) = i;
    count = 1;
    j = i + 1;

    while j <= N && count < nPerGroup
        p = stPosition(j,:);
        d = vecnorm(stPosition(this(1:count),:) - p, 2, 2);
        if all(d >= minDist)
            count = count + 1;
            this(count) = j;
        end
        j = j + 1;
    end

    if count < nPerGroup
        error('Could not form full group %d: only %d/%d points satisfy %.1f m spacing.', ...
              numel(groups)+1, count, nPerGroup, minDist);
    end

    groups{end+1} = this; %#ok<AGROW>
    i = i + nPerGroup;  % move to next starting point
end

index = cell2mat(groups);
end
