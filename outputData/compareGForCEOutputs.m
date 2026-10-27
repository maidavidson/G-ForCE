function report = compareGForCEOutputs(referenceFile, candidateFile, tolerance)
%COMPAREGFORCEOUTPUTS Compare G_ForCE results at matching saved dates.
% report = compareGForCEOutputs('reference.mat','candidate.mat',0)
% Dates present in only one run are reported, not silently ignored.
% Profiles and dissipation are expected as [spatial_points x saved_dates].
if nargin < 3, tolerance = 0; end
assert(isfile(referenceFile) && isfile(candidateFile), 'Both MAT files must exist.');
assert(isnumeric(tolerance) && isscalar(tolerance) && tolerance >= 0, ...
    'Tolerance must be a nonnegative scalar.');
aInfo = whos('-file', referenceFile);
bInfo = whos('-file', candidateFile);
aYears = localYears({aInfo.name});
bYears = localYears({bInfo.name});
years = union(aYears, bYears);
report = struct('year',{},'commonDates',{},'referenceOnlyDates',{}, ...
    'candidateOnlyDates',{},'profileMaxAbsDiff',{},'dissMaxAbsDiff',{}, ...
    'profileEqual',{},'dissEqual',{},'status',{});
for k = 1:numel(years)
    y = years(k);
    report(k).year = y;
    report(k).commonDates = 0;
    report(k).referenceOnlyDates = [];
    report(k).candidateOnlyDates = [];
    report(k).profileMaxAbsDiff = NaN;
    report(k).dissMaxAbsDiff = NaN;
    report(k).profileEqual = false;
    report(k).dissEqual = false;
    report(k).status = 'not compared';
    if ~ismember(y,aYears) || ~ismember(y,bYears)
        report(k).status = 'year absent from one run';
        fprintf('%d: %s\n',y,report(k).status);
        continue
    end
    vars = {sprintf('Profiles_%d',y), sprintf('Dates_%d',y), sprintf('Diss_%d',y)};
    a = load(referenceFile,vars{:});
    b = load(candidateFile,vars{:});
    for j = 1:3
        assert(isfield(a,vars{j}) && isfield(b,vars{j}), ...
            'Missing annual variable %s',vars{j});
    end
    tA = a.(vars{2})(:); tB = b.(vars{2})(:);
    assert(all(isfinite(tA)) && all(isfinite(tB)), 'Nonfinite dates in %d',y);
    assert(numel(unique(tA))==numel(tA) && numel(unique(tB))==numel(tB), ...
        'Duplicate dates in year %d; inspect the outputs.',y);
    pA=a.(vars{1}); pB=b.(vars{1}); dA=a.(vars{3}); dB=b.(vars{3});
    assert(size(pA,2)==numel(tA) && size(pB,2)==numel(tB) && ...
           size(dA,2)==numel(tA) && size(dB,2)==numel(tB), ...
           'Date-to-column mismatch in year %d',y);
    assert(size(pA,1)==size(pB,1) && size(dA,1)==size(dB,1), ...
           'Spatial dimensions differ in year %d',y);
    [tCommon,ia,ib] = intersect(tA,tB,'stable');
    report(k).commonDates = numel(tCommon);
    report(k).referenceOnlyDates = setdiff(tA,tB,'stable');
    report(k).candidateOnlyDates = setdiff(tB,tA,'stable');
    if isempty(tCommon)
        report(k).status = 'no common dates';
        fprintf('%d: no common dates\n',y);
        continue
    end
    [report(k).profileMaxAbsDiff,report(k).profileEqual] = ...
        localCompare(pA(:,ia),pB(:,ib),tolerance);
    [report(k).dissMaxAbsDiff,report(k).dissEqual] = ...
        localCompare(dA(:,ia),dB(:,ib),tolerance);
    if report(k).profileEqual && report(k).dissEqual
        report(k).status='matched values agree';
    else
        report(k).status='NUMERICAL DIFFERENCE';
    end
    fprintf(['%d: common=%d, reference-only=%d, candidate-only=%d, ', ...
        'max profile diff=%g, max diss diff=%g : %s\n'], ...
        y,numel(tCommon),numel(report(k).referenceOnlyDates), ...
        numel(report(k).candidateOnlyDates),report(k).profileMaxAbsDiff, ...
        report(k).dissMaxAbsDiff,report(k).status);
end
r=load(referenceFile,'model'); c=load(candidateFile,'model');
thresholdMatch=true;
if isfield(r,'model') && isfield(c,'model') && ...
   isfield(r.model,'IrThreshold') && isfield(c.model,'IrThreshold')
    thresholdMatch=isequaln(r.model.IrThreshold,c.model.IrThreshold);
    fprintf('IrThreshold: reference=%g, candidate=%g; equal=%d\n', ...
        r.model.IrThreshold,c.model.IrThreshold,thresholdMatch);
else
    fprintf('IrThreshold could not be checked in both files.\n');
end
compared = [report.commonDates] > 0;
allEqual = any(compared) && all([report(compared).profileEqual]) && ...
    all([report(compared).dissEqual]) && thresholdMatch;
if allEqual
    fprintf('PASS: numerical values agree at all %d matched dates (tolerance %g).\n', ...
        sum([report.commonDates]),tolerance);
else
    warning('G_ForCE:ComparisonFailed', ...
        'Numerical mismatch or no comparable dates. Inspect report.');
end
if any(arrayfun(@(s) ~isempty(s.referenceOnlyDates) || ...
        ~isempty(s.candidateOnlyDates),report)) || ~isequal(aYears,bYears)
    fprintf('NOTE: date coverage differs. PASS applies only to matched dates.\n');
end
end

function years=localYears(names)
    tokens=regexp(names,'^Profiles_(\d{4})$','tokens','once');
    tokens=tokens(~cellfun('isempty',tokens));
    years=sort(cellfun(@(x)str2double(x{1}),tokens));
end

function [d,equal]=localCompare(a,b,tolerance)
    if ~isequaln(isnan(a),isnan(b)) || ~isequaln(isinf(a),isinf(b)) || ...
            any(a(isinf(a))~=b(isinf(b)))
        d=Inf; equal=false; return
    end
    mask=isfinite(a)&isfinite(b);
    if any(mask(:)),d=max(abs(a(mask)-b(mask)));else,d=0;end
    equal=d<=tolerance;
end
