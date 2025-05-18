function tbl = plot_ca_rate_by_eye_specific_pref(contra_ipsi_binoc, unique_days)

all_spine_size = [];
% eye_specific_type = {'11', '10', '01'};
% binoc_type = {'001', '000'};
elements = {0:1, 0:1, 0:1}; %cell array with N vectors to combine
eye_specific_type = get_combin_str(elements);


for d = 1:length(unique_days)
    temp_day = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.days, unique_days{d}), :);

    for i = 1:length(eye_specific_type)
        spine_size.(['spine',num2str(eye_specific_type{i})]) = temp_day(strcmp(temp_day.cib_spine, eye_specific_type(i)),:).ca_rate;
    end
%     for i = 1:length(binoc_type)
%         spine_size.(['spine',num2str(binoc_type{i})]) = temp_day(strcmp(temp_day.cib_spine, binoc_type(i)),:).ca_rate;
%     end
    all_spine_size = [all_spine_size, spine_size];
end

%% plot spine size by day and by eye-specific feature
all_eye_by_day = [];
all_type_by_day = [];
all_day = [];
for i = 1:length(all_spine_size)
    eye_specific = fieldnames(all_spine_size(i))
    subplot(1,length(all_spine_size),i)
    spine_size = all_spine_size(i);
    all_eye = [];
    all_type = [];
    for ii = 1:length(eye_specific)
        eye = spine_size.(eye_specific{ii});
        all_eye = [all_eye; eye];
        all_type = [all_type; repmat(eye_specific(ii), size(eye))];
        %all_day = [all_day; repmat(eye_specific(ii), size(eye))];
        
        all_eye_by_day = [all_eye_by_day; eye];
        all_type_by_day = [all_type_by_day; repmat(eye_specific(ii), size(eye))];
        all_day = [all_day; repmat(unique_days(i), size(eye))];
        %bar(ii, nanmean(eye)); hold on;
        %errorbar(ii,nanmean(eye), nanstd(eye)/sqrt(length(eye(~isnan(eye)))), '.k');
    end
    boxchart(categorical(all_type),all_eye);
%     xticks([1:length(eye_specific)])
    xticklabels({'U', 'B', 'C', 'B+C', 'I', 'B+I', 'C+I', 'B+C+I'})
    title(unique_days{i})
    if i == 1
        ylabel('Calcium Event Rate')
    end
    ylim([0,0.2])
end
tbl = table(all_eye_by_day, all_type_by_day, all_day);
end
%%
function result = get_combin_str(elements)
     combinations = cell(1, numel(elements)); %set up the varargout result
     [combinations{:}] = ndgrid(elements{:});
     combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
     temp = [combinations{:}];
     result = [];
     for i = 1:size(temp,1)
         result = [result;{strrep(num2str(temp(i,:)),' ', '')}];

     end
     
end