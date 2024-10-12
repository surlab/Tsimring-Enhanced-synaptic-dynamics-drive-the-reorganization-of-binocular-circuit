function [h] = plot_SEM(time,data,sem,color_data,color_sem)
%plot_SEM plot data with a shaded sem
%time is time vector. data is data vector. sem is sem vector.
%Kyle jenks, 2021-10-22

%Curves are data plus or minus the SEM
curve1 = data' + sem';
curve2 = data' - sem';

%double time to plot a shape (inbetween) of the plus or minus curves
x2 = [time, fliplr(time)];
inbetween = [curve1, fliplr(curve2)];

%Only keep values that are numbers
keepIndex = ~isnan(x2) & ~isnan(inbetween);
x2 = x2(keepIndex);
inbetween = inbetween(keepIndex);

%Plot the filled shape
g = fill(x2, inbetween, color_sem,'facealpha',0.5);
set(g,'EdgeColor','none')
hold on
%plot the data
%h = plot(time, data, 'Color', color_data, 'LineWidth', 2);


end

