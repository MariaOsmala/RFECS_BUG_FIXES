function plot_rfecs(chrmatfilename, chrtxtfilename)

    %read xtx file
    txtfd = fopen(chrtxtfilename);
    txt = textscan(txtfd, '%s %d %f');

    %load mat file
    mean_variable = load(chrmatfilename);

    %x/y values for the bar plot
    barx= (txt{2} / 100)-10;
    bary= txt{3};

    matxy = mean_variable.mean_bg_all;

    %plot barplot
    bar(barx, bary, 'r');
    hold on;
    %plot .mat file values
    plot(matxy, 'b');
    plot([0 max(barx)], [0.5 0.5], 'm')
    lg = legend(chrtxtfilename, chrmatfilename, 'Peak calling threshold');

    set(lg, 'interpreter', 'none'); %turn off TeX interpretation
    axis([6.6e5 8e5 0.01 1.01]) %zoom in to a smaller region

    plot([0 max(bary)], [0.5 0.5], 'm')
    
    zoominpos = [6.82e5 6.86e5 0.001 1.001]; %zoom in further to create an inset plot

    %mark the region that we are creating the inset plot at
    rectangle('position', [zoominpos(1), zoominpos(3), zoominpos(2) - zoominpos(1), zoominpos(4) - zoominpos(3)], 'LineWidth', 4);

    %create the inset plot and plot everything again into the inset plot
%     ax = axes('position', [0.40 0.25 0.4 0.4], 'color', [0 0 0]);
%     bar(barx, bary, 'r');
%     hold on;
%     plot(matxy, '-b');
%     lg = legend(chrtxtfilename, chrmatfilename);
%     set(lg, 'interpreter', 'none');
%     axis(zoominpos);
%     set(ax, 'box','on');
%     set(ax, 'XColor', [1 1 1]); %use white X/Y colors for better visibility
%     set(ax, 'YColor', [1 1 1]);
end