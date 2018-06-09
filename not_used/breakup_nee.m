
for y = 1980:2018; tic
    y
    i = find(date(1,:)==y);
    yearData = nee(:,:,i);
    fname = strcat('./data/nee_year_',num2str(y),'.mat');
    save(fname,'yearData','-v7.3');
    toc
end