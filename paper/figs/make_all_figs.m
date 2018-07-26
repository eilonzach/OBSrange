%% make all
cd 'makeFigures' 
figs_to_make = dir('.');
for i = 1:length(figs_to_make)
    if any(regexp(figs_to_make(i).name,'Figure'))
        fprintf('Making %s ...',figs_to_make(i).name)
        run(figs_to_make(i).name);
        fprintf(' done\n',figs_to_make(i).name)
    end   
end
cd ..
