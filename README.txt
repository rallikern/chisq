README

The production of the Chi-Square Distribution needs a folder with gosia ".out" files, where the parameters are used as the file names (e.g. "out_1.20..21": parmeter = [1.20,0.20]).

The input for the function chisq_plot() and proj_chi() is a list of file names (e.g. outputlist = glob.glob(directory+'total_42_1601/out_*') ) and the parameters, which is only needed for 
the labeling of the axes (e.g. [outputlist,r'$4^+_1||E2||4^+_1$',r'$2^+_1||E2||4^+_1$'])