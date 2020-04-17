'''
Powerful tool to produce different kind of 
distributions from Gosia output files
'''
#
#General stuff
#



##############################################
##############################################

import sys
PATH_TO_TP = "/home/ralli/Physik/hieisolde/off_ana/mathematica_ana/"
sys.path.insert(0, PATH_TO_TP)
import numpy as np
import tran_prob as tp



##############################################
##############################################


#
#Nice fonts and plots
#

import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
plt.rcParams.update({'font.size': 20})
plt.rcParams['lines.linewidth'] = 2.25
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
#from mpl_toolkits.axes_grid1.colorbar import colorbar
from matplotlib.colors import LinearSegmentedColormap
COLS = [(1, 1, 0), (124/255, 252/255, 50/255), (65/255, 105/255, 225/255), (1, 1, 1)]  # R -> G -> B
CM0WN = LinearSegmentedColormap.from_list("hallo", COLS, N=1000)


#
#Produces a matrix of chi squares
#

def chisquare(outputlist, transition):
    '''
    The following is valid if all files are in the same folder, which is recommended
    '''

    outputlist = np.array(outputlist)
    result = []
    file_name = outputlist[0]
    #location and length of mother folder, folder follows right after
    folder_end = file_name.rfind('/')
    folder_start = file_name.rfind('/', 0, folder_end-1)
    #folder name of interest, to correlate the result to different inputs
    folder = file_name[folder_start+1:folder_end]
    result.append(folder)

    #number of chi_square/output files
    lenni_tot = len(outputlist)

    x_axis, fac_x, y_axis, fac_y = parameter_func(outputlist, lenni_tot)

    # range of matrix elements of x axis
    x_min = min(x_axis)
    x_max = max(x_axis)

    # range of matrix elements of y axis
    y_min = min(y_axis)
    y_max = max(y_axis)

    #
    #step size and the resulting number of elements for the x and y axis
    #
    #number of elements determines the size of the chi_matrix and
    #

    step_x = int(tp.stepseeker(x_axis))
    if step_x == 0:
        step_x = 1

    steps_x = int((x_max-x_min)/step_x+1)

    step_y = int(tp.stepseeker(y_axis))
    if step_y == 0:
        step_y = 1

    steps_y = int((y_max-y_min)/step_y+1)

    #
    #Actual filling of the matrix
    #
    x_param = [x_axis, x_min, step_x, steps_x]
    y_param = [y_axis, y_min, step_y, steps_y]
    power_matrix = chisqmatrix(outputlist, lenni_tot, *x_param, *y_param, transition)

    chi_matrix, m1_mx = power_matrix[:2]

    #
    #After the matrix was filled, the Chi plot is prepared
    #

    #Determination of the smallest Chi
    mimi = np.min(chi_matrix[np.nonzero(chi_matrix)])

    #Location of smallest Chi in the matrix can be referred to the linked parameters
    for i in range(len(chi_matrix)):
        for j in range(len(chi_matrix[i])):
            if np.isnan(chi_matrix[i][j]):
                print(i, j, 'YES')
    chi_min_loc = [[(x_min+j*step_x)/fac_x, (y_min+i*step_y)/fac_y]
                   for i in range(steps_y)
                   for j in range(steps_x)
                   if chi_matrix[i][j] <= mimi and chi_matrix[i][j] != 0]

    #
    #cut through the minimum for a 1d projection
    #
    loc_y = int((chi_min_loc[0][0]*fac_x-x_min)/step_x)
    loc_x = int((chi_min_loc[0][1]*fac_y-y_min)/step_y)

    #cut through loc_x and loc_y
    proj_x = [chi_matrix[i][loc_y] for i in range(steps_y)]
    proj_y = [chi_matrix[loc_x][i] for i in range(steps_x)]

    #Determination of the one sigma area (chi_square_min+1)
    s_area_x = [(x_min + j*step_x)/fac_x
                for i in range(steps_y)
                for j in range(steps_x)
                if chi_matrix[i][j] <= mimi+1 and chi_matrix[i][j] != 0]

    s_area_y = [(y_min + i*step_y)/fac_y
                for i in range(steps_y)
                for j in range(steps_x)
                if chi_matrix[i][j] <= mimi+1 and chi_matrix[i][j] != 0]

    result.append(chi_min_loc)
    result.append([[np.min(s_area_x), np.max(s_area_x)], [np.min(s_area_y), np.max(s_area_y)]])

    #M1 minimum
    mini = np.min(m1_mx[np.nonzero(m1_mx)])

    print(result)
    x_out = [x_min, x_max, step_x, fac_x]
    y_out = [y_min, y_max, step_y, fac_y]
    placeh = [m1_mx, mini, result, steps_x, steps_y]
    pros = [proj_x, proj_y]
    return [chi_matrix, mimi, *x_out, *y_out, *placeh, *pros, *power_matrix[2:]]



##############################################
##############################################


def parameter_func(outputlist, lenni_tot):
    '''
    seperate the determination of simple
    parameters from the main function
    '''
    x_axis = []
    y_axis = []
    for i in range(lenni_tot):
        lenni = len(outputlist[i]) #length of filename
        marker1 = outputlist[i].rfind('out_') #location of first parameter
        l_out = 4 #(+4 = length of string 'out_')

        #searches forthe second parameter for different output styles

        marker2 = outputlist[i].rfind('..')# If the second parameter is between 0<par2<1

        if marker2 == -1:
            marker2 = outputlist[i].rfind('.-')# If the second parameter: par2<0

        if marker2 == -1:
            marker2 = outputlist[i].rfind('.0', marker1+10)# If the second parameter: par2=0

        k = 1
        while marker2 == -1:
            marker2 = outputlist[i].rfind('.'+str(k)+'.') #If the second parameter: 1<=par2
            k += 1

        value_x = outputlist[i][marker1+l_out:marker2] # first parameter x
        value_y = outputlist[i][marker2+1:lenni] # second parameter y
        marker3 = outputlist[i].find('.', marker1) # location of '.' in first parameter
        length_x = marker2-marker1+l_out #length of first parameter

        if length_x == 1: # if it is one, the parameter equals 0
            fac_x = 1
        else:
            power = marker2-marker3 # otherwise it is multiplied to get an integer.
            fac_x = 10**power

        if float(value_y) < 0:
            fac_y = 10**(len(value_y) - 2) # factor needed to get an integer also for parameter y
        else:
            fac_y = 10**(len(value_y) - 1) # factor needed to get an integer also for parameter y

        x_temp = round(fac_x*float(value_x)) # x_value of file
        x_axis.append(x_temp)
        y_temp = round(fac_y*float(value_y)) # y_value of file
        y_axis.append(y_temp)

    return x_axis, fac_x, y_axis, fac_y

##############################################
##############################################

def chisqmatrix(outl, lt, x_axis, x_min, step_x, steps_x, y_axis, y_min, step_y, steps_y, trans):
    '''
    Production of the chisq matrices
    '''


    chi_matrix = np.zeros((steps_y, steps_x))

    m1_mx = np.zeros((steps_y, steps_x))

    #
    #Individiual chisq_dist for each experiment
    #

    #number of experiments
    experiment = True
    exp = 1
    while experiment == True:
        string_iter = "EXPERIMENT  "+str(exp)
        try:
            tp.file_checker_max(string_iter, outl[0])
            exp += 1
        except:
            experiment = False
            exp -= 1
    # exp individual chisq_dist for a certain transition
    chi_indi = np.zeros((exp, steps_y, steps_x))

    for i in range(lt):
        #
        #Searching for lowest chisq in the output file (file name = outl[i]):
        #
        #last appearance of "CHISQ=" in the output file
        chi_line = tp.file_checker_max("CHISQ=", outl[i])
        chi = np.genfromtxt(outl[i], skip_header=chi_line-1, max_rows=1, usecols=(2), dtype=None)
        #if there is a problem with finding a min and chisq= NaN, it has to be corrected
        if np.isnan(chi):
            chi = 0
        #
        #Calculation of fit points: UPL + Experimental Yields
        #+ additional observables (branching or delta)
        #
        num_upl = tp.file_checker_num("UPL!", outl[i])#number of "UPL!" in the output file

        #Looking for his line: "k EXPERIMENTAL YIELDS   J MATRIX ELEMENTS TO BE VARIED":
        line_fpar = tp.file_checker_max("MATRIX ELEMENTS TO BE VARIED", outl[i])
        # The first column of this line is the number of experimental yields:
        fparam0 = np.genfromtxt(outl[i], skip_header=line_fpar-1, max_rows=1, usecols=(0), dtype=None)
        #
        #The number of included observables is determined through the number of lines
        #from "MATRIX ELEMENTS TO BE VARIED" (line_free_param) to
        #"*****  RESTART-MATRIX ELEMENTS OVERWRITTEN  *****" (line_free_param_obs_end)
        #
        # It is dependend on the variety of the observables (life times are not included,yet)
        try:
            line_fparam_start = tp.file_checker_min("BRANCHING RATIOS", outl[i])
            free_lines = 25
        except:
            line_fparam_start = tp.file_checker_max("EXPERIMENTAL E2/M1 MIXING RATIOS", outl[i])
            free_lines = 19

        line_fparam_end = tp.file_checker_max("*****  RESTART-MATRIX ELEMENTS OVERWRITTEN  *****", outl[i])
        # many lines of std output and free spaces:
        #(19 from mixing ot restart, 25 from branching to restart)
        free_param1 = line_fparam_end - line_fparam_start - free_lines
        #
        #final number of experimental observables (degrees of freedom)
        #
        dof = fparam0 + free_param1 + num_upl
        #Chi matrix
        #here, x and y is exchanged...
        #the first indice of the me is y and the second is x
        #
        #And the matrix starts withe the highest value [ y_max, y_max-1,...,y_min]
        chi_matrix[int((y_axis[i]-y_min)/step_y), int((x_axis[i]-x_min)/step_x)] = float(chi*(dof))
        #
        #Searching for reduced M1 strength 2+_3 to 2+_1, especially designed for 140Nd
        #
        try:
            lnl = tp.file_checker_min("2.0            2.1400", outl[i]) - 1
        except: #Searching for reduced M1 strength 2+_4 to 2+_1 if 2_3 is not included
            lnl = tp.file_checker_min("2.0            2.3330", outl[i]) - 1
        lvl_nr = str(np.genfromtxt(outl[i], skip_header=lnl, max_rows=1, usecols=(0), dtype=None))
        search_string = '2      '+lvl_nr
        m1_line = tp.file_checker_max(search_string, outl[i])
        m1_str = np.genfromtxt(outl[i], skip_header=m1_line-1, max_rows=1, usecols=(3), dtype=None)
        if np.isnan(m1_str):
            m1_str = 0
        #M1-Matrix
        m1_mx[int((y_axis[i]-y_min)/step_y), int((x_axis[i]-x_min)/step_x)] = abs(float(m1_str))
        #
        #Searching for individual chisq distribution
        #only for Detector 1 for each Experiment
        #

        string_start = ['CALCULATED AND EXPERIMENTAL YIELDS   EXPERIMENT  '+str(int(k))+' DETECTOR  1' for k in range(1,exp+1)]
        start, end = np.zeros(exp), np.zeros(exp)

        for k,j in enumerate(string_start):
            start[k] = tp.file_checker_max(j, outl[i])
            end[k] = tp.file_checker_min_cond("CHISQ SUBTOTAL =", outl[i], start[k])

        #
        #Reading of the deviation of exp to calc of transition yields
        #

        #columns of transition yields tab
        uc = (0, 1, 2, 3, 4, 5, 6, 7, 8)
        max_rows = int(end[k]-start[k]-4)
        yields = np.zeros((exp, max_rows, len(uc)))
        #chisq of a special transition for each exp
        chisq_i = np.zeros(exp)
        num_plots = 1
        for k in range(exp):
            test = np.genfromtxt(outl[i], skip_header=int(start[k]+2), usecols=uc, max_rows=max_rows, dtype=None)
            test = [[test[j][p] for p in range(len(uc))] for j in range(max_rows)]
            yields[k] = np.array(test)
            for j in yields[k]:
                if j[0] == trans[0] and j[1] == trans[1]:
                    chisq_i[k] = j[8]
                    num_plots += 1
                    if np.isnan(chisq_i[k]):
                        chisq_i[k] = 10
                try:
                    chi_indi[k][int((y_axis[i]-y_min)/step_y),
                                int((x_axis[i]-x_min)/step_x)] = abs(float(chisq_i[k]))
                except:
                    num_plots = 1

    indi_min = [np.min(np.abs(chi_indi[k])) for k in range(exp)]
    #chi_indi2[int((y_axis[i]-y_min)/step_y)][int((x_axis[i]-x_min)/step_x)] = abs(float(chisq_i2))

    return chi_matrix, m1_mx, chi_indi, num_plots, indi_min


def matrix_manipulation(cs, val_area, interpolation):
    '''
    !!!!!!Alle eintrage die gleich 0 sind werden
    größer als minimum + val_area gesetzt!!!!!!
    is necessary for this kind of plotting...
    if interpolation = True:
    All zeros between reals Values will be replaced
    by a mean value
    '''

    interpolation = interpolation
    for case in [cs[0], cs[17][0], cs[17][1]]:
        for i in range(cs[14]):
            for j in range(cs[13]):
                #vertical interpolation
                if cs[0][i][j] == 0 and interpolation:
                    kii, pii = i, i
                    while cs[0][kii][j] == 0 and kii > 0:
                        kii -= 1
                    while pii < cs[14]-1 and cs[0][pii][j] == 0:
                        pii += 1
                    if cs[0][kii][j] != 0 and cs[0][pii][j] != 0:
                        cs[0][i][j] = cs[0][kii][j]*(abs(kii-i)/abs(kii-pii))+cs[0][pii][j]*(abs(pii-i)/abs(kii-pii))

        for i in range(cs[14]):
            for j in range(cs[13]):
                #horizontal interpolation
                if cs[0][i][j] == 0 and interpolation:
                    kjj, pjj = j, j
                    while cs[0][i][kjj] == 0 and kjj > 0:
                        kjj -= 1
                    while pjj < cs[13]-1 and cs[0][i][pjj] == 0:
                        pjj += 1
                    if cs[0][i][kjj] != 0 and cs[0][i][pjj] != 0:
                        cs[0][i][j] = cs[0][i][kjj]*(abs(kjj-j)/abs(kjj-pjj))+cs[0][i][pjj]*(abs(pjj-j)/abs(kjj-pjj))

        for i in range(cs[14]):
            for j in range(cs[13]):
                #diagonal interpolation number1
                if cs[0][i][j] == 0 and interpolation:
                    kjj, pjj = j, j
                    kii, pii = i, i
                    while cs[0][kii][kjj] == 0 and kjj > 0 and kii > 0:
                        kjj -= 1
                        kii -= 1
                    while pjj < cs[13]-1 and pii < cs[14]-1 and cs[0][pii][pjj] == 0:
                        pjj += 1
                        pii += 1
                    if cs[0][kii][kjj] != 0 and cs[0][pii][pjj] != 0:
                        cs[0][i][j] = cs[0][kii][kjj]*(abs(kjj-j)/abs(kjj-pjj))+cs[0][pii][pjj]*(abs(pjj-j)/abs(kjj-pjj))

        for i in range(cs[14]):
            for j in range(cs[13]):
                #diagonal interpolation number2
                if case[i][j] == 0 and interpolation:
                    kjj, pjj = j, j
                    kii, pii = i, i
                    while case[kii][kjj] == 0 and kjj > 0 and kii < cs[14]-1:
                        kjj -= 1
                        kii += 1
                    while pjj < cs[13]-1 and pii > 0 and case[pii][pjj] == 0:
                        pjj += 1
                        pii -= 1
                    if case[kii][kjj] != 0 and cs[0][pii][pjj] != 0:
                        case[i][j] = case[kii][kjj]*(abs(kjj-j)/abs(kjj-pjj))+case[pii][pjj]*(abs(pjj-j)/abs(kjj-pjj))
    
    for case in [cs[0], cs[17][0], cs[17][1]]:
        for i in range(cs[14]):
            for j in range(cs[13]):
                if case[i][j] == 0:
                    case[i][j] = cs[1] + 2*val_area

    return cs[0], cs[17]

##############################################
##############################################
#End of matrix manipulation
##############################################
##############################################


#
#Produces the final chi square plot
#

def chisq_plot(out_list, interpolation=False, val_area=1, limitation=[0, 0, 0], transition=[0, 0]):
    '''
    Actual Plotting
    '''
    from PIL import Image, ImageDraw

    for outi in out_list:

        #labeling of x and y axis

        xlabel = r'$2^+_1||E2||4^+_1$' #default
        ylabel = r'$0^+_1||E2||2^+_3$'

        try:
            #works, if x and y label are specified in the input
            cs = chisquare(outi[0], transition)
            xlabel = outi[1]
            ylabel = outi[2]
        except:
            cs = chisquare(outi, transition)     #if x and y label are NOT specified in the input
        val_area = val_area #Plot area (min+val_area), val_area =1 mwans one sigma area
        cs[0], cs[17] = matrix_manipulation(cs, val_area, interpolation)
        #Start of plotting:

        #num_plots = cs[18] = 1 +i if there are 1+ i plots(1 global chisq and i individual),
        #1 if there is only ony
        fig, axe = plt.subplots(cs[18], 1, sharex=True, figsize=(8, cs[18]*5))
        #basewidth = 200
        #
        ##additional option to include another image, i.e. a level scheme
        #
        #image = Image.open('/home/rkern/Physik/hieisolde/pictures/140nd_ls_subset_py.png')
        #wpercent = (basewidth/float(image.size[0]))
        #hsize = int((float(image.size[1])*float(wpercent)))
        #image = image.resize((basewidth,hsize))
        #fig.figimage(image,80,60 ,zorder=5)

        plotc = [cs[0], *cs[17]]
        plotmin = [cs[1], *cs[19]]

        print(plotmin)
        #range of y and x axis
        extent = ((cs[2]-cs[4]/2)/cs[5], (cs[3]+cs[4]/2)/cs[5], (cs[6]-cs[8]/2)/cs[9], (cs[7]+cs[8]/2)/cs[9])
        if cs[18] != 1:
            plot_len = len(axe)
            for i in range(plot_len):
                im1 = axe[i].imshow(plotc[i], cmap=CM0WN, origin='lower', aspect='auto',
                                    Norm=Normalize(vmin=plotmin[i], vmax=plotmin[i]+val_area),
                                    interpolation='nearest', extent=extent, zorder=0)
                cmp = fig.colorbar(im1, ax=axe[i])
                cmp.set_label('χ²', size=24)
                axe[i].set_xlabel(r'$\langle$'+xlabel+r'$\rangle$ (eb)', fontsize=22)
                #ax1.set_xlim(1.01,1.32)
                axe[i].set_ylabel(r'$\langle$'+ylabel+r'$\rangle$ (eb)', fontsize=22)
                #ax1.set_ylim(0.05,0.34)
                axe[i].set_title(cs[12][0])
        else:
            im1 = axe.imshow(plotc[0], cmap=CM0WN, origin='lower', aspect='auto',
                             Norm=Normalize(vmin=cs[1], vmax=cs[1]+val_area),
                             interpolation='nearest', extent=extent, zorder=0)
            cmp = fig.colorbar(im1, ax=axe)
            cmp.set_label('χ²', size=24)
            axe.set_xlabel(r'$\langle$'+xlabel+r'$\rangle$ (eb)', fontsize=22)
            axe.set_xlim(-1.3, 0)
            axe.set_ylabel(r'$\langle$'+ylabel+r'$\rangle$ (eb)', fontsize=22)
            #ax1.set_ylim(0,0.34)
            #axe.set_title(cs[12][0])

        if limitation[0] == 1:
            plt.vlines(limitation[1], extent[2], extent[3], color='black', alpha=0.5)
            plt.vlines(limitation[1]-limitation[2], extent[2], extent[3], linestyles='--', alpha=0.5)
            plt.vlines(limitation[1]+limitation[2], extent[2], extent[3], linestyles='--', alpha=0.5)
        fig.savefig("Plots/"+cs[12][0]+".pdf", bbox_inches='tight')
        plt.show()
    return

##############################################
##############################################

#
#Plots of projections through minimum
#

def proj_chi(out_list):
    '''
    Projection on x- and y-axis through the minimum
    '''
    for outi in out_list:
        try:
            #if label names are in the input, we need the 0th argument of the input
            cs = chisquare(outi[0])
            xlabel = outi[1]
            ylabel = outi[2]
        except:
            cs = chisquare(outi)
            xlabel = '$2^+_3||E2||2^+_3$'
            ylabel = '$2^+_3||E2||0^+_1$'
        #
        #Plot starts
        #
        fig, ax = plt.subplots(2, 1, sharey=False, sharex=False, figsize=(8, 5))
        plt.subplots_adjust(hspace=.5)

        xdata = np.linspace(cs[2]/cs[5], cs[3]/cs[5], cs[13])

        #Interpolation, if a value is zero
        ydata = [(cs[16][i-1]+cs[16][i+1])/2
                 if cs[16][i] == 0 and i+1 < cs[13] and i-1 >= 0
                 else cs[16][i]
                 for i in range(cs[13])]
        fig.suptitle(cs[12][0])

        #ax[0].set_xlim(1.04,1.26)
        #ax[0].set_ylim(0.0001,1.3)
        ax[0].set_ylabel(r'$\langle$'+ylabel+r'$\rangle$', fontsize=22)
        ax[0].set_xlabel(r'$\langle$'+xlabel+r'$\rangle$', fontsize=22)
        ax[0].set_yscale('linear')
        ax[0].plot(xdata, ydata)

        xdata = np.linspace(cs[6]/cs[9], cs[7]/cs[9], cs[14])

        #Interpolation, if a value is zero
        ydata = [(cs[15][i-1]+cs[15][i+1])/2 if cs[15][i] == 0 and i+1 < cs[14] and i-1 >= 0 else cs[15][i] for i in range(cs[14])]

        #ax[1].set_xlim(0.08,0.3)
        #ax[1].set_ylim(0.0001,2.3)
        ax[1].set_xlabel(r'$\langle$'+ ylabel+r'$\rangle$', fontsize=22)
        ax[1].set_ylabel(r'$\chi^2$', fontsize=22)
        ax[1].set_yscale('linear')
        ax[1].plot(xdata, ydata)
        plt.show()
    return None

##############################################
##############################################

#
#Function for visualizing chisqsurface .txt output
#not more, not less
#

def chisqsurface(filename):
    '''
    Import chisqsurface output to python
    '''
    chissi = np.loadtxt(filename)
    lec = len(chissi)
    x_1d = []
    for i in range(lec):
        if chissi[i][0] not in x_1d:
            x_1d.append(chissi[i][0])
    lex = len(x_1d)

    y_1d = []
    for i in range(lec):
        if chissi[i][1] not in y_1d:
            y_1d.append(chissi[i][1])
    ley = len(y_1d)

    # Argumentwerte als 2D Arrays erzeugen
    x_2d, y_2d = np.meshgrid(x_1d, y_1d)

    # Interessante Daten erzeugen
    z_2d = [[chissi[i+j*ley][2] for j in range(lex)] for i in range(ley)]
    z_min = np.min(z_2d)

    # Plotten
    fig = plt.figure(figsize=(8, 8))
    plot = plt.pcolormesh(x_2d, y_2d, z_2d, cmap=CM0WN, norm=Normalize(vmin=z_min, vmax=z_min+1))
    cmp = fig.colorbar(plot)
    cmp.set_label('χ²', size=24)
    #plt.gca().set_aspect("equal") # x- und y-Skala im gleichen Maßstaab
    plt.xlabel(r'$\langle2^+_3||E2||0^+_1\rangle$', fontsize=26)
    plt.ylabel(r'$\langle2^+_1||E2||4^+_1\rangle$', fontsize=26)
    plt.show()
    return
