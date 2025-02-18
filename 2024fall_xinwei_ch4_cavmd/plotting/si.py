import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import PIL.Image as image

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
	x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
	the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

def get_s1():
    ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=12)
    xs, ys = [], []
    data = np.loadtxt(f'./plotting_data/outsize_pac_gaussian.out')
    for e0 in range(0,8,2):
        x, y = data[:,9], data[:,e0]/1e28
        interval = np.where(x<3300)[0]
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0*0.5
        xs.append(x[interval])
        ys.append(y[interval])

    clp.plotone(xs, ys, ax, lw=1.5, 
                xlim=[1000, 3300], ylim=[-0.2,4], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    
    ax.axvline(x=2902.6, linestyle='-.', alpha=0.2)
    ax.axvline(x=1509.7, linestyle='-.', alpha=0.2)
    ax.axvline(x=3004.3, linestyle='-.', alpha=0.2)
    ax.axvline(x=1312.8, linestyle='-.', alpha=0.2)
    ax.text(2875, 4.1, "$v_1$", fontsize=8, color="0.5")
    ax.text(1480, 4.1, "$v_2$", fontsize=8, color="0.5")
    ax.text(2975, 4.1, "$v_3$", fontsize=8, color="0.5")
    ax.text(1285, 4.1, "$v_4$", fontsize=8, color="0.5")

    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./si_figure/s1.png")

def get_s2():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan', 'grey']
    e0list = [1,3,4,5,6]
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    UP_freq = [1362.9, 1409.6, 1459.6, 1511.3, 1561.4, 1619.8, 1668.2, 1719.9, 1769.9, 1818.3] # omega_c=1311.2
    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    kt2au = 3.483e-4
    for i in range(len(e0list)):
        datanoneq = np.loadtxt(f'./plotting_data/photon_110k_6e-3/photon_110k_density_cavity_n400v0_E0_{e0list[i]}e-4.out')
        xxanoneq, xybnoneq = datanoneq[:,6]*ang2au, datanoneq[:,10]*ang2au
        vxanoneq, vybnoneq = datanoneq[:,12], datanoneq[:,16]
        energy = (0.5 * (vxanoneq ** 2 + vybnoneq ** 2) + 0.5 * (xxanoneq ** 2 + xybnoneq ** 2) * (UP_freq[i] * cmn2au) ** 2) / (kt2au*400)
        row, col = np.shape(datanoneq)
        te  = np.linspace(0, 20, row)
        x0, y0 = i // 5, i % 5
        axes[i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[i].set_xticks([0,1,2,3,4,5])
        axes[i].set_yticks([0,3,6,9])
        clp.plotone([te], [energy], axes[i], colors=[color_list[8]], labels=['photon'], lw=1, showlegend=True if i == 4 else False, legendloc=(0.6, 0.5),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if y0 == 0 else None, xlim=[0, 5], ylim=[-0.1,9], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1362.9, 1459.6, 1511.3, 1561.4, 1619.8] # omega_c=1311.2
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(5):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(0.05, y0-0.02, coupling1List[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        
    clp.adjust(savefile=f'./si_figure/s2.png', tight_layout=False)

def get_s3():
    #blue_palette = ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519C', '#08306b']
    #red_palette  = ['#fff5f0', '#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15', '#67000d']
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    labels = [r"$\widetilde{\varepsilon}=%d\times 10^{-4}$ a.u." %n for n in range(1,11)]
    coord_label = [r"S$_{1}$", r"S$_{2a}$", r"S$_{2b}$", r"S$_{3x}$", r"S$_{3y}$", r"S$_{3z}$", r"S$_{4x}$", r"S$_{4y}$", r"S$_{4z}$"]
    e0list = [1,3,4,5,6]
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    for i in range(len(e0list)):
        for size in range(1):
            data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3/noneqcoord_110k_density_cavity_n400v{size}_E0_{e0list[i]}e-4_6e-3.out')
            data2 = np.loadtxt(f'./plotting_data/110k_eqcoord/coord_110k_density_cavity_n400v{size}_E0_{e0list[i]}e-4.out')
            data1 = data1[::4,:]
            data2 = data2[::4,:]
            row, col = np.shape(data2)
            te  = np.linspace(0, 20, row)
            ref = np.mean(data2, axis=0)
            y1  = 1 * (data1[:,0] / ref[0] - 1)
            y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
            y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
            y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
            xs = [te, te, te, te]
            ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
            colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
            labels_local = [r"v$_%d$-$632$ mJ/cm$^2$" %i for i in range(1,5)]
            axes[i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
            axes[i].set_yticks([0,4,8,12,16])
            axes[i].set_xticks([0,1,2,3,4,5])
            clp.plotone(xs, ys, axes[i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 0 else False, legendloc=(0.5, 0.2),
                        xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 2.5], ylim=[0,16], legendFontSize=4.5)

    for i in range(len(e0list)):
        for size in range(1):
            data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_3e-3/noneqcoord_110k_density_cavity_n400v{size}_E0_{e0list[i]}e-4_3e-3.out')
            data2 = np.loadtxt(f'./plotting_data/110k_eqcoord/coord_110k_density_cavity_n400v{size}_E0_{e0list[i]}e-4.out')
            data1 = data1[::4,:]
            data2 = data2[::4,:]
            row, col = np.shape(data2)
            te  = np.linspace(0, 20, row)
            ref = np.mean(data2, axis=0)
            y1  = 4 * 1 * (data1[:,0] / ref[0] - 1)
            y2  = 4 * 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
            y3  = 4 * 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
            y4  = 4 * 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
            xs = [te, te, te, te]
            ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
            colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
            labels_local = [r"4$\times$v$_%d$-$158$ mJ/cm$^2$" %i for i in range(1,5)]
            local_linestyles = ['-','-','-','-','--','--','--','--']
            clp.plotone(xs, ys, axes[i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 0 else False, legendloc=(0.2, 0.175),
                        xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 5], ylim=[0,16], legendFontSize=4, linestyles=local_linestyles, alpha=0.5)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1362.9, 1459.6, 1511.3, 1561.4, 1619.8]
    couplingList = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(len(e0list)):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(x0-0.14, y0-0.02, couplingList[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='right', color="k")
           
    clp.adjust(savefile=f'./si_figure/s3.png', tight_layout=False)

def get_s4():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    for i in range(5):
        mol = 100*(2**i)
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_nmol/noneqcoord_6e-3_{mol}_E0_4e-4.out')
        data2 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_nmol/eqcoord_6e-3_{mol}_E0_4e-4.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        axes[i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[i].set_xticks([0,5,10,15,20])
        axes[i].set_yticks([0,4,8,12])
        clp.plotone(xs, ys, axes[i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 4 else False, legendloc=(0.7, 0.35),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 20], ylim=[-0.1,12], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1511.3, 1511.3, 1511.3, 1511.3, 1511.3] # for density figure
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    labels = [r"$N_{\rm{simu}} = $ %d" %(100*(2**i)) for i in range(5)]
    coupling1List = [coupling1List[i] + "\n" + labels[i] for i in range(5)]
    for i in range(5):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(0.05, y0-0.02, coupling1List[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
       
    clp.adjust(savefile=f'./si_figure/s4.png', tight_layout=False)

def get_s5():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    axes = clp.initialize(2, 5, width=12, height=12/5*0.618*1.1*2, LaTeX=True, fontsize=12, sharey=True)
    for i in range(5):
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3/noneqcoord_110k_density_cavity_n400v{i}_E0_2e-4_6e-3.out')
        data2 = np.loadtxt(f'./plotting_data/110k_eqcoord/coord_110k_density_cavity_n400v{i}_E0_2e-4.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        axes[0,i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[0,i].set_xticks([0,5,10,15,20])
        axes[0,i].set_yticks([0,5,10,15])
        clp.plotone(xs, ys, axes[0,i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 0 else False, legendloc=(0.7, 0.35),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 20], ylim=[-0.1,15], legendFontSize=6)
        
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3/noneqcoord_110k_density_cavity_n400v{i}_E0_4e-4_6e-3.out')
        data2 = np.loadtxt(f'./plotting_data/110k_eqcoord/coord_110k_density_cavity_n400v{i}_E0_4e-4.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        axes[1,i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[1,i].set_xticks([0,5,10,15,20])
        axes[1,i].set_yticks([0,5,10,15])
        clp.plotone(xs, ys, axes[1,i], colors=colors_local, labels=labels_local, lw=1, showlegend=False, legendloc=(0.7, 0.35),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 20], ylim=[-0.1,15], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    label2List = ["(f)", "(g)", "(h)", "(i)", "(j)"]
    up2list = [1511.3, 1511.3, 1511.3, 1511.3, 1511.3]
    up1list = [1409.6, 1409.6, 1409.6, 1409.6, 1409.6]
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in up1list]
    coupling2List = [r"excite UP = %d cm$^{-1}$" %up for up in up2list]
    density = [400/(2.9144683314136138**3*2**n) for n in range(5)]
    labels = [r"$\rho = $ %.2f nm$^{-3}$" % density[i] for i in range(5)]
    coupling1List = [coupling1List[i] + "\n" + labels[i] for i in range(5)]
    coupling2List = [coupling2List[i] + "\n" + labels[i] for i in range(5)]
    for i in range(5):
        axes[0,i].text(x0, y0, label1List[i], transform=axes[0,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[0,i].text(0.05, y0-0.02, coupling1List[i], transform=axes[0,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        axes[1,i].text(x0, y0, label2List[i], transform=axes[1,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[1,i].text(0.05, y0-0.02, coupling2List[i], transform=axes[1,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        
    clp.adjust(savefile=f'./si_figure/s5.png', tight_layout=False)

def get_s6():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    e0list = [1,3,4,5,6]
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    for i in range(len(e0list)):
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian/noeqcoord_gaussian_new_6e-3_size_cavity_400_E0_{e0list[i]}e-4.out')
        data2 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian/eqcoord_size_cavity_400_E0_{e0list[i]}e-4.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        axes[i].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[i].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[i].set_xticks([0,5,10,15,20])
        axes[i].set_yticks([0,2,4,6])
        clp.plotone(xs, ys, axes[i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 4 else False, legendloc=(0.7, 0.35),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 20], ylim=[-0.1,6], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1362.9, 1459.6, 1511.3, 1561.4, 1619.8] # omega_c=1311.2
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(5):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(0.05, y0-0.02, coupling1List[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        
    clp.adjust(savefile=f'./si_figure/s6.png', tight_layout=False)

def get_s7():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan', 'grey']
    e0list = [1,3,4,5,6]
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    
    UP_freq = [1362.9, 1409.6, 1459.6, 1511.3, 1561.4, 1619.8, 1668.2, 1719.9, 1769.9, 1818.3] # omega_c=1311.2
    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    kt2au = 3.483e-4
    for i in range(len(e0list)):
        datanoneq = np.loadtxt(f'./plotting_data/photon_110k_6e-3_loss/noneqout_gaussian_new_6e-3_size_cavity_400_E0_{e0list[i]}e-4.out')
        xxanoneq, xybnoneq = datanoneq[:,6]*ang2au, datanoneq[:,10]*ang2au
        vxanoneq, vybnoneq = datanoneq[:,12], datanoneq[:,16]
        energy = (0.5 * (vxanoneq ** 2 + vybnoneq ** 2) + 0.5 * (xxanoneq ** 2 + xybnoneq ** 2) * (UP_freq[i] * cmn2au) ** 2) / (kt2au*400)
        row, col = np.shape(datanoneq)
        te  = np.linspace(0, 20, row)
        x0, y0 = i // 5, i % 5
        axes[i].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[i].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[i].set_xticks([0,1,2,3,4,5])
        axes[i].set_yticks([0,1,2,3])
        clp.plotone([te], [energy], axes[i], colors=[color_list[8]], labels=['photon'], lw=1, showlegend=True if i == 4 else False, legendloc=(0.6, 0.5),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if y0 == 0 else None, xlim=[0, 5], ylim=[-0.1,3], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1362.9, 1459.6, 1511.3, 1561.4, 1619.8] # omega_c=1311.2
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(5):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(0.05, y0-0.02, coupling1List[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
       
    clp.adjust(savefile=f'./si_figure/s7.png', tight_layout=False)

def get_s8():
    ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=12)
    xs, ys = [], []
    e0list = [0,1,2,3,4,5,6,7,8,9,10]
    for e0 in range(len(e0list)):
        data = np.loadtxt(f'./plotting_data/eqloss_n400v0_1500_110k/eqloss_1500_size_cavity_400_E0_{e0list[e0]}e-4.out')
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        #x, y = data[:,2], data[:,3]/1e28
        interval = np.where(x<3300)[0]
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0
        xs.append(x[interval])
        ys.append(y[interval])

    clp.plotone(xs, ys, ax, lw=1.5, 
                xlim=[750, 3300], ylim=[-0.2,11], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    
    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 11)]
    for e0 in range(11):
        if e0 == 0 : ax.text(2000, 0.3, "outside cavity", fontsize=8, color=cadjust_colors[e0])
        elif e0 == 1 : ax.text(2000, e0+0.3, r"$\widetilde{\varepsilon}=5.0\times 10^{-5}$ a.u.", fontsize=8, color=cadjust_colors[e0])
        else: ax.text(2000, e0+0.3, r"$\widetilde{\varepsilon}=%.1f\times 10^{-4}$ a.u." %(e0/2), fontsize=8, color=cadjust_colors[e0])
    
    ax.axvline(x=2902.6, linestyle='-.', alpha=0.2)
    ax.axvline(x=1509.7, linestyle='-.', alpha=0.2)
    ax.axvline(x=3004.3, linestyle='-.', alpha=0.2)
    ax.axvline(x=1312.8, linestyle='-.', alpha=0.2)
    ax.text(2875, 11.3, "$v_1$", fontsize=8, color="0.5")
    ax.text(1480, 11.3, "$v_2$", fontsize=8, color="0.5")
    ax.text(2975, 11.3, "$v_3$", fontsize=8, color="0.5")
    ax.text(1285, 11.3, "$v_4$", fontsize=8, color="0.5")
    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./si_figure/s8.png")

def get_s9():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    e0list = [1,2,3,4,5]
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    for i in range(len(e0list)):
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian_1500/noeqcoord_gaussian_1500_new_6e-3_size_cavity_400_E0_{e0list[i]}e-4.out')
        data2 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian_1500/eqcoord_gaussian_1500_6e-3_size_cavity_400_E0_{e0list[i]}e-4.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        axes[i].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[i].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[i].set_xticks([0,5,10,15,20])
        axes[i].set_yticks([0,1,2,3])
        clp.plotone(xs, ys, axes[i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 0 else False, legendloc=(0.7, 0.35),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 20], ylim=[-0.1,3], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1509.7, 1534.7, 1566.4, 1604.8, 1643.1] # omega_c=1500
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(5):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(0.05, y0-0.02, coupling1List[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        
    clp.adjust(savefile=f'./si_figure/s9.png', tight_layout=False)

def get_s10():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan', 'grey']
    e0list = [1,2,3,4,5]
    axes = clp.initialize(1, 5, width=12, height=12/5*0.618*1.1, LaTeX=True, fontsize=12, sharey=True)
    
    UP_freq = [1509.7, 1534.7, 1566.4, 1604.8, 1643.1] # omega_c=1500, first 5
    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    kt2au = 3.483e-4
    for i in range(len(e0list)):
        datanoneq = np.loadtxt(f'./plotting_data/photon_110k_6e-3_loss_1500/noneqout_gaussian_1500_new_6e-3_size_cavity_400_E0_{e0list[i]}e-4.out')
        xxanoneq, xybnoneq = datanoneq[:,6]*ang2au, datanoneq[:,10]*ang2au
        vxanoneq, vybnoneq = datanoneq[:,12], datanoneq[:,16]
        energy = (0.5 * (vxanoneq ** 2 + vybnoneq ** 2) + 0.5 * (xxanoneq ** 2 + xybnoneq ** 2) * (UP_freq[i] * cmn2au) ** 2) / (kt2au*400)
        row, col = np.shape(datanoneq)
        te  = np.linspace(0, 20, row)
        x0, y0 = i // 5, i % 5
        axes[i].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[i].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[i].set_xticks([0,1,2,3,4,5])
        axes[i].set_yticks([0,4,8,12])
        clp.plotone([te], [energy], axes[i], colors=[color_list[8]], labels=['photon'], lw=1, showlegend=True if i == 4 else False, legendloc=(0.6, 0.5),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if y0 == 0 else None, xlim=[0, 5], ylim=[-0.1,12], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    uplist = [1509.7, 1534.7, 1566.4, 1604.8, 1643.1] # omega_c=1500
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(5):
        axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[i].text(0.05, y0-0.02, coupling1List[i], transform=axes[i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
         
    clp.adjust(savefile=f'./si_figure/s10.png', tight_layout=False)

def get_s11():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    e0list = [1,2,3,0,4,5,6,7,8,9]
    axes = clp.initialize(2, 5, width=12, height=12/5*0.618*1.1*2, LaTeX=True, fontsize=12, sharey=True)
    for i in range(len(e0list)):
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian_fix/noneq_gaussian_coord_freqidx_{e0list[i]}.out')
        data2 = np.loadtxt(f'./plotting_data/110k_eqcoord_fix/cavity_loss_eqcoord_freqidx_{e0list[i]}.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        x, y = i // 5, i % 5
        axes[x, y].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[x, y].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[x, y].set_xticks([0,5,10,15,20])
        axes[x, y].set_yticks([0,1,2,3])
        clp.plotone(xs, ys, axes[x, y], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 9 else False, legendloc=(0.7, 0.35),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if y == 0 else None, xlim=[0, 20], ylim=[-0.1,3], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    label2List = ["(f)", "(g)", "(h)", "(i)", "(j)"]
    omega_c = np.array([1200, 1250, 1300, 1311.2, 1350, 1400, 1450, 1500, 1550, 1600])
    coupling1List = [r"$\omega_c$ = %d cm$^{-1}$"%up for up in omega_c[:5]]
    coupling2List = [r"$\omega_c$ = %d cm$^{-1}$"%up for up in omega_c[5:]]
    labels = [r"fix UP = 1619 cm$^{-1}$"]*5
    coupling1List = [labels[i] + "\n" + coupling1List[i] for i in range(5)]
    coupling2List = [labels[i] + "\n" + coupling2List[i] for i in range(5)]
    for i in range(5):
        axes[0,i].text(x0, y0, label1List[i], transform=axes[0,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[0,i].text(0.05, y0-0.02, coupling1List[i], transform=axes[0,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        axes[1,i].text(x0, y0, label2List[i], transform=axes[1,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[1,i].text(0.05, y0-0.02, coupling2List[i], transform=axes[1,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
       
    clp.adjust(savefile=f'./si_figure/s11.png', tight_layout=False)

def get_s12():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan', 'grey']
    e0list = [1,2,3,0,4,5,6,7,8,9]
    axes = clp.initialize(2, 5, width=12, height=12/5*0.618*1.1*2, LaTeX=True, fontsize=12, sharey=True)
    
    UP_freq = 1619.8
    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    kt2au = 3.483e-4
    for i in range(len(e0list)):
        datanoneq = np.loadtxt(f'./plotting_data/photon_110k_6e-3_fix/noneq_gaussian_out_fix_freqidx_{e0list[i]}.out')
        xxanoneq, xybnoneq = datanoneq[:,6]*ang2au, datanoneq[:,10]*ang2au
        vxanoneq, vybnoneq = datanoneq[:,12], datanoneq[:,16]
        energy = (0.5 * (vxanoneq ** 2 + vybnoneq ** 2) + 0.5 * (xxanoneq ** 2 + xybnoneq ** 2) * (UP_freq * cmn2au) ** 2) / (kt2au*400)
        row, col = np.shape(datanoneq)
        te  = np.linspace(0, 20, row)
        x, y = i // 5, i % 5
        axes[x, y].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[x, y].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[x, y].set_xticks([0,1,2,3,4,5])
        axes[x, y].set_yticks([0,5,10,15])
        clp.plotone([te], [energy], axes[x, y], colors=[color_list[8]], labels=['photon'], lw=1, showlegend=True if i == 0 else False, legendloc=(0.6, 0.5),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if y == 0 else None, xlim=[0, 5], ylim=[-0.1,15], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    label2List = ["(f)", "(g)", "(h)", "(i)", "(j)"]
    omega_c = np.array([1200, 1250, 1300, 1311.2, 1350, 1400, 1450, 1500, 1550, 1600])
    coupling1List = [r"$\omega_c$ = %d cm$^{-1}$"%up for up in omega_c[:5]]
    coupling2List = [r"$\omega_c$ = %d cm$^{-1}$"%up for up in omega_c[5:]]
    labels = [r"fix UP = 1619 cm$^{-1}$"]*5
    coupling1List = [labels[i] + "\n" + coupling1List[i] for i in range(5)]
    coupling2List = [labels[i] + "\n" + coupling2List[i] for i in range(5)]
    for i in range(5):
        axes[0,i].text(x0, y0, label1List[i], transform=axes[0,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[0,i].text(0.05, y0-0.02, coupling1List[i], transform=axes[0,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        axes[1,i].text(x0, y0, label2List[i], transform=axes[1,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[1,i].text(0.05, y0-0.02, coupling2List[i], transform=axes[1,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        
    clp.adjust(savefile=f'./si_figure/s12.png', tight_layout=False)

def get_s13():
    ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=12)
    xs, ys = [], []
    e0list = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    for e0 in range(len(e0list)):
        data = np.loadtxt(f'./plotting_data/eqloss_n400v0_110k_cd4/cd4_eqdac_400_E0_{e0list[e0]}e-4.out')
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        interval = np.where(x<3300)[0]
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0
        xs.append(x[interval])
        ys.append(y[interval])

    clp.plotone(xs, ys, ax, lw=1.5, 
                xlim=[500, 2500], ylim=[-0.2,11], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    
    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 11)]
    for e0 in range(11):
        if e0 == 0 : ax.text(1350, 0.3, "outside cavity", fontsize=8, color=cadjust_colors[e0])
        elif e0 < 4 : ax.text(1350, e0+0.3, r"$\widetilde{\varepsilon}=%.2f\times 10^{-5}$ a.u."%(e0list[e0]*5), fontsize=8, color=cadjust_colors[e0])
        else: ax.text(1350, e0+0.3, r"$\widetilde{\varepsilon}=%.2f\times 10^{-4}$ a.u." %(e0list[e0]/2), fontsize=8, color=cadjust_colors[e0])
    
    ax.axvline(x=2055.2, linestyle='-.', alpha=0.2)
    ax.axvline(x=1067.6, linestyle='-.', alpha=0.2)
    ax.axvline(x=2220.3, linestyle='-.', alpha=0.2)
    ax.axvline(x=989.2, linestyle='-.', alpha=0.2)
    ax.text(2025, 11.3, "$v_1$", fontsize=8, color="0.5")
    ax.text(1037, 11.3, "$v_2$", fontsize=8, color="0.5")
    ax.text(2190, 11.3, "$v_3$", fontsize=8, color="0.5")
    ax.text(959, 11.3, "$v_4$", fontsize=8, color="0.5")

    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./si_figure/s13.png")

def get_s14():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    axes = clp.initialize(2, 5, width=12, height=12/5*0.618*2.2, LaTeX=True, fontsize=12, sharey=True)
    e0list = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    for i in range(len(e0list)):
        data1 = np.loadtxt(f'./plotting_data/110k_coord_cd4/cd4_noneqcoord_6e-3_400_E0_{e0list[i]}e-4.out')
        data2 = np.loadtxt(f'./plotting_data/110k_coord_cd4/cd4_eqcoord_400_E0_{e0list[i]}e-4.out')
        row, col = np.shape(data2)
        te  = np.linspace(0, 20, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        x = i // 5
        y = i % 5
        axes[x,y].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[x,y].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        #axes[i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[x,y].set_xticks([0,5,10,15,20])
        #axes[i].set_yticks([0,4,8,12,16])
        axes[x,y].set_yticks([0,2,4,6])
        #axes[i].set_yticks([0,4,8,12,16,20])
        clp.plotone(xs, ys, axes[x,y], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 9 else False, legendloc=(0.7, 0.25),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if y == 0 else None, xlim=[0, 20], ylim=[-0.1,6], legendFontSize=6)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    label2List = ["(f)", "(g)", "(h)", "(i)", "(j)"]
    uplist = np.array([1.0159e+03, 1.0393e+03, 1.0643e+03, 1.0893e+03, 1.1143e+03, 1.1377e+03, 1.1627e+03, 1.1861e+03, 1.2111e+03, 1.2344e+03])
    coupling1List = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    amp  = ['\n'+r"$\tilde{\varepsilon}=%.2f \times 10 ^{-5}$ a.u."%(e0*5) for e0 in e0list[:3]]
    amp += ['\n'+r"$\tilde{\varepsilon}=%.2f \times 10 ^{-4}$ a.u."%(e0/2) for e0 in e0list[3:]]
    coupling2List = [coupling1List[5+i] + amp[5+i] for i in range(5)]
    coupling1List = [coupling1List[i] + amp[i] for i in range(5)]
    for i in range(5):
        axes[0,i].text(x0, y0, label1List[i], transform=axes[0,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[0,i].text(0.05, y0-0.02, coupling1List[i], transform=axes[0,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        axes[1,i].text(x0, y0, label2List[i], transform=axes[1,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[1,i].text(0.05, y0-0.02, coupling2List[i], transform=axes[1,i].transAxes, fontsize=9, fontweight='bold', va='top', ha='left', color="k")
        
    clp.adjust(savefile=f'./si_figure/s14.png', tight_layout=False)

def plot_IR():
    ax = clp.initialize(1, 1, width=4.3*0.618*0.618*2, height=4.3*0.618*2.2, LaTeX=True, fontsize=10)
    xs, ys = [], []
    e0list = [0,6]
    for e0 in range(len(e0list)):
        data = np.loadtxt(f'./plotting_data/ml_100_ch4/ml_100_ch4_E0_{e0list[e0]}e-4.out')
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        interval = np.where(x<2000)[0]
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0
        xs.append(x[interval])
        ys.append(y[interval])
    clp.plotone(xs, ys, ax, lw=1.5, 
                xlim=[800, 1800], ylim=[-0.2,2], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    
    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 2)]
    ax.text(1330, 0.3, "outside cavity", fontsize=12, color=cadjust_colors[0])
    ax.text(1030, 1.3, r"$\widetilde{\varepsilon}=6.0\times 10^{-4}$ a.u.", fontsize=12, color=cadjust_colors[1])
    ax.set_yticks([0,1,2])
    ax.axvline(x=1254.4, linestyle='-.', alpha=0.2)
    ax.axvline(x=1453.0, linestyle='-.', alpha=0.2)
    #ax.axvline(x=1581.4, linestyle='-.', alpha=0.2)
    ax.text(1380, 1.9, "$v_2$", fontsize=12, color="0.5")
    ax.text(1180, 1.9, "$v_4$", fontsize=12, color="0.5")
    #ax.text(1495, 1.89, "UP", fontsize=10, color="0.5")
    plt.rcParams["axes.axisbelow"] = False
    ax.text(0.12, 0.98, "(a)", transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    clp.adjust(savefile=f"./si_figure/temp1_s15.png")

def get_temp():
    cmap = plt.colormaps["plasma"]
    cmap = cmap.with_extremes(bad=cmap(0))
    labels = [r"$\widetilde{\varepsilon}=%d\times 10^{-4}$ a.u." %n for n in range(11)]
    axes = clp.initialize(2, 2, width=4.3*2, height=4.3*0.618*2.2, LaTeX=True, fontsize=10)

    te = np.linspace(0, 15, 16)
    ye1 = np.zeros((60,16))
    for time in range(16):

        data0 = np.loadtxt(f'./plotting_data/ml_100_ch4/noneqaac_ml_{time}.out')
        
        x1, y1 = data0[:,2], data0[:,3]/1e28
        large_list = np.where(x1>1150)[0]
        small_list = np.where(x1[large_list]<1550)[0]
        xe = x1[large_list][small_list]
        ye1[:,time] = y1[large_list][small_list]

    Ridx = np.where(xe>1350)[0]
    Iidx = np.where(xe<1350)[0]
    yR1 = np.sum(ye1[Ridx,:],axis=0)
    yI1 = np.sum(ye1[Iidx,:],axis=0) 

    xs = [te]*2
    ys = [yI1, yR1]
    colors_local = ["r-o", "k-o"]
    labels = ["1254.4 cm$^{-1}$ $v_4$", "1453.0 cm$^{-1}$ $v_2$"]
    axes[0,1].set_yticks([0,0.5,1,1.5])
    clp.plotone(xs, ys, axes[0,1], colors=colors_local, labels=labels, lw=1, showlegend=True,
                ylabel="intensity [arb. units]", ylim=[0, 1.5], xlim=[0, 15], legendFontSize=8, legendloc=(0.3, 0.1))

    clp.plotone([], [], axes[0,0], ylabel="frequency [cm$^{-1}$]", showlegend=False)
    pcolor = axes[0,0].pcolormesh(te, xe, ye1, norm=colors.LogNorm(vmin=0.01, vmax=0.19), shading='gouraud', cmap=cmap)
    axes[0,0].set_yticks([1150,1250,1350,1450,1550])
    rect1 = mpatches.Rectangle((0.1,1350),14.8,195, fill=False,color="c",linewidth=2, linestyle='--')
    rect2 = mpatches.Rectangle((0.1,1150),14.8,195, fill=False,color="w",linewidth=2, linestyle='--')
    axes[0,0].add_patch(rect1)
    axes[0,0].add_patch(rect2)

    data1 = np.loadtxt(f'./plotting_data/ml_100_ch4/noneqcoord_6e-4.out')[:7501,:]
    data2 = np.loadtxt(f'./plotting_data/ml_100_ch4/eqcoord_6e-4.out')[:7501,:]
    row, col = np.shape(data2)
    te  = np.linspace(0, 15, row)
    ref = np.mean(data2, axis=0)
    y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
    y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
    xs = [te, te]
    ys = [smooth(y4), smooth(y2)]
    colors_local = ['red', 'black']
    labels_local = labels = ["1254.4 cm$^{-1}$ $v_4$", "1453.0 cm$^{-1}$ $v_2$"]
    axes[1,0].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
    axes[1,0].set_yticks([0,4,8,12,16])
    clp.plotone(xs, ys, axes[1,0], colors=colors_local, labels=labels_local, lw=1, showlegend=True, legendloc=(0.55, 0.45),
                xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]", xlim=[0, 15], ylim=[-0.1,16], legendFontSize=8)
    
    data3 = np.loadtxt(f'./plotting_data/ml_100_ch4/noneqcoord_6e-4_100_cf.out')[:7501,:]
    data4 = np.loadtxt(f'./plotting_data/ml_100_ch4/eqcoord_6e-4_100_cf.out')[:7501,:]
    row, col = np.shape(data2)
    te  = np.linspace(0, 15, row)
    ref = np.mean(data4, axis=0)
    y2  = 2 * ((data3[:,1] + data3[:,2]) / (ref[1] + ref[2]) - 1)
    y4  = 3 * ((data3[:,6] + data3[:,7] + data3[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
    xs = [te, te]
    ys = [smooth(y4), smooth(y2)]
    colors_local = ['red', 'black']
    labels_local = labels = ["1312.8 cm$^{-1}$ $v_4$", "1509.7 cm$^{-1}$ $v_2$"]
    axes[1,1].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
    axes[1,1].set_yticks([0,4,8,12,16])
    clp.plotone(xs, ys, axes[1,1], colors=colors_local, labels=labels_local, lw=1, showlegend=True, legendloc=(0.55, 0.45),
                xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]", xlim=[0, 15], ylim=[-0.1,16], legendFontSize=8)

    # add label of the figures
    x0, y0 = 0.98, 0.95
    label1List = ["(b)", "(c)"]
    label2List = ["(d)", "(e)"]
    epsilonList = [r"$\widetilde{\varepsilon} = 6.0\times 10^{-4}$ a.u."]
    couplingList = [r"excite UP = 1581.4 cm$^{-1}$"]
    text1 = r"$N_{\rm{simu}=100}$ , $\widetilde{\varepsilon} = 6.0\times 10^{-4}$ a.u." + "\n" + r"excite UP = 1581.4 cm$^{-1}$" + "\nML"
    text2 = r"$N_{\rm{simu}=100}$ , $\widetilde{\varepsilon} = 6.0\times 10^{-4}$ a.u." + "\n" + r"excite UP = 1619.8 cm$^{-1}$" + "\nCOMPASS"
    for i in range(2):
        axes[1,i].text(x0, y0, label2List[i], transform=axes[1,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    axes[1,0].text(0.07, y0, text1, transform=axes[1,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
    axes[1,1].text(0.07, y0, text2, transform=axes[1,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
    axes[0,0].text(x0, y0, label1List[0], transform=axes[0,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
    axes[0,1].text(x0, y0, label1List[1], transform=axes[0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")    
    axes[0,0].text(0.07, y0, couplingList[0], transform=axes[0,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
    axes[0,0].text(0.07, 0.12, epsilonList[0], transform=axes[0,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
    
    axes[0,0].set_xticks([])
    axes[0,1].set_xticks([])
    axes[1,0].set_xticks([0,5,10,15])
    axes[1,1].set_xticks([0,5,10,15])

    clp.adjust(savefile=f'./si_figure/temp2_s15.png', tight_layout=False)

def get_s15():

    plot_IR()
    get_temp()
    fig1 = image.open('./si_figure/temp1_s15.png')
    fig2 = image.open('./si_figure/temp2_s15.png')

    fix_size = 923
    new_fig = image.new("RGB", (fix_size+fig2.width, max(fig1.height,fig2.height)), (255,255,255))
    new_fig.paste(fig1, (0, max(fig1.height,fig2.height)-fig1.height))
    new_fig.paste(fig2, (fix_size, max(fig1.height,fig2.height)-fig2.height))
    new_fig.save('./si_figure/s15.png')

if __name__ == "__main__":
    get_s1()
    get_s2()
    get_s3()
    get_s4()
    get_s5()
    get_s6()
    get_s7()
    get_s8()
    get_s9()
    get_s10()
    get_s11()
    get_s12()
    get_s13()
    get_s14()
    get_s15()