import PIL.Image as image
import PIL.ImageDraw as draw
import PIL.ImageFont as font
import matplotlib.pyplot as plt
import numpy as np
import columnplots as clp
import matplotlib.colors as colors
from scipy.optimize import curve_fit as fit
import matplotlib.patches as mpatches

def exp_fit(t,E):
    def func(x, a, b):
      return a * np.exp(-b * x)
    
    t_start = 300
    t_cut = t[t_start:]
    E_cut = E[t_start:]
    popt, pcov = fit(func, t_cut, E_cut)
    print(popt)
    
    E_fit = func(x=t_cut,a=popt[0],b=popt[1])
    E_mean = np.mean(E_cut)
    SSres = np.sum((E_fit - E_cut) ** 2)
    SStot = np.sum((E_cut - E_mean) ** 2)
    print(f'R2 = {1-SSres/SStot}')   
    return t_cut, E_fit, popt[1]

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

def plot_IR():
    ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=12)
    xs, ys = [], []
    e0list = [0,1,2,3,4,5,6,7,8,9,10]
    for e0 in range(len(e0list)):
        data = np.loadtxt(f'./plotting_data/eq_n400v0_110k/eq_size_cavity_400_E0_{e0list[e0]}e-4.out')
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
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
        if e0 == 0 : ax.text(2000, 0.3, "outside cavity", fontsize=10, color=cadjust_colors[e0])
        elif e0 == 1 : ax.text(2000, e0+0.3, r"$\widetilde{\varepsilon}=5.0\times 10^{-5}$ a.u.", fontsize=10, color=cadjust_colors[e0])
        else: ax.text(2000, e0+0.3, r"$\widetilde{\varepsilon}=%.1f\times 10^{-4}$ a.u." %(e0/2), fontsize=10, color=cadjust_colors[e0])

    ax.axvline(x=2902.6, linestyle='-.', alpha=0.2)
    ax.axvline(x=1509.7, linestyle='-.', alpha=0.2)
    ax.axvline(x=3004.3, linestyle='-.', alpha=0.2)
    ax.axvline(x=1312.8, linestyle='-.', alpha=0.2)
    ax.text(2875, 11.3, "$v_1$", fontsize=10, color="0.5")
    ax.text(1480, 11.3, "$v_2$", fontsize=10, color="0.5")
    ax.text(2975, 11.3, "$v_3$", fontsize=10, color="0.5")
    ax.text(1285, 11.3, "$v_4$", fontsize=10, color="0.5")
    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./main_paper_figure/110k_eq_n400v0.png")

def get_polariton_freq():

    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    timecut = 20 # in ps

    freq = [1311.2, 1362.9, 1409.6, 1459.6, 1511.3, 1561.4, 1619.8, 1668.2, 1719.9, 1769.9, 1818.3]
    labels = [r"$F=6.32$ mJ/cm$^2$", r"$F=158$ mJ/cm$^2$", r"$F=632$ mJ/cm$^2$"]
    # labels = [r"$F=632$ mJ/cm$^2$ meaning 6e-3", r"$F=6.32$ mJ/cm$^2$ meaning 6e-4"]
    ex, fx, kx = [], [], []
    for size in range(3):
        E0_list, k_list = [], []
        for i in range(1,11):
            print(f'E0={i}e-4')
            timestepcut = timecut * 500 + 1
            UP_freq = freq[i]
            if size == 0 : data = np.loadtxt(f'./plotting_data/photon_110k_6e-4/photon_110k_density_cavity_n400v0_E0_{i}e-4.out')
            elif size == 1 : data = np.loadtxt(f'./plotting_data/photon_110k_3e-3/photon_110k_density_cavity_n400v0_E0_{i}e-4.out')
            elif size == 2 : data = np.loadtxt(f'./plotting_data/photon_110k_6e-3/photon_110k_density_cavity_n400v0_E0_{i}e-4.out')
            xxa, xya, xxb, xyb = data[:timestepcut,6]*ang2au, data[:timestepcut,7]*ang2au, data[:timestepcut,9]*ang2au, data[:timestepcut,10]*ang2au
            vxa, vya, vxb, vyb = data[:timestepcut,12], data[:timestepcut,13], data[:timestepcut,15], data[:timestepcut,16]
            Ek  = 0.5 * (vxa ** 2 + vyb ** 2)
            Ep  = 0.5 * (xxa ** 2 + xyb ** 2) * (UP_freq * cmn2au) ** 2
            E = Ek + Ep
            t = data[:timestepcut,1]
            t_fit, E_fit, k = exp_fit(t=t,E=Ek + Ep)
            k_list.append(k)
            E0_list.append(1311.2 + 51.7 * i)
        ex.append(np.array(E0_list))
        kx.append(np.array(k_list))
        fx.append(freq[1:])

    figure, ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, return_fig_args=True, fontsize=12, LaTeX=True)    
    clp.plotone(fx, kx, ax, labels=labels, lw=1.5, xlabel=r"UP frequency [cm$^{-1}$]", ylabel="UP decay rate [ps$^{-1}$]",
                xlim=[1280,1850], ylim=[0, 3], legendFontSize=9, legendloc=(0.55,0.65), colors=["k--o", "m--o", "r--o"])

    ax.axvline(x=1311.2, linestyle='-.', color='k', alpha=0.5)
    ax.axvline(x=1509.7, linestyle='-.', color='k', alpha=0.5)
    # add text to indicate these two frequencies
    ax.text(1330, 0.5, "$v_4$\nIR active", fontsize=10, color="0.5")
    ax.text(1530, 0.5, "$v_2$\nIR inactive", fontsize=10, color="0.5")
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)
    ax2 = ax.twiny()
    ax2.set_xticks([0,1,2,3,4,5])
    ax2.set_xlim([-0.3,5.3])
    ax2.tick_params(axis='x',direction='in', which='both')
    ax2.set_xlabel(r"$\widetilde{\varepsilon}$ \ [$\times 10^{-4}$ a.u.]") 
    clp.adjust(savefile=f"./main_paper_figure/photonic_life_n400v0.png")

def get_figure1():

    plot_IR()
    get_polariton_freq()
    figa = image.open('./main_paper_figure/4ch4.png')
    figb = image.open('./main_paper_figure/setup_ch4.png')
    figc = image.open('./main_paper_figure/110k_eq_n400v0.png')
    figd = image.open('./main_paper_figure/photonic_life_n400v0.png')
    fige = image.open('./main_paper_figure/ch4_polariton_mechanism.png')
 
    regiona = figa.crop((60, 136, 888, 836))
    regionb = figb.crop((350, 500, 3750, 2100))
    regiona_image = regiona.resize((1846,1600))
    new_figup = image.new("RGB", (regiona_image.width+regionb.width, regionb.height))
    new_figup.paste(regiona_image, (0, 0))
    new_figup.paste(regionb, (regiona_image.width, 0))
    new_figup = new_figup.resize((2388,728))
    #new_figup.show()

    new_figmiddle = image.new("RGB", (figc.width+figd.width, figd.height), (255,255,255))
    new_figmiddle.paste(figc, (0, 914-839))
    new_figmiddle.paste(figd, (figc.width, 0))
    #new_figdown.show()

    new_figdown = fige.resize((2388,1090))

    new_fig = image.new("RGB", (new_figup.width, new_figup.height+new_figmiddle.height+new_figdown.height), (255,255,255))
    new_fig.paste(new_figup, (0, 0))
    new_fig.paste(new_figmiddle, (0, new_figup.height))
    new_fig.paste(new_figdown, (0, new_figup.height+new_figmiddle.height))

    draw_new_fig = draw.Draw(new_fig)
    font_new_fig = font.truetype('Arial.ttf', size=45)
    draw_new_fig.text((10, 0), "(a)", font=font_new_fig, fill='black')
    draw_new_fig.text((10+regiona_image.width/5246*2388, 0), "(b)", font=font_new_fig, fill='white')
    draw_new_fig.text((10, new_figup.height+20), "(c)", font=font_new_fig, fill='black')
    draw_new_fig.text((10+figc.width, new_figup.height+20), "(d)", font=font_new_fig, fill='black')
    draw_new_fig.text((30, new_figup.height+new_figmiddle.height+20), "(e)", font=font_new_fig, fill='black')
    draw_new_fig.text((20+1250, new_figup.height+new_figmiddle.height+20), "(f)", font=font_new_fig, fill='black')
    new_fig.save('./main_paper_figure/figure1.png')

def get_figure2():
    cmap = plt.colormaps["plasma"]
    cmap = cmap.with_extremes(bad=cmap(0))
    labels = [r"$\widetilde{\varepsilon}=%d\times 10^{-4}$ a.u." %n for n in range(11)]
    e0list = [1,3,4,5,6]
    axes = clp.initialize(3, 5, width=12, height=12/5*0.618*3.3, LaTeX=True, fontsize=10)
    for i in range(len(e0list)):

        te = np.linspace(0, 15, 16)
        ye1 = np.zeros((60,16))
        ye2 = np.zeros((60,16))
        for time in range(16):

            data1 = np.loadtxt(f'./plotting_data/110k_aacvst/noneqaac_110k_density_cavity_n400v0_E0_{e0list[i]}e-4_{time}.out')
            data2 = np.loadtxt(f'./plotting_data/110k_aacvst/eqaac_110k_density_cavity_n400v0_E0_{e0list[i]}e-4_{time}.out')

            x1, y1 = data1[:,2], data1[:,3]/1e28
            x2, y2 = data2[:,2], data2[:,3]/1e28
            large_list = np.where(x1>1200)[0]
            small_list = np.where(x1[large_list]<1600)[0]
            xe = x1[large_list][small_list]
            ye1[:,time] = y1[large_list][small_list]
            ye2[:,time] = y2[large_list][small_list]

        Ridx = np.where(xe>1400)[0]
        Iidx = np.where(xe<1400)[0]
        yR1 = np.sum(ye1[Ridx,:],axis=0)
        yI1 = np.sum(ye1[Iidx,:],axis=0) 

        xs = [te]*2
        ys = [yI1, yR1]
        colors_local = ["r-o", "k-o"]
        labels = ["1313 cm$^{-1}$ $v_4$", "1510 cm$^{-1}$ $v_2$"]
        axes[1,i].set_yticks([0,0.5,1,1.5])
        clp.plotone(xs, ys, axes[1,i], colors=colors_local, labels=labels, lw=1, showlegend=True if i == 4 else False,
                    xlabel="time [ps]", ylabel="intensity [arb. units]" if i == 0 else None, ylim=[0, 1.5], xlim=[0, 15], legendFontSize=8, legendloc=(0.3, 0.02))
    
        clp.plotone([], [], axes[0,i], ylabel="frequency [cm$^{-1}$]" if i == 0 else None, showlegend=False)
        pcolor = axes[0,i].pcolormesh(te, xe, ye1, norm=colors.LogNorm(vmin=0.01, vmax=0.17), shading='gouraud', cmap=cmap)
        axes[0,i].set_yticks([1200,1300,1400,1500,1600])
        rect1 = mpatches.Rectangle((0.1,1400),14.8,195, fill=False,color="c",linewidth=2, linestyle='--')
        rect2 = mpatches.Rectangle((0.1,1200),14.8,195, fill=False,color="w",linewidth=2, linestyle='--')
        axes[0,i].add_patch(rect1)
        axes[0,i].add_patch(rect2)

    for i in range(len(e0list)):
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3/noneqcoord_110k_density_cavity_n400v0_E0_{e0list[i]}e-4_6e-3.out')[:7501,:]
        data2 = np.loadtxt(f'./plotting_data/110k_eqcoord/coord_110k_density_cavity_n400v0_E0_{e0list[i]}e-4.out')[:7501,:]
        row, col = np.shape(data2)
        te  = np.linspace(0, 15, row)
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        xs = [te, te, te, te]
        ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
        colors_local = ['brown', 'black', 'cyan', 'red']
        labels_local = [r"v$_%d$" %i for i in range(1,5)]
        axes[2,i].axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[2,i].set_yticks([0,4,8,12,16])
        clp.plotone(xs, ys, axes[2,i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 4 else False, legendloc=(0.5, 0.4),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 15], ylim=[-0.1,16], legendFontSize=7)

    # add label of the figures
    x0, y0 = 0.98, 0.97
    label1List = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    label2List = ["(f)", "(g)", "(h)", "(i)", "(j)"]
    label3List = ["(k)", "(l)", "(m)", "(n)", "(o)"]
    epsilonList = [r"$\widetilde{\varepsilon} = 5.0\times 10^{-5}$ a.u."]
    epsilonList += [r"$\widetilde{\varepsilon} = %.1f\times 10^{-4}$ a.u." %(e0/2) for e0 in e0list[1:]]
    uplist = [1362.9, 1459.6, 1511.3, 1561.4, 1619.8]
    couplingList = [r"excite UP = %d cm$^{-1}$" %up for up in uplist]
    for i in range(len(e0list)):
        axes[0, i].text(x0, y0, label1List[i], transform=axes[0, i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        axes[1, i].text(x0, y0, label2List[i], transform=axes[1, i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[2, i].text(x0, y0, label3List[i], transform=axes[2, i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[0, i].text(x0-0.15, y0-0.02, couplingList[i], transform=axes[0, i].transAxes, fontsize=9, fontweight='bold', va='top', ha='right', color="w")
        axes[0, i].text(x0-0.12, 0.13, epsilonList[i], transform=axes[0, i].transAxes, fontsize=10, fontweight='bold', va='top', ha='right', color="w")
    
    # remove numbers on the y axis 
    for i in range(len(e0list)):
        if i > 0:
            axes[0, i].set_yticklabels([])
            axes[1, i].set_yticklabels([])
            axes[2, i].set_yticklabels([])
    clp.adjust(savefile=f'./main_paper_figure/figure2.png', tight_layout=False)

def size_denpendence():

    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    timecut = 20 # in ps

    freq = [1311.2, 1362.9, 1409.6, 1459.6, 1511.3, 1561.4, 1619.8, 1668.2, 1719.9, 1769.9, 1818.3]
    mol = [2**n*100 for n in range(5)]
    labels = [r"$N_{\mathrm{simu}}$ = $ %d$" % mol[i] for i in range(5)]

    ex, fx, kx = [], [], []
    for size in range(5):
        E0_list, k_list = [], []
        for i in range(1,11):
            print(f'E0={i}e-4')
            timestepcut = timecut * 500 + 1
            UP_freq = freq[i]
            data = np.loadtxt(f'./plotting_data/photon_size_6e-3/photon_size_cavity_{mol[size]}_E0_{i}e-4.out')
            xxa, xya, xxb, xyb = data[:timestepcut,6]*ang2au, data[:timestepcut,7]*ang2au, data[:timestepcut,9]*ang2au, data[:timestepcut,10]*ang2au
            vxa, vya, vxb, vyb = data[:timestepcut,12], data[:timestepcut,13], data[:timestepcut,15], data[:timestepcut,16]
            Ek  = 0.5 * (vxa ** 2 + vyb ** 2)
            Ep  = 0.5 * (xxa ** 2 + xyb ** 2) * (UP_freq * cmn2au) ** 2
            E = Ek + Ep
            t = data[:timestepcut,1]
            t_fit, E_fit, k = exp_fit(t=t,E=Ek + Ep)
            k_list.append(k)
            E0_list.append(1311.2 + 51.7 * i)
        ex.append(np.array(E0_list))
        kx.append(np.array(k_list))
        fx.append(freq[1:])

    figure, ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, return_fig_args=True, fontsize=12, LaTeX=True)    
    clp.plotone(fx, kx, ax, labels=labels, lw=1.5, xlabel=r"UP frequency [cm$^{-1}$]", ylabel="UP decay rate [ps$^{-1}$]",
                xlim=[1280,1850], ylim=[0, 3], legendFontSize=8, colorMap=plt.cm.twilight, colorMap_startpoint=0.2, colorMap_endpoint=0.8, colors=["-"]*5)
                
    ax.axvline(x=1311.2, linestyle='-.', color='k', alpha=0.5)
    ax.axvline(x=1509.7, linestyle='-.', color='k', alpha=0.5)
    # add text to indicate these two frequencies
    ax.text(1330, 0.5, "$v_4$\nIR active", fontsize=10, color="0.5")
    ax.text(1530, 0.5, "$v_2$\nIR inactive", fontsize=10, color="0.5")
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)
    ax2 = ax.twiny()
    ax2.set_xticks([0,1,2,3,4,5])
    ax2.set_xlim([-0.3,5.3])
    ax2.tick_params(axis='x',direction='in', which='both')
    ax2.set_xlabel(r"$\widetilde{\varepsilon}$ \ [$\sqrt{400/N_{\rm{simu}}}\times 10^{-4}$ a.u.]") 
    clp.adjust(savefile=f"./main_paper_figure/photonic_life_size_6e-3.png")

def density_denpendence():

    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    timecut = 20 # in ps

    freq = [1311.2, 1362.9, 1409.6, 1459.6, 1511.3, 1561.4, 1619.8, 1668.2, 1719.9, 1769.9, 1818.3]
    density = [400/(2.9144683314136138**3*2**n) for n in range(5)]
    labels = [r"$\rho = $ %.2f nm$^{-3}$" % density[i] for i in range(5)]

    ex, fx, kx = [], [], []
    for size in range(5):
        E0_list, k_list = [], []
        for i in range(1,11):
            print(f'E0={i}e-4')
            timestepcut = timecut * 500 + 1
            UP_freq = freq[i]
            data = np.loadtxt(f'./plotting_data/photon_110k_6e-3/photon_110k_density_cavity_n400v{size}_E0_{i}e-4.out')
            xxa, xya, xxb, xyb = data[:timestepcut,6]*ang2au, data[:timestepcut,7]*ang2au, data[:timestepcut,9]*ang2au, data[:timestepcut,10]*ang2au
            vxa, vya, vxb, vyb = data[:timestepcut,12], data[:timestepcut,13], data[:timestepcut,15], data[:timestepcut,16]
            Ek  = 0.5 * (vxa ** 2 + vyb ** 2)
            Ep  = 0.5 * (xxa ** 2 + xyb ** 2) * (UP_freq * cmn2au) ** 2
            E = Ek + Ep
            t = data[:timestepcut,1]
            t_fit, E_fit, k = exp_fit(t=t,E=Ek + Ep)
            k_list.append(k)
            E0_list.append(1311.2 + 51.7 * i)
        ex.append(np.array(E0_list))
        kx.append(np.array(k_list))
        fx.append(freq[1:])

    figure, ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, return_fig_args=True, fontsize=12, LaTeX=True)    
    clp.plotone(fx, kx, ax, labels=labels, lw=1.5, xlabel=r"UP frequency [cm$^{-1}$]", ylabel="UP decay rate [ps$^{-1}$]",
                xlim=[1280,1850], ylim=[0, 3], legendFontSize=8, colorMap=plt.cm.viridis, colorMap_startpoint=0, colorMap_endpoint=0.9, colors=["-"]*5)

    ax.axvline(x=1311.2, linestyle='-.', color='k', alpha=0.5)
    ax.axvline(x=1509.7, linestyle='-.', color='k', alpha=0.5)
    # add text to indicate these two frequencies
    ax.text(1330, 0.5, "$v_4$\nIR active", fontsize=10, color="0.5")
    ax.text(1530, 0.5, "$v_2$\nIR inactive", fontsize=10, color="0.5")
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)
    ax2 = ax.twiny()
    ax2.set_xticks([0,1,2,3,4,5])
    ax2.set_xlim([-0.3,5.3])
    ax2.tick_params(axis='x',direction='in', which='both')
    ax2.set_xlabel(r"$\widetilde{\varepsilon}$ \ [$\times 10^{-4}$ a.u.]") 
    clp.adjust(savefile=f"./main_paper_figure/photonic_life_density_6e-3.png")

def get_figure3():

    size_denpendence()
    density_denpendence()
    fig1 = image.open('./main_paper_figure/photonic_life_size_6e-3.png')
    fig2 = image.open('./main_paper_figure/photonic_life_density_6e-3.png')
    new_fig = image.new("RGB", (fig1.width+fig2.width, max(fig1.height,fig2.height)), (255,255,255))
    new_fig.paste(fig1, (0, max(fig1.height,fig2.height)-fig1.height))
    new_fig.paste(fig2, (fig1.width, max(fig1.height,fig2.height)-fig2.height))

    draw_new_fig = draw.Draw(new_fig)
    font_new_fig = font.truetype('Arial.ttf', size=50)
    draw_new_fig.text((10, 0), "(a)", font=font_new_fig, fill='black')
    draw_new_fig.text((10+fig1.width, 0), "(b)", font=font_new_fig, fill='black')
    new_fig.save('./main_paper_figure/figure3.png')

def get_figure4():
    color_list   = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    axes = clp.initialize(2, 3, width=4.3*3, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12)
    for i in range(3):
        data1 = np.loadtxt(f'./plotting_data/1500_cavity_loss_gaussian/noneqcoord_gaussian_{i}.out')
        data2 = np.loadtxt(f'./plotting_data/1500_cavity_loss_gaussian/eqcoord_gaussian_{i}.out')
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
        axes[0,i].axvspan(xmin=1.75, xmax=2.25, ymin=0, ymax=15, color='orange', alpha=0.4)
        axes[0,i].axvspan(xmin=1.25, xmax=2.75, ymin=0, ymax=15, color='orange', alpha=0.3)
        axes[0,i].set_xticks([0,5,10,15,20])
        axes[0,i].set_yticks([0,2,4,6])
        clp.plotone(xs, ys, axes[0,i], colors=colors_local, labels=labels_local, lw=1, showlegend=True if i == 2 else False, legendloc=(0.7, 0.3),
                    xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]" if i == 0 else None, xlim=[0, 20], ylim=[-0.1,6], legendFontSize=10)

    # add label of the figures
    label1List = ["(a)", "(b)", "(c)"]
    couplingList = [r"excite UP = 1511 cm$^{-1}$", r"excite UP = 1510 cm$^{-1}$", r"excite $v_2$ = 1510 cm$^{-1}$"]
    labels = [r"$\omega_{\mathrm{c}}$ = 1311 cm$^{-1}$", r"$\omega_{\mathrm{c}}$ = 1500 cm$^{-1}$", ' ']
    varepsilon = [r"$\widetilde{\varepsilon}$ = $2.0\times 10^{-4}$ a.u.", r"$\widetilde{\varepsilon}$ = $5.0\times 10^{-5}$ a.u.", 'outside cavity']
    couplingList = [couplingList[i] + "\n" + varepsilon[i] + "\n" + labels[i] for i in range(3)]
    for i in range(3):
        axes[0,i].text(0.01, 0.98, label1List[i], transform=axes[0,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[0,i].text(0.23, 0.96, couplingList[i], transform=axes[0,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
    
    xs, ys = [], []
    for e0 in range(9):
        data = np.loadtxt(f'./plotting_data/eqloss_n400v0_110k_fix/cavity_loss_eqdac_fix_freqidx_{e0+1}.out')
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        #x, y = data[:,2], data[:,3]/1e28
        interval = np.where(x<3300)[0]
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0
        xs.append(x[interval])
        ys.append(y[interval])

    clp.plotone(xs, ys, axes[1,0], lw=1.5, 
                xlim=[750, 3300], ylim=[-0.2,9], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    
    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 11)]

    epsilon = 0.5 * np.array([6.3014, 6.2146, 6.0648, 5.8287, 5.4955, 5.0427, 4.3965, 3.4637, 1.8238])
    for e0 in range(9):
        axes[1,0].annotate('', xy=(1200+50*e0, e0), xytext=(1200+50*e0, e0+0.8), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
        if e0 != 8 : axes[1,0].text(1900, e0+0.3, r"$\widetilde{\varepsilon}=%.2f\times 10^{-4}$ a.u." %epsilon[e0], fontsize=10, color=cadjust_colors[e0])
        else : axes[1,0].text(1900, e0+0.3, r"$\widetilde{\varepsilon}=9.12\times 10^{-5}$ a.u." %epsilon[e0], fontsize=10, color=cadjust_colors[e0])
    axes[1,0].text(1250, 0.3, r"$\omega_{\mathrm{c}}$", fontsize=10, color='blue', alpha=0.8)
    axes[1,0].set_xticks([1000,1500,2000,2500,3000])
    axes[1,0].axvline(x=1619.8, linestyle='-.', alpha=0.2)
    axes[1,0].text(1650, 8.2, "UP", fontsize=10, color="0.5")
    axes[1,0].text(0.01, 0.98, "(d)", transform=axes[1,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
    plt.rcParams["axes.axisbelow"] = False

    UP_freq = [1619.77763771]*9
    ang2au = 1.8897259886
    cmn2au = 7.251632778606085e-7 * 2 * 3.1415926
    kt2au = 3.483e-4
    energy = np.zeros((9,4))
    omega_c = np.array([1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600])
    energy[:,0] = omega_c
    for i in range(len(omega_c)):
        datanoneq = np.loadtxt(f'./plotting_data/photon_110k_6e-3_fix/noneq_gaussian_out_fix_freqidx_{i+1}.out')
        xxanoneq, xybnoneq = datanoneq[:,6]*ang2au, datanoneq[:,10]*ang2au
        vxanoneq, vybnoneq = datanoneq[:,12], datanoneq[:,16]
        energy[i,1] = max(0.5 * (vxanoneq ** 2 + vybnoneq ** 2) + 0.5 * (xxanoneq ** 2 + xybnoneq ** 2) * (UP_freq[i] * cmn2au) ** 2) / (kt2au*400)
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian_fix/noneq_gaussian_coord_freqidx_{i+1}.out')
        data2 = np.loadtxt(f'./plotting_data/110k_eqcoord_fix/cavity_loss_eqcoord_freqidx_{i+1}.out')
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        energy[i,2] = max(y2)
        energy[i,3] = max(y3+y2)

    xs = [energy[:,0]]*1
    ys = [energy[:,1]]
    colors=["orange"]
    linestyles=["-"]*1
    markers=['o']*1
    labels=[r"photon / $N_{\rm{simu}}$"]
    clp.plotone(xs, ys, axes[1,1], colors=colors, linestyles=linestyles, markers=markers, labels=labels, lw=1, showlegend=True, legendloc=(0.26,0.84),
                xlabel=r"$\omega_{\mathrm{c}}$ [cm$^{-1}$]", ylabel=r"energy gain [k$_{\mathrm{B}}$T]", xlim=[1190,1610], ylim=[0,10], legendFontSize=10)
    
    axes[1,1].set_xticks([1200,1300,1400,1500,1600])
    axes[1,1].set_yticks([0,2,4,6,8,10])
    axes[1,1].text(0.01, 0.98, "(e)", transform=axes[1,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
    
    energy = np.zeros((9,4))
    energy[:,0] = omega_c
    for i in range(len(omega_c)):
        datanoneq = np.loadtxt(f'./plotting_data/photon_110k_6e-3_fix/noneq_gaussian_out_fix_freqidx_{i+1}.out')
        xxanoneq, xybnoneq = datanoneq[:,6]*ang2au, datanoneq[:,10]*ang2au
        vxanoneq, vybnoneq = datanoneq[:,12], datanoneq[:,16]
        energy[i,1] = max(0.5 * (vxanoneq ** 2 + vybnoneq ** 2) + 0.5 * (xxanoneq ** 2 + xybnoneq ** 2) * (UP_freq[i] * cmn2au) ** 2) / (kt2au*400)
        data1 = np.loadtxt(f'./plotting_data/110k_noneqcoord_6e-3_gaussian_fix/noneq_gaussian_coord_freqidx_{i+1}.out')
        data2 = np.loadtxt(f'./plotting_data/110k_eqcoord_fix/cavity_loss_eqcoord_freqidx_{i+1}.out')
        ref = np.mean(data2, axis=0)
        y1  = 1 * (data1[:,0] / ref[0] - 1)
        y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
        y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
        y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
        energy[i,2] = max(y2)
        energy[i,3] = max(y3+y2)

    xs = [energy[:,0]]*2
    ys = [energy[:,2],energy[:,3]]
    colors=["k","m"]
    linestyles=["-"]*2
    markers=['o']*2
    labels=[r"v$_2$", r"v$_2$+v$_3$"]
    clp.plotone(xs, ys, axes[1,2], colors=colors, linestyles=linestyles, markers=markers, labels=labels, lw=1, showlegend=True,
                xlabel=r"$\omega_{\mathrm{c}}$ [cm$^{-1}$]", ylabel=r"energy gain [k$_{\mathrm{B}}$T]", xlim=[1190,1610], ylim=[0,5], legendFontSize=10)
    
    axes[1,2].axvline(x=1465.5, linestyle='-.', color="0.5", alpha=0.5)
    axes[1,2].text(1290, 0.7, r"$\omega_{\mathrm{c}}=1466$ cm$^{-1}$"+'\n'+r"max $|X_+^{(\rm{c})}|^4|X_+^{(\rm{B})}|^2$", fontsize=10, color="0.5")
    axes[1,2].set_xticks([1200,1300,1400,1500,1600])
    axes[1,2].set_yticks([0,1,2,3,4,5])
    axes[1,2].text(0.01, 0.98, "(f)", transform=axes[1,2].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
    
    clp.subplots_adjust(hspace=0.3)
    clp.adjust(savefile=f'./main_paper_figure/figure4.png', tight_layout=False)

if __name__ == "__main__" :
    get_figure1()
    get_figure2()
    get_figure3()
    get_figure4()