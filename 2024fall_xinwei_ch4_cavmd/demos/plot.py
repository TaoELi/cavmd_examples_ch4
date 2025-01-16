import numpy as np
import columnplots as clp

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

def plot():
    color_list = ['violet', 'blue', 'green', 'greenyellow', 'gold', 'orange', 'red', 'brown', 'black', 'cyan']
    ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618*1.1, LaTeX=False)

    data1 = np.loadtxt(f'./400/E0_6e-4/noneqcoord/simu_noneq_6e-3_1.xc.xyz.scoord.txt')[:2501,:]
    data2 = np.loadtxt(f'./400/E0_6e-4/coord/simu_1.xc.xyz.scoord.txt')[:2501,:]
    #data1 = np.loadtxt(f'./data/noneqcoord_40/simu_noneq_6e-3_1.xc.xyz.scoord.txt')[:2501,:]
    #data2 = np.loadtxt(f'./data/coord_40/simu_1.xc.xyz.scoord.txt')[:2501,:]
    row, col = np.shape(data2)
    te  = np.linspace(0, 5, row)
    ref = np.mean(data2, axis=0)
    y1  = 1 * (data1[:,0] / ref[0] - 1)
    y2  = 2 * ((data1[:,1] + data1[:,2]) / (ref[1] + ref[2]) - 1)
    y3  = 3 * ((data1[:,3] + data1[:,4] + data1[:,5]) / (ref[3] + ref[4] + ref[5]) - 1)
    y4  = 3 * ((data1[:,6] + data1[:,7] + data1[:,8]) / (ref[6] + ref[7] + ref[8]) - 1)
    xs = [te, te, te, te]
    ys = [smooth(y1), smooth(y2), smooth(y3), smooth(y4)]
    colors_local = [color_list[7], color_list[8], color_list[9], color_list[6]]
    labels_local = [r"v$_%d$" %i for i in range(1,5)]
    ax.axvspan(xmin=0.1, xmax=0.6, ymin=0, ymax=15, color='orange', alpha=0.3)
    ax.set_xticks([0,1,2,3,4,5])
    ax.set_yticks([0,4,8,12,16])
    clp.plotone(xs, ys, ax, colors=colors_local, labels=labels_local, lw=1, showlegend=True, legendloc=(0.7, 0.5),
                xlabel="time [ps]", ylabel=r"energy gain [$\mathrm{k_BT}$]", xlim=[0, 5], ylim=[-0.1,16])
    
    clp.adjust(savefile=f'./figure.png', tight_layout=True)

if __name__ == "__main__":
    plot()
