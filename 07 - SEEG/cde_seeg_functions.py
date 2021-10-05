import matplotlib.pyplot as plt
import scipy
import numpy as np
import copy
import cde_seeg_functions as cs
import cde_seeg_importancesampling as cIS

#=======================================================================
def subject_specify(F):
#=======================================================================    
    import re, os 
    
    # Specify SEEG folders and file locations
    #-------------------------------------------------------------------------------
    names = os.listdir(F['seeg'])
    r     = re.compile('^[A-Z].*[a-z]$')
    folds = list(filter(r.match, names))
    folds.sort()

    Sub = []
    for f in folds:
        edfs = os.listdir(F['seeg'] +os.sep+ f)
        r    = re.compile('^[A-Z].*Baseline.*[EDF|edf]') 
        bl   = list(filter(r.match, edfs))
        for k in range(len(bl)): bl[k] =  F['seeg'] +os.sep+ f +os.sep+ bl[k] 
        bl.sort()

        r    = re.compile('^[A-Z].*Seizure.*[EDF|edf]') 
        sz = list(filter(r.match, edfs))
        for k in range(len(sz)): sz[k] =  F['seeg'] +os.sep+ f +os.sep+ sz[k] 
        sz.sort()

        Sub.append({'Base': bl, 'Seiz': sz})
    
    return Sub

#=======================================================================
def subject_load(sub, PP):
#=======================================================================
    import mne 
    from IPython.display import clear_output
    
    if len(sub['Base']) != len(sub['Seiz']): raise ValueError('Mismatch in the number of Baseline and Seizure segments')
        
    base = []
    seiz = []
    
    for k in range(len(sub['Base'])):

        # Load data, filter and rereference
        #---------------------------------------------------------------------------
        bl = mne.io.read_raw_edf(sub['Base'][k], preload=True)
        bl.filter(PP['Fbp'][0], PP['Fbp'][1])
        bl.set_eeg_reference(ref_channels='average')
        base.append(bl._data)
        
        sz = mne.io.read_raw_edf(sub['Seiz'][k], preload=True)
        sz.filter(PP['Fbp'][0], PP['Fbp'][1])
        sz.set_eeg_reference(ref_channels='average')
        seiz.append(sz._data)

    sub.update({'Base_dat':base, 'Seiz_dat':seiz})
    clear_output()
    return sub


#=======================================================================
def binarise(data, BN, chanmean = [], chanstdv = []):
#=======================================================================
    import scipy 
    import numpy as np
    
    if chanmean == []: chanmean = data.mean(1)
    if chanstdv == []: chanstdv = data.std(1)
    zdat = np.divide(np.subtract(data, chanmean[:,np.newaxis]), chanstdv[:,np.newaxis])
    tdat = np.zeros(zdat.shape)

    for chan in range(zdat.shape[0]):
        pdat = scipy.signal.find_peaks(abs(zdat[chan,:]).astype('float64'), 
                                      height=BN['peak_height'], width=BN['separation_win'])
        for p in pdat[0]: tdat[chan,p] = 1

    drange = (BN['edgewin'], tdat.shape[1]-BN['edgewin'])
    return tdat[:,drange[0]:drange[1]]
        

#=======================================================================        
def subject_binarise(sub, BN, PP = None):    # subject structure; binarisation parameters; preproc params if required
#=======================================================================
    import cde_seeg_functions as cs
    import numpy as np
    import os, re, sys
    
    if len(sub['Base']) != len(sub['Seiz']): raise ValueError('Mismatch in the number of Baseline and Seizure segments')

    # Find out how many channels there are for this subject 
    #-------------------------------------------------------------------------------
    if 'Base_dat' not in sub: sub = cs.subject_load(sub, PP)
    
    for k in range(len(sub['Base_dat'])):
        
        bl = sub['Base_dat'][k]
        sz = sub['Seiz_dat'][k]
        
        if k == 0: 
            Nch    = bl.shape[0]
            Bl_bin = np.ndarray((Nch,0)) 
            bl_id  = np.ndarray(0)
            Sz_bin = np.ndarray((Nch,0))
            sz_id  = np.ndarray(0)
            

        # Segment out the ranges above and keep track of segment boundaries
        #---------------------------------------------------------------------------
        bin_tmp = cs.binarise(bl, BN)
        Bl_bin  = np.append(Bl_bin, bin_tmp, axis=1)
        bl_id   = np.append(bl_id, np.ones((bin_tmp.shape[1]))*k, axis=0)
            
        bin_tmp = cs.binarise(sz, BN, bl.mean(1), bl.std(1))
        Sz_bin = np.append(Sz_bin, bin_tmp, axis=1)
        sz_id  = np.append(sz_id, np.ones((bin_tmp.shape[1]))*k, axis=0)

    sub.update({'bin_base':Bl_bin, 'bl_seg_id':bl_id, 'bin_seiz':Sz_bin, 'sz_seg_id':sz_id})
    return sub

#=======================================================================
def avcount(bintrace, AS):
#=======================================================================
    import numpy as np
    
    # Calculate sums of events within each firing window 
    #-------------------------------------------------------------------------------
    tsums = np.sum(bintrace,0)  
    sums  = np.ndarray(0)
    for s in range(0,len(tsums),AS['dt']):
        sums = np.append(sums, np.sum(tsums[s: (s+AS['dt']-1)]))

    # Concatenate ongoing runs of events
    #-------------------------------------------------------------------------------
    thisid   = 0
    lastseen = -1
    allids   = np.zeros(sums.shape)

    for s in range(len(sums)):
        if sums[s] > 0:
            allids[s] = thisid
            lastseen  = s
        else: 
            if lastseen == s-1: thisid = thisid+1
            allids[s] = 0
    #     print(allids[s])

    # Count number of total events within each cascade
    #-------------------------------------------------------------------------------
    avcounts = np.ndarray(0)
    avtimes  = np.ndarray(0)
    avids    = np.unique(allids)
    avids    = avids[np.where(avids != 0)]
    for a in avids:
        avcounts = np.append(avcounts, np.sum(sums[np.where(allids==a)[0]]))
        avtimes  = np.append(avtimes, len(np.where(allids == a)[0]))
    
    return avcounts, avtimes 

#=======================================================================
def avcalc(Subs, Bands, BN, AS, prange, trange):
#=======================================================================

    # Loop through individual frequency bands
    #-----------------------------------------------------------------------------------
    for band in Bands:          
        Pp    = {}
        Pp['Bands'] = [item for item in Bands if band[0] in item]
        Pp['Fbp']   = (Pp['Bands'][0][1], Pp['Bands'][0][2])

        # Load and preprocess EEG segments
        #===============================================================================
        for s in range(len(Subs)): Subs[s] = cs.subject_load(Subs[s],Pp)

        # Perform parameter sweep for each subject
        #===============================================================================
        for s in range(len(Subs)): 

            # Define parameter ranges
            #-------------------------------------------------------------------------------
            Bn     = copy.deepcopy(BN)
            As     = copy.deepcopy(AS)
            if 'Base_dat' not in Subs[s]: Subs[s] = cs.subject_load(Subs[s], Pp)

            # Perform actual parameter sweep
            #-------------------------------------------------------------------------------
            print('Frequency band ' +str(Pp['Fbp'][0])+' to '+str(Pp['Fbp'][1])+'Hz (' +band[0]+ 'band)')
            print('Working on subject ' +str(s+1)+ ' of ' +str(len(Subs)))
            print('Processing: ', end='')
            for pi in range(len(prange)):
                print('>', end = '')
                for ti in range(len(trange)):
                    Bn['peak_height'] = prange[pi]
                    As['dt']          = trange[ti]
                    sub               = cs.subject_binarise(Subs[s], Bn)

                    bl = sub['bin_base']
                    sz = sub['bin_seiz']
                    blcnt, bltim = cs.avcount(bl, As)
                    szcnt, sztim = cs.avcount(sz, As)
                    if ti == 0: 
                        ava = [(blcnt,szcnt)]
                        avt = [(bltim,sztim)]
                    else:       
                        ava.append((blcnt,szcnt))
                        avt.append((bltim,sztim))
                if pi == 0: 
                    Avas = [ava]
                    Avts = [avt]
                else:       
                    Avas.append(ava)
                    Avts.append(avt)

            sweep = {}
            sweep.update({'Avalanches': Avas,
                          'Peak thresholds': prange,
                          'Time windows': trange,
                          'AvTimes':Avts})
            Subs[s].update({'Sweep':sweep})
            print(' Done')
        print('All done :)')

    return Subs

#=======================================================================
def plot_ava(sub): 
#=======================================================================
    # Plot example avalanche
    #--------------------------------------------------------------------------------------------t
    win   = range(5000,6000)
    offs  = -125
    win2  = range(win[0]+offs, win[-1]+offs)
    chans = 50 
    f, a  = plt.subplots(2,2, figsize=(30,20), sharex = True)

    for k in range(chans):
        zbase = (sub['Base_dat'][0][k] - np.mean(sub['Base_dat'][0][k])) / np.std(sub['Base_dat'][0][k])
        a[0,0].plot(np.linspace(0,0.5,len(win)),zbase[win] - k*4, color='black', linewidth=0.5)
        a[0,0].set_title('Baseline High Gamma Data', fontsize=35)
        a[0,0].set_yticks( () )
        zseiz = (sub['Seiz_dat'][0][k] - np.mean(sub['Base_dat'][0][k])) / np.std(sub['Base_dat'][0][k])
        a[0,1].plot(np.linspace(0,0.5,len(win)),zseiz[win] - k*4, color='black', linewidth=0.5)
        a[0,1].set_title('Seizure High Gamma Data', fontsize=35)
        a[0,1].set_yticks( () )

    a[1,0].imshow([a[win2] for a in sub['bin_base'][:chans]], aspect=0.005,cmap='Greys', extent=(0,0.5,1,k+1))
    a[1,0].set_xlabel('Time in seconds', fontsize=25)
    a[1,0].set_title('Binarised Data', fontsize=25)
    a[1,0].set_ylabel('Channel number', fontsize=25)
    a[1,0].tick_params(axis='both', which='major', labelsize=20)
    a[1,1].imshow([a[win2] for a in sub['bin_seiz'][:chans]], aspect=0.005,cmap='Greys', extent=(0,0.5,1,k+1))
    a[1,1].set_xlabel('Time in seconds', fontsize=25)
    a[1,1].set_ylabel('Channel number', fontsize=25)
    a[1,1].set_title('Binarised Data', fontsize=25)
    a[1,1].tick_params(axis='both', which='major', labelsize=20)
    plt.tight_layout()

#=======================================================================
def plot_ccdf(avc, No_bins=50, ax=None, color=1, plottype = 'scatter', linewidth = 1, linestyle = '-'):
#=======================================================================
    import numpy as np
    import matplotlib.pyplot as plt 
    import matplotlib
    
    if isinstance(color, str):
        tc = color
    else: 
        cmap       = matplotlib.cm.get_cmap('Set1')
        tc         = cmap(color)
    
    if not ax: f, ax = plt.subplots()
        
    # Calcualte logarithmically binned c-cdf 
    #-------------------------------------------------------------------------------
    if plottype == 'distribution':
    # for the distribution plots, we want to track individual fish
        avc_l = avc
        max_length = np.max([np.max(a) for a in avc_l])
        min_length = np.min([np.min(a) for a in avc_l])

        # Convert to matrix of histograms

        bins =  np.linspace(np.log(min_length),np.log(max_length), No_bins)
        binmat = np.array(0) 
        
        for a in avc_l:
            avc_binned = np.histogram(a, bins = np.exp(bins))
            avc_cdf    = len(a) - np.cumsum(avc_binned[0])
            binmat = avc_cdf if np.size(binmat) == 1 else np.vstack((binmat, avc_cdf))
        
    else: 
    # for the other plots we just use one long vector of avalanches 
        bins       = np.linspace(np.log(np.min(avc)),np.log(np.max(avc)), No_bins)
        avc_binned = np.histogram(avc, bins=np.exp(bins))
        avc_cdf    = len(avc) - np.cumsum(avc_binned[0])
    
    # Calcualte logarithmically binned c-cdf 
    #-------------------------------------------------------------------------------
    if plottype == 'scatter':
        ax.scatter(np.log(avc_binned[1][1:]), np.log(avc_cdf), color=tc)
    
    if plottype == 'line':
        ax.plot(np.log(avc_binned[1][1:]), np.log(avc_cdf), color=tc)
    
    if plottype ==  'loglog':            
        ax.loglog(avc_binned[1][1:], avc_cdf, color = tc)
            
    if plottype == 'distribution':
        ste = np.std(binmat,0) / len(avc_l)
        med = np.mean(binmat, 0)

        ax.loglog(avc_binned[1][1:], med, color = tc, linewidth = linewidth, linestyle = linestyle)
        ax.fill_between(avc_binned[1][1:], med+ste, med-ste, color = tc, alpha = 0.5)


    ax.set_ylabel('log observed frequency')
    ax.set_xlabel('log avalanche size')
    
#=======================================================================
def plot_sweep(sweep, param = 'alpha', prange = (0,6), doplot=True, cond='bl'):
#=======================================================================
    import numpy as np
    import powerlaw as pl
    import matplotlib.pyplot as plt
    from IPython.display import clear_output
    
    import logging, sys
    logging.disable(sys.maxsize)
    
    # Calculate alpha
    #-------------------------------------------------------------------   
    if param == 'alpha': 
        alphas = np.ndarray((len(sweep['Peak thresholds']), len(sweep['Time windows'])))
        for pi in range(len(sweep['Peak thresholds'])):
            for ti in range(len(sweep['Time windows'])):
                ci = 0 if cond == 'bl' else 1
                av = sweep['Avalanches'][pi][ti][ci]
                ft = pl.Fit(av, discrete=True, verbose=False); 
                alphas[pi, ti] = ft.alpha
                
        datamap  = alphas
        extentxy = (sweep['Time windows'][0], sweep['Time windows'][-1], 
                   sweep['Peak thresholds'][0], sweep['Peak thresholds'][-1])
    
    if param == 'likelihood':
        bayes_factor = np.ndarray((len(sweep['Peak thresholds']), len(sweep['Time windows'])))
        for pi in range(len(sweep['Peak thresholds'])):
            for ti in range(len(sweep['Time windows'])):
                ci    = 0 if cond == 'bl' else 1 
                av    = sweep['Avalanches'][pi][ti][ci]
                if len(av) > 0:
                    M     = len(av)
                    a     = min(av)                      # define xmin
                    b     = max(av)                      # define xmax
                    npart = 2000                            # number of particles - number of draws from prior distribution
                    ln    = cIS.IS_LN(npart, av, M, a, b) # call lognormal distribution function
                    po    = cIS.IS(npart, av, M, a, b)    # call power law function
                    bayes_factor[pi,ti] = (po[1] - ln[2])
                else: bayes_factor[pi,ti] = np.nan
                
        datamap = bayes_factor
    
    
    clear_output()
    
    # Actual plotting routines
    #-------------------------------------------------------------------           
    if doplot:
        fig, ax = plt.subplots(figsize=(12,12), ncols=1)
        img     = ax.imshow(datamap, vmin = prange[0], vmax = prange[1], origin='lower', 
                            cmap='Spectral_r', aspect='auto', extent = extentxy)
        fig.colorbar(img,ax=ax)
        plt.show
    
    return (datamap)

#=======================================================================
def plot_groupcompare(data, ax, sumline='median', xspread=0.05, connect=False):
#=======================================================================
    oldvec = []
    oldx   = []
    oldsum = []

    for a in range(data.shape[1]):
        vec = data[:,a]
        xvals = np.random.uniform(-xspread, xspread, vec.shape[0])+a+1

        if (connect) & (len(oldvec) > 0):
            for xi in range(len(xvals)):
                if vec[xi] > 1.1 * oldvec[xi]: thiscolor = 'red'
                if (vec[xi] <= 1.1*oldvec[xi]) & (vec[xi] > 0.9*oldvec[xi]): thiscolor='black'
                if vec[xi] < 0.9*oldvec[xi]: thiscolor='blue'
                plt.plot([oldx[xi], xvals[xi]], [oldvec[xi], vec[xi]], color=thiscolor, alpha = 0.1)


        if sumline == 'median':
            sumplot = np.median(vec)
            plt.plot([a+1-1.5*xspread, a+1+1.5*xspread], [sumplot, sumplot], 
                     color='black',linewidth=5)

            if (connect) & (len(oldsum) > 0):
                thiscolor = 'red' if oldsum[0] < sumplot else 'blue' 
                plt.plot([a, a+1], [oldsum[0], sumplot], '--', color=thiscolor, linewidth=3)
                if scipy.stats.ranksums(oldvec, vec).pvalue < 0.05:
                    plt.text(a+0.5, 1.1*np.max([oldvec,vec]), '*', fontsize=24)
        oldvec = vec
        oldx   = xvals
        oldsum = [sumplot]

        plt.scatter(xvals, vec)