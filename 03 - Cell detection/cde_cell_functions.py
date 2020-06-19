def cde_cell_fishspec(Fdata, prefx = ''):
    import os
    import re
    
    pull_planes = 1 if prefx == 'PL' else 2 
    if prefx: prefx = '^' + prefx + '.*'
        
    names = os.listdir(Fdata)
    r     = re.compile(prefx + 'ZFRR.*')
    folds = list(filter(r.match, names))
 
    Zfish = []
    for f in folds:
        cfld = next(os.walk(Fdata + os.sep + f))[1]
        Cond = []
        ci   = 0
        for c in cfld:
            Tpaths = []
            tifs = os.listdir(Fdata + os.sep + f + os.sep + c)
            r    = re.compile('^[0-9A-Z].*[tif|tiff|TIF|TIFF]$')
            tifs = list(filter(r.match, tifs))
            Tpaths = []
            for t in tifs:
                Tpaths.append(Fdata + os.sep + f + os.sep + c + os.sep + t)
            
            # If in single plane tifs, pull them together
            #-----------------------------------------------------------
            if pull_planes == 1:
                plid   = []
                for t in tifs: plid.append(int(t[-6:-4]))
                    
                plst   = set(plid)
                Planes = []
                for p in plst:
                    t = [i for i,x in enumerate(plid) if x==p]
                    PLtifs = []
                    for ti in t:
                        PLtifs.append(str(Tpaths[ti]))
                    Planes.append({'Tpaths':PLtifs})
                        
                        
            Cond.append({'Name':c, 
                         'Path':Fdata + os.sep + f + os.sep + c, 
                         'Tifs':tifs,
                         'Tpaths':Tpaths
                        })
            
            if pull_planes == 1: Cond[ci].update({'Plane':Planes})
            ci = ci + 1
            
        Zfish.append({'Cond':Cond, 'Name':f[len(prefx)-2:]})
    
    return Zfish

#=======================================================================
def cde_cell_planesave(Fdata, Fish, cname = 'all', mxpf = 7500):
#=======================================================================

    from skimage import io
    from PIL import Image
    import os
    import numpy as np
    import datetime

    # Make new folder for plane-divided data
    #--------------------------------------------------------------
    Fplanes = Fdata + os.sep + 'PL_' + Fish["Name"]
    if not os.path.isdir(Fplanes): os.mkdir(Fplanes)

    if cname.lower() == 'all':
        crange = list(range(0,len(Fish["Cond"])))
    else:
        crange = []
        for ci in range(0,len(Fish["Cond"])): 
            if Fish["Cond"][ci]["Name"].lower() == cname.lower():
                crange.append(ci) 

    print('I found ' + str(len(crange)) + ' Conditions')

    # Loop through conditions
    #--------------------------------------------------------------
    for c in crange:

            Fcond = Fplanes + os.sep + Fish["Cond"][c]["Name"]
            if not os.path.isdir(Fcond): os.mkdir(Fcond)

            # Get shape information on the zplanes
            #----------------------------------------------------------
            ttif = Image.open(Fish["Cond"][c]["Tpaths"][0])
            npln = ttif.n_frames
            
#             ttif = io.imread(Fish["Cond"][c]["Tpaths"][0])
#             npln = ttif.shape[0]

            print('Condition ' + Fish["Cond"][c]["Name"])
            print('> There are ' + str(npln) + ' Planes')
            tifs = Fish["Cond"][c]["Tpaths"]

            # Loop through tifs and collate in single numpy array
            #---------------------------------------------------------
            for pl in range(0,npln):
                print('> Processing plane ' + str(pl+1))

                # Load individual tifs in batches
                #----------------------------------------------------
                btch = list(range(0,len(tifs), mxpf))
                if len(tifs) > btch[-1]: btch.append(len(tifs))            
                for bi in range(len(btch)-1):
                    
                    # Load first plane to prime array
                    #-------------------------------------------------
#                     pln = io.imread(tifs[0])[pl,:,:]

                    img = Image.open(tifs[0])
                    img.seek(pl)
                    pln = np.array(img)
                    img.close()
                
                    pln = pln.reshape(1, pln.shape[0], pln.shape[1])
                    for i in range(btch[bi],btch[bi+1]):
                        if i%500 == 0: print('> ', end = ' ') 
                        img = Image.open(tifs[i])
                        img.seek(pl)
                        ldd = np.array(img)
                        img.close()
#                         ldd = io.imread(tifs[i])[pl,:,:]
                        ldd = ldd.reshape(1, ldd.shape[0], ldd.shape[1])
                        pln = np.concatenate((pln, ldd), axis = 0)

                    # Save new tif
                    #---------------------------------------------------
                    ps = str(pl)
                    bs = str(bi)
                    pth  = Fcond + os.sep
                    fnm  = (Fish["Name"] + '_s' + bs.zfill(2) 
                            + '_' + Fish["Cond"][c]["Name"][0] 
                            + '_PL' + ps.zfill(2) + '.tif')
                    currentDT = datetime.datetime.now()
                    print(currentDT.strftime("%Y-%m-%d %H:%M:%S") + '| Saving interim file')
                        
                    io.imsave(pth + fnm, pln)   
         