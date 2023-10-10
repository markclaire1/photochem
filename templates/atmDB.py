def out2DB(templatename,pc,reaction_file,star_file,settings_file):
  
  test=searchDB(templatename,pc,reaction_file,star_file,settings_file,1)
  if (test=='identical'):
    print('Identical case already exists in DB. Not saving...')
    return
  
  import numpy as np
  import h5py # note use of atmDB requires 'conda install h5py' into the photochem environment 
  import getpass
  import os.path
  import yaml
  
  
  x,extras,varnames = mapinput(templatename,pc,reaction_file,star_file,settings_file)
  #retuns x (the vector representation of the atmosphere), the extras.yaml file, and variablebnames
    
  username=getpass.getuser()  #this is as close to platform independent as possible.
                #it apparently works for everything except some Windows edge cases...
  
  username='global'
  
  indexfilename = username+'Index.hdf5'
  DBfilename = username+'AtmDB.hdf5'
  
  
  x[0]=DBfilename  #add DB filename to index (with templatename being the name of the dataset in the DB)
  
  if os.path.isfile(indexfilename):
    with h5py.File(indexfilename, "a") as f:
      nelements=f[varnames[0]].shape[0] #number of atmospheres currently in the DB
      
      tnames=f['templatename']
      if templatename in [tn.decode() for tn in tnames]:
        print('An atmosphere named', templatename,'already exists in the DB.')
        return #although this prevents updating things like Saturn, and could lead to Saturnnew, etc...
        #thinking more, this could cause havoc for the VDS if there was 'Saturn' in global and in locals...
        #so if implementing might need to check names in the VDS file....
      
      for var in range(len(x)):
        f[varnames[var]].resize(nelements+1 , axis=0) #resize dataset for this variable
        f[varnames[var]][nelements] = x[var] #add new elements to datasets
  else: #first time only, create structured hdf5 file
    dTypes=['f8']*32
    dt=h5py.string_dtype(encoding='utf-8')
    dTypes[0]=dt ; dTypes[1]=dt ; dTypes[3]=dt ; dTypes[5]=dt
  
    with h5py.File(indexfilename, "a",libver='latest') as f:
      for var in range(len(x)):
        f.create_dataset(varnames[var],data=[x[var]],dtype=dTypes[var],maxshape=(None,),chunks=True,compression="gzip")

  #now save data to big DB, with header, settings.yaml, and extras.yaml as attributes  
  sol = pc.mole_fraction_dict() #returns a dictionary objects of arrays
  #>>> list(sol)
  #['alt', 'temp', 'pressure', 'density', 'H2SO4aer', 'C2H2aer', 'C2H4aer', 'C2H6aer', 'HCNaer', 'C4H4aer', 'CH3CNaer', 'HCCCNaer', 'N2Oaer', 'NH3aer', 'S8aer', 'HCaer1', 'HCaer2', 'HCaer3', 'He', 'H', 'H2', 'H2O', 'OH', 'O', 'O2', 'CO', 'HCO', 'H2CO', 'C', 'CH', 'CH2', 'CH3', 'CH4', 'S', 'S2', 'S3', 'S4', 'H2S', 'HS', 'SO', 'SO2', 'CS', 'CS2', 'OCS', 'CN', 'HCN', 'N', 'N2', 'NH', 'NH2', 'NH3', 'C2', 'C2H', 'C2H2', 'C2H4', 'NO', 'HNO', 'HO2', 'HNCO', 'N2O', 'H2O2', 'O3', 'NO2', 'NO3', 'SO3', 'HNO2', 'HNO3', 'H2SO4', 'HCl', 'Cl', 'Cl2', 'HOCl', 'N2H4', 'ClO', 'C2H6', 'CH3OH', 'CH2CO', 'CH3CHO', 'C3H4', 'C3H6', 'C4H2', 'C4H4', 'OClO', 'HCS', 'C2H3', 'C2H5', 'NCO', '1CH2', 'HCCO', 'NNH', 'HSO', 'CH3O', 'H2COH', 'H2CN', 'N2H2', 'C4H', 'HCNOH', 'C2H2OH', 'CH3CO', 'CH2CHO', 'C2H3OH', 'C2H4OH', 'CH3O2', 'HSO3', 'N2H3', 'NH2CO', 'HS4', 'C4H3', 'CH2N2', 'CH2CN', 'CH3CN', 'HCCCN']
  

  cols=len(sol.keys())
  rows=len(sol['alt']) #hardcoding to first entry in dictionary ('alt')
  atm=np.zeros((rows,cols))
  keys=list(sol)
  for i in range(cols):
    atm[:,i]=sol[keys[i]]

  #atm=np.asarray(list(sol.values()))  #gets the numeric values only (looses sol.keys() ) #this doesn't work...
  
  with open(settings_file,'r') as file:
      sf=yaml.safe_load(file)
  
  with h5py.File(DBfilename, "a") as f:
      dset = f.create_dataset(templatename, data=atm,chunks=True,compression="gzip")
      dset.attrs['header']=keys
      dset.attrs['settings']=yaml.dump(sf)
      dset.attrs['extras']=yaml.dump(extras)
  
  return

def searchDB(templatename,pc,reaction_file,star_file,settings_file,k):
    import numpy as np
    import getpass, h5py, yaml,math
    
    x,extras,varnames = mapinput(templatename,pc,reaction_file,star_file,settings_file)
    
    #Defining the weights.  Nearest neighbor is closest to 0, so LESS important variables get weighted by say 2-10
    #varnames=['DBfilename','templatename','model','model-version','reactions-file','reactions-version','species','star',
            #'surface-temperature', 'surface-eddy', 'atmosphere-grid','number-of-layers','lower-wavelength','upper-wavelength',
             #'number-of-bins','background-gas','surface-pressure','planet-mass','planet-radius','surface-albedo',
             #'photon-scale-factor','solar-zenith-angle','hydrogen-escape-type','default-gas-lower-boundary',
             #'fix-water-in-troposphere','relative-humidty','gas-rainout','rainfall-rate','tropopause-altitude',
             #'water-condensation','condensation-rate','lightning'] 
  
    weights=[0,0,1,5,1,3,1,1,2,3,3,5,3,3,2,1,1,5,5,3,3,3,1,1,1,3,3,3,4,3,5,2]
    #notes - DBfilename/templatename not included in weights...
    #modelversion (var 3) and reactions-version (var 5) diffs -> vars need to be computed here
  
    #create virtual dataset from all *Index files:
    indexfilename=createVDS(varnames)
        
    #do a k-nearest-neighbors search across the index file, followed by boundary condition search on the k cases...
    with h5py.File(indexfilename, "r") as f:
      #create a (numvariables - 2, nummodels) array
      nummodels=f[varnames[0]].shape[0]
      DB=np.empty((len(x)-2,nummodels))  #DBname/templatename not needed for nearest-neighbor search.
      for var in range(2,len(x)):
        Y=f[varnames[var]][:] #obtain *IndexDB entry for this variable and slice values into array Y
        if var==3 or var==5: # model-version/reactions-versions are strings in the X.X.X mold..
          #unlike the rest of the variables written at time of saving, they are converted to float representations here          
          DB[var-2,:]=versiondiffs(var,nummodels,x,Y,weights,varnames)
        else:
          #print(varnames[var], weights[var]*(x[var]-Y)**2)
          DB[var-2,:]=weights[var]*(x[var]-Y)**2 #squared L2 norm (no need to invoke sqrt if just comparing distances...)
    
    #(with verbose on/off) analyze 2D array
    initialdistances=np.sum(DB,axis=0) #this should sum through the variables for each models
    #so initialdistances is now 1 1D array of dimension (nummodels)
    #find k closest matches - if 1() 
    #k=3 #(setting k to 3 for now...) - could expand based on entries in VDS? - or be a calling parameter?
    note = 'unique'
    idx=np.argpartition(initialdistances,k) #indices of k closest neighbors (not necessarily in order)
    
    if (k>1): print(k,'Nearest-Neighbors to',templatename, 'are:')
    fnames=[]
    tnames=[]
    with h5py.File(indexfilename, "r") as f:
      for i in range(k):
        fnames.append(f['DBfilename'][idx[i]].decode())
        tnames.append(f['templatename'][idx[i]].decode())
        if (k>1): print('  ', '{:2.4e}'.format(initialdistances[idx[i]]),' - ', tnames[i] + ' which is saved in ' + fnames[i] )
        if(DB[4,idx[i]] != 0.0):
          if (k>1): print('  *****WARNING - THIS FILE HAS A DIIFERENT NUMBER OF SPECIES SO REQUIRES TRANSFORMATION*****')
          note='warning'
        if(DB[9,idx[i]] != 0.0):
          if (k>1): print('  *****WARNING - THIS FILE HAS A DIIFERENT NUMBER OF LAYERS SO REQUIRES TRANSFORMATION*****')
          note='warning'
    #at some point (here?) need to check/report on items like
    #different species/vertical heights, etc. which will require modifying the atmosphere file
    #idx are indicies of a vector with nummodels entries, so need to check
    #loop i over idx
    #DB[4,i] - if non-zero report that number of species is different (6th element in varnames, so 4th in DB as constructed)
    #DB[9,i] - report that number of layers is different (11th element in varnames, so 9th in DB as constructed)


    if (k>1): print('Now checking on boundary conditions...')
    #now need to run these k closest matches against boundary conditions to get final distances...


    with open(settings_file,'r') as file:
        new_set=yaml.safe_load(file)
    new_bcs=new_set['boundary-conditions']
  
    BCpenalties=np.empty(k)
    finaldistances=np.empty(k)
    for i in range(k):
      with h5py.File(fnames[i], "r") as f:
        set=yaml.safe_load(f[tnames[i]].attrs['settings'])
        DB_bcs=set['boundary-conditions']
        #print(' Penalty for',tnames[i],'is',checkBoundaryConditions(new_bcs,DB_bcs,verbose=True))
        BCpenalties[i]=checkBoundaryConditions(new_bcs,DB_bcs)
        finaldistances[i]=initialdistances[idx][i]

    if sum(BCpenalties)!=0:  #avoiding div0 in case where all k sets of boundary conditions are equal and 0
      finaldistances=finaldistances+BCpenalties/sum(BCpenalties) #normalizing BC penalties to themselves to get relative values
    finidx=np.where(finaldistances == finaldistances.min())
    bestfitindex=finidx[0][0]
    
    if (k>1): print('Best fit (including boundary conditions) is',tnames[bestfitindex],'in',fnames[bestfitindex])
    if math.isclose(finaldistances[bestfitindex],0.0, abs_tol=1e-10):
      if (k>1): print(' Identical match to an atmosphere already in the database')
      note='identical'

    with h5py.File(fnames[bestfitindex],'r') as f:
      data=f[tnames[bestfitindex]]
      header=data.attrs['header']
      write_atm_txt(header,data) #this writes temp_atm.txt to current working dir.
    
    return(note)  

    

def write_atm_txt(header,data):
  from numpy import format_float_scientific
  
  with open('temp_atm.txt', 'wt', encoding="utf-8") as fnew: #this will overwrite
    fnew.write('   ') #photochem style starts with 3 blank spaces
    for head in header: #followed by headers with toal width 27 characters
      fnew.write(f"{head:<27}") 
    fnew.write("\n")
    #python only reports 2 exponential digits, so using numpy to achieve 'house style' of 3 digit exponents
    for i in range(data.shape[0]): #number of rows (nz)
      fnew.write('   ') #photochem style starts with 3 blank spaces
      for j in range(data.shape[1]):# number of columns (number of vars + species)
        #print(i,j,data[i,j])
        datapoint=format_float_scientific(data[i,j],exp_digits=3,precision=17,unique=False)
        fnew.write(f"{datapoint:<27}")
      fnew.write("\n")
    
    

def checkBoundaryConditions(new_bcs,DB_bcs,verbose=False): #eventually over DBfilename and templatename...
  import h5py, yaml, math

  new_species=[x['name'] for x in new_bcs]
  DB_species=[x['name'] for x in DB_bcs]

  Penalty=0.
  
  #checks are:
  #1) is new species in DB species list?
      #if not - BIG PENALTY
      #if yes, then
  #2 Do the species have the same type of boundary condition?
      #if not - MEDIUM PENALTY
      #if yes, then
  #3 Do the species have the same value for the boundary condition?
      #if not - variable PENALTY
  #4 then check if DB species list has species that are not in new list --> BIG PENALTY    
  
  for i in range(len(new_bcs)):
    try:
      species_match_index=DB_species.index(new_species[i])  #Check 1
      #if verbose: print(new_species[i], 'is index',species_match_index,'in DB_species')  
      
      if verbose: print('  NEW',new_bcs[i])
      if verbose: print('  DB ',DB_bcs[species_match_index])
      #now check on 'type' - Check 2
      try:
        new_bcs[i]['type'] # the only option that has 'type' at top level is short-lived
        try:
          DB_bcs[species_match_index]['type']
          if verbose: print(new_species[i],'is short-lived in both templates. Awesome!')
        except KeyError:
          Penalty +=1
          if verbose: print(new_species[i], 'is short-lived in new template but not in DB. Penalty =',Penalty ) 
      except KeyError: #there are normal upper and lower boundary conditions...
        
        #Lower Boundary Conditions first
        
        newLBCtype=new_bcs[i]['lower-boundary']['type']
        DBLBCtype=DB_bcs[species_match_index]['lower-boundary']['type'] 
        #can be Moses,flux, mix, vdep, vdep + dist flux 
        #Moses vs. anything else should be biggest PENALTY
        if(newLBCtype=='Moses'):
          if (DBLBCtype=='Moses'):
            if verbose: print(new_species[i],'LBC is Moses in both templates. Awesome!')
          else:
            Penalty +=1
            if verbose: print(new_species[i], 'LBC is Moses in new template but',DBLBCtype,'in DB. Penalty =',Penalty)
        else: #LBC is flux, mix, vdep, or vdep + dist flux
          #lets set all combinations other than vdep+flux vs flux to have penalty 0.5
          #but vdep + flux vs flux has penalty 0.2 as potential easier to deal with...
          if (newLBCtype==DBLBCtype): #a match
            if (newLBCtype=='mix'):
              if new_bcs[i]['lower-boundary']['mix']==DB_bcs[species_match_index]['lower-boundary']['mix']:
                if verbose: print(new_species[i],'LBC is identical in both templates. Awesome!')
              else:
                new_mix=float(new_bcs[i]['lower-boundary']['mix'])
                DB_mix=float(DB_bcs[species_match_index]['lower-boundary']['mix'])
                #constructing this so each order of magnitude difference in mix is 0.01 penalty points
                #but weights larger absolute magnitudes mixing ratios by a factor of up to ~6.5
                #weight applies above mixing ratios of 1e-6, increasing to a factor of 6 at 1e-1...
                weight=max(math.log10(max(new_mix,DB_mix))+7,1.0)  #weight the penalty so that larger maximum mixing ratios are weighted a bit more
                Penalty += 0.01*abs(math.log10(new_mix)-math.log10(DB_mix))*weight              #this weight is constructed to stops applying at mixing ratios of 1e-6 
                if verbose: print(new_species[i],'mixing ratio LBC has different values. Penalty =',Penalty)
            if ('flux' in newLBCtype):  #using if/in to catch both flux and flux + vdep
              if new_bcs[i]['lower-boundary']['flux']==DB_bcs[species_match_index]['lower-boundary']['flux']:
                if verbose: print(new_species[i],'flux is identical in both templates. Awesome!')
              else:
                new_flux=float(new_bcs[i]['lower-boundary']['flux'])
                DB_flux=float(DB_bcs[species_match_index]['lower-boundary']['flux'])
                #1 order difference of flux magnitude is a penalty of 0.01, increased as absolute magnittude go higher
                weight=max(math.log10(max(new_flux,DB_flux))-4,1.0)  #1e10 -> factor of 6 increase
                Penalty += 0.01*abs(math.log10(new_flux)-math.log10(DB_flux))*weight 
                if verbose: print(new_species[i],'flux LBC has different values. Penalty =',Penalty) 
            if ('vdep' in newLBCtype): #using if/in to catch both vdep and 'vdep + flux'                
              if new_bcs[i]['lower-boundary']['vdep']==DB_bcs[species_match_index]['lower-boundary']['vdep']:
                if verbose: print(new_species[i],'vdep is identical in both templates. Awesome!')
              else:
                new_vdep=float(new_bcs[i]['lower-boundary']['vdep'])
                DB_vdep=float(DB_bcs[species_match_index]['lower-boundary']['vdep'])
                #1 order difference of vdep magnitude is a penalty of 0.01, increased as absolute magnittude go higher
                #but have to account for case of 0 vdep, which we will recast as vdep=1e-10
                new_vdep=max(new_vdep,1e-10) ; DB_vdep=max(DB_vdep,1e-10)
                weight=max(math.log10(max(new_vdep,DB_vdep))+7,1.0)  #vdep=1 --> factor of 6, 1e-1 --> 5, etc.
                Penalty += 0.01*abs(math.log10(new_vdep)-math.log10(DB_vdep))*weight 
                if verbose: print(new_species[i],'vdep LBC has different values. Penalty =',Penalty)
            if (newLBCtype=='vdep + dist flux'): #only need to check height now..
              if new_bcs[i]['lower-boundary']['height']==DB_bcs[species_match_index]['lower-boundary']['height']:
                if verbose: print(new_species[i],'distributed height is identical in both templates. Awesome!')
              else:
                new_height=float(new_bcs[i]['lower-boundary']['height'])
                DB_height=float(DB_bcs[species_match_index]['lower-boundary']['height'])
                #not super important - go linear - each 1km diff => 0.005
                Penalty += 0.005*abs(new_height-DB_height) 
                if verbose: print(new_species[i],'distributed height has different values. Penalty =',Penalty)
          #reach here if NOT a direct match between types.
          elif ('flux' in newLBCtype and 'flux' in DBLBCtype):
            Penalty += 0.3
            if verbose: print(new_species[i], 'LBC is',newLNCtype,'and',DB_bcs[species_match_index],'is',DBLBCtype,' Penalty =',Penalty)
          else: #reach here for anything other than total match or vdep+flux vs flux
            Penalty += 0.6
            if verbose: print(new_species[i], 'LBC is',newLBCtype,'in new and',DBLBCtype,'in DB. Penalty =',Penalty)
          
        #now UPPER BOUNDARY CONDITIONS
        newUBCtype=new_bcs[i]['upper-boundary']['type']
        DBUBCtype=DB_bcs[species_match_index]['upper-boundary']['type']
        if (newUBCtype == DBUBCtype): #a match
          cond='veff'
          if(newUBCtype==cond):             
            if new_bcs[i]['upper-boundary'][cond]==DB_bcs[species_match_index]['upper-boundary'][cond]:
              if verbose: print(new_species[i],'UBC',cond,'is identical in both templates. Awesome!')
            else:
              new=float(new_bcs[i]['upper-boundary'][cond])
              DB=float(DB_bcs[species_match_index]['upper-boundary'][cond])
              #1 order of magnitude difference of veff magnitude is a penalty of 0.01, increased as absolute magnittude go higher
              #but have to account for case of 0 veff, which we will recast as veff=1e-10
              new=max(new,1e-10) ; DB=max(DB,1e-10)
              weight=max(math.log10(max(new,DB)),1.0)  #1e10 -> factor of 10 increase
              Penalty += 0.01*abs(math.log10(new)-math.log10(DB))*weight
              if verbose: print(new_species[i],cond ,'UBC has different values. Penalty =',Penalty)
          cond='flux'
          if(newUBCtype==cond):             
            if new_bcs[i]['upper-boundary'][cond]==DB_bcs[species_match_index]['upper-boundary'][cond]:
              if verbose: print(new_species[i],'UBC',cond,'is identical in both templates. Awesome!')
            else:
              new=float(new_bcs[i]['upper-boundary'][cond])
              DB=float(DB_bcs[species_match_index]['upper-boundary'][cond])
              #1 order of magnitude difference of flux magnitude is a penalty of 0.01, increased as absolute magnittude go higher
              #but have to account for case of 0 flux (not commonly used), which we will recast as veff=1e-10
              new=max(new,1e-10) ; DB=max(DB,1e-10)
              weight=max(math.log10(max(new,DB)),1.0)  #1e10 -> factor of 10 increase
              Penalty += 0.01*abs(math.log10(new)-math.log10(DB))*weight
              if verbose: print(new_species[i],cond ,'UBC has different values. Penalty =',Penalty)    
        else: #UBC's are different
          Penalty += 0.6
          if verbose: print(new_species[i], 'UBC is',newUBCtype,'in new and',DBUBCtype,'in DB. Penalty =',Penalty)
          
    except ValueError: #penalty for failing check 1
      Penalty += 1  #penalty for missing check 1 (aka new species not in DB species)
      if verbose: print(new_species[i], 'is not in DB_species list -- PENALTY =',Penalty)
      #should I consider a check for vdeps only as perhaps less of a penalty????

  #Check 4
  for i in range(len(DB_bcs)):
    try:
      species_match_index=new_species.index(DB_species[i])  #passes if in the list
    except ValueError:
      Penalty += 1
      if verbose: print(DB_species[i], 'is not in new_species list -- PENALTY =',Penalty)
      
  return(Penalty)
  
def versiondiffs(var,nummodels,x,Y,weights,varnames):
  import numpy as np
  xminusY=np.empty(nummodels)
  for model in range(nummodels):
    version=Y[model].decode() #first - need to convert bytes to unicode
    test=diffver(x[var],version) #returns tuple of (version depth tested, difference)
    if test[1]==0:
      xminusY[model]=0.0  #versions are identical
    else:
      versiondiff,numdiffs=test
      match versiondiff:
        case 1:
          multiplier=1.0  #major version diff
        case 2:
          multiplier=0.1  #minor version diff
        case 3:
          multiplier=0.01 #build version diff
        case _:
          print('Please add a transformation for this version number in atmDB.py/searchDB')
          print(version, versiondiff)
          return
      xminusY[model]=multiplier*numdiffs  #ack this is actually x-Y, so handle differently
  #print(varnames[var], weights[var]*xminusY**2)
  return(weights[var]*xminusY**2)
    
def diffver(v1, v2):
    v1s = map(int, v1.split('.'))
    v2s = map(int, v2.split('.'))

    for ii, (v1r, v2r) in enumerate(zip(v1s, v2s), 1):
        if v1r != v2r:
            return ii, v1r - v2r

    return ii, 0
    #found on stackoverflow - returns a two element tuple (a,b)
    #where a is the element of the version number which differs
    #  (a=1 implies major version number, a=2 is minor version number, a=3 is build number, etc.)
    #and b is the difference in distance.
    #ex >>> print(diffver('2.2.3','1.6.4')) returns (1, 1) - 1 major version number apart
    #   >>> print(diffver('2.2.3','2.0.17')) returns (2, 2) - 2 minor version numbers apart
    #   >>> print(diffver('2.2.3','2.2.3')) returns (3, 0) - 0 version numbers apart
    # since this invovles a comparison of 2, needs to be called from within the nearest neighbor routine.
  


#code for potential devtool to combine *.IndexDB...  
#This copies the data from each dataset in the original file to the new file using the original dataset names. It loops to copy ALL root level datasets. This requires datasets in each file to have different names. The data is not merged into one dataset.

#with h5py.File('table_copy.h5',mode='w') as h5fw:
#    for h5name in glob.glob('file*.h5'):
#        h5fr = h5py.File(h5name,'r') 
#        for obj in h5fr.keys():        
#            h5fr.copy(obj, h5fw)       

def createVDS(varnames):
  import h5py
  import glob
  
  Indexfiles=glob.glob('*Index.hdf5') # a list of index files to concatenate
  shapes=[h5py.File(file)[varnames[0]].shape[0] for file in Indexfiles] #get the shapes (aka number of models) in each dataset
  
  dTypes=['f8']*32
  dt = h5py.string_dtype(encoding='utf-8')
  dTypes[0]=dt ; dTypes[1]=dt ; dTypes[3]=dt ; dTypes[5]=dt

  with h5py.File("VDS.h5","w",libver='latest') as f: #note this is set to 'w' which overwrites each time...
    for j,var in enumerate(varnames): 
      layout=h5py.VirtualLayout(shape=(sum(shapes),),dtype=dTypes[j])
      startingindex=0 #since the shapes are variable, need to keep track of indexing element
      for i in range(len(shapes)):
        vsource=h5py.VirtualSource(Indexfiles[i],var,shape=(shapes[i],),dtype=dTypes[j])
        finalindex=sum(shapes[:i+1]) #python indexes from 0, but layouts are indexed from 1
        layout[startingindex:finalindex]=vsource[:]
        startingindex=sum(shapes[:i+1])
      f.create_virtual_dataset(var, layout)
  return('VDS.h5')

  
    
def mapinput(templatename,pc,reaction_file,star_file,settings_file):
    #describe purpose (file to create mapping - notes available...)
  
    import yaml
    import math
    from photochem import __version__
    from pathlib import Path

    #create extras
    #IMPORTANT NOTE - model photochem is hardcoded here. need to fix before incorporating any ATMOS files...
    
    
    with open(reaction_file,'r') as file:
      rx=yaml.safe_load(file)

    extras={'templatename': templatename, 'model': 'photochem', 'model-version': __version__, 'species': pc.dat.nsp, 'reactions-file': Path(reaction_file).name, 'reactions-version': rx['version'], 'star': Path(star_file).name, 'surface-temperature': pc.var.temperature[0], 'surface-eddy': pc.var.edd[0]}


    x=['DBfilename.DB']  #(0) DB filename - overwritten prior to save
    x.append(templatename) #(1) - add templatename string 
    
    if extras['model'] == 'photochem':  #note currently hardcoded above! fix me!!!!
      normvar=0.0
    elif extras['model'] == 'ATMOS':
      normvar=0.5
    else:
      print('Please add a normalised transform for this PHOTOCHEMICAL MODEL into atmDB.py')
      return
        
    x.append(normvar) #(2) add index for model name
    
    #model-version
    #model versions are strings of integers in "standard" major.minor.build.etc order.
    #this script will break if non-numeric characters are used, so hopefully no one will!
    try:
      normvar=extras['model-version']
    except ValueError:
      print('Model-version is non-numeric which breaks atmDB.py. Please update the model version and re-compile')
      return
    else:
      x.append(normvar) #(3) add model-version string
      #NOTE - No normalization here. This comparison (and weighting) needs to be done in real time in the KNN search

    #reactions-file
    if extras['reactions-file'] == 'zahnle_earth.yaml':
      normvar=0.0
    else:
      print('Please add a normalised transform for this REACTIONS FILE into atmDB.py')
      return
    
    x.append(normvar) #(4) reactions-file

    #reactions-version
    try:
      normvar=extras['reactions-version']
    except ValueError:
      print('Reactions-version is non-numeric which breaks atmDB.py. Please update the reactions version number')
      return
    else:
      x.append(normvar) #(5) add reaction-version string
      #NOTE - No normalization here. This comparison (and weighting) needs to be done in real time in the KNN search  
        
    #species
    #Making the assumption that any additional species groups (e.g. extended Cl, Sorg, etc.)
    # will have a unique 'number of species'.  
    # Probably likely - and if it annoyingly later turns out untrue, could be captured with an if
    # and or captured by boundary condition checks?
    #=111 for the default Wogan case. Assume Min 50 max 250 mean 150 so 
    normvar = (extras['species'] - 150) / (150 - 50)
    x.append(normvar) #(6) species
    #report in NN algorithm if returning different than NN match (as requires transformation of atmosphere file) 

  
      
    #stellar
    #requires human intervention when adding new star - here's a rough guide to assigning normvars:
    #  -1 for A, -0.5 for F, 0 for G2 stars, +0.3 for K, +0.5 for M8, +0.6 for M6, +0.7 for M4, +0.8 for M3, +0.9 for M2, +1.0 for M1
    #  scaling from -0.2 to 0 for sun through time -0.2 is 4Ga

    if extras['star'] =='Sun_now.txt':
	     normvar=0.0
    elif extras['star'] =='Sun_4.0Ga.txt':
	     normvar=-0.2
    else:
	     print('Please add a normalised transform for this STAR into atmDB.py')
	     return
    x.append(normvar) #(7) - add index for star
    #surface-temperature
    # assuming min=100,max=1000, mean=550, mapping min to -1
    normvar= (extras['surface-temperature']-550)/(550-100)
    x.append(normvar)  #(8) - add index for surface temperature.
    
    #surface-eddy-diffusion value
    #map logarithmically between 1e2 and 1e8, min 2, max8 mean 5
    normvar = (math.log10(extras['surface-eddy']) - 5 ) / (5-2)
    x.append(normvar)  #(9)- add index for surface eddy diffusion
    
    ##Now get stuff from the settings file:
    with open(settings_file,'r') as file:
      s=yaml.safe_load(file)
    
    #atmosphere-grid
    #going to use top-bottom as metric here. Max range so far is 1e7 for earth, 1.6e8 for Saturn
    #assuming not much interest in going less than 60km (6e6) if doing photochemistry, lets set loose max at 2000km (2e8).
    #in log space this is 6.7 -> 8.3, mean 7.5
    print(type(s['atmosphere-grid']['top']))
    normvar = (math.log10(float(s['atmosphere-grid']['top']) - float(s['atmosphere-grid']['bottom']))-7.5) / (7.5 - 6.7)
    x.append(normvar) #(10) add index for atmosphere grid height
    #report in NN algorithm if returning different than NN match (as requires transformation of atmosphere file) 
    
    #number-of-layers
    # min=10, max=1000, so mean = 550
    normvar= (float(s['atmosphere-grid']['number-of-layers'])-550)/(550-10)
    x.append(normvar) #11- index for number of layers
    #report in NN algorithm if returning different than NN match (as requires transformation of atmosphere file)    
    
    #photolysis-grid
    #currently code supports only regular-grid, but matching on these three should be future proof as long as lower,upper, and nw are set... 
  
    # assuming lower will probably never be greater than Lyman alpha, but perhaps could go shorter for extended cases?
    #low=1,high=121.6 —> mean 61.3 
    normvar= (float(s['photolysis-grid']['lower-wavelength'])-61.3)/(61.3-1)
    x.append(normvar) #(12) - index for lower wavelength
    
    #while photochemistry becomes boring beyond the visible, some could compute longer to match with climate codes?
    #low:700, high 20000 -> mean 10350 but use log weighting min=2.8, max=4.3, mean=3.55
    normvar=(math.log10(float(s['photolysis-grid']['upper-wavelength']))-3.55)/(3.55-2.8)
    x.append(normvar) #(13) - index for upper wavelength
    
    #assuming never lower than 118 (ATMOS default), but could be huge if anyone ever requires line-by-line
    #low 118, high=100000, mean 50100.0—> log space: min 2.07, max 5, mean 3.54
    normvar= (math.log10(float(s['photolysis-grid']['number-of-bins']))-3.54)/(3.54-2.07)
    x.append(normvar) #(14)- index for number of wavelengths
    #NN note - weight me!
    
    #planet stuff
    
    #background-gas 
    # mapping CO2=-1, H2=-0.5, N2=0 with -2 (??) ,2 (He?) available for future use. Should heavily weight this one
    if s['planet']['background-gas'] =='CO2':
      normvar=-1.0
    elif s['planet']['background-gas']=='H2':
      normvar=1.0
    elif s['planet']['background-gas']=='N2':
	     normvar=0.0
    else:
	     print('Please add a normalised transform for this BACKGROUND GAS into atmDB.py')
       #take
	     return
    x.append(normvar) #(15)- index for background gas
    #NN note - weight me!

    #surface-pressure
    #min=0.001, max 100, so do a logarithmic weighting here. Should heavily weight this one.
    #min=-3, max =2 —> mean —> -0.5
    normvar= (math.log10(float(s['planet']['surface-pressure'])) + 0.5) / ( -0.5 + 3)
    x.append(normvar) #(16) - index for surface pressure
    #NN note - weight me!
    
    #planet-mass
    #min=1e25, max 1e31 —> log so min=25, max=31 —> mean=28
    normvar= (math.log10(float(s['planet']['planet-mass']))-28)/(28 -25)
    x.append(normvar) #(17) - index for planet mass
    
    #planet-radius
    #min=1e8, max 1e10 —> log so min=8, max=10 —> mean=9
    normvar= (math.log10(float(s['planet']['planet-radius']))-9)/(9 - 8)
    x.append(normvar) #(18) - index for planet radius
      
    #surface-albedo
    #min=0, max=1 -> mean=0.5
    normvar= (float(s['planet']['surface-albedo'])-0.5)/(0.5)
    x.append(normvar) #(19) - index for surface albedo
    
    
    #photon-scale-factor - optional parameter, but stored with defaults to 1 if not defined
    #default is 1, range logarithmic - say 1e-3 - 1e3 -> -3,3,0
    
    try:
      psf=s['planet']['photon-scale-factor']
    except KeyError: 
      normvar=0.0 # log(mean) value for photo-scale-factor
    else:
      normvar= (math.log10(float(s['planet']['photon-scale-factor'])))/(-3)
      x.append(normvar) #(20)- index for photon-scale-factor (weighted along with star?)
    #NN note - weight with star somehow?

    #solar-zenith-angle
    #min=0, max=90, mean=45
    normvar= (float(s['planet']['solar-zenith-angle'])-45)/(45)
    x.append(normvar) #(21)
  
    #hydrogen-escape
    if s['planet']['hydrogen-escape']['type'] =='diffusion limited':
      normvar=-1.0
    elif s['planet']['hydrogen-escape']['type']=='zahnle':
      normvar=-0.5
    elif s['planet']['hydrogen-escape']['type']=='none':
      normvar=1.0
    else:
      print('Please add a normalised transform for this HYDROGEN ESCAPE TYPE into atmDB.py')
      return
    x.append(normvar) #(22) - index for hydrogen escape
    
    #default-gas-lower-boundary: Moses (optional - if not set - defaults to deposition velocity (= 0))
    try:
      dlb=s['planet']['default-gas-lower-boundary']
    except KeyError: 
      normvar=1.0 #value for deposition velocity
    else:
      normvar=-1.0 #value for Moses
      if dlb != 'Moses':
        print('Please add a normalised transform for this DEFAULT LOWER BOUNDARY into atmDB.py')
        return
    x.append(normvar) #(23)- index for default lower boundary conditions

    #fix-water-in-troposphere: boolean
    if (s['planet']['water']['fix-water-in-troposphere']):
      normvar=1.0
    else:
      normvar=-1.0
    x.append(normvar) #(24) - index for fixed water in trop...

    #relative-humidity: optional - can be a float between 0 and 1, or ‘manabe’ 
    try:
      RH=float(s['planet']['water']['relative-humidty'])
    except ValueError:
      #value is not a float should imply it is a string and that string is ‘manabe’
      normvar=-0.5
      if s['planet']['water']['relative-humidty'] != 'manabe':
        print('Please add a normalised transform for this RELATIVE HUMIDITY into atmDB.py')
        return
    except KeyError:  #doesn't exist ( if fix-water-in-trop == false)
      normvar=-1.0
    else:  
      #RH should be a float between 0 and 1 - can just map directly!
	    normvar=RH
    x.append(normvar) #(25) - index for relative humidity
    
    
    #gas-rainout - boolean
    if (s['planet']['water']['gas-rainout']):
      normvar=1.0
    else:
      normvar=-1.0
    x.append(normvar) #(26)- index for gas rainout
    
    
    #rainfall-rate: I’m optional
    try:
      rain=s['planet']['water']['rainfall-rate']
    except KeyError:
      normvar=-1.0
    else:
      #rainfall rate is between say 1e-9 and 1e2? -> log scaling min-9, max=2, mean = -3.5
	    normvar=(math.log10(float(rain)) + 3.5) / ( -3.5 + 9)
    x.append(normvar) #27 - index for gas rainout
    
  
        
    #tropopause-altitude - optional
    try:
      trop=s['planet']['water']['tropopause-altitude']
    except KeyError:
      normvar=-1.0
    else:
      # in cm, should be say between 1e5 (1 km) and 1e7 (100km) - log 5,7 mean 6
	    normvar=(math.log10(float(trop))-6) #/(6-5)
    x.append(normvar) #(28) - index for tropopause altitude
    

    #water-condensation: boolean
    if (s['planet']['water']['water-condensation']):
      normvar=1.0
    else:
      normvar=-1.0
    x.append(normvar) #(29) - index for water condensation

    #condensation-rate: optional -  just index on A - range in current templates is 1e-8 - 1e-5
    try:
      cond=s['planet']['water']['condensation-rate']['A']
    except KeyError:
      normvar=-1.0
    else:
      #lets range between 1e-10 and 1e-2 - log space this is -10,-2 mean -6
	    normvar=(math.log10(float(cond)) +6) / (-6 -10)
    x.append(normvar) #(30) - index for condensation rate
    
    
    #lightning
    if extras['model'] == 'photochem':  #no lightning in photochem
      normvar=-1.0
    elif extras['model'] == 'ATMOS': #lightning is an option in ATMOS.  Store rate (mapped from [0 to 1)
      print('Please implement PRONO readin for ATMOS/LIGHTNING into atmDB.py')
      return
      #need  to read in PRONO - scale from say 1e0 to 1e10 in log space and map to [0,1+) min= 0, max = 10, mean = 5
      #first attempt below, but test/incorporate when needed
      PRONO=math.log10(1e8)
      normvar=PRONO/10 #this takes 1e-1 to less than 0, but that's fine...  
      
    else:
      print('Please add a normalised transform for LIGHTNING into atmDB.py')
      return
    
    x.append(normvar) #(31) - index for lightning

#for reference
    varnames=['DBfilename','templatename','model','model-version','reactions-file','reactions-version','species','star',
              'surface-temperature', 'surface-eddy', 'atmosphere-grid','number-of-layers','lower-wavelength','upper-wavelength',
               'number-of-bins','background-gas','surface-pressure','planet-mass','planet-radius','surface-albedo',
               'photon-scale-factor','solar-zenith-angle','hydrogen-escape-type','default-gas-lower-boundary',
               'fix-water-in-troposphere','relative-humidty','gas-rainout','rainfall-rate','tropopause-altitude',
               'water-condensation','condensation-rate','lightning'] 
    
    #hdf5 is picky about datatypes, and in particular does not like python dType = U so be explicit
    #if updating, update also in createVDS()
    dTypes=[float]*32
    dTypes[0]=str ; dTypes[1]=str ; dTypes[3]=str ; dTypes[5]=str

    y=[t(x) for t,x in zip(dTypes,x)]

    
    #if there is a need to extend this file, will also need to add (and fill) a new dataset into the global
    #and local IndexDB.hdf5 files using custom devtools
      
    
    return y,extras,varnames


    
