def RunMod(reaction_file,settings_file,star_file,atmosphere_file):
#this is an attemmpt to automate finding a solution when the input atmosphere 
#is out of equilibrium with the input files.
#returns tuple: success(boolean), atol(float), rtol(float), pc (atmosphere object)

  from photochem import Atmosphere

  pc=Atmosphere(reaction_file,settings_file,star_file,atmosphere_file)
  pc.var.verbose=2
  pc.var.mxsteps=2000 #first successful iteration only, then drops to 500
  
  pc.var.rtol=1e-2
  pc.var.atol=1e-12
  Atolcount=0
  
  while(pc.var.atol > 1e-27):  #need to be sure to not infinite loop when unconverged
  
    print('running with rtol=',pc.var.rtol, 'and atol=',pc.var.atol)
    Atolcount=Atolcount+1
    converged=pc.photochemical_equilibrium()
    
    if (converged):
      pc.out2in()
      pc.var.mxsteps=500
      if (Atolcount < 30):  #nominal case is about 9 runs...
        pc.var.atol=pc.var.atol/50
        print('')
      else:
        print('Photochem is bouncing around but not reaching the end goal')
        return (True,pc.var.atol,pc.var.rtol,pc)  #bouncing around too much...might be salvageable, but needs more attention than this script can provide...
    else:
      print('Photochem did not converge, trying again with higher absolute tollerance')
      pc.var.atol=pc.var.atol*4   #goes x4, x16, x64
      if (pc.var.atol > 1e-6):
        print('Photochem did not converge, and is heading for a runaway. Giving up...')
        #this captures a runaway case....
        return (False,pc.var.atol,pc.var.rtol,pc) #should actually return a boolean, maybe other info?

  print('whoop')
  pc.var.atol=1e-27
  pc.var.rtol=5e-3
  
  while(pc.var.rtol >= 1e-3):
    
    print('running with rtol=',pc.var.rtol, 'and atol=',pc.var.atol)
    
    converged=pc.photochemical_equilibrium()
    if (converged):
      pc.out2in()
      pc.var.rtol=pc.var.rtol/3
      
    else:
      print('Photochem did not converge, trying again with higher relative tollerance')
      pc.var.rtol=pc.var.rtol*2
  
  print('one more for good measure... (needed?)')
  pc.var.rtol=1e-3
  print('running with rtol=',pc.var.rtol, 'and atol=',pc.var.atol)
  
  
  converged=pc.photochemical_equilibrium()
  if (converged):
    pc.out2in()
    print('triple whoop')
  else:
    print('well poopie pants')
    return (False,pc.var.atol,pc.var.rtol,pc) #this shouldn't happen?

  return (True,pc.var.atol,pc.var.rtol,pc)
  
  #pc.out2atmosphere_txt(atmosphere_file,overwrite=True) #update photochem/template folder
  