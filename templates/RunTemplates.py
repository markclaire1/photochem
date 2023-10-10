#run all templates and updates
#also using as (temporary?) creater of global DB files from existing templates

#need to re-write me to include a pc.var.atol ramp, as well as some temp saves/are you sures?
#should also re-run MondernEarth (and others that are dicking around...)


from photochem import Atmosphere
import atmDB
reaction_file = "../photochem/data/reaction_mechanisms/zahnle_earth.yaml"

#setup if working directory is ~/photochem/templates/
templates=['Mars','ModernEarth','Titan','Saturn','Jupiter']

#'Hadean',, - temp removing becuase broken?

for t in templates:
  settings_file = t+"/settings_"+t+".yaml"
  star_file = "ModernEarth/Sun_now.txt" #keeping this fixed for now...
  if t=='Hadean':
    star_file = "Hadean/Sun_4.0Ga.txt"
  atmosphere_file = t+"/atmosphere_"+t+".txt"
  
  print(t,settings_file,star_file,atmosphere_file) #kill me eventuallly
  print('Running '+t)
  print('')
  pc=Atmosphere(reaction_file,settings_file,star_file,atmosphere_file)
  pc.var.verbose=2
  #pc.var.atol=1e-12
  pc.photochemical_equilibrium()
  pc.out2in()  # this and next not technically neccessary, but just a force of habit - redox conservation can get better on a re-run
  print('Re-Running '+t)
  print('')
  pc.photochemical_equilibrium()
  print('')
  print('Saving '+t)
  print('')
  #pc.out2atmosphere_txt(atmosphere_file,overwrite=True) #update photochem/template folder
  atmDB.out2DB(t,pc,reaction_file,star_file,settings_file) #save to DB (end eventually update?) 
