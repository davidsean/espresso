import numpy as np
import matplotlib.pyplot as plt




name="cyl_ext"
name="cyl_int"
name="Wall"
name="Rhom"

shapes= ["Wall"]
shapes= ["Rhom_int","Rhom_ext"]
shapes=["cyl_int", "cyl_ext"]
shapes=["Rhom_int", "Rhom_ext"]
shapes=["Sphere_int", "Sphere_ext"]
shapes=["Cone_int", "Cone_ext"]
shapes= ["Line"]

shapes=["Triangle"]
for shape in shapes:
  map=np.loadtxt('phi%s.dat' %(shape))
  max=np.max(map)
  len= int(np.sqrt(np.shape(map)))
  
  nx=np.sqrt(np.shape(map))
  extent = (0, len*0.01, 0, len*0.01)

  map=map.reshape((len,len))
  map=map.transpose()

  fig=plt.figure()
  ax=fig.add_subplot(111)

  cmap=ax.imshow(map,origin='lower', clim=(0, 3.0), extent=extent)#, interpolation='none')
  cb=plt.colorbar(cmap) 
  cb.set_label(r'Potential ($k_\mathrm{B}T$)')

  ax.set_xlabel(r'$x$ pos $[\sigma]$')
  ax.set_ylabel(r'$y$ pos $[\sigma]$')

  fig.savefig("potential_%s.png" %(shape))
  plt.close(fig)


