from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scienceplots

# Import netCDF file
ncfile = './test.nc'
data = Dataset(ncfile)
var = data.variables

# Prepare Data to Plot
x = var['x'][:]
y = var['y'][:]  

# Plot params
pparam = dict(
    xlabel = r'$x$',
    ylabel = r'$y$',
    xscale = 'linear',
    yscale = 'linear',
)

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    ax.plot(x, y, label=r'$grf$')
    ax.legend()
    fig.savefig('plot.png', dpi=300, bbox_inches='tight')
