from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scienceplots

# Import netCDF file
ncfile = './test.nc'
data = Dataset(ncfile)
var = data.variables

# Prepare Data to Plot
x = var['x'][:]
y0 = var['y0'][:]  
y1 = var['y1'][:]
y2 = var['y2'][:]
y3 = var['y3'][:]
y4 = var['y4'][:]

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
    ax.plot(x, y0, label=r'$y_0$')
    ax.plot(x, y1, label=r'$y_1$')
    ax.plot(x, y2, label=r'$y_2$')
    ax.plot(x, y3, label=r'$y_3$')
    ax.plot(x, y4, label=r'$y_4$')
    ax.legend()
    fig.savefig('plot.png', dpi=300, bbox_inches='tight')
