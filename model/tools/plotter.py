import glob
import h5py
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from matplotlib.ticker import ScalarFormatter
from numpy.polynomial.legendre import legfit
from matplotlib.colors import Normalize

def model_info(verbose=False):
    """
    List out the avaliable models in the directory
    """
    dir = os.getcwd()
    file_pattern = dir+'/*Central*/outfile'
    if verbose:
        print(f"{file_pattern} : ")
        print("Avaliable models : ")
    # get all files
    files = sorted(glob.glob(file_pattern))
    if verbose:
        for file in files:
            print(f"{file}")
    return files

class HDF5Plotter:
    def __init__(self, file_pattern="rkiter-*.hdf5", files_range=None, exclude_last_n_files=0):
        """
        Initialize the HDF5Plotter class.

        Parameters:
        - file_pattern: str, pattern to match HDF5 files
        - files_range: tuple, (start_index, end_index) to select a subset of files
        """
        # Constants for unit conversions
        self.rhocgs2code = 1.6190e-18
        self.timecgs2code = 2.03017e5
        self.lencgs2code = 6.77193e-06
        self.gausscgs2code = 1.1974017397874592e-20
        self.prscgs2code = 1.80139e-39
        self.epscgs2code = 5.59429e-55

        # Variable mapping with unit conversion factors or functions
        # Each entry is a dictionary with keys:
        # - group: the HDF5 group name
        # - index: index within the group dataset (None if not applicable)
        # - unit_conversion: function or lambda for unit conversion (None if no conversion)
        self.variable_map = {
            'rho': {
                'group': 'primitive',
                'index': 0,
                'unit_conversion': lambda data: data/self.rhocgs2code,
            },
            'vx': {
                'group': 'primitive',
                'index': 1,
                'unit_conversion': lambda data: data/self.lencgs2code*self.timecgs2code,
            },
            'vy': {
                'group': 'primitive',
                'index': 2,
                'unit_conversion': lambda data: data/self.lencgs2code*self.timecgs2code,
            },
            'vz': {
                'group': 'primitive',
                'index': 3,
                'unit_conversion': lambda data: data/self.lencgs2code*self.timecgs2code,
            },
            'p': {
                'group': 'primitive',
                'index': 4,
                'unit_conversion': lambda data: data/self.prscgs2code,
            },
            'ye': {
                'group': 'primitive',
                'index': 5,
                'unit_conversion': None,
            },
            'bx': {
                'group': 'bfield',
                'index': 0,
                'unit_conversion': lambda data: data/self.gausscgs2code,
            },
            'by': {
                'group': 'bfield',
                'index': 1,
                'unit_conversion': lambda data: data/self.gausscgs2code,
            },
            'bz': {
                'group': 'bfield',
                'index': 2,
                'unit_conversion': lambda data: data/self.gausscgs2code,
            },
            'epsilon': {
                'group': 'epsilon',
                'index': None,
                'unit_conversion': lambda data: data/self.epscgs2code,
            },
            'temp': {
                'group': 'temp',
                'index': None,
                'unit_conversion': None,
            },
            'phi': {
                'group': 'phi',
                'index': None,
                'unit_conversion': None,
            },
            'omega': {
                'group': 'primitive',
                'index': 3,  # Using 'vz' for 'omega'
                'unit_conversion': lambda data, r: (data/r*self.timecgs2code),
            },
            'face': {
                'group': 'face',
                'index': None,
                'unit_conversion': None,
            },
        }

        self.default_plot_params = {
            'xscale': 'linear',
            'yscale': 'linear',
            'xlabel': None,
            'ylabel': None,
            'title': None,
            'grid': True,
            'label': None,
            'legend': True,
            'formatter': False,
            'scientific': False,
            'useOffset': False,
            'xlim': None,
            'ylim': None
        }

        self.var_axlabels = {
            'rho': r'$\rho$ (g/cm$^3$)',
            'vx': r'$v_x$ (cm/s)',
            'vy': r'$v_y$ (cm/s)',
            'vz': r'$v_z$ (cm/s)',
            'p': r'$P$ (dyn/cm$^2$)',
            'ye': None,
            'bx': r'$B_x$ (G)',
            'by': r'$B_y$ (G)',
            'bz': r'$B_z$ (G)',
            'epsilon': r'$\epsilon$ (erg/cm$^3$)',
            'temp': r'$T$ (MeV)',
            'phi': r'$\phi/c^2$',
            'omega': r'$\omega$ (s$^{-1}$)',
            'face': None,
            'linestyle': 'solid',
            'color': None
        }

        self.file_pattern = file_pattern
        self.files_range = files_range
        self.exclude_last_n_files = exclude_last_n_files
        self.filenames = self._find_files()

    def _find_files(self):
        # Find all relevant HDF5 files based on pattern and files_range
        matched_files = glob.glob(self.file_pattern)
        # Extract iteration numbers and sort based on them
        files_with_numbers = []
        for filename in matched_files:
            match = re.search(r'rkiter-(\d+)-nm\.hdf5', filename)
            if match:
                iter_num = int(match.group(1))
                files_with_numbers.append((iter_num, filename))
        # Sort the list based on iteration numbers
        sorted_files_with_numbers = sorted(files_with_numbers, key=lambda x: x[0])
        # Extract filenames in order
        sorted_filenames = [filename for (iter_num, filename) in sorted_files_with_numbers]
        # Exclude the last n files if specified
        if self.exclude_last_n_files > 0:
            sorted_filenames = sorted_filenames[:-self.exclude_last_n_files]
        # Now apply files_range if provided
        if self.files_range:
            start, end = self.files_range
            sorted_filenames = sorted_filenames[start:end]
        return sorted_filenames

    def _get_variable_data(self, f, var_name):
        mapping = self.variable_map[var_name]
        group = mapping['group']
        index = mapping['index']
        # Read the data from the HDF5 file
        if group == 'bfield':
            rface = f['x-interface']
            thface = f['y-interface']
            dset = f['bfield']
            bfield = dset[:]
            bfield = bfield.T
            bfacex = bfield[0,:,:,:]
            bfacey = bfield[1,:,:,:]
            bfacez = bfield[2,:,:,:]
            if var_name == 'bx':
                data = np.zeros((len(rface)-3, len(thface)-3,1))
                for i in range(0, len(rface)-3):
                    data[i,:,:] = 1.5*(rface[i+2]-rface[i+1])/(rface[i+2]**3-rface[i+1]**3)\
                            * (rface[i+2]**2*bfacex[i+2,2:-1,2:-1] + rface[i+1]**2*bfacex[i+1,2:-1,2:-1])
                data = data.T[0]
            elif var_name == 'by':
                data = np.zeros((len(rface)-3, len(thface)-3,1))
                for i in range(0, len(thface)-3):
                    data[:,i,:] = 0.5*(thface[i+2]-thface[i+1])\
                            * (np.sin(thface[i+2])*bfacey[3:,i+2,2:-1] + np.sin(thface[i+1])*bfacey[3:,i+1,2:-1])\
                            / np.abs(np.cos(thface[i+2])-np.cos(thface[i+1]))
                data = data.T[0]
            elif var_name == 'bz':
                data = f['bfield'][:].T[2].T[2][2:-1,2:-1]
        else:
            if index is not None:
                data = f[group][:].T[index].T[0][1:-1,1:-1]
            else:
                data = f[group][0][1:-1,1:-1]
        # Apply unit conversion if needed
        unit_conversion = mapping.get('unit_conversion')
        if unit_conversion:
            if var_name == 'omega':
                r = self._get_radius_data(f)*self.lencgs2code
                data = unit_conversion(data, r)
            else:
                data = unit_conversion(data)
        return data

    def _get_radius_data(self, f):
        # Adjust r based on group
        r = (f['x-interface'][2:-1] + f['x-interface'][1:-2])/2
        r = r / self.lencgs2code
        return r
    
    def _get_dr_data(self, f):
        # Adjust r based on group
        dr = f['x-interface'][2:-1] - f['x-interface'][1:-2]
        dr = dr / self.lencgs2code
        return dr
    
    def _get_dV_data(self, f):
        r = self._get_radius_data(f)
        th = self._get_theta_data(f)
        dr = self._get_dr_data(f)
        dth = self._get_dth_data(f)
        r_grid, th_grid = np.meshgrid(r, th)
        dr_grid, dth_grid = np.meshgrid(dr, dth)
        
        dV = r_grid**2*np.sin(th_grid)*dr_grid*dth_grid
        return dV
    
    def _get_theta_data(self, f):
        th = f['y-interface'][:]
        th = (f['y-interface'][2:-1] + f['y-interface'][1:-2])/2
        return th
    
    def _get_dth_data(self, f):
        dth = f['y-interface'][2:-1] - f['y-interface'][1:-2]
        return dth
    
    def _get_latitude_data(self, f):
        y_interface = f['y-interface'][:]
        lat = (y_interface[1:] + y_interface[:-1]) / 2
        lat = np.rad2deg(lat) - 90  # Convert from radians to degrees and adjust
        return lat

    def _select_data_slice(self, data):
        # Select the data slice (e.g., index 10 as in original code)
        if data.ndim > 1:
            return data[20]
        else:
            return data

    def plotevolveprofile(self, var_name, bounce=False,plot_params=None, custom=None, close_up=False, ax=None):
        """
        Plot the specified variable from the HDF5 files.

        Parameters:
        - var_name: str, the variable name to plot
        - use_log_axes: dict, e.g., {'x': True, 'y': True} (default None)
        - xlim: tuple, e.g., (xmin, xmax)
        - ylim: tuple, e.g., (ymin, ymax)
        - close_up: bool, if True, adjust xlim for a close-up plot
        """
        if var_name not in self.variable_map:
            raise ValueError(f"Invalid variable name. Must be one of {list(self.variable_map.keys())}")

        # Get default plot parameters
        if plot_params is not None:
            plot_settings = {**self.default_plot_params, **plot_params}
        else:
            plot_settings = self.default_plot_params
        
        plot_settings['xlabel'] = plot_settings['xlabel'] or r'$r_{eq}$ (km)'
        plot_settings['ylabel'] = plot_settings['ylabel'] or self.var_axlabels[var_name]


        # Get the colormap
        cmap = plt.get_cmap('viridis')

        if ax is None:
            fig, ax = plt.subplots()

        # Set up plot
        plt.title(f"{var_name}")
        if bounce:
            self.filenames = self.filenames[60:][::4]
        num_files = len(self.filenames)
        time = np.zeros(num_files)
        for i, filename in enumerate(self.filenames):
            with h5py.File(filename, 'r') as f:
                r = self._get_radius_data(f)/1e5  # Convert to km
                # Select appropriate data slice
                if var_name == 'face':
                    if custom == 'effU':
                        vphi = self._get_variable_data(f, 'vz')*self.lencgs2code/self.timecgs2code
                        phi = self._get_variable_data(f, 'phi')
                        rho = self._get_variable_data(f, 'rho')*self.rhocgs2code
                        prs = self._get_variable_data(f, 'p')*self.prscgs2code
                        data = 0.5*vphi**2 + phi + prs/rho
                        data = self._select_data_slice(data)
                    if custom == 'j':
                        vphi = self._get_variable_data(f, 'vz')*self.lencgs2code/self.timecgs2code
                        data = vphi*self._get_radius_data(f)
                        data = self._select_data_slice(data)
                elif var_name == 'bz':
                    # For 'bz', sum over absolute values along axis=0
                    data = np.sum(np.abs(data), axis=0)
                else:
                    data = self._get_variable_data(f, var_name)
                    data = self._select_data_slice(data)
                # Prepare color
                time[i] = f['time'][()] / self.timecgs2code  # Convert time to physical units (e.g., milliseconds)
                color = cmap(i / num_files)
                # Plot data
                ax.plot(r, data, label=filename, color=color)

        # print(vphi[20])
        ax.set_xscale(plot_settings['xscale'])
        ax.set_yscale(plot_settings['yscale'])
        ax.set_xlabel(plot_settings['xlabel'])
        ax.set_ylabel(plot_settings['ylabel'])
        ax.set_title(plot_settings['title'])
        ax.grid(plot_settings['grid'])
        ax.legend().remove()
        
        # Apply axis limits
        if plot_settings['xlim']:
            ax.set_xlim(plot_settings['xlim'])
        if plot_settings['ylim']:
            ax.set_ylim(plot_settings['ylim'])

        # Add a colorbar to indicate the temporal order
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=time[0], vmax=time[-1]))
        cbar = plt.colorbar(sm)
        cbar.set_label('Time (s)')
        # plt.show()
        return ax


    def plotevolveline(self, var_name, location='center', method='value', plot_params=None, ax=None):
        """
        Plot the time evolution of a variable at a specific location or using a specific method.

        Parameters:
        - var_name (str): The variable name to plot. Must be one of the keys in self.variable_map.
        - location (str or int, optional): Specifies the location to plot the variable. Can be:
            - 'center': The center of the domain.
            - 'max': The maximum value in the domain.
            - 'min': The minimum value in the domain.
            - int: An index specifying a specific location.
            Default is 'center'.
        - method (str, optional): Specifies how to process data over the domain. Can be:
            - 'value': Plot the value at the specified location.
            - 'sum': Sum the values over the domain.
            - 'average': Average the values over the domain.
            Default is 'value'.
        - plot_params (dict, optional): A dictionary containing plot parameters to override the default settings. 
            Example keys include:
            - 'xscale': 'linear' or 'log' (default is 'linear')
            - 'yscale': 'linear' or 'log' (default is 'linear')
            - 'xlim': tuple (xmin, xmax) to set the x-axis limits
            - 'ylim': tuple (ymin, ymax) to set the y-axis limits
            - 'xlabel': str to set the x-axis label
            - 'ylabel': str to set the y-axis label
            - 'title': str to set the plot title
            - 'grid': bool to enable/disable grid (default is True)
            - 'legend': bool to enable/disable legend (default is True)
        - ax (matplotlib.axes.Axes, optional): An existing matplotlib axes object to plot on. If None, a new figure and axes are created.

        Raises:
        - ValueError: If var_name is not in self.variable_map.

        Returns:
        - matplotlib.axes.Axes: The axes object with the plot.

        Example usage:
        plotter = HDF5Plotter()
        
        # Using default settings
        fig, ax = plotter.plotevolveline('rho')

        # Overriding default settings with plot_params dictionary
        fig, ax = plotter.plotevolveline(
            'rho',
            location='max',
            method='average',
            plot_params={
                'xscale': 'log',
                'yscale': 'linear',
                'xlim': (1e-1, 1e3),
                'ylim': (1e-2, 1e2),
                'title': 'Custom Title',
                'xlabel': 'Custom X Label',
                'ylabel': 'Custom Y Label'
            }
        )
        """
        if var_name not in self.variable_map:
            raise ValueError(f"Invalid variable name. Must be one of {list(self.variable_map.keys())}")

        # Get default plot parameters
        if plot_params is not None:
            plot_settings = {**self.default_plot_params, **plot_params}
        else:
            plot_settings = self.default_plot_params

        plot_settings['xlabel'] = plot_settings['xlabel'] or 'Time (s)'
        plot_settings['ylabel'] = plot_settings['ylabel'] or self.var_axlabels[var_name]

        # Initialize arrays to store time and variable values
        num_files = len(self.filenames)
        time = np.zeros(num_files)
        var_values = np.zeros(num_files)

        for i, filename in enumerate(self.filenames):
            with h5py.File(filename, 'r') as f:
                if var_name == 'face':
                    pass
                # Process data based on location and method
                elif var_name == 'bz':
                    # For 'bz', sum over absolute values along axis=0
                    data = np.sum(np.abs(data), axis=0)
                elif var_name == 'omega':
                    # For 'omega', compute omega = vz / r
                    r = self._get_radius_data(f)
                    lat = self._get_latitude_data(f)
                    vz_data = data
                    data = self._select_data_slice(vz_data)/r/np.cos(np.deg2rad(lat))
                else:
                    data = self._get_variable_data(f, var_name)
                    data = self._select_data_slice(data)

                # Determine the value based on the specified location or method
                if var_name == 'face':
                    pass
                elif location == 'center':
                    value = data[0]
                elif location == 'max':
                    value = np.max(data)
                elif location == 'min':
                    value = np.min(data)
                elif isinstance(location, int):
                    value = data[location]
                else:
                    raise ValueError("Invalid location. Must be 'center', 'max', 'min', or an integer index.")

                if var_name == 'face':
                    pass
                elif method == 'value':
                    var_values[i] = value
                elif method == 'sum':
                    var_values[i] = np.sum(data)
                elif method == 'average':
                    var_values[i] = np.mean(data)
                elif method == 'mass':
                    data = self._get_variable_data(f, var_name)
                else:
                    raise ValueError("Invalid method. Must be 'value', 'sum', or 'average'.")
        
                if var_name == 'rho' and method == 'mass':
                    # Compute the total mass
                    dV = self._get_dV_data(f)
                    mask = data > 1e11
                    var_values[i] = np.sum(data[mask]*dV[mask])*2*np.pi/1.989e33

                if var_name == 'face' and method == 'btotenergy':
                    # Compute the total magnetic energy
                    bx = self._get_variable_data(f, 'bx')
                    by = self._get_variable_data(f, 'by')
                    bz = self._get_variable_data(f, 'bz')
                    dV = self._get_dV_data(f)
                    benergy = 0.5*(bx**2 + by**2 + bz**2)
                    var_values[i] = np.sum(benergy*dV)*2*np.pi
                elif var_name == 'face' and method == 'btorenergy':
                    # Compute the toroidal magnetic energy
                    bz = self._get_variable_data(f, 'bz')
                    dV = self._get_dV_data(f)
                    benergy = 0.5*(bz**2)
                    var_values[i] = np.sum(benergy*dV)*2*np.pi
                elif var_name == 'face' and method == 'bpolenergy':
                    # Compute the poloidal magnetic energy
                    bx = self._get_variable_data(f, 'bx')
                    by = self._get_variable_data(f, 'by')
                    dV = self._get_dV_data(f)
                    benergy = 0.5*(bx**2 + by**2)
                    var_values[i] = np.sum(benergy*dV)*2*np.pi
                elif var_name == 'face' and method == 'gravenergy':
                    # Compute the gravitational potential energy
                    phi = self._get_variable_data(f, 'phi')/self.lencgs2code*self.timecgs2code**2
                    dV = self._get_dV_data(f)
                    var_values[i] = abs(np.sum(phi*dV)*2*np.pi)
                elif var_name == 'face' and method == 'bratio':
                    # Compute the ratio of toroidal to poloidal magnetic energy
                    bx = self._get_variable_data(f, 'bx')
                    by = self._get_variable_data(f, 'by')
                    bz = self._get_variable_data(f, 'bz')
                    dV = self._get_dV_data(f)
                    bpol = 0.5*(bx**2 + by**2)
                    btor = 0.5*(bz**2)
                    var_values[i] = np.sum(btor*dV)/np.sum(bpol*dV)
                elif var_name == 'face' and method == 'maxbnorm':
                    # Compute the normal magnetic field
                    bx = self._get_variable_data(f, 'bx')
                    by = self._get_variable_data(f, 'by')
                    bz = self._get_variable_data(f, 'bz')
                    bnorm = np.sqrt(bx**2 + by**2 + bz**2)
                    var_values[i] = np.max(bnorm)
                elif var_name == 'face' and method == 'rotenergy':
                    # Compute the rotational kinetic energy
                    rho = self._get_variable_data(f, 'rho')
                    vz = self._get_variable_data(f, 'vz')
                    dV = self._get_dV_data(f)
                    energy = 0.5*rho*(vz**2)
                    var_values[i] = np.sum(energy*dV)*2*np.pi
                elif var_name == 'face' and method == 'radenergy':
                    # Compute the radial kinetic energy
                    rho = self._get_variable_data(f, 'rho')
                    vx = self._get_variable_data(f, 'vx')
                    dV = self._get_dV_data(f)
                    energy = 0.5*rho*(vx**2)
                    var_values[i] = np.sum(energy*dV)*2*np.pi
                elif var_name == 'face' and method == 'angenergy':
                    # Compute the angular momentum
                    rho = self._get_variable_data(f, 'rho')
                    vy = self._get_variable_data(f, 'vy')
                    dV = self._get_dV_data(f)
                    energy = 0.5*rho*(vy**2)
                    var_values[i] = np.sum(energy*dV)*2*np.pi
                elif var_name == 'face' and method == 'Lz':
                    # Compute the angular momentum
                    r = self._get_radius_data(f)
                    theta = self._get_theta_data(f)
                    r_grid, theta_grid = np.meshgrid(r, theta)
                    R = r_grid*np.sin(theta_grid)
                    rho = self._get_variable_data(f, 'rho')
                    vz = self._get_variable_data(f, 'vz')
                    dV = self._get_dV_data(f)
                    Lz = R*rho*vz
                    var_values[i] = np.sum(Lz*dV)*2*np.pi
                    

                time[i] = f['time'][()] / self.timecgs2code  # Convert time to physical units (e.g., milliseconds)

        # Plot the time evolution
        if ax is None:
            fig, ax = plt.subplots()
            alpha = 1
        else:
            alpha = 0.8
        
        plot_args = {
            'label': plot_settings['label'],
            'linestyle': plot_settings['linestyle'],
            'alpha': alpha
        }

        if 'color' in plot_settings:
            plot_args['color'] = plot_settings['color']

        ax.plot(time, var_values, **plot_args)
        ax.set_xscale(plot_settings['xscale'])
        ax.set_yscale(plot_settings['yscale'])
        ax.set_xlabel(plot_settings['xlabel'])
        ax.set_ylabel(plot_settings['ylabel'])
        ax.set_title(plot_settings['title'])
        ax.grid(plot_settings['grid'])
        if plot_settings['legend']:
            ax.legend()
        
        # Apply axis limits
        if plot_settings['xlim']:
            ax.set_xlim(plot_settings['xlim'])
        if plot_settings['ylim']:
            ax.set_ylim(plot_settings['ylim'])
        return ax
    
    def plot_data_slice(self, var_name, slice_index, plot_params=None):
        """
        Plot a slice of the specified variable from the HDF5 files.

        Parameters:
        - var_name: str, the variable name to plot
        - slice_index: int, the index of the slice to plot
        - plot_params: dict, plot parameters to override the default settings

        Returns:
        - matplotlib.axes.Axes: The axes object with the plot
        """
        if var_name not in self.variable_map:
            raise ValueError(f"Invalid variable name. Must be one of {list(self.variable_map.keys())}")

        # Get default plot parameters
        if plot_params is not None:
            plot_settings = {**self.default_plot_params, **plot_params}
        else:
            plot_settings = self.default_plot_params

        # Initialize arrays to store time and variable values
        num_files = len(self.filenames)
        time = np.zeros(num_files)
        var_values = np.zeros(num_files)

        filename = self.filenames[slice_index]
        with h5py.File(filename, 'r') as f:
            r = self._get_radius_data(f) / 1e5
            # Select appropriate data slice
            data = self._get_variable_data(f, var_name)
            data = data[slice_index]

    def plot_map(self, var_name, method=None, plot_params=None, ax=None):
        """
        This function generates a 2D map of coordinates.
        Example case is for looking at the trajectory of the location of maximum field strength.
        We find the coordinates of the maximum field strength in each time step and plot a line on the map.

        Parameters:
        - var_name: str, the variable name to plot
        - plot_params: dict, plot parameters to override the default settings
        
        Returns:
        - matplotlib.axes.Axes: The axes object with the plot
        """
        if var_name not in self.variable_map:
            raise ValueError(f"Invalid variable name. Must be one of {list(self.variable_map.keys())}")

        # Get default plot parameters
        if plot_params is not None:
            plot_settings = {**self.default_plot_params, **plot_params}
        else:
            plot_settings = self.default_plot_params
        
        plot_settings['xlabel'] = plot_settings['xlabel'] or r'$r_{eq}$ (km)'
        plot_settings['ylabel'] = plot_settings['ylabel'] or self.var_axlabels[var_name]

        # Initialize arrays to store time and coordinates of maximum field strength
        num_files = len(self.filenames)
        time = np.zeros(num_files)
        max_r = np.zeros(num_files)
        max_th = np.zeros(num_files)

        for i, filename in enumerate(self.filenames):
            with h5py.File(filename, 'r') as f:
                r = self._get_radius_data(f) / 1e5  # Convert to km
                th = self._get_theta_data(f)
                # print(data.shape)
                if var_name == 'face':
                    if method == 'bnorm':
                        bx = self._get_variable_data(f, 'bx')
                        by = self._get_variable_data(f, 'by')
                        bz = self._get_variable_data(f, 'bz')
                        data = np.sqrt(bx**2 + by**2 + bz**2)
                    elif method == 'btor':
                        bz = self._get_variable_data(f, 'bz')
                        data = np.abs(bz)
                    elif method == 'bpol':
                        bx = self._get_variable_data(f, 'bx')
                        by = self._get_variable_data(f, 'by')
                        data = np.sqrt(bx**2 + by**2)
                else:
                    data = self._get_variable_data(f, var_name)
                max_idx = np.unravel_index(np.argmax(data, axis=None), data.shape)
                max_r[i] = r[max_idx[1]]
                max_th[i] = th[max_idx[0]]
                # Store the time
                time[i] = f['time'][()] / self.timecgs2code  # Convert time to physical units (e.g., milliseconds)

        if ax is None:
            fig, ax = plt.subplots()
        # Create the plot
        sc = ax.plot(max_r, max_th)
        sc = ax.scatter(max_r, max_th, c=time, cmap='viridis', norm=Normalize(vmin=time.min(), vmax=time.max()))
        plt.colorbar(sc, label='Time (ms)')
        ax.set_xlabel('Radius (km)')
        ax.set_ylabel('Theta (radians)')
        ax.set_title(f'Trajectory of Maximum {var_name}')
        ax.grid(True)

        plot_args = {
            'label': plot_settings['label'],
        }

        if 'color' in plot_settings:
            plot_args['color'] = plot_settings['color']

        ax.set_xscale(plot_settings['xscale'])
        ax.set_yscale(plot_settings['yscale'])
        ax.set_xlabel(plot_settings['xlabel'])
        ax.set_ylabel(plot_settings['ylabel'])
        ax.set_title(plot_settings['title'])
        ax.grid(plot_settings['grid'])
        if plot_settings['legend']:
            ax.legend()
        
        # Apply axis limits
        if plot_settings['xlim']:
            ax.set_xlim(plot_settings['xlim'])
        if plot_settings['ylim']:
            ax.set_ylim(plot_settings['ylim'])
        return ax
        
    def plot_variable_vs_lat_time(
        self,
        var_name,
        lat_axis_index=2,
        slice_indices=None,
        legendre_order=None,
        use_divergent_colormap=False,
        cmap=None
    ):
        """
        Plot a variable as a function of latitude and time.

        Parameters:
        - var_name: str, the variable name to plot
        - lat_axis_index: int, the axis index corresponding to latitude (default 2)
        - slice_indices: dict, indices to slice other axes (e.g., {'radius': slice(5, 10)})
        - legendre_order: int, order of Legendre decomposition (if None, no decomposition)
        - use_divergent_colormap: bool, whether to use a divergent colormap
        - cmap: str or Colormap, specify the colormap to use
        """
        if var_name not in self.variable_map:
            raise ValueError(f"Invalid variable name. Must be one of {list(self.variable_map.keys())}")

        # Initialize arrays to store time and data
        num_files = len(self.filenames)
        time = np.zeros(num_files)
        data_list = []

        # Collect latitude data from the first file
        with h5py.File(self.filenames[0], 'r') as f:
            lat = self._get_latitude_data(f)[1:-1]

        for i, filename in enumerate(self.filenames):
            with h5py.File(filename, 'r') as f:
                # Read time data
                time_data = f['time'][()] / self.timecgs2code  # Convert time to physical units
                time[i] = time_data

                # Read variable data using existing methods
                data = self._get_variable_data(f, var_name)[2:-1,3:]

                # Apply slicing if specified
                if slice_indices:
                    for axis_name, sl in slice_indices.items():
                        axis = {'radius': 0, 'theta': 1, 'phi': 2}.get(axis_name)
                        if axis is not None and axis < data.ndim:
                            data = np.take(data, indices=sl, axis=axis)
                        else:
                            raise ValueError(f"Invalid axis name '{axis_name}' in slice_indices.")

                # For Legendre decomposition
                if legendre_order is not None:
                    # Perform Legendre decomposition along the latitude axis
                    # data shape: (radius, theta, phi) or similar
                    # Here, we assume latitude corresponds to theta (axis=lat_axis_index)
                    cos_lat = np.cos(np.deg2rad(lat))
                    coeffs = legfit(cos_lat, data, deg=legendre_order)[1]
                    # Use the specified Legendre coefficient (e.g., order 1)
                    # For example, take the first order coefficient
                    data = coeffs[1]

                else:
                    # Sum over other axes if necessary
                    data = np.mean(data, axis=tuple(range(data.ndim))[:lat_axis_index] + tuple(range(data.ndim))[lat_axis_index+1:])

                data_list.append(data)

        # Prepare data for plotting
        time_grid, lat_grid = np.meshgrid(time, lat)

        data_array = np.array(data_list).T  # Shape: (lat, time)

        # Plotting
        plt.figure()
        if cmap is None:
            if use_divergent_colormap:
                cmap = 'seismic'
            else:
                cmap = 'plasma'
        plt.pcolormesh(lat_grid, time_grid, data_array, cmap=cmap)
        plt.colorbar(label=var_name)
        plt.xlabel('Latitude (degrees)')
        plt.ylabel('Time (ms)')
        plt.title(f'{var_name} vs. Latitude and Time')
        plt.show()

    def plot_variable_vs_radius_time(
        self,
        var_name,
        rad_axis_index=0,
        slice_indices=None,
        use_divergent_colormap=False,
        cmap=None
    ):
        """
        Plot a variable as a function of radius and time.

        Parameters:
        - var_name: str, the variable name to plot
        - rad_axis_index: int, the axis index corresponding to radius (default 0)
        - slice_indices: dict, indices to slice other axes (e.g., {'theta': slice(5, 10)})
        - use_divergent_colormap: bool, whether to use a divergent colormap
        - cmap: str or Colormap, specify the colormap to use
        """
        if var_name not in self.variable_map:
            raise ValueError(f"Invalid variable name. Must be one of {list(self.variable_map.keys())}")

        # Initialize arrays to store time and data
        num_files = len(self.filenames)
        time = np.zeros(num_files)
        data_list = []

        # Collect radius data from the first file
        with h5py.File(self.filenames[0], 'r') as f:
            radius = self._get_radius_data(f, var_name)

        for i, filename in enumerate(self.filenames):
            with h5py.File(filename, 'r') as f:
                # Read time data
                time_data = f['time'][()] / self.timecgs2code  # Convert time to physical units
                time[i] = time_data

                # Read variable data using existing methods
                data = self._get_variable_data(f, var_name)

                # Apply slicing if specified
                if slice_indices:
                    for axis_name, sl in slice_indices.items():
                        axis = {'radius': 0, 'theta': 1, 'phi': 2}.get(axis_name)
                        if axis is not None and axis < data.ndim:
                            data = np.take(data, indices=sl, axis=axis)
                        else:
                            raise ValueError(f"Invalid axis name '{axis_name}' in slice_indices.")

                # Reduce data to 1D along radius
                data = np.mean(data, axis=tuple(range(data.ndim))[rad_axis_index+1:])

                data_list.append(data)

        # Prepare data for plotting
        time_grid, radius_grid = np.meshgrid(time, radius)

        data_array = np.array(data_list).T  # Shape: (radius, time)

        # Plotting
        plt.figure()
        if cmap is None:
            if use_divergent_colormap:
                cmap = 'seismic'
            else:
                cmap = 'plasma'
        plt.pcolormesh(radius_grid, time_grid, data_array, cmap=cmap)
        plt.colorbar(label=var_name)
        plt.xlabel('Radius (km)')
        plt.ylabel('Time (ms)')
        plt.title(f'{var_name} vs. Radius and Time')
        plt.show()