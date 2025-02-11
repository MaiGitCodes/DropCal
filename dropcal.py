# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:15:04 2025

@author: Maria Teresa Alameda

Dropcal library is designed to calibrate fluid dispensing systems under
various conditions. It uses a known calibration volume (mL) - pressure (kPa)
curve (for example water) to extrapolate the volume dispensed of other liquid
with a known viscosity (for example collagen) at a certain pressure.
"""

from numpy import linspace, ndarray
from numpy import array as nparray
from pandas import read_csv
from pandas import errors
# import tkinter as tk
from tkinter import W, X, LEFT, RIGHT, BOTH
from tkinter import ttk, filedialog, messagebox
from matplotlib.pyplot import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.interpolate import interp1d


# Main tkinter window
class CalibrationApp:

    def __init__(self, root):
        self.root = root
        self.root.title("Fluid Dispensation Calibration")
        
        # Frame for inputs
        input_frame = ttk.Frame(root)
        input_frame.pack(padx=10, pady=10, fill=X)
        
        ttk.Label(input_frame, text="Pressure Array (kPa):").grid(row=0, column=0, sticky=W)
        self.pressure_entry = ttk.Entry(input_frame, width=50)
        self.pressure_entry.grid(row=0, column=1, padx=5)

        ttk.Label(input_frame, text="Volume Array (mL):").grid(row=1, column=0, sticky=W)
        self.volume_entry = ttk.Entry(input_frame, width=50)
        self.volume_entry.grid(row=1, column=1, padx=5)

        ttk.Label(input_frame, text="Known Viscosity (Pa·s):").grid(row=2, column=0, sticky=W)
        self.known_viscosity_entry = ttk.Entry(input_frame)
        self.known_viscosity_entry.grid(row=2, column=1, padx=5)

        ttk.Label(input_frame, text="Fluid Viscosity (Pa·s):").grid(row=3, column=0, sticky=W)
        self.target_viscosity_entry = ttk.Entry(input_frame)
        self.target_viscosity_entry.grid(row=3, column=1, padx=5)

        ttk.Label(input_frame, text="Target Pressure (kPa):").grid(row=4, column=0, sticky=W)
        self.target_pressure_entry = ttk.Entry(input_frame)
        self.target_pressure_entry.grid(row=4, column=1, padx=5)
        
        ttk.Label(input_frame, text="Fluid Name:").grid(row=5, column=0, sticky=W)
        self.fluid_name_entry = ttk.Entry(input_frame)
        self.fluid_name_entry.grid(row=5, column=1, padx=5)

        # Buttons
        button_frame = ttk.Frame(root)
        button_frame.pack(pady=10, fill=X)

        ttk.Button(button_frame, text="Load from File", command=self.load_csv).pack(side=LEFT, padx=5)
        ttk.Button(button_frame, text="Calculate", command=self.update).pack(side=RIGHT, padx=5)

        # Result Display
        self.result_label = ttk.Label(root, text="", font=('Helvetica', 12))
        self.result_label.pack(pady=10)

        # Canvas for plot
        self.figure = Figure(dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.get_tk_widget().pack(expand=True, fill=BOTH)

    def load_csv(self):
        pressure, volume, metadata = _load_calibration_data()
        if pressure is not None and volume is not None:
            self.pressure_entry.insert(0, ";".join(map(str, pressure)))
            self.volume_entry.insert(0, ";".join(map(str, volume)))
            self.known_viscosity_entry.insert(0, float(metadata['Known Viscosity (Pa · s)']))
            self.target_viscosity_entry.insert(0, float(metadata['Fluid Viscosity (Pa · s)']))
            self.target_pressure_entry.insert(0, float(metadata['Target Pressure (kPa)']))
            self.fluid_name_entry.insert(0, metadata['Fluid Name'])
            
    def update(self):
        # Capture inputs
        try:
            pressure_array = _array_checker(self.pressure_entry.get())
            volume_array = _array_checker(self.volume_entry.get())
            known_viscosity = float(self.known_viscosity_entry.get())
            collagen_viscosity = float(self.target_viscosity_entry.get())
            target_pressure = float(self.target_pressure_entry.get())
            fluid_name = self.fluid_name_entry.get()

        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input: {str(e)}")
            return

        # Calculate interpolations and results
        interpolation_water = interp1d(pressure_array, volume_array, kind='linear', fill_value='extrapolate')

        def adjust_volume(pressure_kPa, known_viscosity, new_viscosity):
            known_volume_mL = interpolation_water(pressure_kPa)
            known_flow_rate_mL_s = known_volume_mL
            new_flow_rate_mL_s = known_flow_rate_mL_s * (known_viscosity / new_viscosity)
            new_volume_mL = new_flow_rate_mL_s
            return new_volume_mL

        volume_fluid_mL = adjust_volume(target_pressure, known_viscosity, collagen_viscosity)
        message = f"Dispensed volume for {fluid_name} at {target_pressure} kPa: {volume_fluid_mL:.4f} mL"
        # message = 'Dispensed volume for %s at %s kPa: %.4f mL'%(
        #             fluid_name, target_pressure, volume_fluid_mL)
        
        self.result_label.config(text=message)

        new_pressure_kPa = linspace(min(pressure_array), max(pressure_array), 100)
        new_volume_mL = adjust_volume(new_pressure_kPa, known_viscosity, collagen_viscosity)

        # Plot calibration curves
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.scatter(pressure_array, volume_array, marker='o', s=100,
                   color='royalblue', label='Calibration curve (water)')
        ax.plot(new_pressure_kPa, new_volume_mL, '-',
                color='tomato', label='Extrapolated curve (%s)'%fluid_name)
        ax.scatter(target_pressure, volume_fluid_mL, marker='x', s=100,
                   color='black', label='Target volume (%s)'%fluid_name)
        ax.set_xlabel('Pressure (kPa)')
        ax.set_ylabel('Volume (mL)')
        ax.legend()
        ax.set_title(message)
        # self.figure.savefig('C:/Users/Default/Desktop/Calibration.png', dpi=100)

        self.canvas.draw()
        
#%% AUXILIAR FUNCTIONS

def _array_checker(array: str) -> ndarray:

    """
    AUXILIAR FUNCTION
    Prepares the array values given in a string format into a numpy array of floats.
    
    Expects the input string to separate pressure values by either newline ('\\n') or comma (';').
    Handles both dot ('.') and comma (',') as decimal separators within numeric values.
    
    Args:
        array (str): A string containing values separated by newline or comma.
                     Values may use either a comma or a dot as a decimal separator.
    
    Returns:
        np.ndarray: A NumPy array containing the converted float values.
    
    Raises:
        ValueError: If the input string contains non-numeric values that cannot be converted to float.
    """
    
    try:
        if '\n' in array:
            new_array = nparray([float(p.strip().replace(',', '.')) for p in array.split('\n') if p.strip()])
        else:
            print(array)
            print(array[0].strip().replace(',', '.'))
            new_array = nparray([float(p.strip().replace(',', '.')) for p in array.split(';') if p.strip()])
        return new_array
    except ValueError as e:
        raise ValueError("Input string contains non-numeric values.") from e

# Function to load CSV data
def _load_calibration_data():
    """
    AUXILIAR FUNCTION
    Opens a file dialog for the user to select a CSV or text file containing calibration data,
    then reads and returns pressure, volume values, and metadata from the file.
    
    The function relies on the user selecting a file with columns labeled 'Pressure' and 'Volume'
    and a metadata header starting with '#'.
    
    Returns:
        tuple: A tuple containing three elements:
            - A numpy array of 'Pressure' values in kPa.
            - A numpy array of 'Volume' values in mL.
            - A dictionary containing metadata (e.g., viscosity, target pressure, fluid name).
        If no file is selected, returns (None, None, None).
            
    Raises:
        FileNotFoundError: Raised if the selected file does not exist or cannot be read.
        pd.errors.ParserError: Raised if the file cannot be parsed as a CSV.
        KeyError: Raised if the file does not contain 'Pressure' and 'Volume' columns.
    """
    
    file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ])
    if file_path:
        try:
            # Read the file line by line to extract metadata
            metadata = {}
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            # Extract metadata (lines starting with #)
            for line in lines:
                if line.startswith('#'):
                    key_value = line.strip('#').strip().split(':')
                    if len(key_value) == 2:
                        key = key_value[0].strip()
                        value = key_value[1].strip()
                        metadata[key] = value
            
            # Convert Target Pressure to kPa if necessary
            if 'Target Pressure' in metadata:
                target_pressure = metadata['Target Pressure']
                if '(Pa)' in target_pressure:
                    metadata['Target Pressure (kPa)'] = float(target_pressure.replace('(Pa)', '').strip()) * 1e-3
                elif '(hPa)' in target_pressure:
                    metadata['Target Pressure (kPa)'] = float(target_pressure.replace('(hPa)', '').strip()) * 0.1
                elif '(kPa)' in target_pressure:
                    metadata['Target Pressure (kPa)'] = float(target_pressure.replace('(kPa)', '').strip())
                else:
                    metadata['Target Pressure (kPa)'] = float(target_pressure)  # Assume kPa if no unit is specified
            
            # Read the data table
            data = read_csv(file_path, sep='\t', comment='#')
            data_dict = {}
            keys = data.keys()
            
            for key in keys:
                # Handle Pressure units
                if 'Pressure' in key or 'pressure' in key:
                    if '(Pa)' in key:  # Convert Pa to kPa
                        data_dict['Pressure'] = nparray(data[key].values) * 1e-3
                    elif '(hPa)' in key:  # Convert hPa to kPa
                        data_dict['Pressure'] = nparray(data[key].values) * 0.1
                    elif '(kPa)' in key:  # Already in kPa
                        data_dict['Pressure'] = nparray(data[key].values)
                    else:  # Assume kPa if no unit is specified
                        data_dict['Pressure'] = nparray(data[key].values)
                
                # Handle Volume units
                if 'Volume' in key or 'volume' in key:
                    if '(L)' in key:  # Convert L to mL
                        data_dict['Volume'] = nparray(data[key].values) * 1e3
                    elif '(dL)' in key:  # Convert dL to mL
                        data_dict['Volume'] = nparray(data[key].values) * 100
                    elif '(mL)' in key:  # Already in mL
                        data_dict['Volume'] = nparray(data[key].values)
                    else:  # Assume mL if no unit is specified
                        data_dict['Volume'] = nparray(data[key].values)
            
            # Ensure both 'Pressure' and 'Volume' are present
            if 'Pressure' not in data_dict or 'Volume' not in data_dict:
                raise KeyError("File must contain 'Pressure' and 'Volume' columns.")
            print(metadata)
            return data_dict['Pressure'], data_dict['Volume'], metadata
        
        except FileNotFoundError:
            messagebox.showerror("Error", "File not found or cannot be read.")
            return None, None, None
        except errors.ParserError:
            messagebox.showerror("Error", "Unable to parse the file as a CSV.")
            return None, None, None
        except KeyError as e:
            messagebox.showerror("Error", str(e))
            return None, None, None
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")
            return None, None, None
    else:
        messagebox.showerror("Error", "No file selected.")
        return None, None, None
