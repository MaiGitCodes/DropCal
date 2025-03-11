# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:15:04 2025

@author: Maria Teresa Alameda

Dropcal library is designed to calibrate fluid dispensing systems under
various conditions. It uses a known calibration volume (mL) - pressure (kPa)
curve (for example water) to extrapolate the volume dispensed of other liquid
with a known viscosity (for example collagen) at a certain pressure.
"""

# -*- coding: utf-8 -*-
from numpy import linspace, ndarray
from numpy import array as nparray
from pandas import read_csv
from pandas import errors
import customtkinter as ctk
from tkinter import messagebox  # Import messagebox from tkinter
from matplotlib.pyplot import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.interpolate import interp1d
import pyscreenshot as ImageGrab
import pyperclip

# Main tkinter window
class CalibrationApp:

    def __init__(self, root):
        self.root = root
        self.root.title("Fluid Dispensation Calibration")
        
        # Frame for inputs
        input_frame = ctk.CTkFrame(root)
        input_frame.pack(padx=10, pady=10, fill="both", expand=True)
        
        ctk.CTkLabel(input_frame, text="Pressure Array (kPa):").grid(row=0, column=0, sticky="w")
        self.pressure_entry = ctk.CTkEntry(input_frame, width=300)
        self.pressure_entry.grid(row=0, column=1, padx=5, pady=5)

        ctk.CTkLabel(input_frame, text="Volume Array (mL):").grid(row=1, column=0, sticky="w")
        self.volume_entry = ctk.CTkEntry(input_frame, width=300)
        self.volume_entry.grid(row=1, column=1, padx=5, pady=5)

        ctk.CTkLabel(input_frame, text="Known Viscosity (Pa·s):").grid(row=2, column=0, sticky="w")
        self.known_viscosity_entry = ctk.CTkEntry(input_frame, width=300)
        self.known_viscosity_entry.grid(row=2, column=1, padx=5, pady=5)
        
        # Frame for New Fluid Parameters
        new_fluid_frame = ctk.CTkFrame(root, border_width=2, border_color="#ADD8E6")
        new_fluid_frame.pack(padx=10, pady=10, fill="both", expand=True)
        
        ctk.CTkLabel(new_fluid_frame, text="New Fluid Name:").grid(row=0, column=0, sticky="w")
        self.fluid_name_entry = ctk.CTkEntry(new_fluid_frame, width=300)
        self.fluid_name_entry.grid(row=0, column=1, padx=5, pady=5)
        self.fluid_name_entry.insert(0, 'set target fluid name here')
        
        ctk.CTkLabel(new_fluid_frame, text="New Fluid Viscosity (Pa·s):").grid(row=1, column=0, sticky="w")
        self.target_viscosity_entry = ctk.CTkEntry(new_fluid_frame, width=300)
        self.target_viscosity_entry.grid(row=1, column=1, padx=5, pady=5)
        self.target_viscosity_entry.insert(0, 'set new fluid name here')
        
        ctk.CTkLabel(new_fluid_frame, text="Target Volume (mL):").grid(row=2, column=0, sticky="w")
        self.target_volume_entry = ctk.CTkEntry(new_fluid_frame, width=300)
        self.target_volume_entry.grid(row=2, column=1, padx=5, pady=5)
        self.target_volume_entry.insert(0, 'set target volume here')

        # Buttons
        button_frame = ctk.CTkFrame(root)
        button_frame.pack(pady=10, fill="both", expand=True)

        ctk.CTkButton(button_frame, text="Load from File", command=self.load_csv).pack(side="left", padx=5, pady=5)
        ctk.CTkButton(button_frame, text="Calculate", command=self.update).pack(side="right", padx=5, pady=5)

        # Result Display
        self.result_label = ctk.CTkLabel(root, text="", font=('Helvetica', 12))
        self.result_label.pack(pady=10)

        # Canvas for plot
        self.figure = Figure(dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.get_tk_widget().pack(expand=True, fill="both")

        # Button to copy plot to clipboard
        ctk.CTkButton(root, text="Copy Plot to Clipboard", command=self.copy_plot_to_clipboard).pack(pady=10)

    def load_csv(self):
        pressure, volume, metadata = _load_calibration_data()
        if pressure is not None and volume is not None:
            self.pressure_entry.insert(0, ";".join(map(str, pressure)))
            self.volume_entry.insert(0, ";".join(map(str, volume)))
            
            print(metadata)
            self.known_viscosity_entry.insert(0, float(metadata['Known Viscosity (Pa s)']))           
        
    def update(self):
        # Capture inputs
        try:
            pressure_array = _array_checker(self.pressure_entry.get())
            volume_array = _array_checker(self.volume_entry.get())
            known_viscosity = float(self.known_viscosity_entry.get())
            target_viscosity = float(self.target_viscosity_entry.get())
            target_volume = float(self.target_volume_entry.get())
            fluid_name = self.fluid_name_entry.get()
    
        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input: {str(e)}")  # Use tkinter.messagebox
            return
    
        # Calculate interpolations and results
        interpolation_water = interp1d(pressure_array, volume_array, kind='linear', fill_value='extrapolate')
    
        def adjust_volume(pressure_kPa, known_viscosity, new_viscosity):
            known_volume_mL = interpolation_water(pressure_kPa)
            known_flow_rate_mL_s = known_volume_mL
            new_flow_rate_mL_s = known_flow_rate_mL_s * (known_viscosity / new_viscosity)
            new_volume_mL = new_flow_rate_mL_s
            return new_volume_mL
    
           
        # Invert the calculation: find pressure for the target volume
        def find_pressure(target_volume, known_viscosity, target_viscosity):
            # Create a function to calculate the difference between target volume and calculated volume
            def volume_difference(pressure):
                calculated_volume = adjust_volume(pressure, known_viscosity, target_viscosity)
                return calculated_volume - target_volume
    
            # Use a root-finding method to solve for pressure
            from scipy.optimize import fsolve
            initial_guess = pressure_array.mean()  # Use the mean pressure as an initial guess
            required_pressure = fsolve(volume_difference, initial_guess)[0]
            return required_pressure

        # Calculate the required pressure for the target volume
        required_pressure = find_pressure(target_volume, known_viscosity, target_viscosity)
        message = f"Required pressure for {fluid_name} to dispense {target_volume:.4f} mL: {required_pressure:.4f} kPa"
        self.result_label.configure(text=message)
    
        # Generate new pressure values for plotting
        new_pressure_kPa = linspace(min(pressure_array), max(pressure_array), 100)
        new_volume_mL = adjust_volume(new_pressure_kPa, known_viscosity, target_viscosity)
    
        # Plot calibration curves
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.scatter(pressure_array, volume_array, marker='o', s=100,
                   color='steelblue', label='Calibration curve (water)', zorder = 1)
        ax.plot(new_pressure_kPa, new_volume_mL, '-',
                color='salmon', lw = 2.0, label='Extrapolated curve (%s)' % fluid_name,
                zorder = 0)
        ax.scatter(required_pressure, target_volume, marker='x', s=100,
                   color='black', label='Target volume (%s)' % fluid_name,
                   linewidths = 3.0, zorder = 2)
        ax.set_xlabel('Pressure (kPa)')
        ax.set_ylabel('Volume (mL)')
        ax.legend()
        ax.set_title(message)
        self.canvas.draw()

    def copy_plot_to_clipboard(self):
        
        import os
        import tempfile
        import subprocess
        import platform
        
        # Get the plot as an image
        self.canvas.draw()
        bbox = (
            self.canvas.get_tk_widget().winfo_rootx(),
            self.canvas.get_tk_widget().winfo_rooty(),
            self.canvas.get_tk_widget().winfo_rootx() + self.canvas.get_tk_widget().winfo_width(),
            self.canvas.get_tk_widget().winfo_rooty() + self.canvas.get_tk_widget().winfo_height()
        )
        image = ImageGrab.grab(bbox=bbox)
        
        # Save the image to a temporary file
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp_file:
            image.save(temp_file.name, format="PNG")
            temp_file_path = temp_file.name
        
        # Copy the image to the clipboard based on the operating system
        try:
            system_platform = platform.system()
            if system_platform == "Linux":
                if os.environ.get("WAYLAND_DISPLAY"):
                    # Use wl-copy for Wayland
                    subprocess.run(["wl-copy", "--type", "image/png", "<", temp_file_path], shell=True, check=True)
                else:
                    # Use xclip for X11
                    subprocess.run(["xclip", "-selection", "clipboard", "-t", "image/png", "-i", temp_file_path], check=True)
            elif system_platform == "Windows":
                # Use the clipboard library for Windows
                import clipboard
                with open(temp_file_path, "rb") as f:
                    image_data = f.read()
                clipboard.copy(image_data)
            elif system_platform == "Darwin":  # macOS
                # Use osascript for macOS
                subprocess.run(["osascript", "-e", f'set the clipboard to (read (POSIX file "{temp_file_path}") as PNG picture)'], check=True)
            else:
                raise OSError(f"Unsupported platform: {system_platform}")
            
            # Notify the user
            messagebox.showinfo("Success", "Plot copied to clipboard!")  # Use tkinter.messagebox
        except Exception as e:
            messagebox.showerror("Error", f"Failed to copy image to clipboard: {str(e)}")  # Use tkinter.messagebox
        finally:
            # Clean up the temporary file
            os.remove(temp_file_path)


# Auxiliary functions (unchanged)
def _array_checker(array: str) -> ndarray:
    try:
        if '\n' in array:
            new_array = nparray([float(p.strip().replace(',', '.')) for p in array.split('\n') if p.strip()])
        else:
            new_array = nparray([float(p.strip().replace(',', '.')) for p in array.split(';') if p.strip()])
        return new_array
    except ValueError as e:
        raise ValueError("Input string contains non-numeric values.") from e


def _load_calibration_data():
    file_path = ctk.filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ])
    if file_path:
        try:
            metadata = {}
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            for line in lines:
                if line.startswith('#'):
                    key_value = line.strip('#').strip().split(':')
                    if len(key_value) == 2:
                        key = key_value[0].strip()
                        value = key_value[1].strip()
                        metadata[key] = value
            
            data = read_csv(file_path, sep='\t', comment='#')
            data_dict = {}
            keys = data.keys()
            
            for key in keys:
                if 'Pressure' in key or 'pressure' in key:
                    if '(Pa)' in key:
                        data_dict['Pressure'] = nparray(data[key].values) * 1e-3
                    elif '(hPa)' in key:
                        data_dict['Pressure'] = nparray(data[key].values) * 0.1
                    elif '(kPa)' in key:
                        data_dict['Pressure'] = nparray(data[key].values)
                    else:
                        data_dict['Pressure'] = nparray(data[key].values)
                
                if 'Volume' in key or 'volume' in key:
                    if '(L)' in key:
                        data_dict['Volume'] = nparray(data[key].values) * 1e3
                    elif '(dL)' in key:
                        data_dict['Volume'] = nparray(data[key].values) * 100
                    elif '(mL)' in key:
                        data_dict['Volume'] = nparray(data[key].values)
                    else:
                        data_dict['Volume'] = nparray(data[key].values)
            
            if 'Pressure' not in data_dict or 'Volume' not in data_dict:
                raise KeyError("File must contain 'Pressure' and 'Volume' columns.")
            print(metadata)
            return data_dict['Pressure'], data_dict['Volume'], metadata
        
        except FileNotFoundError:
            messagebox.showerror("Error", "File not found or cannot be read.")  # Use tkinter.messagebox
            return None, None, None
        except errors.ParserError:
            messagebox.showerror("Error", "Unable to parse the file as a CSV.")  # Use tkinter.messagebox
            return None, None, None
        except KeyError as e:
            messagebox.showerror("Error", str(e))  # Use tkinter.messagebox
            return None, None, None
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")  # Use tkinter.messagebox
            return None, None, None
    else:
        messagebox.showerror("Error", "No file selected.")  # Use tkinter.messagebox
        return None, None, None

# Run the application
if __name__ == "__main__":
    root = ctk.CTk()
    app = CalibrationApp(root)
    root.mainloop()
