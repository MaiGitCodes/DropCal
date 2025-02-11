# Dropcal

Overview

The Dropcal library is designed to calibrate fluid dispensing systems under various conditions. It uses a known calibration curve (volume vs. pressure) for a reference fluid (e.g., water) to extrapolate the volume dispensed for another fluid with a known viscosity (e.g., collagen) at a specific pressure. The library provides a graphical user interface (GUI) for easy interaction and visualization of results.
Features

    Calibration Curve Extrapolation: Uses a reference calibration curve to calculate dispensed volumes for fluids with different viscosities.

    GUI Support: Provides a user-friendly interface for loading calibration data, inputting parameters, and visualizing results.

    File Support: Reads calibration data and metadata from .txt or .csv files.

    Unit Conversion: Automatically handles pressure and volume units (e.g., kPa, Pa, hPa, mL, L, dL).

Installation

    Ensure you have Python 3.x installed.

    Install the required dependencies:
    bash
    Copy

    pip install numpy pandas matplotlib scipy tkinter

    Download or clone the Dropcal library files (dropcal.py and main.py).

Usage
Running the Application

    Place the dropcal.py and main.py files in the same directory.

    Run the main.py script:
    bash
    Copy

    python main.py

    The GUI will open, allowing you to:

        Load calibration data from a file or input the different values by hand.

        Input known viscosity, target viscosity, target pressure, and fluid name.

        Calculate and visualize the dispensed volume for the target fluid.

Calibration File Format
    The calibration file should be a .txt or .csv file with the following structure:
    Metadata Header: Lines starting with # contain metadata (e.g., viscosity, target pressure, fluid name).
    Data Table: A tab-separated table with columns for Pressure (kPa) and Volume (mL).
    Pressure and Volume can be given in other units but they should be indicated in the column header. 

Example (calibration_file.txt):

        # Known Viscosity (Pa · s): 1.0
        # Fluid Viscosity (Pa · s): 2.0
        # Target Pressure (kPa): 100
        # Fluid Name: collagen

        Pressure (kPa)	Volume (mL)
        5	0,0008
        10	0,006433333333
        25	0,0164
        50	0,0281
        75	0,03696666667
        100	0,04343333333
        125	0,04916666667
        150	0,05486666667
        175	0,0605
        200	0,0676

Example Workflow

    Load Calibration Data: Click the "Load from File" button to load a calibration file.
    Input Parameters: Enter the known viscosity, target viscosity, target pressure, and fluid name.
    Calculate: Click the "Calculate" button to compute the dispensed volume for the target fluid.
    Visualize Results: The GUI displays the calibration curve, extrapolated curve, and target volume on a plot.

Code Structure

    dropcal.py: Contains the CalibrationApp class and auxiliary functions for handling calibration data and calculations.
    main.py: The entry point for running the application.

Dependencies

    numpy: For numerical operations.
    pandas: For reading and processing calibration data.
    matplotlib: For plotting calibration curves.
    scipy: For interpolation.
    tkinter: For the graphical user interface.

License
This project is open-source under GNU GPL v3.

by Maria Teresa Alameda Felgueiras.
