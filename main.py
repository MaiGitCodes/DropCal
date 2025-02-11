#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:00:13 2025

@author: Maria Teresa Alameda

This program is designed to calibrate fluid dispensing systems under
various conditions. It uses a known calibration volume (mL) - pressure (kPa)
curve (for example water) to extrapolate the volume dispensed of other liquid
with a known viscosity (for example collagen) at a certain pressure.

"""
from tkinter import Tk
from dropcal import CalibrationApp

# Run the application
if __name__ == "__main__":
    root = Tk()
    app = CalibrationApp(root)
    root.mainloop()
