import os
import sys
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import argparse



parser = argparse.ArgumentParser(
                    prog='Force Plotter V1.0',
                    description='Plots the integral force development in your case',
                    epilog='Read the documentation for more information on how to use this')


args = parser.parse_args()

#getting case information
casePath = os.getcwd() #path of directory this script is being run in
caseName = casePath.split('/')[-1] #directory name
