"""
Diagnostic OHC
================

This module contains several python tools used for evaluation of
Oean Heat Content from model runs.

"""
import logging
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def folder(name):
    """
    Take a string or a list of strings, convert it to a directory style,
    then make the folder and the string.
    Returns folder string and final character is always os.sep. ('/')
    Arguments
    ---------
    name: list or string
        A list of nested directories, or a path to a directory.

    Returns
    ---------
    str
        Returns a string of a full (potentially new) path of the directory.
    """
    sep = os.sep
    if isinstance(name, list):
        name = os.sep.join(name)
    if name[-1] != sep:
        name = name + sep
    if os.path.exists(name) is False:
        os.makedirs(name)
        logger.info('Making new directory:\t%s', str(name))
    return name

def get_input_files(cfg, index=''):
    """
    Load input configuration file as a Dictionairy.


    Arguments
    ---------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    index: int
        the index of the file in the cfg file.

    Returns
    ---------
    dict
        A dictionairy of the input files and their linked details.
    """
    if isinstance(index, int):
        metadata_file = cfg['input_files'][index]
        with open(metadata_file) as input_file:
            metadata = yaml.safe_load(input_file)
        return metadata
    return _get_input_data_files(cfg)


def errinfo(Nlen,Aref,Bmod):
    """
    Calculate the error matrix for bias, rmse, and correlation between Aref and Bmod array

    return the bias, rmse, and the correlation for SEACLIM.

    """
    # data control
    if Nlen != len(Aref) or Nlen != len(Bmod):
       print("Check the array dimensions!")
       print("broken in diag_errinfo")
       quit()

    # Convert inputs to numpy arrays for vectorization
    A, B = np.array(Aref), np.array(Bmod)

    # Calculate the bias:
    bias=np.nanmean(B-A)

    # Calculate the squared differences
    bsqared=np.square(B-A)
    rmse=np.sqrt(np.nanmean(bsqared))
    
    #Calculate the correlation matrix
    #print("debug correlation ..")
    #print(B.shape)
    #print(A.shape)
    AA = A[~np.isnan(A)]
    BB = B[~np.isnan(B)]
    corr = np.corrcoef(BB, AA)
    # The correlation coefficient between the two variables is the off-diagonal value
    coef = corr[0, 1]

    return bias,rmse,coef


def write_onefld(filename,var,varname,varunit,varlong,varhist):
    import scipy.io.netcdf as sionet
    import netCDF4
    idm=var.shape[1]
    jdm=var.shape[0]
    nc=sionet.netcdf_file(filename,"w")
    nc.createDimension("x",idm)
    nc.createDimension("y",jdm)
    nc.createVariable(varname,"double",("y","x",))
    nc.variables[varname].units      =varunit
    nc.variables[varname].long_name  =varlong
    nc.variables[varname].history    =varhist
    nc.variables[varname][:]         =np.array(var)
    nc.close()

