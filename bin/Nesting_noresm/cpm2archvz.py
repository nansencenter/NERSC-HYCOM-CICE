#!/usr/bin/env python3
import argparse
import abfile.abfile as abf
import logging
import re
import os.path
from netCDF4 import Dataset, num2date
import numpy as np

# Set up logger
_loglevel = logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate = False

def maplev(a, lpp=1):
    """Gap-filling and smoothing routine via 2D Laplacian land-filler."""
    jm, im = a.shape
    J, I = np.where(~np.isnan(a))
    with np.errstate(invalid='ignore'):
        if len(I) > 0 and len(J) > 0:
            av = np.nansum(a[J, I]) / (len(I) * len(J))
        else:
            av = 0.0
            
    J, I = np.where(np.isnan(a))
    a[J, I] = av
    b = a
    
    cc = np.zeros(a.shape)
    for k in range(lpp):
        cc[1:-2, 1:-2] = b[1:-2, 1:-2] + .5/4 * (
            b[1:-2, 2:-1] + b[0:-3, 1:-2] + b[1:-2, 0:-3] + b[2:-1, 1:-2] - 4.*b[1:-2, 1:-2]
        )
        cc[:, 0] = cc[:, 1]
        cc[:, im-1] = cc[:, im-2]
        cc[jm-1, :] = cc[jm-2, :]
        b[J, I] = cc[J, I]

    a[J, I] = cc[J, I]
    J, I = np.where(np.isnan(a))
    a[J, I] = 0.
    return a

def u_to_hycom_u(field2d):
    """Shift Arakawa C-grid U points to HYCOM indexing."""
    return np.roll(field2d, 1, axis=1)

def v_to_hycom_v(field2d):
    """Shift Arakawa C-grid V points to HYCOM indexing."""
    myfield = np.copy(field2d)
    myfield[1:, :] = myfield[:-1, :]
    return myfield

def read_mesh(filemesh):
    """Reads longitude, latitude, and depth from the companion mesh files."""
    filepath = os.path.dirname(filemesh)
    filename = os.path.basename(filemesh)
    meshfile_depth = os.path.join(filepath, "deptho" + filename[9:])

    with Dataset(filemesh, "r") as ncida, Dataset(meshfile_depth, "r") as ncidd:
        depth = ncidd.variables["deptho"][:, :] 
        plon = ncida.variables["longitude"][:, :] 
        plat = ncida.variables["latitude"][:, :]  
    
    return depth, plon, plat

def main(filemesh, merged_nc_files, bio_file=None, iexpt=1, iversn=22, yrflag=3):
    depth, plon, plat = read_mesh(filemesh)
    ip = (depth == 0)  # Land mask: True where depth is 0

    for nc_file in merged_nc_files:
        if not os.path.exists(nc_file):
            logger.warning(f"File not found, skipping: {nc_file}")
            continue

        logger.info(f"Processing unified merged file: {os.path.basename(nc_file)}")
        
        with Dataset(nc_file, "r") as ncid:

            # Open optional biogeochemical variables if provided
            bio = Dataset(bio_file, "r") if bio_file is not None else None

            # Time parsing
            time = ncid.variables["time"][0]
            tunit = ncid.variables["time"].units
            t_cal = ncid.variables["time"].calendar if "calendar" in ncid.variables["time"].ncattrs() else "standard"
            date = num2date(time, units=tunit, calendar=t_cal)

            # Read 3D variables from single file
            u = np.squeeze(ncid.variables["uo"][:])
            v = np.squeeze(ncid.variables["vo"][:])
            t = np.squeeze(ncid.variables["thetao"][:])
            s = np.squeeze(ncid.variables["so"][:])

            if bio is not None:
                no3 = np.squeeze(bio.variables["no3"][:])
                po4 = np.squeeze(bio.variables["po4"][:])
                si  = np.squeeze(bio.variables["si"][:])
                o2  = np.squeeze(bio.variables["o2"][:])

                # Same conversions as esm2archvz.py
                no3 *= 6.625 * 12.01 * 1000.0
                si   *= 6.625 * 12.01 * 1000.0
                po4 *= 106.0 * 12.01 * 1000.0
                o2  *= 1000.0

            # Read Pre-Calculated 2D Barotropic currents and SSH
            ubaro_raw = np.squeeze(ncid.variables["ubaro_netcdf"][:])
            vbaro_raw = np.squeeze(ncid.variables["vbaro_netcdf"][:])
            ssh = np.squeeze(ncid.variables["zos"][:])

            # Vertical grid layout
            lev_bnds = ncid.variables["depth_bnds"][:, :]
            z_levels = ncid.variables["depth"][:]
            nlev = np.size(z_levels)
            dz = lev_bnds[:, 1] - lev_bnds[:, 0]

            if bio is not None:
                bio.close()

        # Clean fill values
        u = np.where(np.abs(u) < 1e10, u, 0.)
        v = np.where(np.abs(v) < 1e10, v, 0.)
        ubaro_raw = np.where(np.abs(ubaro_raw) < 1e10, ubaro_raw, 0.)
        vbaro_raw = np.where(np.abs(vbaro_raw) < 1e10, vbaro_raw, 0.)
        ssh = np.where(np.abs(ssh) < 1e10, 0., ssh * 9.81) # Convert to geopotential
        montg1 = np.zeros(ssh.shape)

        # Apply Arakawa C-grid shifts to align velocities with HYCOM
        uu_x = np.zeros(u.shape)
        vv_y = np.zeros(v.shape)
        for k in range(nlev):
            uu_x[k, :, :] = u_to_hycom_u(u[k, :, :])  
            vv_y[k, :, :] = v_to_hycom_v(v[k, :, :])  

        ubaro = u_to_hycom_u(ubaro_raw)
        vbaro = v_to_hycom_v(vbaro_raw)
        del u, v

        # Open intermediate file writer
        oname = date.strftime("archv.%Y_%j_00")
        header1 = f"Converted CPM Physical Data\n"
        header2 = "Archive files tailored for nesting interpolation\n"
        header3 = "Physical parameters only\n"
        
        out_path = os.path.join("./data", oname)
        outfile = abf.ABFileArchv(out_path, "w", iexpt=iexpt, iversn=iversn, yrflag=yrflag, 
                                  cline1=header1, cline2=header2, cline3=header3)

        logger.info("Writing 2D parameters...")
        outfile.write_field(montg1, ip, "montg1", 0, time, 1, 0)
        outfile.write_field(ssh, ip, "srfhgt", 0, time, 0, 0)
        outfile.write_field(np.zeros(ssh.shape), ip, "surflx", 0, time, 0, 0) 
        outfile.write_field(np.zeros(ssh.shape), ip, "salflx", 0, time, 0, 0) 
        outfile.write_field(np.zeros(ssh.shape), ip, "bl_dpth", 0, time, 0, 0) 
        outfile.write_field(np.zeros(ssh.shape), ip, "mix_dpth", 0, time, 0, 0) 
        outfile.write_field(ubaro, ip, "u_btrop", 0, time, 0, 0) 
        outfile.write_field(vbaro, ip, "v_btrop", 0, time, 0, 0) 

        logger.info("Processing and writing 3D fields...")
        for k in range(nlev):
            ul = np.squeeze(uu_x[k, :, :]) - ubaro 
            vl = np.squeeze(vv_y[k, :, :]) - vbaro 

            dzl = np.zeros(ul.shape)
            if k < nlev - 1:
                active_mask = (depth > lev_bnds[k, 0])
                dzl = np.where(active_mask, dz[k], dzl)
                bottom_cell_mask = (depth <= lev_bnds[k, 1]) & (depth > lev_bnds[k, 0])
                dzl = np.where(bottom_cell_mask, depth - lev_bnds[k, 0], dzl)
            else:
                bottom_cell_mask = (depth > lev_bnds[k, 0])
                dzl = np.where(bottom_cell_mask, depth - lev_bnds[k, 0], dzl)
            
            sl = np.squeeze(s[k, :, :])
            tl = np.squeeze(t[k, :, :])
            
            sl = np.where(sl < 1e2, sl, np.nan)
            sl = np.minimum(np.maximum(maplev(sl), 25.0), 80.0)
            tl = np.where(tl <= 5e2, tl, np.nan)
            tl = np.minimum(np.maximum(maplev(tl), -5.0), 50.0)

            if bio_file is not None:
                no3l = np.squeeze(no3[k])
                po4l = np.squeeze(po4[k])
                sil  = np.squeeze(si[k])
                o2l  = np.squeeze(o2[k])

                no3l = np.where(no3l < 1e8, no3l, np.nan)
                no3l = np.minimum(np.maximum(maplev(no3l), 0.0), 1.0e8)

                po4l = np.where(po4l < 1e8, po4l, np.nan)
                po4l = np.minimum(np.maximum(maplev(po4l), 0.0), 1.0e8)

                sil = np.where(sil < 1e8, sil, np.nan)
                sil = np.minimum(np.maximum(maplev(sil), 0.0), 1.0e8)

                o2l = np.where(o2l < 1e8, o2l, np.nan)
                o2l = np.minimum(np.maximum(maplev(o2l), 0.0), 1.0e8)

            if k > 0:
                empty_layer_mask = (dzl < 1e-4)
                tl[empty_layer_mask] = tl_above[empty_layer_mask]
                sl[empty_layer_mask] = sl_above[empty_layer_mask]

                if bio_file is not None:
                    no3l[empty_layer_mask] = no3_above[empty_layer_mask]
                    po4l[empty_layer_mask] = po4_above[empty_layer_mask]
                    sil[empty_layer_mask]  = si_above[empty_layer_mask]
                    o2l[empty_layer_mask]  = o2_above[empty_layer_mask]
            
            tl_above = np.copy(tl)
            sl_above = np.copy(sl)

            if bio_file is not None:
                no3_above = np.copy(no3l)
                po4_above = np.copy(po4l)
                si_above  = np.copy(sil)
                o2_above  = np.copy(o2l)

            onem = 9806.0
            outfile.write_field(ul, ip, "u-vel.", 0, time, k+1, 0) 
            outfile.write_field(vl, ip, "v-vel.", 0, time, k+1, 0) 
            outfile.write_field(dzl * onem, ip, "thknss", 0, time, k+1, 0)
            outfile.write_field(tl, ip, "temp", 0, time, k+1, 0)
            outfile.write_field(sl, ip, "salin", 0, time, k+1, 0)

            if bio_file is not None:
                outfile.write_field(no3l, ip, "ECO_no3", 0, time, k+1, 0)
                outfile.write_field(po4l, ip, "ECO_pho", 0, time, k+1, 0)
                outfile.write_field(sil,  ip, "ECO_sil", 0, time, k+1, 0)
                outfile.write_field(o2l,  ip, "ECO_oxy", 0, time, k+1, 0)

        outfile.close()
        logger.info(f"Finalized: {oname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert preprocessed unified NorCPM NetCDF files to HYCOM format.')
    parser.add_argument('meshfile', type=str, help="Mesh footprint file path")
    parser.add_argument('merged_nc_files', type=str, nargs="+", help="Target .merged_YYYY-MM.nc files")
    parser.add_argument('--bio-file', type=str, default=None, help="Optional NetCDF file containing biogeochemical variables.")
    parser.add_argument('--iexpt', type=int, default=1)
    parser.add_argument('--iversn', type=int, default=22)
    parser.add_argument('--yrflag', type=int, default=3)

    args = parser.parse_args()
    main(args.meshfile, args.merged_nc_files, bio_file=args.bio_file, iexpt=args.iexpt, iversn=args.iversn, yrflag=args.yrflag)
