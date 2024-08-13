import numpy as np
import xarray as xr
import metpy
from metpy.units import units
from metpy.constants import water_heat_vaporization, dry_air_gas_constant, earth_gravity
from metpy.interpolate import interpolate_1d
import glob


def extract_everything(datestring, dx=2.5, hrs_to_shift=0, hrs_to_simulate=24):
    h = 30.0  # surface elevation of ENA site in m
    lat = 39.0916  # deg N
    lon = -28.0257  # deg E
    forc_dir = "/ccsopen/home/h1x/scratch/ena_forcing_check/forcing"
    sonde_dir = "/ccsopen/home/h1x/scratch/ena_forcing_check/obs/enasondewnpnC1"

    dz = 25.0  # m
    # ztop = 6500.0
    ztop = 10750.0
    ztopo = ztop + 300.0
    zout = np.arange(dz * 0.5, ztop + dz * 0.6, dz)
    dz_sndout = 10 # m
    z_snd_out2 = np.insert(np.arange(dz_sndout * 0.5, ztop + dz_sndout * 0.6, dz_sndout), 0, 1.)

    # all the dates!
    casedate = np.datetime64(datestring)
    datem1 = (casedate - np.timedelta64(1, "D")).astype(object).strftime("%Y%m%d")
    date0 = casedate.astype(object).strftime("%Y%m%d")
    datep1 = (casedate + np.timedelta64(1, "D")).astype(object).strftime("%Y%m%d")
    ts = casedate + np.timedelta64(hrs_to_shift, "h")  # time to extract sounding
    ds = casedate.astype(object).timetuple().tm_yday + hrs_to_shift / 24.0
    # start/end time to extract forcing in hours before/after ts
    h1, h2 = -6, hrs_to_simulate + 6
    # start/end time to extract forcing in hours before/after ts
    h1_snd, h2_snd = 0, hrs_to_simulate + 12
    # start/end time to extract forcing
    t1 = ts + np.timedelta64(h1, "h")
    t2 = ts + np.timedelta64(h2, "h")
    # start/end time to extract sounding
    t1_snd = ts + np.timedelta64(h1_snd, "h")
    t2_snd = ts + np.timedelta64(h2_snd, "h")
    d1 = casedate.astype(object).timetuple().tm_yday + (hrs_to_shift + h1) / 24.0
    d2 = casedate.astype(object).timetuple().tm_yday + (hrs_to_shift + h2) / 24.0
    d1_snd = (
        casedate.astype(object).timetuple().tm_yday + (hrs_to_shift + h1_snd) / 24.0
    )
    d2_snd = (
        casedate.astype(object).timetuple().tm_yday + (hrs_to_shift + h2_snd) / 24.0
    )

    print(f"Simulation starts (init. sounding) at {ts}, doy={ds}")
    print(f"Forcing extraction starts at {t1}, doy={d1}")
    print(f"Forcing extraction ends at {t2} doy={d2}")
    print(f"Sounding extraction starts at {t1_snd}, doy={d1_snd}")
    print(f"Sounding extraction ends at {t2_snd} doy={d2_snd}")

    # extract all the available ARM soundings for nudging
    # searching for soundings from t1_snd to t2_snd
    # with a given time interval (12 h for now)
    snd_times = []
    snd_doys = []
    snd_vars = []
    for isnd, snd_time in enumerate(np.arange(t1_snd, t2_snd, np.timedelta64(12, "h"))):
        snd_doy = d1_snd + 12 / 24.0 * isnd
        snd_time_adj = snd_time - np.timedelta64(1, "h")
        snd_datestring = snd_time_adj.astype(object).strftime("%Y%m%d")
        snd_hour = snd_time_adj.astype(object).hour
        snd_files = glob.glob(
            f"{sonde_dir}/enasondewnpnC1.b1.{snd_datestring}.{snd_hour:02d}*.cdf"
        )
        if len(snd_files) != 1:
            print(f"Cannot find the sounding file for {snd_time}.")
        else:
            print(f"ARM sounding file found: {snd_files[0]} for {snd_time}")
            snd_times.append(snd_time)
            snd_doys.append(snd_doy)
            snd_vars.append(extract_sounding_sonde(snd_files[0]))

    # output sounding
    with open(f"snd.{date0}", "w") as snd:

        snd.write("z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]\n")

        for snd_var, snd_doy in zip(snd_vars, snd_doys):
            ps_snd, z_snd, tp_snd, rv_snd, u_snd, v_snd = snd_var
            tp_snd_out = np.interp(z_snd_out2, z_snd-h, tp_snd.magnitude)
            rv_snd_out = np.interp(z_snd_out2, z_snd-h, rv_snd.magnitude)
            rv_snd_out = rv_snd_out * 1.0e3  # convert to g/kg
            u_snd_out = np.interp(z_snd_out2, z_snd-h, u_snd)
            v_snd_out = np.interp(z_snd_out2, z_snd-h, v_snd)
            nz = z_snd_out2.shape[0]
            snd.write(f"{snd_doy},  {nz},  {ps_snd:8.3f}  day,levels,pres0\n")
            for z, t, r, u, v in zip(
                z_snd_out2,
                tp_snd_out,
                rv_snd_out,
                u_snd_out,
                v_snd_out,
            ):
                snd.write(
                    f"{z:8.2f}  -999.9  {t:8.4f}  {r:8.4f}  {u:9.4f}  {v:9.4f} \n"
                )

    # extract ERA5 surface fluxes and large-scale forcings
    era5_sfc = xr.open_dataset(f"{forc_dir}/era5/data/ERA5-{datem1}-{datep1}-sfc.nc")
    era5_l = xr.open_dataset(f"{forc_dir}/era5/data/ERA5-{datem1}-{datep1}-ml.nc")
    print("Extracting ERA5 sfc data ...")
    era5_sfc = era5_sfc.sortby(era5_sfc.time)
    toy_era5_sfc, sst_era5, shf_era5, lhf_era5, _ = extract_sfc_fluxes_era5(
        era5_sfc, t1, t2, lat, lon, dx
    )
    print("Extracting ERA5 vertical levels data ...")
    toy_era5, zs, ps, z, p, w, omega, u, v, ta, ma, ug, vg = extract_forcing_era5(
        era5_l, t1, t2, lat, lon, dx
    )
    (
        w_era5,
        u_era5,
        v_era5,
        ta_era5,
        ma_era5,
        ug_era5,
        vg_era5,
    ) = interpolate_forcing_era5(zout, z, w, u, v, ta, ma, ug, vg, zs)
    ps_era5 = ps.mean(axis=(1, 2))
    del zs, ps, z, p, w, omega, u, v, ta, ma, ug, vg
    del era5_l, era5_sfc

    # output sounding with surface pressure from forcing dataset
    _, z_snd, tp_snd, rv_snd, u_snd, v_snd = snd_vars[0]
    nt_snd = np.nonzero(toy_era5 == ds)[0][0]
    print(f"Writing out sounding data with ERA5 P_surf at {toy_era5[nt_snd]}...")
    lim = z_snd < ztopo
    z_snd_out = np.compress(lim, z_snd)
    tp_snd_out = np.compress(lim, tp_snd)
    rv_snd_out = np.compress(lim, rv_snd)
    rv_snd_out = rv_snd_out * 1.0e3  # convert to g/kg
    u_snd_out = np.compress(lim, u_snd)
    v_snd_out = np.compress(lim, v_snd)
    nz = z_snd_out.shape[0]
    with open(f"snd_era5ps.{date0}", "w") as snd:
        snd.write("z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]\n")
        snd.write(
            f"{ds},  {nz},  {ps_era5.values[nt_snd]/100.0:8.3f}  day,levels,pres0\n"
        )
        for z, t, r, u, v in zip(
            z_snd_out, tp_snd_out.magnitude, rv_snd_out.magnitude, u_snd_out, v_snd_out
        ):
            snd.write(
                f"{z - h:8.2f}  -999.9  {t:8.4f}  {r:8.4f}  {u:9.4f}  {v:9.4f} \n"
            )
        snd.write(
            f"{d2},  {nz},  {ps_era5.values[nt_snd]/100.0:8.3f}  day,levels,pres0\n"
        )
        for z, t, r, u, v in zip(
            z_snd_out, tp_snd_out.magnitude, rv_snd_out.magnitude, u_snd_out, v_snd_out
        ):
            snd.write(
                f"{z - h:8.2f}  -999.9  {t:8.4f}  {r:8.4f}  {u:9.4f}  {v:9.4f} \n"
            )

    # output surface fluxes
    print(f"Writing out ERA5 sfc data from {toy_era5_sfc[0]} to {toy_era5_sfc[-1]}...")
    with open(f"sfc_era5.{date0}", "w") as sfc:
        sfc.write("day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2) \n")
        for nt in range(len(toy_era5_sfc)):
            sfc.write(
                f"{toy_era5_sfc[nt]:13.9f}"
                f"  {sst_era5[nt]:8.4f}"
                f"  {-shf_era5[nt]:9.4f}"
                f"  {-lhf_era5.values[nt]:9.4f}"
                "  -999.9 \n"
            )

    # output large-scale forcing
    print(f"Writing out ERA5 lsf data from {toy_era5[0]} to {toy_era5[-1]}...")
    nz = zout.shape[0]
    with open(f"lsf_era5.{date0}", "w") as lsf:
        lsf.write(" z[m] p[mb]  tls[K/s] qls[kg/kg/s] uls  vls  wls[m/s] \n")
        for nt in range(len(toy_era5)):
            lsf.write(
                f"{toy_era5[nt]:13.9f},  {nz},  {ps_era5.values[nt]/100.0:8.3f}  day,levels,pres0 \n"
            )
            for zl in range(nz):
                lsf.write(
                    f"{zout[zl]:9.2f}  -999.9"
                    f"  {ta_era5[nt, zl]:12.4e}"
                    f"  {ma_era5[nt, zl]:12.4e}"
                    f"  {ug_era5[nt, zl]:9.4f}"
                    f"  {vg_era5[nt, zl]:9.4f}"
                    f"  {w_era5[nt, zl]:12.4e} \n"
                )

    # extract MERRA2 surface fluxes and large-scale forcings
    merra2_l = xr.open_mfdataset(
        [
            f"{forc_dir}/merra2/data/MERRA2_400.inst3_3d_asm_Nv.{datem1}.nc4",
            f"{forc_dir}/merra2/data/MERRA2_400.inst3_3d_asm_Nv.{date0}.nc4",
            f"{forc_dir}/merra2/data/MERRA2_400.inst3_3d_asm_Nv.{datep1}.nc4",
        ]
    )
    merra2_sfc = xr.open_mfdataset(
        [
            f"{forc_dir}/merra2/data/MERRA2_400.tavg1_2d_flx_Nx.{datem1}.nc4",
            f"{forc_dir}/merra2/data/MERRA2_400.tavg1_2d_flx_Nx.{date0}.nc4",
            f"{forc_dir}/merra2/data/MERRA2_400.tavg1_2d_flx_Nx.{datep1}.nc4",
        ]
    )
    merra2_sfc2 = xr.open_mfdataset(
        [
            f"{forc_dir}/merra2/data/MERRA2_400.tavg1_2d_slv_Nx.{datem1}.nc4",
            f"{forc_dir}/merra2/data/MERRA2_400.tavg1_2d_slv_Nx.{date0}.nc4",
            f"{forc_dir}/merra2/data/MERRA2_400.tavg1_2d_slv_Nx.{datep1}.nc4",
        ]
    )
    print("Extracting MERRA2 sfc data ...")
    toy_merra2_sfc, sst_merra2, shf_merra2, lhf_merra2 = extract_sfc_fluxes_merra2(
        merra2_sfc, merra2_sfc2, t1, t2, lat, lon, dx
    )
    print("Extracting MERRA2 vertical levels data ...")
    toy_merra2, zs, ps, z, p, w, omega, u, v, ta, ma, ug, vg = extract_forcing_merra2(
        merra2_l, t1, t2, lat, lon, dx
    )
    (
        w_merra2,
        u_merra2,
        v_merra2,
        ta_merra2,
        ma_merra2,
        ug_merra2,
        vg_merra2,
    ) = interpolate_forcing_merra2(zout, z, w, u, v, ta, ma, ug, vg, zs)
    ps_merra2 = ps.mean(axis=(1, 2))
    del zs, ps, z, p, w, omega, u, v, ta, ma, ug, vg
    del merra2_l, merra2_sfc, merra2_sfc2

    # output sounding with surface fluxes from forcing datasets
    _, z_snd, tp_snd, rv_snd, u_snd, v_snd = snd_vars[0]
    nt_snd = np.nonzero(toy_merra2 == ds)[0][0]
    print(f"Writing out sounding data with MERRA2 P_surf at {toy_merra2[nt_snd]}...")
    lim = z_snd < ztopo
    z_snd_out = np.compress(lim, z_snd)
    tp_snd_out = np.compress(lim, tp_snd)
    rv_snd_out = np.compress(lim, rv_snd)
    rv_snd_out = rv_snd_out * 1.0e3  # convert to g/kg
    u_snd_out = np.compress(lim, u_snd)
    v_snd_out = np.compress(lim, v_snd)
    nz = z_snd_out.shape[0]
    with open(f"snd_merra2ps.{date0}", "w") as snd:
        snd.write("z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]\n")
        snd.write(
            f"{ds},  {nz},  {ps_merra2.values[nt_snd]/100.0:8.3f}  day,levels,pres0\n"
        )
        for z, t, r, u, v in zip(
            z_snd_out, tp_snd_out.magnitude, rv_snd_out.magnitude, u_snd_out, v_snd_out
        ):
            snd.write(
                f"{z - h:8.2f}  -999.9  {t:8.4f}  {r:8.4f}  {u:9.4f}  {v:9.4f} \n"
            )
        snd.write(
            f"{d2},  {nz},  {ps_merra2.values[nt_snd]/100.0:8.3f}  day,levels,pres0\n"
        )
        for z, t, r, u, v in zip(
            z_snd_out, tp_snd_out.magnitude, rv_snd_out.magnitude, u_snd_out, v_snd_out
        ):
            snd.write(
                f"{z - h:8.2f}  -999.9  {t:8.4f}  {r:8.4f}  {u:9.4f}  {v:9.4f} \n"
            )

    # output surface fluxes

    print(
        f"Writing out MERRA2 sfc data from {toy_merra2_sfc[0]} to {toy_merra2_sfc[-1]}..."
    )
    with open(f"sfc_merra2.{date0}", "w") as sfc:
        sfc.write("day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2) \n")
        for nt in range(len(toy_merra2_sfc)):
            sfc.write(
                f"{toy_merra2_sfc[nt]:13.9f}"
                f"  {sst_merra2.values[nt]:8.4f}"
                f"  {shf_merra2.values[nt]:9.4f}"
                f"  {lhf_merra2.values[nt]:9.4f}"
                "  -999.9 \n"
            )

    # output large-scale forcings
    print(f"Writing out MERRA2 lsf data from {toy_merra2[0]} to {toy_merra2[-1]}...")
    nz = zout.shape[0]
    with open(f"lsf_merra2.{date0}", "w") as lsf:
        lsf.write(" z[m] p[mb]  tls[K/s] qls[kg/kg/s] uls  vls  wls[m/s] \n")
        for nt in range(len(toy_merra2)):
            lsf.write(
                f"{toy_merra2[nt]:13.9f},  {nz},  {ps_merra2.values[nt]/100.0:8.3f}  day,levels,pres0 \n"
            )
            for zl in range(nz):
                lsf.write(
                    f"{zout[zl]:9.2f}  -999.9"
                    f"  {ta_merra2[nt, zl]:12.4e}"
                    f"  {ma_merra2[nt, zl]:12.4e}"
                    f"  {ug_merra2[nt, zl]:9.4f}"
                    f"  {vg_merra2[nt, zl]:9.4f}"
                    f"  {w_merra2[nt, zl]:12.4e} \n"
                )

    return


def extract_sounding_sonde(sonde_file):
    sonde = xr.open_dataset(sonde_file)
    ps = sonde.pres.values[0]
    tp = metpy.calc.potential_temperature(
        sonde.pres.metpy.quantify(),
        sonde.tdry.metpy.magnitude * units("degC"),
    )
    rv = metpy.calc.mixing_ratio_from_relative_humidity(
        sonde.pres.metpy.quantify(),
        sonde.tdry.metpy.magnitude * units("degC"),
        sonde.rh.metpy.quantify(),
    )
    z = sonde.alt
    u = sonde.u_wind
    v = sonde.v_wind
    return [ps, z, tp, rv, u, v]


def extract_sounding_sonde_pin(sonde_file):
    sonde = xr.open_dataset(sonde_file)
    ps = sonde.pres.values[0] * 100.
    t  = sonde.tdry + 273.15 # deg C to K
    p = sonde.pres * 100. 
    rv = metpy.calc.mixing_ratio_from_relative_humidity(
        sonde.pres.metpy.quantify(),
        sonde.tdry.metpy.magnitude * units("degC"),
        sonde.rh.metpy.quantify(),
    )
    z = sonde.alt
    u = sonde.u_wind
    v = sonde.v_wind
    return [ps, z, p, t, rv, u, v]


def extract_sfc_fluxes_merra2(d1, d2, t1, t2, lat, lon, dx):
    shf = d1.HFLUX.loc[t1:t2, lat - dx : lat + dx, lon - dx : lon + dx].mean(
        axis=(1, 2)
    )
    lhf = d1.EFLUX.loc[t1:t2, lat - dx : lat + dx, lon - dx : lon + dx].mean(
        axis=(1, 2)
    )
    sst = d2.TS.loc[t1:t2, lat - dx : lat + dx, lon - dx : lon + dx].mean(axis=(1, 2))
    print(
        "Latitude points in the averageing box: "
        + repr(d1.lat.loc[lat - dx : lat + dx].values.tolist())
    )
    print(
        "Longitude points in the averageing box: "
        + repr(d1.lon.loc[lon - dx : lon + dx].values.tolist())
    )
    # Strange usage of np.datetime64 to get fractional day of the year
    toy = [
        np.timedelta64(
            t - np.datetime64(d1.time.sel(time=t1, method="nearest").values, "Y"), "m"
        )
        / np.timedelta64(1, "D")
        + 1.0
        for t in d1.time.loc[t1:t2].values
    ]
    return np.asarray(toy), sst, shf, lhf


def extract_sfc_fluxes_era5(d, t1, t2, lat, lon, dx):
    g = earth_gravity
    shf = d.ishf.loc[t1:t2, lat + dx : lat - dx, lon - dx : lon + dx].mean(axis=(1, 2))
    lhf = (
        (water_heat_vaporization * d.ie.metpy.quantify())
        .loc[t1:t2, lat + dx : lat - dx, lon - dx : lon + dx]
        .mean(axis=(1, 2))
    )
    sst = d.sst.loc[t1:t2, lat + dx : lat - dx, lon - dx : lon + dx].mean(axis=(1, 2))
    lsm = d.lsm.loc[t1, lat + dx : lat - dx, lon - dx : lon + dx].max()
    print(
        "Latitude points in the averageing box: "
        + repr(d.latitude.loc[lat + dx : lat - dx].values.tolist())
    )
    print(
        "Longitude points in the averageing box: "
        + repr(d.longitude.loc[lon - dx : lon + dx].values.tolist())
    )
    # Strange usage of np.datetime64 to get fractional day of the year
    toy = [
        np.timedelta64(
            t - np.datetime64(d.time.sel(time=t1, method="nearest").values, "Y"), "m"
        )
        / np.timedelta64(1, "D")
        + 1.0
        for t in d.time.loc[t1:t2].values
    ]
    zs = (
        d.z.loc[t1:t2, lat + dx : lat - dx, lon - dx : lon + dx]
        .metpy.quantify()
        .mean(axis=(1, 2))
    )
    # print(d.z.loc[t1, lat+dx:lat-dx, lon-dx:lon+dx].values)
    if lsm > 0.5:
        # According to ECMWF documentation, lsm > 0.5 is considered ocean/water.
        # But actually none of the points around ENA has lsm = 0.0 either.
        print("Averaging area contains over-land grid points!")
    return np.asarray(toy), sst, shf, lhf, zs.metpy.quantify() / g


def extract_forcing_era5(d, t1, t2, lat, lon, dx_in):
    """
    Extract large-scale subsidence, horizontal winds, geostrophic horizontal winds, large-scale horizontal
    advection tendencies of temperature and moisture from the forcing dataset.
    """
    # use a larger box to calculate things
    dx = dx_in * 5.0

    """
    interpolation to pressure levels and calculate geopotential heights following:
    https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    """
    rd = dry_air_gas_constant
    g = earth_gravity
    ab = np.genfromtxt(
        "era5_table.csv", delimiter=",", skip_header=1, missing_values="-"
    )
    a = ab[:, 1]
    b = ab[:, 2]

    t = d.t.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    qv = d.q.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    u = d.u.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    v = d.v.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    omega = d.w.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    # note here that the second index (1) here means getting the first level data, 
    # not the actual array index.Also, only first level data is valid
    lnsp = d.lnsp.loc[t1:t2, 1, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    # the surface geopotential looks so noisy because of the spectral decomposition/representation used in IFS
    sgp = d.z.loc[t1:t2, 1, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx]
    rv = metpy.calc.mixing_ratio_from_specific_humidity(qv)

    nt, nz, ny, nx = t.shape

    pi = np.zeros((nt, nz + 1, ny, nx))
    ps = np.exp(lnsp)
    pi[:] = ps.values[:, np.newaxis, :, :]
    pi = (
        a[np.newaxis, :, np.newaxis, np.newaxis]
        + pi * b[np.newaxis, :, np.newaxis, np.newaxis]
    )
    p = (pi[:, 1:, :, :] + pi[:, :-1, :, :]) * 0.5
    pi[:, 0, :, :] = 0.1
    dpi = pi[:, 1:, :, :] - pi[:, :-1, :, :]
    dlnpi = np.log(pi[:, 1:, :, :] / pi[:, :-1, :, :])

    # have not got time to derive alpha, this is just what is given in the ERA5 documentation
    alpha = 1.0 - dlnpi * pi[:, :-1, :, :] / dpi
    alpha[:, 0, :, :] = np.log(2.0)

    tm = t.metpy.quantify() * (1.0 + 0.609133 * rv.metpy.quantify())
    # tm = t*(1.0+0.609133*qv)
    dphi = rd.magnitude * tm.values * dlnpi
    phi = np.zeros((nt, nz + 1, ny, nx))
    phi[:, :-1, :, :] = np.flip(np.cumsum(dphi[:, ::-1, :, :], axis=1), 1)
    phi[:] = phi[:] + sgp.values[:, np.newaxis, :, :]
    ph = phi[:, 1:, :, :] + rd.magnitude * tm.values * alpha

    ph = t.copy(deep=True, data=ph)
    ph.attrs["units"] = "m**2/s**2"
    del ph.attrs["long_name"]
    del ph.attrs["standard_name"]
    ph.metpy.quantify()

    p = t.copy(deep=True, data=p)
    p.attrs["units"] = "Pa"
    del p.attrs["long_name"]
    del p.attrs["standard_name"]
    p.metpy.quantify()

    z = ph.metpy.quantify() / g
    z = t.copy(deep=True, data=z)
    z.attrs["units"] = "m"
    del z.attrs["long_name"]
    del z.attrs["standard_name"]

    w = -omega.metpy.quantify() * rd * tm.metpy.quantify() / p.metpy.quantify() / g
    # print(omega.metpy.units)
    # print(tm.metpy.units)
    # print(rd.units)
    # print(p.metpy.units)
    # print(w.metpy.units)

    ph = ph.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    )
    p = p.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    )
    z = z.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    )

    # calculate geostrophic winds
    lnp = np.log(p)
    lnp.attrs["units"] = ""
    lnp.metpy.quantify()
    dpx1, dpy1 = metpy.calc.geospatial_gradient(ph)
    dpx2, dpy2 = metpy.calc.geospatial_gradient(lnp)
    dpx = dpx1 + rd * tm.metpy.quantify() * dpx2
    dpy = dpy1 + rd * tm.metpy.quantify() * dpy2
    f = metpy.calc.coriolis_parameter(lat * units("degrees"))
    vg = dpx / f
    ug = -dpy / f

    # calculate horizontal advection tendencies
    # tp = calc_tp_era5(d)
    tp = metpy.calc.potential_temperature(p, t)
    ma = metpy.calc.advection(rv, u, v)
    ta = metpy.calc.advection(tp, u, v)

    # Strange usage of np.datetime64 to get fractional day of the year
    toy = [
        np.timedelta64(
            t - np.datetime64(d.time.sel(time=t1, method="nearest").values, "Y"), "m"
        )
        / np.timedelta64(1, "D")
        + 1.0
        for t in d.time.loc[t1:t2].values
    ]

    # output just points within the specified box
    dx = dx_in
    print(
        "Latitude points in the averageing box: "
        + repr(d.latitude.loc[lat + dx : lat - dx].values.tolist())
    )
    print(
        "Longitude points in the averageing box: "
        + repr(d.longitude.loc[lon + 360.0 - dx : lon + 360.0 + dx].values.tolist())
    )
    return (
        np.asarray(toy),
        d.z.loc[
            t1:t2, 1, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx
        ].metpy.quantify()
        / g,
        ps.loc[t1:t2, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        z.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        p.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        w.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        d.w.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        d.u.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        d.v.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        ta.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        ma.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        ug.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
        vg.loc[t1:t2, :, lat + dx : lat - dx, lon + 360 - dx : lon + 360 + dx],
    )


def extract_forcing_merra2(d, t1, t2, lat, lon, dx_in):
    rd = dry_air_gas_constant
    g = earth_gravity

    dx = 5.0 * dx_in

    p = d.PL.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx].metpy.quantify()
    z = d.H.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx].metpy.quantify()
    t = d.T.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx].metpy.quantify()
    qv = d.QV.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx].metpy.quantify()
    rv = metpy.calc.mixing_ratio_from_specific_humidity(qv)
    u = d.U.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx].metpy.quantify()
    v = d.V.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx].metpy.quantify()
    omega = d.OMEGA.loc[
        t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx
    ].metpy.quantify()
    tm = t.metpy.quantify() * (1.0 + 0.609133 * rv.metpy.quantify())
    p = p.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    )
    z = z.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    )

    # calculate geostrophic winds
    lnp = np.log(p.metpy.dequantify())
    lnp.attrs["units"] = ""
    dpx1, dpy1 = metpy.calc.geospatial_gradient(z * g)
    dpx2, dpy2 = metpy.calc.geospatial_gradient(lnp)
    dpx = dpx1 + rd * tm.metpy.quantify() * dpx2
    dpy = dpy1 + rd * tm.metpy.quantify() * dpy2
    f = metpy.calc.coriolis_parameter(lat * units("degrees"))
    vg = dpx / f
    ug = -dpy / f

    w = -omega.metpy.quantify() * rd * tm.metpy.quantify() / p.metpy.quantify() / g
    # print(omega.metpy.units)
    # print(tm.metpy.units)
    # print(rd.units)
    # print(p.metpy.units)
    # print(w.metpy.units)

    # calculate horizontal advection tendencies
    tp = metpy.calc.potential_temperature(p, t)
    ma = metpy.calc.advection(rv, u, v)
    ta = metpy.calc.advection(tp, u, v)

    # Strange usage of np.datetime64 to get fractional day of the year
    toy = [
        np.timedelta64(
            t - np.datetime64(d.time.sel(time=t1, method="nearest").values, "Y"), "m"
        )
        / np.timedelta64(1, "D")
        + 1.0
        for t in d.time.loc[t1:t2].values
    ]

    dx = dx_in
    print(
        "Latitude points in the averageing box: "
        + repr(d.lat.loc[lat - dx : lat + dx].values.tolist())
    )
    print(
        "Longitude points in the averageing box: "
        + repr(d.lon.loc[lon - dx : lon + dx].values.tolist())
    )
    return (
        np.asarray(toy),
        (
            d.PHIS.loc[t1:t2, lat - dx : lat + dx, lon - dx : lon + dx]
            * units("m**2/s**2")
            / g
        ).metpy.dequantify(),
        d.PS.loc[t1:t2, lat - dx : lat + dx, lon - dx : lon + dx],
        z.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        p.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        w.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        d.OMEGA.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        u.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        v.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        ta.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        ma.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        ug.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
        vg.loc[t1:t2, :, lat - dx : lat + dx, lon - dx : lon + dx],
    )


def interpolate_forcing_era5(zout, z, w, u, v, ta, ma, ug, vg, zs):
    """
    interpolate all the forcing profiles to SAM's vertical grid given by `zout`
    """

    nt, nz, ny, nx = z.shape
    zi = np.zeros((nt, nz + 1, ny, nx))
    """
    This is just one way to take into account surface elevation.
    The other way is to interpolate using z (heights above MSL) directly,
    then add a mean surface elevation to zout. Using the latter method, the averaging
    will be done on heights above MSL rather than heights above surface.
    But given the small surface elevations we have around ENA, it shouldn't make much difference. 
    """
    zi[:, :-1, :, :] = z.values - zs.values[:, np.newaxis, :, :]
    # zi[:, :-1, :, :] = z.values
    # zout = zout + zs.values.mean()
    out = []
    sfc = [0, 0, 0, 1, 1, 1, 1]
    for v, s in zip([w, u, v, ta, ma, ug, vg], sfc):
        vi = np.zeros_like(zi)
        vi[:, :-1, :, :] = v.values
        # surface winds are all zero
        if s == 1:
            vi[:, -1, :, :] = v.values[:, -1, :, :]
        # vi[:, -1, :, :] = v.values[:, -1, :, :]
        vo = interpolate_1d(zout, zi, vi, axis=1)
        out.append(vo.mean(axis=(2, 3)))
    return out


def interpolate_forcing_merra2(zout, z, w, u, v, ta, ma, ug, vg, zs):
    """
    interpolate all the forcing profiles to SAM's vertical grid given by `zout`
    """

    nt, nz, ny, nx = z.shape
    zi = np.zeros((nt, nz + 1, ny, nx))
    """
    This is just one way to take into account surface elevation.
    The other way is to interpolate using z (heights above MSL) directly,
    then add a mean surface elevation to zout. Using the latter method, the averaging
    will be done on heights above MSL rather than heights above surface.
    But given the small surface elevations we have around ENA, it shouldn't make much difference. 
    """
    zi[:, :-1, :, :] = z.values - zs.values[:, np.newaxis, :, :]
    # zi[:, :-1, :, :] = z.values
    # zout = zout + zs.values.mean()
    out = []
    sfc = [0, 0, 0, 1, 1, 1, 1]
    for v, s in zip([w, u, v, ta, ma, ug, vg], sfc):
        vi = np.zeros_like(zi)
        vi[:, :-1, :, :] = v.values
        # surface winds are all zero
        if s == 1:
            vi[:, -1, :, :] = v.values[:, -1, :, :]
        # vi[:, -1, :, :] = v.values[:, -1, :, :]
        vo = interpolate_1d(zout, zi, vi, axis=1)
        out.append(vo.mean(axis=(2, 3)))
    return out


def calc_tp_era5(data):
    v = np.genfromtxt(
        "era5_table.csv", delimiter=",", skip_header=1, missing_values="-"
    )
    a = v[:, 1]
    b = v[:, 2]
    t = data.t
    lnsp = data.lnsp.values[:, 0, :, :]
    nt, nz, ny, nx = t.shape
    pi = np.zeros((nt, nz + 1, ny, nx))
    pi[:] = np.exp(lnsp[:, np.newaxis, :, :])
    pi = (
        a[np.newaxis, :, np.newaxis, np.newaxis]
        + pi * b[np.newaxis, :, np.newaxis, np.newaxis]
    )
    p = (pi[:, 1:, :, :] + pi[:, :-1, :, :]) * 0.5
    tpv = metpy.calc.potential_temperature(p * units.pascal, t)
    tp = t.copy(deep=True, data=tpv)
    tp.rename("Air Potential Temperature")
    return tp