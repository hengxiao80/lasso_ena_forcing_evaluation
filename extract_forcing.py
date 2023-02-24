import numpy as np
import xarray as xr
import metpy
from metpy.units import units
from metpy.constants import water_heat_vaporization, dry_air_gas_constant, earth_gravity
from metpy.interpolate import interpolate_1d


def extract_sounding_sonde(sonde_file):
    sonde = xr.open_dataset(sonde_file)
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
    return z, tp, rv, u, v


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
    print(omega.metpy.units)
    print(tm.metpy.units)
    print(rd.units)
    print(p.metpy.units)
    print(w.metpy.units)

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
    print(omega.metpy.units)
    print(tm.metpy.units)
    print(rd.units)
    print(p.metpy.units)
    print(w.metpy.units)

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
