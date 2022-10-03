import numpy as np
import datetime
import argparse

from typing import Union
from netCDF4 import Dataset


def parse_args():

    args_parser = argparse.ArgumentParser()

    args_parser.description = "This script returns an height-integrated refractive indice for the atmosphere over a " \
                              "specified geological event of interest. The integrated refractive indice is used to " \
                              "calculate the phase shift cause by the troposphere for inSAR imagery correction."
    args_parser.add_argument("--dataset", "-d",
                             help="Path to the .nc file containing atmospheric data to be parsed.",
                             required=True)
    # args_parser.add_argument("--timestep", "-ts",
    #                          help="Time step between each data point of the atmospheric dataset.",
    #                          required=True)
    args_parser.add_argument("--timestamp", "-t",
                             help="The date and time of the geological event of interest.",
                             required=True)
    args_parser.add_argument("--time_format", "-f",
                             help="The time format of the timestamp. See "
                                  "https://www.geeksforgeeks.org/python-datetime-strptime-function/"
                                  "for alternative time formats",
                             default="%Y-%m-%d:%H:%M:%S")
    args_parser.add_argument("--latitude", "-lt",
                             help="The latitude of the geological event of interest.",
                             type=float,
                             required=True)
    args_parser.add_argument("--start_height", "-sh",
                             help="The starting height at which the refractive indice integration is computed.",
                             type=float,
                             required=True)
    args_parser.add_argument("--longitude", "-lg",
                             help="The longitude of the geological event of interest.",
                             type=float,
                             required=True)

    return args_parser.parse_args()


def millibar_to_mmHg(pressure: Union[float, np.array]) -> Union[float, np.array]:
    """
    Converts presssure in millibar units to mmHg units
    :param pressure:    Atmospheric pressure in millibar
    :return:            Pressure converted to mmHg
    """
    return np.multiply(pressure, 0.75006)


def mmHg_to_millibar(pressure: Union[float, np.array]) -> Union[float, np.array]:
    """
    Converts presssure from mmHg to millibar units
    :param pressure:    Atmospheric pressure in mmHg
    :return:            Pressure converted to millibar
    """
    return np.multiply(pressure, 1/0.75006)


def get_closest_coordinates(data: dict, position: Union[tuple, np.array]):
    """
    :param data:        The parsed atmospheric data from a nc file.
    :param position:    Latitude/longitude coordinates of a geological event of interest.
    :return:            The indexes of the coordinates of atmospheric data closest to the coordinates of the geological
                        event of interest.
    """
    latitutdes = np.array(data['latitude'])
    longitudes = np.array(data['longitude'])

    diffs = np.subtract(latitutdes, position[0])
    lat_ix = np.where(abs(diffs) == min(abs(diffs)))[0][0]
    # clat = latitutdes[lat_ix]

    diffs = np.subtract(longitudes, position[1])
    long_ix = np.where(abs(diffs) == min(abs(diffs)))[0][0]
    # clong = longitudes[long_ix]

    # return clat, clong
    return lat_ix, long_ix


def get_closest_time(data: dict, timestamp: datetime.datetime) -> tuple:
    """
    This function calculates the time of the atmospheric data closest to the time of the geological event of interest.
    :param data:        The parsed atmospheric data.
    :param timestamp:   The time and date of the geological event of interest.
    :return:            The index of the time from the atmospheric dataset closest to the time of the geological event
                        of interest.
    """
    min_diff = datetime.timedelta(days=1e8)
    start = data.variables['time'].units.split()[2] + ":" + data.variables['time'].units.split()[3][:-2:]
    start_format = "%Y-%m-%d:%H:%M:%S"
    start = datetime.datetime.strptime(start, start_format)
    for ix in range(data.variables['time'].shape[0]):
        diff = abs(timestamp - (start + datetime.timedelta(days=data.variables['time'][ix] / 24)))
        if diff < min_diff:
            min_diff = diff
        else:
            break
    # return data.variables['time'][ix], start + datetime.timedelta(days=data.variables['time'][ix]/24)
    return ix


def get_height_from_acceleration(acceleration: Union[float, np.array]) -> Union[float, np.array]:
    """
    This function calculate the distance from the center of the earth in function to an gravitational acceleration from
    the Newton's equation F = GM1M2 / r ^ 2.
    :param acceleration:    An accelaratin in m/s^2.
    :return:                A height in meters.
    """
    # gravitational constant 6.67430 x 10^−11 N⋅m2⋅kg−2 from https://en.wikipedia.org/wiki/Gravitational_constant
    # earth mass 5.9722 x 10^24 kg https://en.wikipedia.org/wiki/Earth_mass
    # earth radius  6.3781 x 10^6 https://en.wikipedia.org/wiki/Earth_radius
    M1 = 1  # kg *only for units homogeneity
    M2 = 5.9722e24  # kg
    G = 6.6743e-11
    r = 6.3781e6
    R = np.sqrt(np.divide(np.multiply(G, M1, M2), acceleration))

    return np.subtract(R, r)


def get_closest_height_ix(heights_vector: np.array, height) -> float:
    diffs = abs(np.subtract(heights_vector, height))
    min_height_ix = np.where(diffs == min(diffs))[0][0]
    return min_height_ix


def geopotential_to_height(geopotential: Union[float, np.array]) -> Union[float, np.array]:
    avg_r = 6.3781e6
    M2 = 5.9722e24
    G = 6.6743e-11

    height = np.divide(np.multiply(np.multiply(avg_r, G), M2), np.subtract(np.multiply(G, M2), np.multiply(avg_r, geopotential))) - avg_r
    return np.array(height)


def heigh_to_geopotential(height: Union[float, np.array]) -> Union[float, np.array]:
    avg_r = 6.3781e6
    G = 6.6743e-11
    M2 = 5.9722e24

    geopotential = np.multiply(np.multiply(G, M2), np.subtract(1/avg_r, np.divide(1, np.add(avg_r, height))))
    return geopotential


def get_refractive_indice(data: dict, timestamp: datetime.datetime, position: tuple) -> float:
    """
    Returns the refractive indice in function the time.
    :param timestamp: Time, in hours by slice of 6 hours since 1900-01-01.
    :return: Refraction indice
    """

    # air refractive equation constants
    K1 = 0.776  # K * Pa^-1
    K2 = 0.716  # K * Pa^-1
    K3 = 3.75e3  # K**2 * Pa^-1
    K4 = 1.45e3  # m^3 * kg^-1

    # air density equation constants
    Rd = 287.058
    Rv = 461.495

    # Antoine's equation constants
    A = 8.07131
    B = 1730.63
    C = 233.426

    # get space-time coordinates
    time_ix = get_closest_time(data=data, timestamp=timestamp)
    lat_ix, long_ix = get_closest_coordinates(data=data, position=position)

    # height_ix = get_height_from_acceleration(data.variables['z'])
    # heights = geopotential_to_height(data.variables['z'][time_ix, :, lat_ix, long_ix])
    # geopotentials = heigh_to_geopotential(heights)

    # dimensions order: [time, level, latitude, longitude]
    rel_humidity = np.divide(data['r'][time_ix, :, lat_ix, long_ix], 100)  # relative humidity
    specific_humidity = data['q'][time_ix, :, lat_ix, long_ix]  # q: specific humidity
    air_temperature_K = data['t'][time_ix, :, lat_ix, long_ix]  # air temperature

    air_temperature_C = np.subtract(data['t'][time_ix, :, lat_ix, long_ix], 273.15)  # Kelvin to Celsius

    # Antoin's equation: log10(P) = A - B/(C+T) => P = 10 ^ (A - (B/(C+T)) where P is water vapor pressur at equilibrium
    eq_water_vapor_pressur = np.power(10, np.subtract(A, (np.divide(B, (np.add(C, air_temperature_C))))))

    # relative humidity = water vapor pressure / water vapor pressure at equilibrium =>
    # water vapor pressure = relative humidity * water vapor pressure at equilibrium
    water_vapor_pressur = mmHg_to_millibar(eq_water_vapor_pressur) * rel_humidity

    # Dry air pressure: Dp = 0.622 * water_vapor_pressur/Specific humidity
    dry_air_pressure = 0.622 * water_vapor_pressur / specific_humidity

    # moist air density = dry_air_pressure/Rd*T + water_vapor_pressur/Rv*T
    air_density = dry_air_pressure / (Rd * air_temperature_K) + water_vapor_pressur / (Rv * air_temperature_K)

    # Specific clouds water content is estimated as the sum of specifics cloud: ice content (ciwc),
    # liquid water content (clwc), rain water content (crwc), snow water content (cswc)
    specific_cloud_wc = data.variables['ciwc'][time_ix, :, lat_ix, long_ix] + \
                        data.variables['clwc'][time_ix, :, lat_ix, long_ix] + \
                        data.variables['crwc'][time_ix, :, lat_ix, long_ix] + \
                        data.variables['cswc'][time_ix, :, lat_ix, long_ix]
    cloud_water_content = specific_cloud_wc * air_density

    ref_indice = np.divide(np.multiply(K1, dry_air_pressure), air_temperature_K) + \
                 np.divide(np.multiply(K2, water_vapor_pressur), air_temperature_K) + \
                 np.divide(np.multiply(K3, water_vapor_pressur), np.power(air_temperature_K, 2)) + \
                 np.multiply(K4, cloud_water_content)

    # sanity check
    if ref_indice.any() < 1:
        raise ValueError("You would exceed the speed of light! Something is wrong!!!")

    return ref_indice  # "integration" of the refractive indices for the given height slices


def get_optical_path(data: dict,
                     timestamp: datetime.datetime,
                     input_height_max: float = None,
                     input_height_min: float = None) -> float:

    if input_height_min is None:
        input_height_min = np.inf
    if input_height_max is None:
        input_height_max = -np.inf

    # get relevant space-time coordinates
    time_ix = get_closest_time(data=data, timestamp=timestamp)
    # lat_ix, long_ix = get_closest_coordinates(data=data, position=position)

    ref_ixs = np.zeros(shape=(data.variables['t'].shape[2], data.variables['t'].shape[3], data.variables['t'].shape[1]))

    # extract the max and min geopotential from atmospheric data
    for i in range(data.variables['t'].shape[2]):
        for j in range(data.variables['t'].shape[3]):
            heights = np.sort(data.variables['z'][time_ix, :, i, j])
            heights = geopotential_to_height(heights)
            # max_data_geopotential = data.variables['z'][time_ix, 0, lat_ix, long_ix]
            # min_data_geopotential = data.variables['z'][time_ix, data.variables['z'].shape[1]-1, lat_ix, long_ix]

            # compute max and min height from geopotential data
            # max_data_height = geopotential_to_height(max_data_geopotential)
            # min_data_height = geopotential_to_height(min_data_geopotential)

            # compute respectives maxs and mins between atmospheric data and interferometric data
            # h_max = max(input_height_max, max_data_height)
            # h_min = min(input_height_min, min_data_height)

            # compute refractive indice
            ref_indice = get_refractive_indice(data=data,
                                               timestamp=timestamp,
                                               position=(data.variables['latitude'][i], data.variables['latitude'][j]))

            # return ref_indice

            # compute optical path
            # np.multiply(ref_indice, heights)
            # height_ix = get_closest_height_ix(heights_vector=heights, height=position[2])

            #
            # np.multiply(ref_indice, heights)
            ref_ixs[i, j, :] = np.cumsum(np.multiply(ref_indice, heights))
    # return np.cumsum(np.multiply(ref_indice, heights))
    return ref_ixs
            # return np.trapz()


if __name__ == "__main__":

    args = parse_args()

    data = Dataset(args.dataset)
    dummy_position = (args.latitude, args.longitude, args.start_height)
    ts_tamp = args.timestamp
    timstamp_format = args.time_format
    try:
        timestamp = datetime.datetime.strptime(ts_tamp, timstamp_format)
    except Exception:
        raise ValueError(ts_tamp + " is not recognized with regard to the time format " + timstamp_format)

    ref_indice = get_refractive_indice(data=data, timestamp=timestamp, position=dummy_position)
    print("refractive indice is: " + str(ref_indice))

    o_path = get_optical_path(data=data, timestamp=timestamp, position=dummy_position)
    print("optical path is: " + str(o_path))
