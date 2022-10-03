import cdsapi

if __name__ == "__main__":
    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'pressure_level': [
                '1', '2', '3',
                '5', '7', '10',
                '20', '30', '50',
                '70', '100', '125',
                '150', '175', '200',
                '225', '250', '300',
                '350', '400', '450',
                '500', '550', '600',
                '650', '700', '750',
                '775', '800', '825',
                '850', '875', '900',
                '925', '950', '975',
                '1000',
            ],
            'year': '2019',
            'month': '07',
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [
                36.25, -117, 35.25,
                -116,
            ],
            'variable': [
                'fraction_of_cloud_cover', 'geopotential', 'relative_humidity',
                'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity',
                'specific_rain_water_content', 'specific_snow_water_content', 'temperature',
            ],
            'day': [
                '02', '03', '04',
                '05', '06',
            ],
        },
        'RidgecrestEQ2019July2-6.nc')

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'pressure_level': [
                '1', '2', '3',
                '5', '7', '10',
                '20', '30', '50',
                '70', '100', '125',
                '150', '175', '200',
                '225', '250', '300',
                '350', '400', '450',
                '500', '550', '600',
                '650', '700', '750',
                '775', '800', '825',
                '850', '875', '900',
                '925', '950', '975',
                '1000',
            ],
            'year': '2016',
            'month': '04',
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [
                34, 129, 32,
                131,
            ],
            'variable': [
                'fraction_of_cloud_cover', 'geopotential', 'relative_humidity',
                'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity',
                'specific_rain_water_content', 'specific_snow_water_content', 'temperature',
            ],
            'day': [
                '13', '14', '15',
                '16', '17',
            ],
        },
        'Kumamoto2016April13-17.nc')

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'pressure_level': [
                '1', '2', '3',
                '5', '7', '10',
                '20', '30', '50',
                '70', '100', '125',
                '150', '175', '200',
                '225', '250', '300',
                '350', '400', '450',
                '500', '550', '600',
                '650', '700', '750',
                '775', '800', '825',
                '850', '875', '900',
                '925', '950', '975',
                '1000',
            ],
            'year': '2018',
            'month': '07',
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [
                20, -155.5, 19,
                -154.5,
            ],
            'variable': [
                'fraction_of_cloud_cover', 'geopotential', 'relative_humidity',
                'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity',
                'specific_rain_water_content', 'specific_snow_water_content', 'temperature',
            ],
            'day': [
                '02', '03', '04',
                '05', '06',
            ],
        },
        'KilaueaEruption2018May2-6.nc')
