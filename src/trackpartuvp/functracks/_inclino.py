# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Handle inclino data
2022
@author: Manon Laget
------------------------------
"""


import pandas as pd
import numpy as np
from datetime import timedelta


def add_inclination(inclino_params, df_particles):
    """
    Add trap inclination for each particle.

    Parameters
    ----------
    inclino_params: pandas dataframe
        Dataframe containing parameters to use to prepare inclino data.
    df_particles: pandas dataframe
        Dataframe containing particle data.

    Returns
    -------
    None.

    """

    # Get parameters
    # Path to inclino file
    path = inclino_params['path'].iloc[0]
    # Axis to use to correct angle
    axis = inclino_params['axis'].iloc[0]
    # Whether to add or substract the value
    sign = inclino_params['sign'].iloc[0]

    # Reading file
    df_inclino = pd.read_csv(path, sep=r'[\t;]')

    # Converting string datetime to datetime format
    df_inclino['Date-Time'] = pd.to_datetime(
        df_inclino['Date-Time'], dayfirst=True)
    arr_date = df_inclino['Date-Time'].dt.to_pydatetime()
    df_inclino['Date-Time'] = pd.Series(arr_date, dtype=object)

    if df_inclino is None:
        return df_particles

    # Select good column

    if axis == 'x':
        df_inclino = df_inclino[['Date-Time', 'Tilt-X']]
        df_inclino['Tilt-X'] = 90 - df_inclino['Tilt-X']
        name = 'Tilt-X'

    elif axis == 'y':
        df_inclino = df_inclino[['Date-Time', 'Tilt-Y']]
        name = 'Tilt-Y'

    elif axis == 'z':
        df_inclino = df_inclino[['Date-Time', 'Tilt-Z']]
        name = 'Tilt-Z'

    elif axis == 'roll':
        df_inclino = df_inclino[['Date-Time', 'roll']]
        name = 'roll'

    else:
        print('Wrong axis selected, inclinometer data ignored')
        return df_particles

    # Get unique datetime first

    df_particles['datetime2'] = pd.to_datetime(
        df_particles['datetime']).dt.to_pydatetime()

    date_times = np.unique(
        pd.to_datetime(df_particles['datetime']).dt.to_pydatetime())

    for date_time in date_times:

        start = date_time - timedelta(minutes=1)
        end = date_time + timedelta(minutes=1)

        # Get the timestamps inside the bin
        rows = df_inclino.loc[df_inclino['Date-Time'] >= start]
        rows = rows.loc[rows['Date-Time'] < end]

        # Get mean value
        if rows.empty:
            mean_val = None
        elif sign == '-':
            mean_val = -rows[name].mean()
        else:
            mean_val = rows[name].mean()

        # Add this value to df particles
        df_particles.loc[
            df_particles['datetime2'] == date_time, 'roll'] = mean_val

    return df_particles
