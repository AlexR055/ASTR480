import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body
import astropy.coordinates
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def object_visibility(obj_name, obs_lat, obs_lon, obs_alt, date):
    location = EarthLocation(lat=obs_lat * u.deg, lon=obs_lon * u.deg, height=obs_alt * u.m)  # obsevers location, altitude in m above sea level
    times = Time(date) + np.linspace(0, 24, 100) * u.hour
    altaz_frame = AltAz(obstime=times, location=location)

    obj = get_body(obj_name, times)  # fetches object of interests position
    obj_altaz = obj.transform_to(altaz_frame)

    sun_altaz = get_body("sun", times).transform_to(altaz_frame)  # sun altitude
    moon_altaz = get_body("moon", times).transform_to(altaz_frame)  # moon altitude

    plt.figure(figsize=(10, 5))
    plt.plot(times.datetime, obj_altaz.alt, label=obj_name, color='b')
    plt.plot(times.datetime, sun_altaz.alt, label='Sun', color='r', linestyle='dashed')
    plt.plot(times.datetime, moon_altaz.alt, label='Moon', color='gray', linestyle='dotted')

    plt.axhline(0, color='black', linestyle='dashed', linewidth=0.7)
    plt.xlabel("Time (UTC)")
    plt.ylabel("Altitude (°)")
    plt.title(f"Visibility of {obj_name} on {date}\nfrom lat={obs_lat}, lon={obs_lon}")
    plt.legend()
    plt.xticks(rotation=30, ha="right")  # rotated for interpretability
    plt.gca().xaxis.set_major_locator(MaxNLocator(6))  # limits to 6 major ticks with maxnlocator
    
    plt.grid()
    plt.show()

# Call the function
object_visibility("Jupiter", 47.6062, -122.3321, 50, "2025-03-01") #format: object_visibility("obj_name", observer lat, observer long, observer alt, date)
