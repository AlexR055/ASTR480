import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_moon, SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt

def object_visibility(obj_name, observer_lat, observer_lon, observer_alt, date):
    location = EarthLocation(lat=obs_lat * u.deg, lon=obs_lon * u.deg, height=obs_alt * u.m) #input user/obs position
    times = Time(date) + np.linspace(0, 24, 100) * u.hour
    altaz_frame = AltAz(obstime=times, location=location)
    
    obj = SkyCoord.from_name(obj_name) # these two lines define object of interests position
    obj_altaz = obj.transform_to(altaz_frame)
    
    sun_altaz = get_sun(times).transform_to(altaz_frame) # gets sun altitude
    moon_altaz = get_moon(times).transform_to(altaz_frame) # gets moon altitude
    
    #creates plot
    plt.figure(figsize=(10, 5))
    plt.plot(times.datetime, obj_altaz.alt, label=obj_name, color='b')
    plt.plot(times.datetime, sun_altaz.alt, label='Sun', color='r', linestyle='dashed')
    plt.plot(times.datetime, moon_altaz.alt, label='Moon', color='gray', linestyle='dotted')
    plt.axhline(0, color='black', linestyle='dashed', linewidth=0.7)
    plt.xlabel("Time (UTC)")
    plt.ylabel("Altitude (degrees)")
    plt.title(f"Visibility of {obj_name} on {date} from lat={obs_lat}, lon={obs_lon}")
    plt.legend()
    plt.xticks(rotation=45)
    plt.grid()
    plt.show()

# call the function
object_visibility("Jupiter", 47.6062, -122.3321, 50, "2025-03-01")
