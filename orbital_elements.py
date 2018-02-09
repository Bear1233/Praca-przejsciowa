"""
Created: 2018-02-04
Author: Filip Perczyński

"""
import requests
import re
import astropy.units as u
import numpy as np

def orbitalements(planet):
    tle_url = "http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt"
    orb_elem = requests.get(tle_url).text
    orb_file = open('Orbital_elements.txt','w')
    orb_file.write(orb_elem)
    orb_file.close()
    orbital_elem = open('Orbital_elements.txt', 'r').readlines()
    
    orbital_elem = [a for a in orbital_elem if a != '\n']
    orbital_elem = ''.join(orbital_elem)
    lines = orbital_elem.split("\n")
    
    lines = [re.split("\s{2,}", line) for line in lines[12:]]

    for i in range(0, len(lines)-1, 2):
        name = lines[i][0]
        elements = lines[i][1:]
 
        if name == "EM Bary":
            name = "Earth"
            
        (semi_major_axis, eccentricity, inclination, mean_longitude,
            longitude_of_periapsis, longitude_of_ascending_node) = elements    
        
        #body = bodies.setdefault(name, {})
        if name == planet:
            body = {"name": name,
                "semi_major_axis": float(semi_major_axis) ,
                "eccentricity": eccentricity,
                "inclination": float(inclination),
                "longitude_of_periapsis": longitude_of_periapsis,
                "longitude_of_ascending_node": longitude_of_ascending_node,
            }
        
    a = body["semi_major_axis"]
    a = a 
    inc = body['inclination']    
    inc = np.deg2rad(inc)
        
    return(a,inc)