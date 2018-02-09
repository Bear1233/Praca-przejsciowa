"""
Created: 2018-01-05
Author: Filip Perczyński
------------------------------------------------------------------------------
Funkcja wprowadzania danych wejsciowych:
- funkcja mission_date() - wprowadzenie daty startu i rozpoczęcia misji
    w formacie [rrrr]-[mm]-[dd], sprawdzenie poprawnosci danych np. formatu,
    liczby dni, miesiecy, lat przestepnych
- funkcja planets() - zwraca nazwę planety w formacie odczytywanym przez poliastro oraz
    wartosc minimalnej i maksymalnej liczby dni potrzebnych na podroz do planety
- funkcja inputf() - wprowadzenie danych przez użytkownika oraz sprawdzenie ich 
    poprawnosci, funkcja zwraca następujące dane:
        - nr - wybór opcji obliczen, 
        - m_spacecraft - masa statku kosmicznego, 
        - I_sp - impuls wlasciwy, 
        - H - wysokosc orbity poczatkowej, 
        - planet - wybrana planeta, 
        - date_launch - data startu, 
        - date_arrival - data końca misji, 
        - transit_min - minimalna liczba dni podrozy, 
        - transit_max - maksymalna liczba dni podrozy.
------------------------------------------------------------------------------
"""
import astropy.units as u
from astropy import time

def mission_date():
    check = True
    while check:
        print()
        print('Data w granicach 1900 - 2100 rok')
        date_ymd =input(('Rok-Miesiac-Dzien [rrrr]-[mm]-[dd]: '))
        
        A=[]
        A=date_ymd.split("-")
        if "" in A:
            A.remove("")
        else:
            A = A
        
        if date_ymd.count("-") != 2:
            print()
            print (('Zly format daty'))
        elif len(A) != 3:
            print()
            print (('Zly format daty'))
        else:
            y,m,d = date_ymd.split("-") 
            
            if 1899 < int(y) < 2101:
                if 0 < int(m) < 13:
                    if int(m) in [1, 3, 5, 7, 8, 10, 12]:
                        if 0 < int(d) < 32:
                            check = False
                        else:
                            print()
                            print(('Podany dzien jest bledny'))
                    elif int(m) in [4, 6, 9, 11]:
                        if 0 < int(d) < 31:
                            check = False
                        else:
                            print()
                            print(('Podany dzien jest bledny'))
                    elif int(m) == 2:
                        if (int(y) % 4) == 0:
                            if (int(y) % 100) == 0:
                                if (int(y) % 400) == 0:
                                    if 0 < int(d) < 30:
                                        check = False
                                    else:
                                        print()
                                        print(('Podany dzien jest bledny'))
                                else:
                                    if 0 < int(d) < 29:
                                        check = False
                                    else:
                                        print()
                                        print(('Podany dzien jest bledny'))
                            else:
                                if 0 < int(d) < 30:
                                    check = False
                                else:
                                    print()
                                    print(('Podany dzien jest bledny'))
                        else:
                             if 0 < int(d) < 29:
                                 check = False
                             else:
                                 print()
                                 print(('Podany dzien jest bledny'))
                else:
                    print()
                    print(('Podany miesiac jest bledny')) 
            else:
                print()
                print(('Podany rok jest poza zakresem'))
              
    date = str(date_ymd) + ' 12:00'   
    date = time.Time(date, format='iso', scale='utc')
    return (date)


def planets(planet):
    if planet == 'Merkury':
        planet = 'mercury'
        as1 = 0
        as2 = 0
        as3 = 0
        trans_min = 100
        trans_min_as1 = 0
        trans_min_as2 = 0
        trans_min_as3 = 0
        trans_min_v = 0
        trans_max = 200
        trans_max_as1 = 0
        trans_max_as2 = 0
        trans_max_as3 = 0
        trans_max_v = 0
    if planet == 'Wenus':
        planet = 'venus'
        as1 = 0
        as2 = 0
        as3 = 0
        trans_min = 50
        trans_min_as1 = 0
        trans_min_as2 = 0
        trans_min_as3 = 0
        trans_min_v = 0
        trans_max = 100
        trans_max_as1 = 0
        trans_max_as2 = 0
        trans_max_as3 = 0
        trans_max_v = 0
    if planet == 'Mars':
        planet = 'mars'
        as1 = 'venus'
        as2 = 0
        as3 = 0
        trans_min = 100
        trans_min_as1 = 50
        trans_min_as2 = 0
        trans_min_as3 = 0
        trans_min_v = 100
        trans_max = 150
        trans_max_as1 = 100
        trans_max_as2 = 0
        trans_max_as3 = 0
        trans_max_v = 300
    if planet == 'Jowisz':
        planet = 'jupiter'
        as1 = 'venus'
        as2 = 0
        as3 = 0
        trans_min = 300
        trans_min_as1 = 50
        trans_min_as2 = 200
        trans_min_as3 = 0
        trans_min_v = 300
        trans_max = 500
        trans_max_as1 = 100
        trans_max_as2 = 600
        trans_max_as3 = 0
        trans_max_v = 500
    if planet == 'Saturn':
        planet = 'saturn'
        as1 = 'venus' 
        as2 = 'jupiter'
        as3 = 0
        trans_min = 1000
        trans_min_as1 = 50
        trans_min_as2 = 700
        trans_min_as3 = 700
        trans_min_v = 1000
        trans_max = 2000
        trans_max_as1 = 100
        trans_max_as2 = 1000
        trans_max_as3 = 1000
        trans_max_v = 2000
    if planet == 'Uran':
        planet = 'uranus'
        as1 = 'venus'
        as2 = 'jupiter'
        as3 = 'saturn'
        trans_min = 4000
        trans_min_as1 = 50
        trans_min_as2 = 2000
        trans_min_as3 = 2000
        trans_min_v = 3000
        trans_max = 5000
        trans_max_as1 = 100
        trans_max_as2 = 2500
        trans_max_as3 = 2500
        trans_max_v = 3500
    if planet == 'Neptun':
        planet = 'neptune'
        as1 = 'jupiter'
        as2 = 'saturn'
        as3 = 'uranus'
        trans_min = 5000
        trans_min_as1 = 1000
        trans_min_as2 = 2000
        trans_min_as3 = 3000
        trans_min_v = 3000
        trans_max = 6000
        trans_max_as1 = 1500
        trans_max_as2 = 2500
        trans_max_as3 = 3500
        trans_max_v = 3500
    
    transit_min = [trans_min, trans_min_as1, trans_min_as2, trans_min_as3, trans_min_v]* u.day
    transit_max = [trans_max, trans_max_as1, trans_max_as2, trans_max_as3, trans_max_v]* u.day
    return planet, as1, as2, as3, transit_min, transit_max

    
def inputf():
    
    print ('Możliwe są dwie rodzaje analiz:')
    print ('1 - trajektoria lotu bezposredniego do wybranej planety na podstawie daty startu i ladowania')
    print ('2 - optymalizacja daty startu i lądowania ze względu na koszt transferu dla lotu bezporedniego')
    print ('3 - optymalizacja daty startu i lądowania ze względu na koszt transferu dla lotu z wykorzystaniem asysty grawitacyjnej dla zewnętrznych planet Układu Słonecznego')
    print ('4 - obliczenia transferu z małym ciągiem')
    opt = input(('Wybierz odpowiedni numer:  '))
    nr = int(opt)
    print()
    print('-'*80)
    print (('Wprowadź parametry statku kosmicznego'))
    print('-'*80)
    
    if nr in [2,3,4]:
        check = True        
        while check:
            print()
            print ('Masa własna statku kosmicznego w granicach 0-150 000 kg')
            m_spacecraft = int(input(('M [kg]: ')))
            
            if 0 < m_spacecraft < 100000:
                check = False
            else:
                print(('Podana masa jest poza zakresem'))
        
        check = True        
        while check:
            print()
            print ('Impuls własciwy statku kosmicznego w granicach 0-150 000 kg')
            I_sp = int(input(('Isp [m/s]: ')))
            
            if 0 < I_sp < 100000:
                check = False
            else:
                print(('Podany impuls jest poza zakresem'))
        
        check = True
        while check:
            print()
            print('Wysokosc poczatkowej orbity kolowej w granicach 100 - 10 000 km')
            H = int(input(('H [km]: ')))
    
            if 100 < H < 10000:
                check = False
            else:
                print(('Podana wysokosc jest poza zakresem'))
    
    check = True
    while check:
        print()
        print(('Wybierz cel misji - jedna z planet Ukladu Slonecznego'))
        print(('Merkury, Wenus, Mars, Jowisz, Saturn, Uran, Neptun'))
        planet = input(('Cel: '))

        if planet in ['Merkury', 'Wenus', 'Mars', 'Jowisz', 'Saturn', 'Uran', 'Neptun']:
            check = False
            planet, as1, as2, as3, transit_min, transit_max = planets(planet)
            if nr == 3:
                print()
                if as2 == 0:
                    print('Lot z wykorzystaniem asysyty grawitacyjnej:', as1.title())
                elif as3 ==0:
                    print('Lot z wykorzystaniem asysyty grawitacyjnej:', as1.title(), 'i', as2.title())
                else:
                    print('Lot z wykorzystaniem asysyty grawitacyjnej:', as1.title(), ',', as2.title(), 'i', as3.title())
            
        else:
            print(('Podana planeta jest bledna'))
                
    print()
    print(('Wprowadź datę startu misji'))
    date_launch = mission_date()
    print()
    print(('Wprowadź datę zakończenia misji'))
    date_arrival = mission_date()
    
# Parametry niewykorzystywane w analize 1
    if nr == 1:
        m_spacecraft = 0
        H = 0
        I_sp = 0
# nadanie wprowadzonym parametrom jednostek    
    m_spacecraft = m_spacecraft * u.kg
    H = H * u.km
    I_sp = I_sp/1000 * u.km / u.s
    date_launch = time.Time(date_launch.jd, format='jd', scale='utc')
    date_arrival = time.Time(date_arrival.jd, format='jd', scale='utc')

    return nr, m_spacecraft, I_sp, H, planet, date_launch, date_arrival, as1, as2, as3, transit_min, transit_max
