import numpy as np

import mpmath as mp

mp.mp.dps = 50

pr_swiatla = 299792458
c = 1


class Czastka():
    predkosc_ul = None
    predkosc_ms = None
    masa_spoczynkowa = None
    masa_relatywistyczna = None
    przyspieszenie_ul = None
    przyspieszenie_ms = None
    wektor_predkosci = None

    def __init__(self, pr, mas_spo, przysp=0, typ=0, trojwymiar=0, pr_x=0, pr_y=0, pr_z=0, a_x=0, a_y=0, a_z=0):
        if trojwymiar==0:
            if typ == 0:
                while pr >= 1:
                    pr = float(input('Podaj poprawna predkosc (mniej niz 1)'))
                self.predkosc_ul = np.array([pr])
                self.predkosc_ms = np.array([pr * pr_swiatla])
                self.przyspieszenie_ul = np.array([przysp])
                self.przyspieszenie_ms = np.array([przysp * pr_swiatla])
            else:
                while pr >= pr_swiatla:
                    pr = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
                self.predkosc_ul = np.array([pr / pr_swiatla])
                self.predkosc_ms = np.array([pr])
                self.przyspieszenie_ms = np.array([przysp])
                self.przyspieszenie_ul = np.array([przysp / pr_swiatla])
        else:
            if typ==0:
                while abs(pr_x) >= 1:
                    pr_x = float(input('Podaj poprawna predkosc (mniej niz 1)'))
                while abs(pr_y) >= 1 or ((pr_x**2)+(pr_y**2))**(0.5)>=1:
                    pr_y = float(input('Podaj poprawna predkosc (mniej niz 1)'))
                while abs(pr_z) >= 1 or ((pr_x**2)+(pr_y**2)+(pr_z**2))**(0.5)>=1:
                    pr_z = float(input('Podaj poprawna predkosc (mniej niz 1)'))
                self.predkosc_ul=((pr_x**2)+(pr_y**2)+(pr_z**2))**(0.5)
                self.predkosc_ms = ((pr_x**2)+(pr_y**2)+(pr_z**2)) * pr_swiatla
                self.przyspieszenie_ul = ((a_x**2)+(a_y**2)+(a_z**2))**(0.5)
                self.przyspieszenie_ms = ((a_x**2)+(a_y**2)+(a_z**2))**(0.5) * pr_swiatla
                wektor=np.array([pr_x,pr_y,pr_z])
                self.wektor_predkosci=wektor/wektor.max()
            else:
                while abs(pr_x) >= pr_swiatla:
                    pr_x = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
                while abs(pr_y) >= pr_swiatla or ((pr_x**2)+(pr_y**2))**(0.5)>=pr_swiatla:
                    pr_y = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
                while abs(pr_z) >= pr_swiatla or ((pr_x**2)+(pr_y**2)+(pr_z**2))**(0.5)>=pr_swiatla:
                    pr_z = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
                self.predkosc_ms = ((pr_x**2)+(pr_y**2)+(pr_z**2))**(0.5)
                self.predkosc_ul = ((pr_x**2)+(pr_y**2)+(pr_z**2))**(0.5)/pr_swiatla
                self.przyspieszenie_ms = ((a_x**2)+(a_y**2)+(a_z**2))**(0.5)
                self.przyspieszenie_ul = ((a_x**2)+(a_y**2)+(a_z**2))**(0.5) / pr_swiatla
                wektor=np.array([pr_x,pr_y,pr_z])
                self.wektor_predkosci=wektor/wektor.max()
        self.czynnik_lor = np.array([(1 - self.predkosc_ul ** 2) **(-0.5)])
        self.masa_spoczynkowa = mas_spo
        # self.masa_relatywistyczna = self.masa_spoczynkowa * self.predkosc_ms
        # self.energia_kinetyczna=(1/2)*self.masa*self.predkosc**2

    def droga(self, t):
        if self.przyspieszenie_ms == 0:
            t_dylatacja = t/self.czynnik_lor
            droga = self.predkosc_ms * t_dylatacja
            droga_new=self.predkosc_ms*t
        else:
            droga = np.array([((np.sqrt(1+(((self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor)**2) \
                                  / pr_swiatla**2)))- self.czynnik_lor)* (pr_swiatla ** 2)/self.przyspieszenie_ms])
            droga_new=self.predkosc_ms*t+(1/2)*self.przyspieszenie_ms*t**2
        return droga,droga_new

dane = np.array([10000, 1, 300])
czasteczka = Czastka(mp.mpf(100), 0, przysp = mp.mpf(3),typ=1)


print(czasteczka.droga(1))