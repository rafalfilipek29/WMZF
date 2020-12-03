import numpy as np

pr_swiatla = 299792458
c = 1


class Czastka():
    predkosc_ul = None
    predkosc_ms = None
    masa_spoczynkowa = None
    masa_relatywistyczna = None
    przyspieszenie_ul = None
    przyspieszenie_ms = None

    def __init__(self, pr, mas_spo, przysp=0, typ=0):
        if typ == 0:
            while pr >= 1:
                pr = float(input('Podaj poprawna predkosc (mniej niz 1)'))
            self.predkosc_ul = pr
            self.predkosc_ms = pr * pr_swiatla
            self.przyspieszenie_ul = przysp
            self.przyspieszenie_ms = przysp * pr_swiatla
        else:
            while pr >= pr_swiatla:
                pr = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
            self.predkosc_ul = pr / pr_swiatla
            self.predkosc_ms = pr
            self.przyspieszenie_ms = przysp
            self.przyspieszenie_ul = przysp / pr_swiatla
        self.czynnik_lor = np.power(1 - self.predkosc_ul ** 2, -0.5)
        self.masa_spoczynkowa = mas_spo
        # self.masa_relatywistyczna = self.masa_spoczynkowa * self.predkosc_ms
        # self.energia_kinetyczna=(1/2)*self.masa*self.predkosc**2

    def droga(self, t):
        if self.przyspieszenie_ms == 0:
            t_dylatacja = t/self.czynnik_lor
            droga = self.predkosc_ms * t_dylatacja
            droga_new=self.predkosc_ms*t
        else:
            droga = ((np.sqrt(1+(((self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor)**2) \
                                  / pr_swiatla**2)))- self.czynnik_lor)* (pr_swiatla ** 2)/self.przyspieszenie_ms
            droga_new=self.predkosc_ms*t+(1/2)*self.przyspieszenie_ms*t**2
        return droga,droga_new


czasteczka = Czastka(10000000, 1, przysp = 3,typ=1)

print(czasteczka.droga(1))

