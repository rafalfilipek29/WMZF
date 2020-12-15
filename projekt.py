import numpy as np

import mpmath as mp

import matplotlib.pyplot as plt

mp.mp.dps = 50

pr_swiatla = mp.mpf(299792458)
c = 1


class Czastka():
    predkosc_ul = None
    predkosc_ms = None
    masa_spoczynkowa = None
    masa_relatywistyczna = None
    przyspieszenie_ul = None
    przyspieszenie_ms = None
    wektor_predkosci = None

    def __init__(self, pr, mas_spo, przysp=0, typ=0, trojwymiar=0, pr_x=0, pr_y=0, pr_z=0):
        if trojwymiar == 0:
            if typ == 0:
                while pr >= 1:
                    pr = float(input('Podaj poprawna predkosc (mniej niz 1)'))
                self.predkosc_ul = mp.mpf(pr)
                self.predkosc_ms = mp.mpf(pr * pr_swiatla)
                self.przyspieszenie_ul = mp.mpf(przysp)
                self.przyspieszenie_ms = mp.mpf(przysp * pr_swiatla)
            else:
                while pr >= pr_swiatla:
                    pr = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
                self.predkosc_ul = mp.mpf(pr / pr_swiatla)
                self.predkosc_ms = mp.mpf(pr)
                self.przyspieszenie_ms = mp.mpf(przysp)
                self.przyspieszenie_ul = mp.mpf(przysp / pr_swiatla)
        else:
            if typ == 0:
                z = True
                while z:
                    self.predkosc_ul = mp.mpf(((pr_x ** 2) + (pr_y ** 2) + (pr_z ** 2)) ** (0.5))
                    if self.predkosc_ul >= 1:
                        pr_x, pr_y, pr_z = float(
                            input("Zbyt duza predkosc, podaj wartosci wspolrzednych ponownie (mniej niz 1)"))
                    else:
                        z = False
                self.predkosc_ms = mp.mpf(self.predkosc_ul * pr_swiatla)
            else:
                z = True
                while z:
                    self.predkosc_ms = mp.mpf(((pr_x ** 2) + (pr_y ** 2) + (pr_z ** 2)) ** (0.5))
                    if self.predkosc_ul >= pr_swiatla:
                        pr_x, pr_y, pr_z = float(input(
                            "Zbyt duza predkosc, podaj wartosci wspolrzednych ponownie (mniej niz 299 792 458 m/s)"))
                    else:
                        z = False
                self.predkosc_ul = mp.mpf(self.predkosc_ms / pr_swiatla)
            self.wektor_predkosci = np.array([mp.mpf(pr_x), mp.mpf(pr_y), mp.mpf(pr_z)])
        self.czynnik_lor = mp.mpf((1 - self.predkosc_ul ** 2) ** (-0.5))
        self.masa_spoczynkowa = mas_spo
        # self.masa_relatywistyczna = self.masa_spoczynkowa * self.predkosc_ms

    def droga(self, t):
        if self.przyspieszenie_ms == 0:
            t_dylatacja = mp.mpf(t / self.czynnik_lor)
            droga = self.predkosc_ms * t_dylatacja
            droga_new = self.predkosc_ms * t
        else:
            droga = ((np.sqrt(1 + (((self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor) ** 2) \
                                   / pr_swiatla ** 2))) - self.czynnik_lor) * (pr_swiatla ** 2) / self.przyspieszenie_ms
            droga_new = self.predkosc_ms * t + (1 / 2) * self.przyspieszenie_ms * t ** 2
        return droga, droga_new

    def predkosc_przyspieszanego_obiektu(self, t):
        # t_dylatacja = mp.mpf(t/self.czynnik_lor)
        predkosc = (self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor) * (
                    1 + (((self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor) ** 2) \
                         / pr_swiatla ** 2)) ** (-0.5)
        return predkosc

    def czas(self):
        t = (((299000000 * pr_swiatla) ** 2 / (pr_swiatla ** 2 - 299000000 ** 2)) ** (
            0.5) - self.predkosc_ms * self.czynnik_lor) / self.przyspieszenie_ms
        return t
    # def przysp_dla_0


czasteczka = Czastka(1000, 0, przysp=1000, typ=1)

print(czasteczka.predkosc_przyspieszanego_obiektu(10000))

print(czasteczka.droga(1))


x = np.linspace(1, float(czasteczka.czas()), 1000)
y = czasteczka.droga(x)[0]
z = pr_swiatla * x - pr_swiatla ** 2 / czasteczka.przyspieszenie_ms
w = czasteczka.droga(x)[1]
plt.plot(x, y, color='g', lw=1, ls='-.', label='sin(x)')
plt.plot(x, z)
plt.plot(x, w)
plt.show()
