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

    def __init__(self, pr=0, mas_spo=0, przysp=0, typ=0, trojwymiar=0, pr_x=0, pr_y=0, pr_z=0):
        if trojwymiar == 0:
            if typ == 0:
                while pr >= 1:
                    pr = float(input('Podaj poprawna predkosc (mniej niz 1)\n'))
                self.predkosc_ul = mp.mpf(pr)
                self.predkosc_ms = mp.mpf(pr * pr_swiatla)
                self.przyspieszenie_ul = mp.mpf(przysp)
                self.przyspieszenie_ms = mp.mpf(przysp * pr_swiatla)
            else:
                while pr >= pr_swiatla:
                    pr = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)\n'))
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
                        pr_x = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi x ponownie (mniej niz 1)\n"))
                        pr_y = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi y ponownie (mniej niz 1)\n"))
                        pr_z = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi z ponownie (mniej niz 1)\n"))
                    else:
                        z = False
                self.predkosc_ms = mp.mpf(self.predkosc_ul * pr_swiatla)
                self.przysp_3D(typ=0)
            else:
                z = True
                while z:
                    self.predkosc_ms = mp.mpf(((pr_x ** 2) + (pr_y ** 2) + (pr_z ** 2)) ** (0.5))
                    if self.predkosc_ms >= pr_swiatla:
                        pr_x = float(input(
                        "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi x ponownie (mniej niz 299 792 458 m/s)\n"))
                        pr_y = float(input(
                        "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi y ponownie (mniej niz 299 792 458 m/s)\n"))
                        pr_z = float(input(
                        "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi z ponownie (mniej niz 299 792 458 m/s)\n"))

                    else:
                        z = False
                self.predkosc_ul = mp.mpf(self.predkosc_ms / pr_swiatla)
                self.przysp_3D(typ=1)
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
    def przysp_3D(self, typ=0):
        if self.predkosc_ul==0:
            a_x = float(input('podaj przyspiesznie wzgledem osi x\n'))
            a_y = float(input('podaj przyspiesznie wzgledem osi y\n'))
            a_z = float(input('podaj przyspiesznie wzgledem osi z\n'))
            if typ==0:
                self.przyspieszenie_ul=mp.mpf(((a_x ** 2) + (a_y ** 2) + (a_z ** 2)) ** (0.5))
                self.przyspieszenie_ms=self.przyspieszenie_ul*pr_swiatla
            else:
                self.przyspieszenie_ms = mp.mpf(((a_x ** 2) + (a_y ** 2) + (a_z ** 2)) ** (0.5))
                self.przyspieszenie_ul = self.przyspieszenie_ms / pr_swiatla
        else:
            przysp=float(input('podaj przyspieszenie czastki\n'))
            if typ==0:
                self.przyspieszenie_ul=mp.mpf(przysp)
                self.przyspieszenie_ms=self.przyspieszenie_ul*pr_swiatla
            else:
                self.przyspieszenie_ms = mp.mpf(przysp)
                self.przyspieszenie_ul = self.przyspieszenie_ms * pr_swiatla


def menu():

    print('to program liczacy droge pokonana przez czastki, z wykorzystaniem wzorow relatywistycznych\nczy działać w trójwymiarze? nie(0)/(1)tak')
    z=True
    while z:
        troj=int(input())
        if troj==1 or troj==0:
            z=False
        else:
            print('podaj poprawna wartosc nie(0)/(1)tak ')
    z=True
    print('czy wartosci beda podawane jako ulamek predkosci swiatla (0) czy w metrach na sekunde (1)')
    while z:
        typdanych=int(input())
        if typdanych==1 or typdanych==0:
            z=False
        else:
            print('nie poprawna wartosc, podaj ponownie, wartosci jako ulamek predkosci swiatla (0) czy w metrach na sekunde (1)')
    if troj==1:
        v_x = float(input('podaj predkosc wzgledem osi x\n'))
        v_y = float(input('podaj predkosc wzgledem osi y\n'))
        v_z = float(input('podaj predkosc wzgledem osi z\n'))
        czasteczka = Czastka(typ=typdanych, trojwymiar=1, pr_x=v_x, pr_y=v_y, pr_z=v_z)
    else:
        v=float(input('podaj predkosc czastki\n'))
        a=float(input('podaj przyspieszenie czastki\n'))
        czasteczka = Czastka(pr=v, przysp=a, typ=typdanych)
    t=float(input('w jakim czasie czastka wykonuje ruch?\n'))
    rel, new=czasteczka.droga(t)
    print('droga obliczona ze wzorow mechaniki klasycznej wynosi\n',new,'\ndroga obliczona ze wzorow relatywistycznych wynosi\n',rel,'\ni jest o ',(new-rel)*100/new,'% mniejsza')
    while True:
        dalej = int(input('\n\nopcje:\n(1)-narysuj wykres\n(2)-zacznij ponownie \n(3)-zakończ działanie programu\n'))
        if dalej ==1:
            x = np.linspace(1, float(czasteczka.czas()), 1000)
            y = czasteczka.droga(x)[0]
            z = pr_swiatla * x - pr_swiatla ** 2 / czasteczka.przyspieszenie_ms
            w = czasteczka.droga(x)[1]
            plt.plot(x, y, color='g', lw=1, ls='-.', label='droga rel.')
            plt.plot(x, z)
            plt.plot(x, w, label='droga new.')
            plt.show()
        elif dalej==2:
            czasteczka=None
            menu()
        elif dalej==3:
            break





menu()

#czasteczka = Czastka(0, 0, typ=1, trojwymiar=1,pr_x=30, pr_y=30, pr_z=30)

#print(czasteczka.predkosc_przyspieszanego_obiektu(10000))

#print(czasteczka.droga(1)[0], czasteczka.droga(1)[1])


#x = np.linspace(1, float(czasteczka.czas()), 1000)
#y = czasteczka.droga(x)[0]
#z = pr_swiatla * x - pr_swiatla ** 2 / czasteczka.przyspieszenie_ms
#w = czasteczka.droga(x)[1]
#plt.plot(x, y, color='g', lw=1, ls='-.', label='droga rel.')
#plt.plot(x, z)
#plt.plot(x, w, label='droga new.')
#plt.show()