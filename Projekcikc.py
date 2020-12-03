import numpy as np

pr_swiatla = 299792458
c = 1


class Czastka():
    predkosc_ul = None
    predkosc_ms = None
    masa_spoczynkowa = None
    masa_relatywistyczna = None

    def __init__(self, pr, mas_spo, typ=0):
        if typ == 0:
            while pr >= 1:
                pr = float(input('Podaj poprawna predkosc (mniej niz 1)'))
            self.predkosc_ul = pr
            self.predkosc_ms=pr*pr_swiatla
        else:
            while pr >= pr_swiatla:
                pr = float(input('Podaj poprawna predkosc (mniej niz 299 792 458 m/s)'))
            self.predkosc_ul = pr / pr_swiatla
            self.predkosc_ms=pr
        self.czynnik_lor=np.power(1-self.predkosc_ul**2,-0.5)
        self.masa_spoczynkowa = mas_spo
        #self.masa_relatywistyczna = self.masa_spoczynkowa * self.predkosc_ms
        # self.energia_kinetyczna=(1/2)*self.masa*self.predkosc**2




    def droga(self,t):
        t_dylatacja=self.czynnik_lor*t
        droga=self.predkosc_ms*t_dylatacja
        xk=self.predkosc_ms*t
        x_rel=xk*self.czynnik_lor
        x_lor=xk*np.power(1-self.predkosc_ul**2,0.5)
        return droga,x_rel,x_lor



czasteczka = Czastka(0.9, 5)


#komentarz do sprawdzenia
#komentarz2



print(czasteczka.droga(5))