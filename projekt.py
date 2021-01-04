# Biblioteka numpy, sluzaca do optymalizowania projektu
import numpy as np

"Biblioteka mpmath, sluzaca do wykonywania precyzyjnych obliczen matematycznych"
import mpmath as mp

"Biblioteka matplotlib, sluzaca do tworzenia wykresow"
import matplotlib.pyplot as plt

"Okreslenie precyzji wykonywania obliczen"
mp.mp.dps = 50

"Zdefiniowane zmienne, zawierajace predkosc swiatla w prozni"
pr_swiatla = mp.mpf(299792458)
c = 1

"Definiowanie klasy, przedstawiajace pojedyncza czastke"


class Czastka:
    # Parametry klasy
    predkosc_ul = None
    predkosc_ms = None
    masa_spoczynkowa = None
    przyspieszenie_ul = None
    przyspieszenie_ms = None
    wektor_predkosci = None
    wektor_przyspieszenia = None
    "Zapisanie parametrow funkcji"

    def __init__(self, pr_x=0.0, pr_y=0.0, pr_z=0.0, przysp_x=0.0, przysp_y=0.0, przysp_z=0.0, mas_spo=0.0, typ=0):
        # parametr "typ" odpowiada za zaznaczenie przez uzytkownika, czy chce podac predkosc oraz przyspieszenie
        # jako ulamek pr. swiatla czy jednak w jednostkach m/s oraz m/s^2. Jezeli typ = 0, to predkosc i przysp sa podawane
        # jako ulamek pr. swiatla.
        if typ == 0:
            z = True
            while z:
                # Definiowanie predkosci jako ulamek pr. swiatla, z dodatkowym warunkiem na wartosc
                self.predkosc_ul = mp.mpf(((pr_x ** 2) + (pr_y ** 2) + (pr_z ** 2)) ** 0.5)
                if self.predkosc_ul >= 1:
                    pr_x = float(input(
                        "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi x ponownie (mniej niz 1)\n"))
                    if pr_y != 0:
                        pr_y = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi y ponownie (mniej niz 1)\n"))
                    if pr_z != 0:
                        pr_z = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi z ponownie (mniej niz 1)\n"))
                else:
                    z = False
            "Definiowanie reszty parametrow"
            self.predkosc_ms = mp.mpf(self.predkosc_ul * pr_swiatla)
            self.przyspieszenie_ul = mp.mpf(((przysp_x ** 2) + (przysp_y ** 2) + (przysp_z ** 2)) ** 0.5)
            self.przyspieszenie_ms = mp.mpf(self.przyspieszenie_ul * pr_swiatla)
            # self.przysp_3D(typ=0)
        else:
            z = True
            while z:
                self.predkosc_ms = mp.mpf(((pr_x ** 2) + (pr_y ** 2) + (pr_z ** 2)) ** 0.5)
                if self.predkosc_ms >= pr_swiatla:
                    pr_x = float(input(
                        "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi x ponownie (mniej niz 299 792 458 "
                        "m/s)\n"))
                    if pr_y != 0:
                        pr_y = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi y ponownie (mniej niz 299 792 "
                            "458 m/s)\n"))
                    if pr_z != 0:
                        pr_z = float(input(
                            "Zbyt duza predkosc, podaj wartosc predkosci wzgledem osi z ponownie (mniej niz 299 792 "
                            "458 m/s)\n"))
                else:
                    z = False
            self.predkosc_ul = mp.mpf(self.predkosc_ms / pr_swiatla)
            self.przyspieszenie_ms = mp.mpf(((przysp_x ** 2) + (przysp_y ** 2) + (przysp_z ** 2)) ** 0.5)
            self.przyspieszenie_ul = mp.mpf(self.przyspieszenie_ms / pr_swiatla)
            # self.przysp_3D(typ=1)
        self.wektor_predkosci = np.array([mp.mpf(pr_x), mp.mpf(pr_y), mp.mpf(pr_z)])
        self.wektor_przyspieszenia = np.array([mp.mpf(przysp_x), mp.mpf(przysp_y), mp.mpf(przysp_z)])
        self.czynnik_lor = mp.mpf((1 - self.predkosc_ul ** 2) ** (-0.5))
        self.masa_spoczynkowa = mas_spo

    "Funkcja liczaca droge przebyta przez otoczenie, wzgledem czastki poruszajacej sie jednostajnie"
    "lub droge przebyta przez obiekt przyspieszajacy jednostajnie, wzgledem obserwatora. Liczy rowniez w "
    "obu przypadkach droge w sensie newtonowskim, dla porownania"

    def droga(self, t):
        if self.przyspieszenie_ms == 0:
            t_dylatacja = t / self.czynnik_lor
            droga = self.predkosc_ms * t_dylatacja
            droga_new = self.predkosc_ms * t
        else:
            droga = ((np.sqrt(1 + (((self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor) ** 2)
                                   / pr_swiatla ** 2))) - self.czynnik_lor) * (pr_swiatla ** 2) / self.przyspieszenie_ms
            droga_new = self.predkosc_ms * t + (1 / 2) * self.przyspieszenie_ms * t ** 2
        return droga, droga_new

    "Funkcja liczaca predkosc obiektu przyspieszanego."

    def predkosc_przyspieszanego_obiektu(self, t):
        predkosc = (self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor) * (
                1 + (((self.przyspieszenie_ms * t + self.predkosc_ms * self.czynnik_lor) ** 2)
                     / pr_swiatla ** 2)) ** (-0.5)
        return predkosc

    "Funkcja liczaca czas, po ktorym predkosc przyspieszanego obiektu jest bardzo bliska pr_swiatla"

    def czas(self):
        t = (((299000000 * pr_swiatla) ** 2 / (pr_swiatla ** 2 - 299000000 ** 2)) ** (
            0.5) - self.predkosc_ms * self.czynnik_lor) / self.przyspieszenie_ms
        return t

    "Funkcja liczaca predkosc jednej czastki, poruszajacej sie ruchem jednostajnym prostoliniowym, wzgledem innej"
    "czastki poruszajacej ruchem jednostajnym prostoliniowym, lecz z innym wektorem predkosci"

    def wzgledna_predkosc(self, czasteczka2):
        pr_2_wzgl_1 = (
            (1 / self.czynnik_lor * (
                    1 - mp.fdot(self.wektor_predkosci, czasteczka2.wektor_predkosci) / pr_swiatla ** 2) *
             (czasteczka2.wektor_predkosci - self.wektor_predkosci + self.wektor_predkosci * (self.czynnik_lor - 1) *
              ((mp.fdot(self.wektor_predkosci, czasteczka2.wektor_predkosci) / self.predkosc_ms ** 2) - 1))))
        szyb_2_wzgle_1 = (pr_2_wzgl_1[0] ** 2 + pr_2_wzgl_1[1] ** 2 + pr_2_wzgl_1[2] ** 2) ** 0.5
        return pr_2_wzgl_1, szyb_2_wzgle_1

    "Funkcja obliczajaca mase relatywistyczna obiektu"

    def masa_relatywistyczna(self, v):
        m_relat = self.masa_spoczynkowa * ((1 - (v / pr_swiatla) ** 2) ** (-0.5))
        return m_relat


"Funkcja bedaca menu"


def menu():
    # Zmienna troj sluzy do przechowania informacji, czy uzytkownik chce rozpatrywac ruch w trojwymiarze czy nie
    # Zmienna typdanych sluzy do przechowania informacji, czy uzytkownik chce miec predkosci i przyspieszenia
    # jako ulamek pr. swiatla, czy jednak jednostki odpowiednio m/s i m/s^2
    global troj, typdanych

    p = True
    while p:
        print(
            'to program liczacy droge pokonana przez czastki, z wykorzystaniem wzorow relatywistycznych.\nCzy działać '
            'w trójwymiarze? nie(0)/(1)tak')
        q = True
        while q:
            try:
                troj = int(input())
            except ValueError:
                troj = 2
            if troj == 1 or troj == 0:
                q = False
            else:
                print('podaj poprawna wartosc nie(0)/(1)tak ')
        q = True
        print('czy wartosci beda podawane jako ulamek predkosci swiatla (0) czy w metrach na sekunde (1)')
        while q:
            try:
                typdanych = int(input())
            except ValueError:
                typdanych = 2
            if typdanych == 1 or typdanych == 0:
                q = False
            else:
                print(
                    'nie poprawna wartosc, podaj ponownie, wartosci jako ulamek predkosci swiatla (0) czy w metrach '
                    'na sekunde (1)\n')
        wybor = int(input('Chcesz policzyc droge dla jednej poruszajacej sie czastki, czy dwoch? Jedna(1)/dwie(2)\n'))
        if wybor == 1:
            m = int(input("podaj mase czastki\n"))
            if troj == 1:
                v_x = float(input('podaj predkosc wzgledem osi x\n'))
                v_y = float(input('podaj predkosc wzgledem osi y\n'))
                v_z = float(input('podaj predkosc wzgledem osi z\n'))
                if v_x == 0 and v_y == 0 and v_z == 0:
                    a_x = float(input('podaj przyspieszenie względem osi x\n'))
                    a_y = float(input('podaj przyspieszenie względem osi y\n'))
                    a_z = float(input('podaj przyspieszenie względem osi z\n'))
                    a = ((a_x ** 2) + (a_y ** 2) + (a_z ** 2)) ** 0.5
                else:
                    a = float(input('podaj przyspieszenie w kierunku predkosci\n'))
                czasteczka = Czastka(typ=typdanych, pr_x=v_x, pr_y=v_y, pr_z=v_z, przysp_x=a, mas_spo=m)
            else:
                v = float(input('podaj predkosc czastki\n'))
                a = float(input('podaj przyspieszenie czastki\n'))
                czasteczka = Czastka(pr_x=v, przysp_x=a, typ=typdanych, mas_spo=m)

            t = float(input('w jakim czasie czastka wykonuje ruch?\n'))
            rel, new = czasteczka.droga(t)
            print('droga obliczona ze wzorow mechaniki klasycznej wynosi\n', new,
                  '\ndroga obliczona ze wzorow relatywistycznych wynosi\n', rel, '\ni jest o ', (new - rel) * 100 / new,
                  '% mniejsza')
            if czasteczka.przyspieszenie_ms == 0:
                print('masa czasteczki poruszajacej sie z taka predkoscia wynosi',
                      czasteczka.masa_relatywistyczna(czasteczka.predkosc_ms))

            if czasteczka.przyspieszenie_ms == 0:
                l = 0
            else:
                l = 1
            while l == 0:
                dalej = int(
                    input('\n\nopcje:\n(1)-narysuj wykres\n(2)-zacznij ponownie \n'
                          '(3)-zakończ działanie programu\n'))
                if dalej == 1:
                    x = np.linspace(1, t, 1000)
                    y = czasteczka.droga(x)[0]
                    plt.plot(x, y, color='g', label='droga rel.')
                    plt.legend('droga rel.', loc='upper center', shadow=True)
                    plt.xlabel('czas')
                    plt.ylabel('droga')
                    plt.show()
                if dalej == 2:
                    l = False
                    break
                if dalej == 3:
                    p = False
                    break
            while l == 1:
                dalej = int(
                    input('\n\nopcje:\n(1)-narysuj wykres\n(2)-zacznij ponownie \n'
                          '(3)-zakończ działanie programu\n'))
                if dalej == 1:
                    rodzaj = int(input('\n\nna wykresie ma byc\n(1) - droga relatywistyczna i newtonowska\n '
                                       '(2) - droga relatywistyczna i droga fotonu\n'
                                       '(3) droga relatywistyczna, droga newtonowska i droga fotonu\n (4) masy od '
                                       'predkosci\n'))
                    x = np.linspace(1, t, 1000)
                    y = czasteczka.droga(x)[0]
                    z = pr_swiatla * x - pr_swiatla ** 2 / czasteczka.przyspieszenie_ms
                    w = czasteczka.droga(x)[1]
                    n = czasteczka.predkosc_przyspieszanego_obiektu(x)
                    m = czasteczka.masa_relatywistyczna(n)
                    if rodzaj == 1:
                        plt.plot(x, y, color='g', label='droga rel.')
                        plt.plot(x, w, color='r', label='droga new.')
                        plt.legend(('droga rel.', 'droga new.'), loc='upper center', shadow=True)
                        plt.xlabel('czas')
                        plt.ylabel('droga')
                        plt.show()

                    if rodzaj == 2:
                        plt.plot(x, y, color='g', label='droga rel.')
                        plt.plot(x, z, color='b', label='asymptota (przesunieta droga, ktora pokanalby foton')
                        plt.legend(('droga rel.', 'asymptota (przesunieta droga, ktora pokanalby foton'),
                                   loc='upper center', shadow=True)
                        plt.xlabel('czas')
                        plt.ylabel('droga')
                        plt.show()
                    if rodzaj == 3:
                        plt.plot(x, y, color='g', lw=1, ls='-.', label='droga rel.')
                        plt.plot(x, z, color='b', label='asymptota (przesunieta droga, ktora pokanalby foton')
                        plt.plot(x, w, label='droga new.')
                        plt.legend(('droga rel.', 'asymptota (przesunieta droga, ktora pokanalby foton', 'droga new.'),
                                   loc='upper center', shadow=True)
                        plt.xlabel('czas')
                        plt.ylabel('droga')
                        plt.show()
                    if rodzaj == 4:
                        plt.plot(n, m, color='g', label='masa od predkosci')
                        plt.legend('masa od predkosci',
                                   loc='upper center', shadow=True)
                        plt.xlabel('predkosc')
                        plt.ylabel('masa')
                        plt.show()
                elif dalej == 2:
                    l = False
                    break
                elif dalej == 3:
                    p = False
                    break
        elif wybor == 2:
            if troj == 1:
                v1_x = float(input('podaj predkosc pierwszej czastki wzgledem osi x\n'))
                v1_y = float(input('podaj predkosc pierwszej czastki wzgledem osi y\n'))
                v1_z = float(input('podaj predkosc pierwszej czastki wzgledem osi z\n'))
                m1 = float(input('podaj podaj mase spoczynkowa pierwszej czasteczki\n'))
                czasteczka1 = Czastka(typ=typdanych, pr_x=v1_x, pr_y=v1_y, pr_z=v1_z, mas_spo=m1)
                v2_x = float(input('podaj predkosc drugiej czastki wzgledem osi x\n'))
                v2_y = float(input('podaj predkosc drugiej czastki wzgledem osi y\n'))
                v2_z = float(input('podaj predkosc drugiej czastki wzgledem osi z\n'))
                m2 = float(input('podaj podaj mase spoczynkowa pierwszej czasteczki\n'))
                czasteczka2 = Czastka(typ=typdanych, pr_x=v2_x, pr_y=v2_y, pr_z=v2_z, mas_spo=m2)
                wzgl_predkosc, wzgl_szybkosc = czasteczka1.wzgledna_predkosc(czasteczka2)
                print("wzgledna predkosc czastki 1 wzgledem czastki 2 to: ", wzgl_predkosc)
                print("natomiast wzgledna szybkosc to: ", wzgl_szybkosc)
                print("masa czastki 1 wzgledem czastki 2 wynosi: ",
                      czasteczka1.masa_relatywistyczna(czasteczka2.wzgledna_predkosc(czasteczka1)[1]))
                print("masa czastki 2 wzgledem czastki 1 wynosi: ",
                      czasteczka2.masa_relatywistyczna(czasteczka1.wzgledna_predkosc(czasteczka2)[1]))
            else:
                v1 = float(input('podaj predkosc pierwszej czastki czastki\n'))
                m1 = float(input('podaj podaj mase spoczynkowa pierwszej czasteczki\n'))
                czasteczka1 = Czastka(pr_x=v1, typ=typdanych, mas_spo=m1)
                v2 = float(input('podaj predkosc drugiej czastki czastki\n'))
                m2 = float(input('podaj podaj mase spoczynkowa drugiej czasteczki\n'))
                czasteczka2 = Czastka(pr_x=v2, typ=typdanych, mas_spo=m2)
                wzgl_predkosc, wzgl_szybkosc = czasteczka1.wzgledna_predkosc(czasteczka2)
                print("wzgledna szybkosc to: ", wzgl_szybkosc)
                print("masa czastki 1 wzgledem czastki 2 wynosi: ",
                      czasteczka1.masa_relatywistyczna(czasteczka2.wzgledna_predkosc(czasteczka1)[1]))
                print("masa czastki 2 wzgledem czastki 1 wynosi: ",
                      czasteczka2.masa_relatywistyczna(czasteczka1.wzgledna_predkosc(czasteczka2)[1]))
            while True:
                dalej = int(input("Rozpoczac ponownie (1) czy zakonczyc dzialanie programu (2)\n"))
                if dalej == 1:
                    break
                elif dalej == 2:
                    p = False
                    break


menu()
