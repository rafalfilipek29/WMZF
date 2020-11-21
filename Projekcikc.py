
class Czastka():
    predkosc=None
    masa_spoczynkowa=None
    masa_relatywistyczna=None


    def __init__(self,pr,mas_spo):
        self.predkosc=pr
        self.masa_spoczynkowa=mas_spo
        self.masa_relatywistyczna=self.masa_spoczynkowa*self.predkosc
        #self.energia_kinetyczna=(1/2)*self.masa*self.predkosc**2

    


