from dsmacc import model

class dynenv(model):
    def custom_updateenv(self):
        return

mod = dynenv()
mod.run(2000172.25, 8, 180, conc_ppb = dict(O3 = 30., C5H8 = 850., NO = 128.), globvar = dict(LAT = 30., LON = 0, TEMP = 298., PRESS = 101325.))
