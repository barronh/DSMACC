from . import model
if __name__ == '__main__':
    mod = model()
    # Use default environment
    mod.run(2012184.0, 8., 180, conc_ppb = dict(O3 = 30, C5H8 = 850., NO = 128), globvar = dict(LAT = 30., LON = 0))
