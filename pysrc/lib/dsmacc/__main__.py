from . import model
import argparse


class dynenv(model):
    def custom_updateenv(self):
        return


mod = model()
parser = argparse.ArgumentParser()
parser.add_argument(
    'JDAY', type=float,
    help='YYYYJJJ.FFFFFF where FFF is fractional day in GMT'
)

parser.add_argument('RUNHOURS', type=int, help='hours to run')
parser.add_argument('DT', type=int, help='seconds between rconst updates')
parser.add_argument(
    '--globvar', type=eval,
    default=dict(LAT=30., LON=0, TEMP=298., PRESS=101325.),
    help=(
        'python dictionary of global variables (default='
        + '"dict(LAT=30., LON=0, TEMP=298., PRESS=101325.)")'
    )
)
parser.add_argument(
    'conc_ppb', type=eval,
    help=(
        'python dictionary of initial concentrations (e.g., "dict(O3=30., '
        + 'C5H8=850., NO=128.)")'
    )
)

if __name__ == '__main__':
    args = parser.parse_args()
    # Use default environment
    mod.run(
        args.JDAY, args.RUNHOURS, args.DT, conc_ppb=args.conc_ppb,
        globvar=args.globvar
    )
