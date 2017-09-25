#!/usr/bin/env python

"""
Provides plotting functionality for DSMACC output.

run with -h for help

$ python plot.py -h

or

$ ./plot.py -h

an example usage:

python plot.py jingqiu.Spec_1.dat.timeseries -v O3 -v "(NO+NO2)*10;color='r'" --factor="1.e9/CFACTOR" --unit='ppb' --axes="[.1,.2,.8,.7]" --plotall --title="O3 and NOx: on a similar scale"

"""

from argparse import ArgumentParser
import numpy as np
parser = ArgumentParser(description = 'DSMACC Plotting Script to make individual species or combined species plots')

parser.add_argument('ifile', type = str, help='path to a file formatted as type -f')
parser.add_argument("-v", "--vardesc", dest = "variables", default = [], action = "append", metavar = '-v vardesc1[ -v vardesc[... -v vardesc]',
                        help = "Variable description; made of an expression and optional keywords expr[;key1=val1,[key2=val2]]; e.g., -v NO or -v \"NO+NO2;label=NOx,color='g'\"")


parser.add_argument("-d", "--sdate", dest = "startdate", type = str, default = '2000-01-01T00:00:00', metavar = 'YYYY-mm-ddTHH:MM:SS', help = "Reference date for data start.")

parser.add_argument("--backend", dest = "backend", type = str, default = 'Agg', metavar = 'BACKEND', help = "Agg; TkAgg; or PDF")

parser.add_argument("--format", dest = "format", type = str, default = 'png', metavar = 'FMT', help = "png; pdf; ps; etc")

parser.add_argument("--prefix", dest = "prefix", type = str, default = None, metavar = 'PREFIX', help = "Figure output start.")

parser.add_argument("-s", "--slice", dest = "slice", type = str, default = ':', metavar = 'start[,stop[,step]]',
                    help = "Variables have only one dimension at this time (time,), which can be subset using dim,start,stop,stride (e.g., --slice=0,47,5 would sample every fifth layer starting at 0)")

parser.add_argument("--title", dest = "axtitle", type = str, default = "%(name)s: min/max=%(min).1e/%(max).1e", metavar = "TEMPLATE", help = 'String template for title')

parser.add_argument("--varlabel", dest = "varlabel", type = str, default = "%(name)s", metavar = "TEMPLATE", help = 'String template for var label')

parser.add_argument("--unit", dest = "unit", type = str, default = "ppb", metavar = "UNIT", help = 'Unit name')

parser.add_argument("--xlabel", dest = "xlabel", type = str, default = "time", metavar = "STR", help = 'X Label')

parser.add_argument("--factor", dest = "factor", type = str, default = "1.e9/CFACTOR", help = '1 or CFACTOR; var * 1 or var * CFACTOR')

parser.add_argument("--axes", dest = "axes", type = str, default = "[.1, .15, .8, .75]", help = 'See pylab')

parser.add_argument("--xtickformat", dest = "xtickformat", type = str, default = "%jT%H", help = 'See pylab')

parser.add_argument("--xtickrotation", dest = "xtickrotation", type = float, default = 45, help = 'See pylab')

parser.add_argument("--plotall", dest = "plotall", default = False, action = "store_true", help = 'Each variable should have its own plot')

parser.add_argument("--figprops", dest = "figprops", default = "", type = str, help = 'Key words string for figure')

parser.add_argument("--axprops", dest = "axprops", default = "", type = str, help = 'Key words string for axes')

args = parser.parse_args()
data = np.recfromtxt(args.ifile, delimiter = '!', dtype = 'd', names = True)
data = dict([(k, data[k]) for k in data.dtype.names])

from matplotlib import use
use(args.backend)
from pylab import *
if args.prefix is None:
    args.prefix = args.ifile + '_'

from datetime import datetime, timedelta
if not 'JDAY_GMT' in data.keys():
    sdate = datetime.strptime(args.startdate, '%Y-%m-%dT%H:%M:%S')
    dates = [sdate + timedelta(seconds = float(x)) for x in data['TIME']]
else:
    dates = [datetime.strptime('%.0f' % (j//1), '%Y%j') + timedelta(hours = (j % 1) * 24) for j in data['JDAY_GMT']]
    
CFACTOR = data['CFACTOR']
if len(args.variables) == 0:
    args.variables.extend([k for k in data.keys() if k in ('O3', 'NO2')])
if len(args.variables) == 0:
    args.variables.extend([k for k in data.keys() if k != 'TIME'])


if args.plotall:
    fig = figure()
    figprops = eval('dict(%s)' % args.figprops)
    if len(figprops) > 0: setp(fig, **figprops)
    exec('ax = fig.add_axes(%s)' % args.axes)
    axprops = eval('dict(%s)' % args.axprops)
    if len(axprops) > 0: setp(ax, **axprops)
for vardesc in args.variables:
    pieces = vardesc.split(';')
    name = pieces[0]
    if len(pieces) == 1:
        properties = {}
    elif len(pieces) == 2:
        properties = eval('dict(%s)' % pieces[1])
    else:
        raise ValueError('vardesc: expr[;key1=val1;key2=val2]')
    properties.setdefault('label', args.varlabel % globals())
    var = eval(name, None, data)
    plot_date = eval('dates[%s]' % args.slice)
    plot_var = eval('var[%s] * %s' % (args.slice, args.factor))
    min = plot_var.min()
    max = plot_var.max()
    std = plot_var.std()
    if not args.plotall:
        fig = figure()
        figprops = eval('dict(%s)' % args.figprops)
        if len(figprops) > 0: setp(fig, **figprops)
        exec('ax = fig.add_axes(%s)' % args.axes)
        axprops = eval('dict(%s)' % args.axprops)
        if len(axprops) > 0: setp(ax, **axprops)
    ax.set_title(args.axtitle % globals())
    ax.set_ylabel(args.unit)
    ax.set_xlabel(args.xlabel)
    ax.xaxis.set_major_formatter(DateFormatter(args.xtickformat))
    setp(ax.get_xticklabels(), rotation = args.xtickrotation)
    ax.plot(plot_date, plot_var, **properties)
    if not args.plotall:
        fig.savefig(args.prefix + name + '.' + args.format)
        close(fig)
if args.plotall:
    legend()
    fig.savefig(args.prefix + 'ALL' + '.' + args.format)
   