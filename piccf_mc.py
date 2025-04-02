#!/usr/bin/env python

"""
A interface to piccf_mc.c. It's used for PICCF calculation and
Monte-Carlo simulations.
"""

import sys
import math
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing

LIB_PATH = "./libpiccf_mc.so"

def piccf(jdc, fc, jdl, fl, nt, tbeg, tend):
    """
    Calculate PICCF of two light curves.  jdc, fc are the arrays of continuum
    lc.  jdl, fl are the arrays of line lc. nt is the length of output time lag
    array, tbeg is the beginning of time lag array and tend is the end of that
    array.
    """
    cnc = ctypes.c_int(len(jdc))
    cjdc = (ctypes.c_double * len(jdc))()
    cfc = (ctypes.c_double * len(jdc))()

    for i in range(len(jdc)):
        cjdc[i] = jdc[i]
        cfc[i] = fc[i]

    cnl = ctypes.c_int(len(jdl))
    cjdl = (ctypes.c_double * len(jdl))()
    cfl = (ctypes.c_double * len(jdl))()

    for i in range(len(jdl)):
        cjdl[i] = jdl[i]
        cfl[i] = fl[i]

    cnt = ctypes.c_int(nt)
    ctbeg = ctypes.c_double(tbeg)
    ctend = ctypes.c_double(tend)
    ct = (ctypes.c_double * nt)()
    cr = (ctypes.c_double * nt)()
    crmax = ctypes.c_double(0.0)
    ctau_cent = ctypes.c_double(0.0)
    ctau_peak = ctypes.c_double(0.0)

    cdll = ctypes.CDLL(LIB_PATH)
    cdll.piccf(cnc, cjdc, cfc, cnl, cjdl, cfl, cnt, ctbeg, 
            ctend, ct, cr, ctypes.byref(crmax), ctypes.byref(ctau_cent), ctypes.byref(ctau_peak))

    t = []
    r = []
    for i in range(nt):
        t.append(ct[i])
        r.append(cr[i])
    rmax = crmax.value
    tau_cent = ctau_cent.value
    tau_peak = ctau_peak.value

    (t, r) = np.array((t, r))

    return t, r, rmax, tau_cent, tau_peak

def piccf_mp_single(p):
    t, r, rmax, tau_cent, tau_peak = piccf(p[0], p[1], p[2], p[3], p[4], p[5], p[6])
    return t, r, rmax, tau_cent, tau_peak

def piccf_mp(p, num_process):
    """p is an array with multiple combination of [jdc, fc, jdl, fl, nt, tbeg, tend]"""

    process_pool = multiprocessing.Pool(num_process)
    result = process_pool.map(piccf_mp_single, p)
    process_pool.close()
    process_pool.join()
    process_pool.terminate()
    return result

def piccf_mc_mp_single(p):
    tau_cent_tot, tau_peak_tot = piccf_mc(p[0], p[1], p[2], p[3], p[4], p[5], p[6], \
                                              p[7], p[8], p[9])
    return tau_cent_tot, tau_peak_tot

def piccf_mc_mp(p, num_process):
    """p is an array with multiple combination of [jdc, fc, efc, jdl, fl, efl, nt, tbeg, tend, num_mc]"""

    process_pool = multiprocessing.Pool(num_process)
    result = process_pool.map(piccf_mc_mp_single, p)
    process_pool.close()
    process_pool.join()
    process_pool.terminate()
    return result

def piccf_mc(jdc, fc, efc, jdl, fl, efl, nt, tbeg, tend, num_mc):
    """
    Calculate PICCF of two light curves and apply Monte-Carlo simulation to
    obtain the error bar of time lag. fc and efc are the array of continuum lc
    and the error bars, fl and efl are the array of line lc and the error bars.
    nt is the length of output time lag array, tbeg is the beginning of time
    lag array and tend is the end of that array. num_mc is the number of
    Monte-Carlo simulation.
    """
    cnc = ctypes.c_int(len(jdc))
    cjdc = (ctypes.c_double * len(jdc))()
    cfc = (ctypes.c_double * len(jdc))()
    cefc = (ctypes.c_double * len(jdc))()

    for i in range(len(jdc)):
        cjdc[i] = jdc[i]
        cfc[i] = fc[i]
        cefc[i] = efc[i]

    cnl = ctypes.c_int(len(jdl))
    cjdl = (ctypes.c_double * len(jdl))()
    cfl = (ctypes.c_double * len(jdl))()
    cefl = (ctypes.c_double * len(jdl))()

    for i in range(len(jdl)):
        cjdl[i] = jdl[i]
        cfl[i] = fl[i]
        cefl[i] = efl[i]

    cnt = ctypes.c_int(nt)
    ctbeg = ctypes.c_double(tbeg)
    ctend = ctypes.c_double(tend)
    cnum_mc = ctypes.c_int(num_mc)
    ctau_cent_tot = (ctypes.c_double * num_mc)()
    ctau_peak_tot = (ctypes.c_double * num_mc)()

    cdll = ctypes.CDLL(LIB_PATH)
    cdll.piccf_mc(cnc, cjdc, cfc, cefc, cnl, cjdl, cfl, cefl, cnt, 
            ctbeg, ctend, ctypes.byref(cnum_mc), ctau_cent_tot, ctau_peak_tot)

    tau_cent_tot = []
    tau_peak_tot = []
    for i in range(num_mc):
        tau_cent_tot.append(ctau_cent_tot[i])
        tau_peak_tot.append(ctau_peak_tot[i])

    (tau_cent_tot, tau_peak_tot) = np.array((tau_cent_tot, tau_peak_tot))

    return tau_cent_tot, tau_peak_tot

def main():

    if (len(sys.argv) == 2):
        if sys.argv[1] == "-h":
            print("usage: python3 piccf_mc.py con_name line_name tbeg tend nt nbin num_mc")
            print("       or")
            print("       python3 piccf_mc.py param")
        else:
            paramname = sys.argv[1]
            f = open(sys.argv[1])
            l = f.readlines()
            continuum = l[0].split()[1]
            line = l[1].split()[1]
            tbeg = float(l[2].split()[1])
            tend = float(l[3].split()[1])
            nt = int(l[4].split()[1])
            nbin = int(l[5].split()[1])
            num_mc = int(l[6].split()[1])
    elif len(sys.argv) == 8:
        continuum = sys.argv[1]
        line = sys.argv[2]
        tbeg = float(sys.argv[3])
        tend = float(sys.argv[4])
        nt = int(sys.argv[5])
        nbin = int(sys.argv[6])
        num_mc = int(sys.argv[7])
    else:
        print("usage: python3 piccf_mc.py con_name line_name tbeg tend nt nbin num_mc")
        print("       or")
        print("       python3 piccf_mc.py param")
        sys.exit()
        

    f = open(continuum)
    l = f.readlines()
    f.close()
    jdc = [float(i.split()[0]) for i in l]
    fc = [float(i.split()[1]) for i in l]
    efc = [float(i.split()[2]) for i in l]
    (jdc, fc, efc) = np.array((jdc, fc, efc))

    f = open(line)
    l = f.readlines()
    f.close()
    jdl = [float(i.split()[0]) for i in l]
    fl = [float(i.split()[1]) for i in l]
    efl = [float(i.split()[2]) for i in l]
    (jdl, fl, efl) = np.array((jdl, fl, efl))

    unitc = 10 ** int(round(np.log10(np.mean(fc))))
    unitl = 10 ** int(round(np.log10(np.mean(fl))))

    fc = fc / unitc
    efc = efc / unitc

    fl = fl / unitl
    efl = efl / unitl

    jd0 = int(round(np.min([np.min(jdc), np.min(jdl)])))
    
    jdc = jdc - jd0
    jdl = jdl - jd0

    at, ar, armax, atau_cent, atau_peak = piccf(jdc, fc, jdc, fc, nt, tbeg, tend)
    print('ACF rmax:%.6f tau_cent:%.6f tau_peak:%.6f' % (armax, atau_cent, atau_peak))

    t, r, rmax, tau_cent, tau_peak = piccf(jdc, fc, jdl, fl, nt, tbeg, tend)
    print('CCF rmax:%.6f tau_cent:%.6f tau_peak:%.6f' % (rmax, tau_cent, tau_peak))

    output = open('pypiccf_ccf_out.txt', 'w')
    output.write('t  r\n')
    for i in range(len(t)):
        text = '%.6f %.6f\n' % (t[i], r[i])
        output.write(text)
    output.close()

    tau_cent_mc, tau_peak_mc = piccf_mc(jdc, fc, efc, jdl, fl, efl, nt, tbeg, tend, num_mc)
    print('CCCD from FR/RSS:%.6f %.6f %.6f' % (np.percentile(tau_cent_mc, 15.85), 
            np.percentile(tau_cent_mc, 50.0), np.percentile(tau_cent_mc, 84.15)))
    print('CCPD from FR/RSS:%.6f %.6f %.6f' % (np.percentile(tau_peak_mc, 15.85), 
            np.percentile(tau_peak_mc, 50.0), np.percentile(tau_peak_mc, 84.15)))

    print('time lag (centroid): %.6f - %.6f + %.6f' % (tau_cent, tau_cent - np.percentile(tau_cent_mc, 15.85), 
            np.percentile(tau_cent_mc, 84.15) - tau_cent))
    print('time lag (peak): %.6f - %.6f + %.6f' % (tau_peak, tau_peak - np.percentile(tau_peak_mc, 15.85), 
            np.percentile(tau_peak_mc, 84.15) - tau_peak))

    output = open('pypiccf_mc_out.txt', 'w')
    output.write('tau_cent  tau_peak\n')
    for i in range(len(tau_cent_mc)):
        text = '%.6f %.6f\n' % (tau_cent_mc[i], tau_peak_mc[i])
        output.write(text)
    output.close()

    # ==== plot

    mpl.rcParams['text.usetex'] = 'True'
    mpl.rcParams['axes.linewidth'] = 1.1
    mpl.rcParams['font.size'] = 12.0
    mpl.rcParams['xtick.major.pad'] = 6
    mpl.rcParams['xtick.minor.pad'] = 6
    mpl.rcParams['ytick.major.pad'] = 4
    mpl.rcParams['ytick.minor.pad'] = 4

    fig = plt.figure(figsize = (10, 6))
    
    ax1 = fig.add_axes([0.1, 0.55, 0.5, 0.4])
    ax1.errorbar(jdc, fc, efc, fmt = 'o')
    ax1.set_ylabel(r'$F_{\mathrm{con}}\;\mathrm{(\times 10^{%i})}$' % (np.log10(unitc)))
    [i.set_visible(False) for i in ax1.get_xticklabels()]
    ax1.minorticks_on()

    ax2 = fig.add_axes([0.1, 0.12, 0.5, 0.4])
    ax2.errorbar(jdl, fl, efl, fmt = 'o')
    ax2.set_ylabel(r'$F_{\mathrm{line}}\;\mathrm{(\times 10^{%i})}$' % (np.log10(unitl)))
    ax2.set_xlabel(r'$\mathrm{JD\ -\ %i\ (days)}$' % (jd0))
    ax2.minorticks_on()

    ax3 = fig.add_axes([0.62, 0.55, 0.22, 0.4])
    ax3.errorbar(at, ar)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position('right')
    ax3.minorticks_on()
    ax3.set_xlim(at[0], at[-1])
    [i.set_visible(False) for i in ax3.get_xticklabels()]
    ax3.set_ylabel(r'$\mathrm{ACF}$')

    ax4 = fig.add_axes([0.62, 0.12, 0.22, 0.4], sharex = ax3)
    ax4.errorbar(t, r)
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position('right')
    ax4.minorticks_on()
    ax4.set_xlim(t[0], t[-1])
    ax4.set_ylabel(r'$\mathrm{CCF}$')
    ax4.set_xlabel(r'$\mathrm{Time\ lag\ (days)}$')

    ylim1, ylim2 = ax4.get_ylim()

    counts, bins = np.histogram(tau_cent_mc, bins = nbin, range = [tbeg, tend])
    ax4.stairs(counts / np.max(counts) * (ylim2 - ylim1) * 0.9 + ylim1, bins, fill = True, baseline = ylim1, alpha = 0.8, label = r"$\mathrm{CCCD}$")

    counts, bins = np.histogram(tau_peak_mc, bins = nbin, range = [tbeg, tend])
    ax4.stairs(counts / np.max(counts) * (ylim2 - ylim1) * 0.9 + ylim1, bins, fill = True, baseline = ylim1, alpha = 0.8, label = r"$\mathrm{CCPD}$")

    ax4.legend(loc = "upper right")

    plt.show()
    fig.savefig("CCF_result.pdf", format = "pdf", bbox_inches = "tight")

    

if __name__ == "__main__":
    main()
