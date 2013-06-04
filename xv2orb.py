import numpy as np
import struct
import pylab
import math
import os

def xv2orb(x,y,z,vx,vy,vz,mcent):
    """this program converts xv coordinates to orbital elements"""
    tiny=4e-15
    dr = 180.0/math.pi

    #calculate angular momenta
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h2 = hx * hx + hy * hy + hz * hz
    h = np.sqrt(h2)
    inc = np.arccos(hz / h)

    u = np.empty(np.size(x))
    capom = np.empty(np.size(x))
    omega = np.empty(np.size(x))
    capm = np.empty(np.size(x))
    e = np.empty(np.size(x))
    a = np.empty(np.size(x))

    #calculating long. of ascending node
    fac = np.sqrt(hx**2 + hy**2) / h
    sel = np.where(fac < tiny)[0]
    capom[sel] = 0.0
    u[sel] = np.arctan2(y[sel], x[sel])
    sel2 = np.where(np.fabs(inc[sel] - math.pi) < 10.0 * tiny)[0]
    u[sel[sel2]] = -u[sel[sel2]]

    sel = np.where(fac >= tiny)[0]
    capom[sel] = np.arctan2(hx[sel], -hy[sel])
    u[sel] = (np.arctan2((z[sel] / np.sin(inc[sel])), (x[sel] *
        np.cos(capom[sel]) + y[sel] * np.sin(capom[sel]))))

    sel = np.where(capom < 0)[0]
    capom[sel] = capom[sel] + 2 * math.pi
    sel = np.where(u < 0)[0]
    u[sel] = u[sel] + 2 * math.pi

    #calculating energy, etc.
    r = np.sqrt(x * x + y * y + z * z)
    v2 = vx * vx + vy * vy + vz * vz
    v = np.sqrt(v2)
    vdotr = x * vx + y * vy + z * vz
    energy = np.asarray(v2 / 2 - mcent / r)

    #determing conic section and label it with ialpha
    ialpha = np.empty(np.size(x), dtype = int)
    sel = np.where(energy < 0)[0]
    ialpha[sel] = -1
    sel = np.where(energy > 0)[0]
    ialpha[sel] = 1
    sel = np.where(np.fabs(energy * r / mcent) < np.sqrt(tiny))[0]
    ialpha[sel] = 0

    #ellipse
    ellipse = np.where(ialpha == -1)[0]
    a[ellipse] = -mcent[ellipse] / energy[ellipse] / 2
    fac = 1 - h2[ellipse] / mcent[ellipse] / a[ellipse]
    w = np.empty(np.size(fac))
    cape = np.empty(np.size(fac))

    sel = np.where(fac > tiny)[0]
    e[ellipse[sel]] = np.sqrt(fac[sel])
    face = ((a[ellipse[sel]] - r[ellipse[sel]]) / (a[ellipse[sel]] *
        e[ellipse[sel]]))

    sel2 = np.where(face > 1)[0]
    cape[sel[sel2]] = 0.0
    sel2 = np.where((face > -1) & (face <= 1))[0]
    cape[sel[sel2]] = np.arccos(face[sel2])
    sel2 = np.where(face <= -1)[0]
    cape[sel[sel2]] = math.pi

    sel2 = np.where(vdotr[ellipse[sel]] < 0)[0]
    cape[sel[sel2]] = 2 * math.pi - cape[sel[sel2]]

    cw = np.empty(np.size(sel))
    sw = np.empty(np.size(sel))
    cw = ((np.cos(cape[sel]) - e[ellipse[sel]])/ (1 - e[ellipse[sel]] *
        np.cos(cape[sel])))
    sw = (np.sqrt(1 - e[ellipse[sel]]**2) * np.sin(cape[sel]) / (1 -
        e[ellipse[sel]] * np.cos(cape[sel])))
    w[sel] = np.arctan2(sw, cw)
    sel2 = np.where(w[sel] < 0)[0]
    w[sel[sel2]] = w[sel[sel2]] + 2 * math.pi

    #taking care of cases of almost perfectly circular orbits
    sel = np.where(fac <= tiny)[0]
    e[ellipse[sel]] = 0.0
    w[sel] = u[ellipse[sel]]
    cape[sel] = u[ellipse[sel]]

    capm[ellipse] = cape - e[ellipse] * np.sin(cape)
    omega[ellipse] = u[ellipse] - w

    sel = np.where(omega[ellipse] < 0)[0]
    omega[ellipse[sel]] = omega[ellipse[sel]] + 2 * math.pi

    omega[ellipse] = (omega[ellipse] - np.floor(omega[ellipse] / 2 / math.pi) *
        2 * math.pi)

    #hyperbola
    hyper = np.where(ialpha == 1)[0]
    a[hyper] = mcent[hyper] / energy[hyper] / 2
    fac = h2[hyper] / (mcent[hyper] * a[hyper])
    w = np.empty(np.size(fac))
    tmpf = np.empty(np.size(fac))
    capf = np.empty(np.size(fac))

    sel = np.where(fac > tiny)[0]
    e[hyper[sel]] = np.sqrt(1 + fac[sel])
    tmpf[sel] = ((a[hyper[sel]] + r[hyper[sel]]) / (a[hyper[sel]] *
        e[hyper[sel]]))

    sel2 = np.where(tmpf[sel] < 1)[0]
    tmpf[sel[sel2]] = 1

    capf[sel] = np.log10(tmpf[sel] + np.sqrt(tmpf[sel]**2 - 1))

    sel2 = np.where(vdotr[hyper[sel]] < 0)[0]
    capf[sel[sel2]] = -capf[sel[sel2]]

    cw = ((e[hyper[sel]] - np.cosh(capf[sel])) / (e[hyper[sel]] *
        np.cosh(capf[sel])))
    sw = ((e[hyper[sel]]**2 - 1) * np.sinh(capf[sel]) / (e[hyper[sel]] *
        np.cosh(capf[sel])))
    w[sel] = np.arctan2(sw, cw)

    sel2 = np.where(w[sel] < 0)[0]
    w[sel[sel2]] = w[sel[sel2]] + 2 * math.pi

    #taking care of near-parabola cases
    sel = np.where(fac <= tiny)[0]
    e[hyper[sel]] = 1
    tmpf[sel] = h2[hyper[sel]] / mcent[hyper[sel]] / 2
    w[sel] = np.arccos(2 * tmpf[sel] / r[hyper[sel]] -1)

    sel2 = np.where(vdotr[hyper[sel]] < 0)[0]
    w[sel[sel2]] = 2 * math.pi - w[sel[sel2]]

    tmpf[sel] = ((a[hyper[sel]] + r[hyper[sel]]) / (a[hyper[sel]] *
        e[hyper[sel]]))
    capf[sel] = np.log10(tmpf[sel] + np.sqrt(tmpf[sel]**2 - 1))

    capm[hyper] = e[hyper] * np.sinh(capf) - capf
    omega[hyper] = u[hyper] - w

    sel = np.where(omega[hyper] < 0)[0]
    omega[hyper[sel]] = omega[hyper[sel]] + 2 * math.pi
    omega[hyper] = (omega[hyper] - np.floor(omega[hyper] / 2 / math.pi) * 2
        * math.pi)

    #parabola
    para = np.where(ialpha == 0)[0]
    a[para] = h2[para] / mcent[para] / 2
    e[para] = 1
    w = np.arccos(2 * a[para] / r[para] - 1)

    sel = np.where(vdotr[para] < 0)[0]
    w[sel] = 2 * math.pi - w[sel]

    tmpf = np.tan(w / 2)
    capm[para] = tmpf * (1 + tmpf**2 / 3)
    omega[para] = u[para] - w

    sel = np.where(omega[para] < 0)[0]
    omega[para[sel]] = omega[para[sel]] + 2 * math.pi

    omega[para] = (omega[para] - np.floor(omega[para] / 2 / math.pi) * 2
        * math.pi)

    return(a,e,inc,capom,omega,capm)
