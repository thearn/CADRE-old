import pickle
import numpy as np
import pylab
import pygmaps

"""
Must install pygamps-extended first: https://github.com/thearn/pygmaps-extended
"""

data = pickle.load(open("src/CADRE/test/data1346.pkl"))
savedir = "docs/maps"


def val2hex(v, scale=1.0):
    if not v:
        return "#000044"
    j = pylab.get_cmap("jet")
    v = v / scale
    nums = [int(255 * i) for i in j(v)][:3]
    return ''.join(["00" if i == 0 else hex(i)[2:] for i in nums])


def calc_lat_lon(r_e2b_I, O_IE):
    r2d = 180 / np.pi
    rgs = 6378.137
    lats, lons = [], []
    n = r_e2b_I.shape[1]
    for i in xrange(n):

        r_e2g_I = r_e2b_I[:3, i]
        r_e2g_I = r_e2g_I / np.linalg.norm(r_e2g_I) * rgs
        r_e2g_E = np.dot(O_IE[:, :, i].T, r_e2g_I)

        lat = np.arcsin(r_e2g_E[2] / rgs) * r2d
        lon = np.arctan2(r_e2g_E[1], r_e2g_E[0]) * r2d

        lats.append(lat), lons.append(lon)

    return lats, lons

mxdata = max([max(data["%s:Dr" % str(i)]) for i in xrange(6)])

gmap_all = pygmaps.gmap(41, -88, 2)

#[0,0,1] = 90,0
#[0,1,0] = 0,90
#[1,0,0] = 0,0
#[-1, 0, 0] = 0,180

for i in xrange(6):
    si = str(i)
    dr = data[si + ":Dr"]
    p = data[si + ":P_comm"]
    g = data[si + ":gamma"]
    s = data[si + ":SOC"]
    pylab.figure()
    pylab.suptitle("CADRE Design Point " + si)
    pylab.subplot(411)
    pylab.title("Dr")
    pylab.plot(dr)
    pylab.gca().get_xaxis().set_visible(False)

    pylab.subplot(412)
    pylab.title("P_comm")
    pylab.plot(p)
    pylab.gca().get_xaxis().set_visible(False)

    pylab.subplot(413)
    pylab.title("gamma")
    pylab.plot(g)
    pylab.gca().get_xaxis().set_visible(False)

    pylab.subplot(414)
    pylab.title("SOC")
    pylab.plot(s[0])

    pylab.gcf().savefig(savedir + "/" + si + '.png', bbox_inches='tight')

    O_IE = data["%s:O_IE" % si]

    gmap = pygmaps.gmap(41, -88, 2)
    r_e2b_I = data["%s:r_e2b_I" % si]

    lats, lons = calc_lat_lon(r_e2b_I, O_IE)

    path = zip(lats, lons)
    gmap.add_weighted_path(path, dr, scale=mxdata)
    gmap_all.add_weighted_path(path, dr, scale=mxdata)
    gmap.draw(savedir + "/" + si + '_data.html')
gmap_all.draw(savedir + '/all_data.html')
