from astropy.wcs import WCS
import casa_cube as casa
from matplotlib.patches import Ellipse
import matplotlib.colors as colors
from matplotlib.patches import Circle



cont = casa.Cube("/Users/cpinte/Downloads/Elias24_continuum.fits")
line = casa.Cube("/Users/cpinte/Downloads/Elias24_CO.fits")

#line = casa.Cube("/Users/cpinte/Observations/DSHARP/ELias24/Elias24_COcube_contsub_robust-1.0_uvtaper0.05.fits")
#line = casa.Cube("/Users/cpinte/Observations/DSHARP/ELias24/Elias24_COcube_contsub_robust0.5_uvtaper0.05.fits")
#line = casa.Cube("/Users/cpinte/Observations/DSHARP/ELias24/Elias24_COcube_contsub_robust-1.5_uvtaper0.07.fits")

#print(np.sqrt(line.beam[0] * line.beam[1])) # 0.07462732423289271


# Dec, RA index in that order !!!
#imax = np.unravel_index(cont.image.argmax(), cont.image.shape)

#wcs_cont = WCS(cont.header)
#sky  = wcs.pixel_to_world(imax[0], imax[1], 0, 0)

#wcs_line = WCS(line.header)
#wcs_line[:,:,0,0].world_to_array_index(sky)


# star in the lines in at pixel 740,5, 712.5 in x, y
# image is 1500 by 1500
# blob is at 773, 721

star = np.array([740.5,712.5])
blob = np.array([773,721])


PA = 360 - (90 - np.arctan2((blob-star)[1],(blob-star)[0]) *180./np.pi)
# -75  == 285
pixel_scale = 0.01


dist = np.hypot((star-blob)[0],(star-blob)[1]) * pixel_scale

star_shift_dx = 0.105
star_shift_dy = -0.383

limit = 1.0


planet_dist = 0.411
planet_PA = (302.1 + 90) / 180. * np.pi

x_planet = -planet_dist * np.cos(planet_PA)
y_planet = planet_dist * np.sin(planet_PA)

x_planet= -0.335
y_planet = 0.105



plt.plot(x_planet, y_planet ,"o",color="cyan",ms=2)


#--------- plot continuum
plt.figure(1)
ax = plt.gca()
cont.plot(limit=limit,Tb=True,fmax=90,fmin=0,shift_dx=star_shift_dx,shift_dy = star_shift_dy, ax=ax)
ax.plot(0,0,"*",color="white",ms=2)
ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)



#--------- plot lines
win=2
fig_xfactor = 5
fig_yfactor = 5

n_channels = 6
iv0=18

fig2, axes2 = plt.subplots(nrows=1, ncols=n_channels, figsize=(fig_xfactor*n_channels,fig_yfactor*1), num=win)

for i in range(n_channels):
    ax = axes2[i]
    line.plot(limit=limit,Tb=True,fmax=55,fmin=5,shift_dx=star_shift_dx,shift_dy = star_shift_dy, iv=iv0 + i, ax=ax)
    ax.plot(0,0,"*",color="white",ms=4)
    ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)


#--- moment map

v0 = 0.5
dv = 0.8
n_moments = 2

win=3
fig3, axes3 = plt.subplots(nrows=1, ncols=n_moments+1, figsize=((1.4*fig_xfactor)*(n_moments+1),fig_yfactor*1), num=win)


ax = axes3[0]
cont.plot(limit=limit,Tb=True,fmax=90,fmin=0,shift_dx=star_shift_dx,shift_dy = star_shift_dy, ax=ax, cmap="viridis")
ax.plot(0,0,"*",color="white",ms=2)
ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)


ax = axes3[1]
line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], no_ylabel = True)

ax.plot(0,0,"*",color="white",ms=4)
#ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
circle = Circle((x_planet, y_planet), 0.1, clip_on=False, zorder=10, linewidth=2, edgecolor='white', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
ax.add_artist(circle)

ax = axes3[2]
line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=9, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], M0_threshold=5e-3, fmin=-dv/2+0.1,fmax=dv/2+0.1, no_ylabel = True)

ax.plot(0,0,"*",color="black",ms=4)
#ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
circle = Circle((x_planet, y_planet), 0.1, clip_on=False, zorder=10, linewidth=2, edgecolor='black', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
ax.add_artist(circle)

plt.tight_layout()

plt.show()


#---- PV diagram

from pvextractor import Path


path = Path([(691, 765), (800., 660.)], width=3)

path = Path([ (705, 755), (800, 655)], width=3)

center = (739, 712)
center = (745, 710)
center = (740.5, 712.5)
disk_PA = np.deg2rad(48)
dx=150
path = Path([ (center[0] - dx * np.cos(disk_PA), center[1] + dx * np.sin(disk_PA)), (center[0] + dx * np.cos(disk_PA), center[1] - dx * np. sin(disk_PA))], width=3)


slice = extract_pv_slice("/Users/cpinte/Downloads/Elias24_CO.fits", path)


plt.figure(13)
plt.clf()
halfsize = slice.shape[1] / 2 * line.pixelscale
extent = [-halfsize,halfsize, np.min(line.velocity), np.max(line.velocity)]

norm = colors.Normalize(vmin=0., vmax=np.max(slice.data), clip=True)

plt.imshow(slice.data, origin='lower', extent=extent, aspect="auto", norm=norm, cmap="inferno")

def v_kepler(Mstar,r):
    from scipy import constants as sc
    GxMsun   = 1.3271244e20

    return np.sign(r) * np.sqrt(GxMsun * Mstar  / (np.maximum(np.abs(r),1e-3) * sc.au))

distance = 136

inc = np.deg2rad(28.5)
v_syst = 2.85

Mstar = 1.5
dv = 0.25

r = np.linspace(0.01,1,100)
plt.plot(r , -v_kepler(Mstar,r*distance)*1e-3 * np.sin(inc) + v_syst, color="white")
plt.plot(r , -v_kepler(Mstar+dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")
plt.plot(r , -v_kepler(Mstar-dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")

r = -r
plt.plot(r , -v_kepler(Mstar,r*distance)*1e-3 * np.sin(inc) + v_syst, color="white")
plt.plot(r , -v_kepler(Mstar+dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")
plt.plot(r , -v_kepler(Mstar-dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")

plt.xlim(-0.6,0.6)
plt.ylim(-4,12)
plt.xlabel(r'Distance along the midplane ["]')
plt.ylabel(r"Velocity [km/s]")

#plt.plot([np.max(np.abs(r)),-np.max(np.abs(r))], [v_syst,v_syst], linestyle=":", color="white")

beam_position=(0.125, 0.125)
dx, dy = beam_position
bpa = 0.
bmaj = 0.7 # km/s
bmin = line.bmaj
ax = plt.gca()
beam = Ellipse(
    ax.transLimits.inverted().transform((dx, dy)),
    width=bmin,
    height=bmaj,
    angle=-bpa,
    fill=False,
    color="white",
)
ax.add_patch(beam)




model = mcfost.Line("/Users/cpinte/mcfost+phantom/Elias24/5Mj/data_CO")

win=5
fig4, axes4 = plt.subplots(nrows=1, ncols=n_moments, figsize=((1.4*fig_xfactor)*n_moments,fig_yfactor*1), num=win)

iv_support = [107,108,109,110,111]

ax = axes4[0]
model.plot_map(iv_support=iv_support,psf_FWHM=0.075,plot_stars=True,Tb=True,moment=0,ax=ax)

ax = axes4[1]
model.plot_map(iv_support=iv_support, psf_FWHM=0.075,plot_stars=True,moment=1,ax=ax, fmin=1,fmax=1.5)
