from astropy.wcs import WCS
import casa_cube as casa
from matplotlib.patches import Ellipse
import matplotlib.colors as colors
from matplotlib.patches import Circle
import dynamite

import cmasher as cmr

cmap = cmr.arctic
NACO_cmap = cmr.dusk
#NACO_cmap = cmr.amber


dir = "/Users/cpinte/Observations/DSHARP/Elias24/"

cont = casa.Cube(dir+"Elias24_continuum.fits")
line = casa.Cube(dir+"Elias24_CO.fits")

levels = np.arange(30,90,13)

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

#x_planet= -0.335
#y_planet = 0.105

x_planet= -0.33
y_planet = 0.17

x_planet_CO= -0.33
y_planet_CO = 0.12

# linewidth plot ? for instance comparison of 2 pixels at same radius
#plt.plot(line.velocity, line.image[:,720,772]) # planet
#plt.plot(line.velocity, line.image[:,710,771]) # nearby



fwhm = fits.getdata('NACO/fwhm.fits')[0]
NACO_pixelscale = 0.027208 # arcsec/pixel
#NACO = casa.Cube("NACO/final_PCA-ADI_ann_npc1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19-20.fits",pixelscale=NACO_pixelscale, unit="adu") # wrong rotation
#NACO = casa.Cube("NACO/simple_cADI.fits",pixelscale=NACO_pixelscale, unit="adu")
NACO = casa.Cube("NACO/PCA-annular_npc1-50_delta_rot1.fits",pixelscale=NACO_pixelscale, unit="adu")

#NACO = casa.Cube("NACO/pca_annular_res_npc19.fits",pixelscale=NACO_pixelscale, unit="adu")
#NACO.image[np.abs(NACO.image) < 1e-6] = np.nan


NACO.bmaj = fwhm * NACO_pixelscale
NACO.bmin = fwhm * NACO_pixelscale

NACO_empty = casa.Cube('NACO/pca_ann_empty_rerun.fits',pixelscale=NACO_pixelscale, unit="adu")
NACO_empty.bmaj = fwhm * NACO_pixelscale
NACO_empty.bmin = fwhm * NACO_pixelscale

r_mask = 2*fwhm * NACO_pixelscale


mcmc = fits.getdata('NACO/mcmc_results_gauss.fits')
r = mcmc[0][0]
theta = mcmc[1][0]
flux = mcmc[2][0]

#x_planet= - r * NACO_pixelscale * np.cos(np.deg2rad(theta))
#y_planet = r * NACO_pixelscale * np.sin(np.deg2rad(theta))

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



#--------- plot continuum
if False:
    plt.figure(1)
    ax = plt.gca()
    cont.plot(limit=limit,Tb=True,fmax=90,fmin=0,shift_dx=star_shift_dx,shift_dy = star_shift_dy, ax=ax)
    ax.plot(0,0,"*",color="white",ms=2)
    ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)

# -0.35, 0.11

#--------- plot channel maps on 1 line
#win=2
#fig_xfactor = 3
#fig_yfactor = 3
#
#n_channels = 6
#iv0=18
#
#fig2, axes2 = plt.subplots(nrows=1, ncols=n_channels, figsize=(fig_xfactor*n_channels,fig_yfactor*1), num=win, sharex='all', sharey='all')
#
#limit = 0.5
#
#for i in range(n_channels):
#    ax = axes2[i]
#
#    if i==n_channels-1:
#        cbar = True
#    else:
#        cbar = False
#
#    im = line.plot(limit=limit,Tb=True,shift_dx=star_shift_dx,shift_dy = star_shift_dy, iv=iv0 + i, no_ylabel=i>0, ax=ax, cmap=cmap, fmin=5, colorbar=True, colorbar_label=cbar)
#
#    ax.plot(0,0,"*",color="white",ms=4)
#    ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
#
#plt.tight_layout()


#---- Channel maps v2 : channel 2x3 + moment 1
win=13

fig = plt.figure(figsize=(15.5, 5), constrained_layout=False, num=win)
#fig = plt.figure(constrained_layout=False)
gs1 = fig.add_gridspec(nrows=2, ncols=3, left=0.05, right=0.6, wspace=0.0, hspace=0.05)
ax0 = fig.add_subplot(gs1[0, 0])
ax1 = fig.add_subplot(gs1[0, 1])
ax2 = fig.add_subplot(gs1[0, 2])
ax3 = fig.add_subplot(gs1[1, 0])
ax4 = fig.add_subplot(gs1[1, 1])
ax5 = fig.add_subplot(gs1[1, 2])

axes3 = [ax0,ax1,ax2,ax3,ax4,ax5]


gs2 = fig.add_gridspec(nrows=1, ncols=1, left=0.6, right=0.98, hspace=0.05)
ax_m1 = fig.add_subplot(gs2[0,0])

n_channels = 6
iv0=18

limit = 0.5
limits = [0.25,-0.6,-0.6,0.25]


Vsyst = 2.5
Mstar = 0.8
toy = dynamite.toy_model(Mstar=Mstar, dist=140, inc=28.5, cube=line, r0=100, z0=30, beta=1., PA=48)

for i in range(n_channels):
    ax = axes3[i]

    if i<3:
        ax.axes.xaxis.set_ticklabels([])


    if i%3==0:
        no_ylabel = False
    else:
        no_ylabel = True

    if i%3==2:
        cbar=True
    else:
        cbar=False

    if i>2:
        no_xlabel = False
    else:
        no_xlabel = True

    im = line.plot(limits=limits,Tb=True,shift_dx=star_shift_dx,shift_dy = star_shift_dy, iv=iv0 + i, ax=ax, cmap=cmap, fmin=5, colorbar=True, colorbar_label=cbar,no_xlabel=no_xlabel,no_ylabel=no_ylabel)

    ax.plot(0,0,"*",color="white",ms=4)
    ax.plot(x_planet, y_planet ,"o",color="darksalmon",ms=2)


    #--- Adding toy model
    isov_color= "white"
    alpha = 0.3

    #sin_i = np.sin(inc[i]/180 * np.pi)
    #toy.plot_isovelocity_curve(v= (line.velocity[iv0+i] - Vsyst),ax=ax,colors=isov_color,flip_v = False, linestyles="-", linewidths=1.25, rmax=1400, alpha=alpha) # axes[k,2]
#    toy.plot_isovelocity_curve(v= (line.velocity[civ0+i - Vsyst) - 0.2 * v_Kep * sin_i,ax=axes[k,2],colors=isov_color,flip_v = flip_v[i], linestyles="-",linewidths=1.25, rmax=1400, alpha=alpha) # axes[k,2]


v0 = 0.5
dv = 0.8
ax = ax_m1
line.plot(limits=limits,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=1, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], M0_threshold=5e-3, fmin=-dv/2+0.1+v0,fmax=dv/2+0.1+v0, no_ylabel = True, interpolation="nearest")

ax.plot(0,0,"*",color="black",ms=4)
#ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
circle = Circle((x_planet_CO, y_planet_CO), 0.1, clip_on=False, zorder=10, linewidth=2, edgecolor='black', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
ax.add_artist(circle)

ax.annotate('planet\ncandidate', xy=(x_planet, y_planet), xytext=(0.0, 0.16),
            arrowprops=dict(arrowstyle="->", color="0.3",
                            shrinkA=40,
                            shrinkB=35,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0"),
             color='0.3', fontsize=12, va='center', ha='center', rotation=0.0)

ax.annotate('inner planet\nwake?', xy=(-0.07, -0.33), xytext=(0.12,0.03), #xytext=(0.34,-0.07),
            arrowprops=dict(arrowstyle="->", color="0.3",
                            shrinkA=20,
                            #shrinkB=5,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0.3"),
            color='0.3', fontsize=12, va='center', ha='center', rotation=0.0)


ax.annotate('gap ?', xy=(-0.18, -0.43), xytext=(0.14,-0.4), #xytext=(0.34,-0.07),
            arrowprops=dict(arrowstyle="->", color="0.3",
                            shrinkA=20,
                            #shrinkB=5,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0."),
            color='0.3', fontsize=12, va='center', ha='center', rotation=0.0)


#--- plot moment maps with toy model to find stellar mass

win=15

fig, axes = plt.subplots(nrows=2, ncols=5, num=win, figsize=(16.5, 5))

V = [-3.4,-3,-2.6,-2.2,-1.8,1.8,2.2,2.6,3.,3.4]

Vsyst = 3.0
Mstar = 1.5

Vsyst = 2.0
Mstar = 0.8

toy = dynamite.toy_model(Mstar=Mstar, dist=140, inc=28.5, cube=line, r0=100, z0=30, beta=1., PA=48)

for i in range(len(V)):
    ax = axes.ravel()[i]

    im = line.plot(limits=np.array(limits)*2.5,Tb=True,shift_dx=star_shift_dx,shift_dy = star_shift_dy, v = V[i]+Vsyst, ax=ax, cmap=cmap, fmin=5, colorbar=True, colorbar_label=cbar,no_xlabel=no_xlabel,no_ylabel=no_ylabel)

    ax.plot(0,0,"*",color="white",ms=4)
    ax.plot(x_planet, y_planet ,"o",color="darksalmon",ms=2)


    #--- Adding toy model
    isov_color= "white"
    alpha = 0.3

    #sin_i = np.sin(inc[i]/180 * np.pi)
    toy.plot_isovelocity_curve(v= V[i],ax=ax,colors=isov_color,flip_v = False, linestyles="-", linewidths=1.25, rmax=1400, alpha=alpha) # axes[k,2]

#--- moment map v1

#if False:
#    v0 = 0.5
#    dv = 0.8
#    n_moments = 2
#    limit = 1.0
#
#    fig_xfactor = 4
#    fig_yfactor = 4
#
#    win=3
#    fig3, axes3 = plt.subplots(nrows=1, ncols=n_moments+1, figsize=((1.4*fig_xfactor)*(n_moments),fig_yfactor*1), num=win)
#
#
#    ax = axes3[0]
#    cont.plot(limit=limit,Tb=True,fmax=90,fmin=0,shift_dx=star_shift_dx,shift_dy = star_shift_dy, ax=ax, cmap="inferno")
#    ax.plot(0,0,"*",color="white",ms=2)
#    ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
#
#
#    ax = axes3[1]
#    line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], no_ylabel = True, cmap=cmap, Tb=True)
#
#    ax.plot(0,0,"*",color="white",ms=4)
#    #ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
#    circle = Circle((x_planet, y_planet), 0.1, clip_on=False, zorder=10, linewidth=2, edgecolor='white', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
#    ax.add_artist(circle)
#
#    ax = axes3[2]
#    line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=1, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], M0_threshold=5e-3, fmin=-dv/2+0.1+v0,fmax=dv/2+0.1+v0, no_ylabel = True)
#
#    ax.plot(0,0,"*",color="black",ms=4)
#    #ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
#    circle = Circle((x_planet, y_planet), 0.1, clip_on=False, zorder=10, linewidth=2, edgecolor='black', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
#    ax.add_artist(circle)
#
#    plt.tight_layout()
#

#---------------- Figure 1 : ALMA detcetion + NACO

v0 = 0.5
dv = 0.8
n_moments = 2
limit = 0.99

fig_xfactor = 4
fig_yfactor = 4

win=4
fig4, axes4 = plt.subplots(nrows=2, ncols=n_moments, figsize=(1.4*(fig_xfactor)*(n_moments),fig_yfactor*2),  constrained_layout=False, num=win)
plt.subplots_adjust(hspace=0.05)
axes4[0,1].axes.xaxis.set_ticklabels([])
axes4[0,0].axes.xaxis.set_ticklabels([])

# ALMA moment 0
ax = axes4[0,1]
line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], cmap=cmap, Tb=True, no_xlabel=True, no_ylabel=True)


ax.plot(0,0,"*",color="white",ms=4)
#ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
circle = Circle((x_planet_CO, y_planet_CO), 0.1, clip_on=False, zorder=10, linewidth=1, edgecolor='white', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
ax.add_artist(circle)

ax.annotate('planet\ncandidate', xy=(x_planet, y_planet), xytext=(-0.74, 0.5),
            arrowprops=dict(arrowstyle="->", color="white",
                            shrinkA=20,
                            shrinkB=10,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0"),
             color="white", fontsize=8, va='center', ha='center', rotation=0.0)

inner_wake = (-0.35, -0.18)
ax.annotate('inner planet\nwake?', xy=inner_wake, xytext=(-0.51, -0.8),
            arrowprops=dict(arrowstyle="->", color="white",
                            shrinkA=20,
                            #shrinkB=5,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0.3"),
            color="white", fontsize=8, va='center', ha='center', rotation=0.0)
ax.text(0.9,0.9,"Integrated $^{12}$CO J=2-1",color="white",fontsize=10,va="center")
ax.text(0.9,0.78,"(v=0.5$\pm$0.8$\,$km$\,$s$^{-1}$)",color="white",fontsize=10,va="center")


# NACO image
ax = axes4[0,0]
NACO.plot(ax=ax,iv=14,interpolation="bicubic",stellar_mask=r_mask,fmax=15,cmap=NACO_cmap,limit=limit,fmin=-10, no_ylabel=False, no_xlabel=True) # 20th PCA component
#NACO.plot(ax=ax,interpolation="bicubic",stellar_mask=r_mask,cmap=NACO_cmap,limit=limit, no_ylabel=True, no_xlabel=True, color_scale="lin", fmin=0, fmax=15)

#ax.annotate('inner planet\nwake?', xy=inner_wake, xytext=(-0.51, -0.8),
#            arrowprops=dict(arrowstyle="->", color="white",
#                            shrinkA=20,
#                            #shrinkB=5,
#                            patchA=None,
#                            #patchB=None,
#                            connectionstyle="arc3,rad=0.3"),
#            color="white", fontsize=8, va='center', ha='center', rotation=0.0)

#ax.annotate('outer planet\n wake?', xy=(-0.23, 0.47), xytext=(0.6, 0.75),
#            arrowprops=dict(arrowstyle="->", color="white",
#                            shrinkA=35,
#                            #shrinkB=5,
#                            patchA=None,
#                            #patchB=None,
#                            connectionstyle="arc3,rad=-0.3"),
#            color="white", fontsize=8, va='top', ha='center', rotation=0.0)

ax.plot(0,0,"*",color="white",ms=4)
#ax.plot(x_planet, y_planet ,"o",color="black",ms=2)

ax.annotate('planet\ncandidate', xy=(x_planet, y_planet), xytext=(-0.74, 0.5),
            arrowprops=dict(arrowstyle="->", color="white",
                            shrinkA=20,
                            shrinkB=5,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0"),
             color="white", fontsize=8, va='center', ha='center', rotation=0.0)
ax.text(0.9,0.9,"NACO L' PCA-ADI",color="white",fontsize=10,va="center")
ax.plot(x_planet, y_planet ,".",color="black",ms=2)

# continuum overlay
ax = axes4[1,1]
cont.plot(limit=limit,Tb=True,fmax=90,fmin=3,shift_dx=star_shift_dx,shift_dy = star_shift_dy, ax=ax, cmap="inferno", color_scale="sqrt")
ax.plot(0,0,"*",color="white",ms=2)
#ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)

line.plot(ax=ax,limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0,v0=v0, v_minmax=[v0-dv,v0+dv], no_ylabel = False, cmap="viridis", Tb=True,plot_type="contourf",levels=levels,alpha=0.75, colorbar=False, no_xlabel=True,zorder=10, plot_beam=False)

#ax.plot(0,0,"*",color="black",ms=4)
ax.text(0.9,0.9,"",color="white",fontsize=10,va="center")
ax.text(0.9,0.9,"ALMA 1.3mm continuum + $^{12}$CO",color="white",fontsize=10,va="center")



# NACO overlay
ax = axes4[1,0]
NACO.plot(ax=ax,iv=14,interpolation="bicubic",stellar_mask=r_mask,fmax=15,cmap=NACO_cmap,limit=limit,fmin=-10, no_ylabel=False, no_xlabel=False) # 20th PCA component
#NACO.plot(ax=ax,interpolation="bicubic",stellar_mask=r_mask,cmap=NACO_cmap,limit=limit, no_ylabel=True, no_xlabel=True, color_scale="lin", fmin=0, fmax=15)

line.plot(ax=ax,limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0,v0=v0, v_minmax=[v0-dv,v0+dv], no_ylabel = True, Tb=True,plot_type="contour",levels=levels,alpha=1, colorbar=False, no_xlabel=True,zorder=10,linewidths=0.6, colors="white") #cmap=cmr.arctic)
ax.plot(0,0,"*",color="white",ms=4)
ax.text(0.9,0.9,"NACO + $^{12}$CO",color="white",fontsize=10,va="center")


# moment 1
#ax = axes4[1,0]
#line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=1, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], M0_threshold=5e-3, fmin=-dv/2+0.1+v0,fmax=dv/2+0.1+v0)

#plt.tight_layout()

#---------------- Moment map v3

if False:
    win=5
    fig5, axes5 = plt.subplots(nrows=1, ncols=3, figsize=(1*(fig_xfactor)*3,fig_yfactor*1), num=win, constrained_layout=True)
    plt.subplots_adjust(hspace=0.3)

    # moment 0
    ax = axes5[0]
    line.plot(limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0, ax=ax,v0=v0, v_minmax=[v0-dv,v0+dv], cmap=cmap, Tb=True, no_xlabel=True)

    ax.plot(0,0,"*",color="white",ms=4)
    #ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
    circle = Circle((x_planet, y_planet), 0.1, clip_on=False, zorder=10, linewidth=2, edgecolor='white', linestyle=":", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
    ax.add_artist(circle)


    # NACO
    ax = axes5[1]
    NACO.plot(ax=ax,iv=19,interpolation="bicubic",stellar_mask=r_mask,fmax=15,cmap="viridis",limit=limit,fmin=-10, no_ylabel=True) # 20th PCA component

    ax.annotate('inner planet\nwake?', xy=(-0.35, -0.18), xytext=(-0.51, -0.48),
                arrowprops=dict(arrowstyle="->", color="white",
                                shrinkA=10,
                                #shrinkB=5,
                                patchA=None,
                                #patchB=None,
                                connectionstyle="arc3,rad=0.3"),
                color="white", fontsize=10, va='center', ha='center', rotation=0.0)

    ax.annotate('outer planet\n wake?', xy=(-0.23, 0.47), xytext=(0.6, 0.75),
                arrowprops=dict(arrowstyle="->", color="white",
                                shrinkA=10,
                                #shrinkB=5,
                                patchA=None,
                                #patchB=None,
                                connectionstyle="arc3,rad=-0.3"),
                color="white", fontsize=10, va='top', ha='left', rotation=0.0)

    # NACO overlay
    ax = axes5[2]
    NACO.plot(ax=ax,iv=19,interpolation="bicubic",stellar_mask=r_mask,fmax=15,cmap="viridis",limit=limit,fmin=-10, no_ylabel=True, colorbar=False) # 20th PCA component

    line.plot(ax=ax,limit=limit,shift_dx=star_shift_dx,shift_dy = star_shift_dy, moment=0,v0=v0, v_minmax=[v0-dv,v0+dv], no_ylabel = True, cmap=cmap, Tb=True,plot_type="contourf",levels=levels,alpha=0.5, colorbar=False,zorder=10)

#plt.tight_layout()


# NACO subtraction
#win=6
#vmin = np.percentile(NACO.image[19,:,:],0.1)
#vmax = np.percentile(NACO.image[19,:,:],99.5)
#
#
#fig6, axes6 = plt.subplots(nrows=1, ncols=2, figsize=(1*(fig_xfactor)*2,fig_yfactor*1), num=win, constrained_layout=True)
#
#limit=0.99
#ax = axes6[0]
#NACO.plot(ax=axes6[0],iv=19,interpolation="bicubic",stellar_mask=r_mask,fmax=vmax,cmap="viridis",limit=limit,fmin=vmin, colorbar=False, color_scale="lin") # 20th PCA compon
#plt.subplots_adjust(wspace=0.01)
#
#
#x_planet= - r * NACO_pixelscale * np.cos(np.deg2rad(theta))
#y_planet = r * NACO_pixelscale * np.sin(np.deg2rad(theta))
#
#
#ax.annotate('best fit', xy=(x_planet, y_planet), xytext=(-0.72, 0.45),
#            arrowprops=dict(arrowstyle="->", color="1.0",
#                            shrinkA=15,
#                            shrinkB=3,
#                            patchA=None,
#                            #patchB=None,
#                            connectionstyle="arc3,rad=0"),
#             color='1.0', fontsize=12, va='center', ha='center', rotation=0.0)
#
#ax.plot(x_planet, y_planet ,".",color="black",ms=2)
#ax.text(0.9,0.9,"NACO L' PCA-ADI",color="white",fontsize=10,va="center")
#
#ax = axes6[1]
#ax.axes.yaxis.set_visible(False)
#NACO_empty.plot(ax=axes6[1],interpolation="bicubic",stellar_mask=r_mask,fmax=vmax,cmap="viridis",limit=limit,fmin=vmin, no_ylabel=True, colorbar=False, color_scale="lin") # 20th PCA compon
#ax.text(0.9,0.9,"Planet subtracted",color="white",fontsize=10,va="center")


#---- PV diagram : to check
#
#from pvextractor import Path
#from pvextractor import extract_pv_slice
#
#
#path = Path([(691, 765), (800., 660.)], width=3)
#
#path = Path([ (705, 755), (800, 655)], width=3)
#
#center = (739, 712)
#center = (745, 710)
#center = (740.5, 712.5)
#disk_PA = np.deg2rad(48)
#dx=150
#path = Path([ (center[0] - dx * np.cos(disk_PA), center[1] + dx * np.sin(disk_PA)), (center[0] + dx * np.cos(disk_PA), center[1] - dx * np. sin(disk_PA))], width=3)
#
#
#slice = extract_pv_slice(dir+"Elias24_CO.fits", path)
#
#
#plt.figure(13)
#plt.clf()
#halfsize = slice.shape[1] / 2 * line.pixelscale
#extent = [-halfsize,halfsize, np.min(line.velocity), np.max(line.velocity)]
#
#norm = colors.Normalize(vmin=0., vmax=np.max(slice.data), clip=True)
#
#plt.imshow(slice.data, origin='lower', extent=extent, aspect="auto", norm=norm, cmap=cmap)
#
#def v_kepler(Mstar,r):
#    from scipy import constants as sc
#    GxMsun   = 1.3271244e20
#
#    return np.sign(r) * np.sqrt(GxMsun * Mstar  / (np.maximum(np.abs(r),1e-3) * sc.au))
#
#distance = 136
#
#inc = np.deg2rad(28.5)
#v_syst = 2.85
#
#Mstar = 1.5
#dv = 0.25
#
#r = np.linspace(0.01,1,100)
#plt.plot(r , -v_kepler(Mstar,r*distance)*1e-3 * np.sin(inc) + v_syst, color="white")
#plt.plot(r , -v_kepler(Mstar+dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")
#plt.plot(r , -v_kepler(Mstar-dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")
#
#r = -r
#plt.plot(r , -v_kepler(Mstar,r*distance)*1e-3 * np.sin(inc) + v_syst, color="white")
#plt.plot(r , -v_kepler(Mstar+dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")
#plt.plot(r , -v_kepler(Mstar-dv,r*distance)*1e-3 * np.sin(inc) + v_syst, linestyle="--", color="white")
#
#plt.xlim(-0.6,0.6)
#plt.ylim(-4,12)
#plt.xlabel(r'Distance along the midplane ["]')
#plt.ylabel(r"Velocity [km/s]")
#
##plt.plot([np.max(np.abs(r)),-np.max(np.abs(r))], [v_syst,v_syst], linestyle=":", color="white")
#
#beam_position=(0.125, 0.125)
#dx, dy = beam_position
#bpa = 0.
#bmaj = 0.7 # km/s
#bmin = line.bmaj
#ax = plt.gca()
#beam = Ellipse(
#    ax.transLimits.inverted().transform((dx, dy)),
#    width=bmin,
#    height=bmaj,
#    angle=-bpa,
#    fill=False,
#    color="white",
#)
#ax.add_patch(beam)
#
#
#
#
#model = mcfost.Line("/Users/cpinte/mcfost+phantom/Elias24/5Mj_hr0.1/no_Lacc/data_CO")
#
#win=5
#fig4, axes4 = plt.subplots(nrows=1, ncols=n_moments, figsize=((1.4*fig_xfactor)*n_moments,fig_yfactor*1), num=win)
#
#iv_support = [107,108,109,110,111]
#
#ax = axes4[0]
#model.plot_map(iv_support=iv_support,psf_FWHM=0.075,plot_stars=True,Tb=True,moment=0,ax=ax)
#
#ax = axes4[1]
#model.plot_map(iv_support=iv_support, psf_FWHM=0.075,plot_stars=True,moment=1,ax=ax, fmin=1,fmax=1.5)


plt.show()
