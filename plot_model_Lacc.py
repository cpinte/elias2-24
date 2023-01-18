# %%

from matplotlib.patches import Circle
import casa_cube as casa
import cmasher as cmr

x_planet= -0.33
y_planet = 0.12

#------------------
dir = "/Users/cpinte/Observations/DSHARP/Elias24/"

# %%

cont = casa.Cube(dir+"Elias24_continuum.fits")
CO = casa.Cube(dir+"Elias24_CO.fits")

dx = 0.105
dy = -0.383

Vsyst = 3.0


# NACO :
fwhm = fits.getdata('NACO/fwhm.fits')[0]
NACO_pixelscale = 0.027208 # arcsec/pixel
#NACO = casa.Cube("NACO/final_PCA-ADI_ann_npc1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19-20.fits",pixelscale=NACO_pixelscale, unit="adu")
NACO = casa.Cube("NACO/PCA-annular_npc1-50_delta_rot1.fits",pixelscale=NACO_pixelscale, unit="adu")
NACO.bmaj = fwhm * NACO_pixelscale
NACO.bmin = fwhm * NACO_pixelscale

r_mask = 2*fwhm * NACO_pixelscale

#todo :
#try to change PA a tiny bit

# 1Mj --> 3Mj  visible at 1e-7 : cont way too bright
# 3Mj --> 5.5Mj visible at 1e-8
# 5Mj --> 8Mj  visible at 1e-9


# fluffy : gap in CO visible if fluffy closer to 1, not sure why, CPD less visible with higher fluffy parameter too
# 1 Mj  --> 0.3
# 3Mj --> 0.1
# 5Mj --> 0.03



#------------------
Mplanet = ["10","5","5","3","3"]
Lacc = [0.,0.,1e-7,0.,1e-7]
Lacc_name = ["","","$\dot{\mathrm{M}} = 10^{-7}$ M$_\odot$","","$\dot{\mathrm{M}} = 10^{-7}$ M$_\odot$"]


Mplanet = ["5","5","3","1"]
Lacc = [1e-9,0.0,1e-8,1.e-7]
Lacc_name = ["$\dot{\mathrm{M}} = 3\,10^{-8}$ M$_\odot$","","$\dot{\mathrm{M}} = 10^{-8}$ M$_\odot$","",""]

#Mplanet = ["1","3","5","7"]
Mplanet = ["1","1","5","5"]
Lacc = [1e-7,0.0,3e-9,3e-9]
#Lacc = [1e-7,0.0,0.0,0.0]
#Lacc_name = ["$\dot{\mathrm{M}} = 3\,10^{-8}$ M$_\odot$","","$\dot{\mathrm{M}} = 10^{-8}$ M$_\odot$","",""]
Lacc_name = ["$10^{-7}$M$_\mathrm{\odot}$/yr","0","0","0","0"]
Lacc_name = ["$10^{-4}$M$_\mathrm{Jup}$/yr","0","0","0","0"]

fluffy = [1,1,1,0.01]

#------------------


n_models = len(Mplanet)

V = [-1.5,-1.8,-2.1, -2.4]
V = [-1.35,-1.8,-2.15, -2.5, -2.85]
V = [-1.8,-2.15, -2.5]
nv = len(V)


mod_dir_basename = "/Users/cpinte/mcfost+phantom/Elias24/" ; mod_dir_ext= "Mj/" ;


win=16 ; factor = 3.0
fig, axes = plt.subplots(nrows=n_models+1, ncols=nv+2, figsize=((nv+1) * factor,(n_models) * factor), sharex='all', sharey='all', num=win)
plt.subplots_adjust(wspace=0.03, hspace=0.0)

limits = [0.99,-0.99,-0.99,0.99]


star_shift_dx = 0.105
star_shift_dy = -0.383


cmap_cont = "inferno"
NACO_cmap = "viridis"
NACO_cmap = "gist_earth" #
NACO_cmap = cmr.dusk



# We plot the observations on the first row
im = cont.plot(colorbar=False, ax=axes[0,0], no_xlabel=True, limits=limits, cmap=cmap_cont, color_scale="lin", shift_dx=star_shift_dx, shift_dy = star_shift_dy, Tb=True, fmin=5, fmax=35)
#colorbar2(im)

axes[0,0].text(0.92,0.85,"ALMA 1.25mm continuum",color="white",fontsize=10,va="center")


win=2
fig_xfactor = 5
fig_yfactor = 5

n_channels = 6
iv0=18

for i in range(nv):
    ax = axes[0,1+i]
    no_ylabel = True #i>0
    iv = np.abs(CO.velocity - (Vsyst + V[i])).argmin()
    im = CO.plot(iv=iv ,colorbar=False,Tb=True, v0=Vsyst, ax=ax, no_xlabel=True, no_ylabel=no_ylabel, limits=limits, shift_dx=star_shift_dx, shift_dy = star_shift_dy, cmap=cmr.arctic)

    #line.plot(limit=limit,Tb=True,fmax=55,fmin=5,shift_dx=star_shift_dx,shift_dy = star_shift_dy, iv=iv0 + i, ax=ax)
    ax.plot(0,0,"*",color="white",ms=4)
    #    ax.plot(x_planet, y_planet ,"o",color="cyan",ms=2)
    circle = Circle((x_planet, y_planet), 0.2, clip_on=False, zorder=10, linewidth=1, edgecolor='lightgrey', linestyle="-", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
    ax.add_artist(circle)

axes[0,1].text(0.92,0.85,"$^{12}$CO J=2-1",color="white",fontsize=10,va="center")


ax = axes[0,nv+1]
NACO.plot(ax=ax,iv=14,interpolation="bicubic",stellar_mask=r_mask,fmax=15,cmap=NACO_cmap,limits=limits,fmin=-10, no_ylabel=True, no_xlabel=True, colorbar=False) # 20th PCA component

#inner_wake = (-0.35, -0.18)
#ax.annotate('inner planet\nwake?', xy=inner_wake, xytext=(-0.51, -0.8),
#            arrowprops=dict(arrowstyle="->", color="white",
#                            shrinkA=20,
#                            #shrinkB=5,
#                            patchA=None,
#                            #patchB=None,
#                            connectionstyle="arc3,rad=0.3"),
#            color="white", fontsize=8, va='center', ha='center', rotation=0.0)
#
#ax.annotate('outer planet\n wake?', xy=(-0.23, 0.47), xytext=(0.55, 0.75),
#            arrowprops=dict(arrowstyle="->", color="white",
#                            shrinkA=35,
#                            #shrinkB=5,
#                            patchA=None,
#                            #patchB=None,
#                            connectionstyle="arc3,rad=-0.3"),
#            color="white", fontsize=8, va='top', ha='center', rotation=0.0)
#
ax.plot(0,0,"*",color="white",ms=4)
#ax.plot(x_planet, y_planet ,"o",color="black",ms=2)

ax.annotate('planet\ncandidate', xy=(x_planet-0.01, y_planet+0.03), xytext=(-0.68, 0.5),
            arrowprops=dict(arrowstyle="->", color="white",
                            shrinkA=12,
                            shrinkB=5,
                            patchA=None,
                            #patchB=None,
                            connectionstyle="arc3,rad=0"),
             color="white", fontsize=8, va='center', ha='center', rotation=0.0)
ax.text(0.92,0.85,"NACO L' PCA-ADI",color="white",fontsize=10,va="center")



#------------ Models
for k, mplanet in enumerate(Mplanet):
    mod_dir_cont = mod_dir_basename+mplanet+mod_dir_ext+"/phantom_Lacc/"
    mod_dir = mod_dir_basename+mplanet+mod_dir_ext+"/01001_Mdot_"+str(Lacc[k])+"_fluffy_"+str(fluffy[k])+"/"
    mod_dir_cont = mod_dir

    print(mod_dir)

    modcont = mcfost.Image(mod_dir_cont+"data_1300/")
    mod_NACO = mcfost.Image(mod_dir_cont+"data_3.8/")
    modCO = mcfost.Line(mod_dir+"data_CO/")


    if (Lacc[k] > 1e-10):
        #mod_dir = mod_dir_basename+mplanet+mod_dir_ext+"Lacc/01001_Mdot_"+str(Lacc[k])+"/"
        title = Mplanet[k]+"M$_\mathrm{jup}$, $\dot{\mathrm{M}}="+str(Lacc[k])+", f="+str(fluffy[k])
    else:
        #mod_dir = mod_dir_cont
        title = Mplanet[k]+"M$_\mathrm{jup}$, f="+str(fluffy[k])

    #title = "M="+"{:.1f}".format(modcont.star_properties[1,0] * 1047.57)+"M$_\mathrm{jup}$\n"+"\'M="+str(Lacc[k])+"M$_\mathrm{\odot}$/yr\n"+"p="+str(fluffy[k])
    title = "M="+"{:.1f}".format(modcont.star_properties[1,0] * 1047.57)+"M$_\mathrm{jup}$\n"+"$\dot{\mathrm{M}}$="+Lacc_name[k]+"\n"+"p="+str(1-fluffy[k])


    no_xlabel = k<n_models-1

    #modcont.image += 100 * modcont.image[4,0,0,:,:]
    modcont.plot(ax = axes[k+1,0], colorbar=False, bmaj=cont.bmaj, bmin=cont.bmin, bpa=cont.bpa, no_xlabel = no_xlabel, limits=limits, cmap = cmap_cont, scale="lin", Tb=True, vmin=5, vmax=35, interpolation="bicubic")
    #modcont.plot(ax = axes[k+1,0], colorbar=False, no_xlabel = no_xlabel, limits=limits, cmap = cmap_cont, scale="lin", Tb=True)

    Delta_v = 0.225 * 0.
    Delta_v = 0.7
    for i in range(nv):
        no_ylabel = True # i>0
        ax = axes[k+1,1+i]
        modCO.plot_map(v=V[i], ax = ax, colorbar=False, bmaj=CO.bmaj, bmin=CO.bmin, bpa=CO.bpa, no_ylabel=no_ylabel, no_xlabel = no_xlabel, limits=limits, Tb=True, Delta_v=Delta_v, subtract_cont=True, cmap=cmr.arctic, no_vlabel=True, rms=3, fmin=0, interpolation="lanczos")

        circle = Circle((x_planet, y_planet), 0.2, clip_on=False, zorder=10, linewidth=1, edgecolor='lightgrey', linestyle="-", facecolor=(0, 0, 0, .0125), alpha=0.5, fill=False)
        ax.add_artist(circle)

    mod_NACO.plot(ax = axes[k+1,nv+1], colorbar=False, no_xlabel = no_xlabel, limits=limits, cmap = NACO_cmap, scale="log", vmax=1e-16, interpolation="bicubic")

    xc = mod_NACO.nx//2
    dx = int(np.rint(-mod_NACO.star_positions[0,0,0,1]/mod_NACO.pixelscale))
    dy = int(np.rint(mod_NACO.star_positions[1,0,0,1]/mod_NACO.pixelscale))
    Delta_m = -2.5 * np.log10(mod_NACO.last_image[xc+dy,xc+dx] / mod_NACO.last_image[xc,xc])

    #print(Delta_m)
    NACO_title = "$m_\mathrm{p} - m_*$ = "+"{:.2f}".format(Delta_m)


    #, psf_FWHM=mcfost.telescope_beam(3.8*1e-6,8.2))


    # label
    ax = axes[k+1,0]
    ax.text(0.1,0.65,title,horizontalalignment='left',color="white",transform=ax.transAxes,fontsize=10)

    ax = axes[k+1,-1]
    ax.text(0.1,0.85,NACO_title,horizontalalignment='left',color="white",transform=ax.transAxes,fontsize=10)


#colorbar2(im, ax=axes[0,1:], trim_left=0.08)


plt.show()
