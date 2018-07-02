import numpy as np
import matplotlib.pyplot as plt
import matplotlib  as matplotlib
from scipy import ndimage

# create a power per unit volume image in a plane 
# returning x,y,hist to plot with pcolormesh
#   dx is resolution of sum of power, pixel width 
#   dz is height limit for points from  midplane
#   passing xrot,yrot,zarr coordinate arrays
#   power array is parr
#   summing powers within distance of midplane
#   each pixel is energy per time summed in box
#   power is divided by 3d pixel volume  dx^2*dz

pltxmax = 1.0   # set plot limits
def heatxy(xrot,yrot,zarr,parr,dx,dz):
	pixvolume =  dx*dx*dz;  # divide by this if you want energy density/time
	n = np.size(xrot);
	xmax = pltxmax; xmin = -xmax; 
	ymax = pltxmax; ymin = -ymax; 
	nx = np.int((xmax-xmin)/dx) + 1;
	yarrh,xarrh = np.meshgrid(np.arange(ymin,ymax+dx,dx),\
		np.arange(xmin,xmax+dx,dx) )  #nontrivial order here, is flipped 
	#xarrh[i,*] gives x depending on i and is all the same for second index
	mhist = xarrh*0.0;
	poff = 0.0;  # pixel offset
	for i in range (0,n):
		if (np.abs(zarr[i]) < dz):
			xi = np.int((xrot[i] - xmin)/dx + poff);
			yi = np.int((yrot[i] - xmin)/dx + poff);
			if ((xi >=0) and (xi<nx)):
				if ((yi >=0) and (yi<nx)):
					mhist[xi,yi]= mhist[xi,yi]+parr[i]/pixvolume
					#checked order!
					
	return xarrh,yarrh,mhist   # note order 


# create a power per unit volume image in a plane 
# returning x,y,mhist to plot with pcolormesh
#   dx is resolution of sum of power, pixel width 
#   dz is height limit for points from  midplane
#   passing xrot,yrot,zarr coordinate arrays
#   power array is parr
#   summing powers within distance of midplane
#   each pixel is energy per time summed in box
#   power is divided by 3d pixel volume  dx^2*dz
# use Gaussian weighting as a function of distance from plane
def heatxysig(xrot,yrot,zarr,parr,dx,sigz,pmax):
	pixvolume =  dx*dx*sigz;  # divide by this if you want energy density/time
	# mistake here? w0 has sigz in it
	n = np.size(xrot);
	xmax = pmax; xmin = -xmax; 
	ymax = pmax; ymin = -ymax; 
	nx = np.int((xmax-xmin)/dx) + 1;
	yarrh,xarrh = np.meshgrid(np.arange(ymin,ymax+dx,dx),\
		np.arange(xmin,xmax+dx,dx) )  #nontrivial order here, is flipped 
	#xarrh[i,*] gives x depending on i and is all the same for second index
	mhist = xarrh*0.0;
	poff = 0.0;  # pixel offset
	w0 = 1.0/np.sqrt(2.0*sigz*sigz*np.pi)
	for i in range (0,n):
		if (np.abs(zarr[i]) < 3*sigz):
			zg = 0.5*(zarr[i]/sigz)**2
			weight = w0*np.exp(-zg) # gaussian weight
			xi = np.int((xrot[i] - xmin)/dx + poff);
			yi = np.int((yrot[i] - xmin)/dx + poff);
			if ((xi >=0) and (xi<nx)):
				if ((yi >=0) and (yi<nx)):
					mhist[xi,yi]= mhist[xi,yi]+weight*parr[i]/pixvolume
					#checked order!
					
	return xarrh,yarrh,mhist   # note order 


# make 3x1 images plot given arrays
# save figure in file fname
def plt3_hist(x1,y1,hist1,x2,y2,hist2,x3,y3,hist3,fname):
	zmin=0
	z1=hist1.max();
	z2=hist2.max();
	z3=hist3.max();
	zmax = max(z1,z2,z3)
	matplotlib.rcParams.update({'font.size': 20})
	f, axarr = plt.subplots(1,3,sharex=True,sharey=True,figsize=(8.5,3));
	f.subplots_adjust(hspace=0);
	f.subplots_adjust(wspace=0.05);
	f.subplots_adjust(left=0.07)
	f.subplots_adjust(right=0.82)
	plt.axis('tight');
	plt.axis([x1.min(), x1.max(), x1.min(), x1.max()]);
	plt.setp([a.get_xticklabels() for a in f.axes], visible=False) ;
	plt.setp([a.get_yticklabels() for a in f.axes], visible=False) ;
	plt.setp(axarr[0].get_xticklabels(), visible=True) ;
	plt.setp(axarr[0].get_yticklabels(), visible=True) ;
	axarr[0].set_aspect(1.0);
	axarr[1].set_aspect(1.0);
	axarr[2].set_aspect(1.0);
	axarr[0].set_adjustable('box-forced')
	axarr[1].set_adjustable('box-forced')
	axarr[2].set_adjustable('box-forced')
	axarr[0].pcolormesh(x1,y1,hist1,cmap='spectral',vmin=zmin,vmax=zmax);
	axarr[1].pcolormesh(x3,y3,hist3,cmap='spectral',vmin=zmin,vmax=zmax); #order changed
	xx=axarr[2].pcolormesh(x2,y2,hist2,cmap='spectral',vmin=zmin,vmax=zmax);
	cbar_ax = f.add_axes([0.830, 0.17, 0.03, 0.68])
	f.colorbar(xx, cax=cbar_ax,format='%.1e')
	d1=0.07; v1=0.05;
	d2=0.07; v2=0.09;
	pm = x1.max();
	axarr[0].text( pm-d1,-pm+v1,'x',ha='center',va='bottom',fontsize=16,color='white');
	axarr[0].text(-pm+d2, pm-v2,'y',ha='left',va='center',fontsize=16,color='white');
	axarr[1].text( pm-d1,-pm+v1,'x',ha='center',va='bottom',fontsize=16,color='white');
	axarr[1].text(-pm+d2, pm-v2,'z',ha='left',va='center',fontsize=16,color='white');
	axarr[2].text( pm-d1,-pm+v1,'y',ha='center',va='bottom',fontsize=16,color='white');
	axarr[2].text(-pm+d2, pm-v2,'z',ha='left',va='center',fontsize=16,color='white');
	plt.savefig(fname);


# return latitute and longitude angles in degrees from x,y,z vectors
def latlon(xarr,yarr,zarr):
	n = np.size(xarr);
	rarr = np.sqrt(xarr*xarr + yarr*yarr + zarr*zarr)
	phi = np.arctan2(yarr,xarr)*180.0/np.pi;  # [-pi,pi]
	theta = np.arcsin(zarr/rarr)*180.0/np.pi; # [-pi/2,pi/2]
	return rarr,theta,phi


# make a heat flux in lat lon space
# x is longitude
# y is latitude
def latlonhist(tfile,ifile,lbin):
	iarr,xarr,yarr,zarr,parr,xrot,yrot = np.loadtxt(tfile, skiprows=2, unpack=True)
	iarr0,xarr0,yarr0,zarr0,parr0,xrot0,yrot0 = np.loadtxt(ifile, skiprows=2, unpack=True)
	rarr,theta,phi = latlon(xarr0,yarr0,zarr0); # lat lon vectors
	n = np.size(xarr);
	ymin = -90.0; 
	ymax =  90.0 ;  # latitude
	xmin = -180.0;  
	xmax = 180.0;  # longitude
	dx = dy = lbin;
	latarrh,lonarrh = np.meshgrid(np.arange(ymin,ymax+dy,dy),\
		np.arange(xmin,xmax+dx,dx) )  #nontrivial order here, is flipped 
	mhist = latarrh*0.0;
	angfac = np.pi/180.0
	pixvolume = 1.0;  
	# dx*dy*angfac*angfac; # divide by if you want energy per surface area? 
	nx = np.int((xmax-xmin)/dx) + 1;
	ny = np.int((ymax-ymin)/dx) + 1;
	poff = 0.0;  # pixel offset
	for i in range (0,n):
		xi = np.int((phi[i] - xmin)/dx + poff);   # longitude
		yi = np.int((theta[i] - ymin)/dy + poff); #  latitude
		if ((xi >=0) and (xi<nx)):
			if ((yi >=0) and (yi<ny)):
				if rarr[i] > 0.5:
					ymid = ((yi + 0.5)*dy + ymin)*np.pi/180.0
					cymid = np.cos(ymid)
					mhist[xi,yi]= mhist[xi,yi]+parr[i]/pixvolume/cymid
					#checked order!

	#hist = hist/np.cos(latarrh*np.pi/180.0); # correct for area dependence on lat
	plt.figure()
	plt.axis('tight');
	plt.axis([lonarrh.min(), lonarrh.max(), latarrh.min(), latarrh.max()]);
	smo=2
	img = ndimage.gaussian_filter(mhist, smo)
	plt.pcolormesh(lonarrh,latarrh,img)
	plt.colorbar()
	return latarrh,lonarrh,mhist

					


# make a 3x1 plot of different heats from a single file
# using positions from an initial file!
# bin is pixel size (dx)
# pwidth is distance from midplane for each pixel (dz)
# reading in file and calling all above routines
# fname is output figure name
# ifile is t=0 file, tfile is at time we want to see heat
# smo is a smoothing disp. set to zero if you don't want smoothing
def plt3heat_init(tfile,ifile,bin,pwidth,fname,smo,xmax):
	#file = 'a1_000100_heat.txt'
	iarr,xarr,yarr,zarr,parr,xrot,yrot = np.loadtxt(tfile, skiprows=2, unpack=True)
	iarr0,xarr0,yarr0,zarr0,parr0,xrot0,yrot0 = np.loadtxt(ifile, skiprows=2, unpack=True)
	# add up energys in width
	x1,y1,histxy=heatxysig(xarr0,yarr0,zarr0,parr,bin,pwidth,xmax);
	y2,z2,histyz=heatxysig(yarr0,zarr0,xarr0,parr,bin,pwidth,xmax);
	x3,z3,histxz=heatxysig(xarr0,zarr0,yarr0,parr,bin,pwidth,xmax);
	histxy_s = ndimage.gaussian_filter(histxy,smo)
	histyz_s = ndimage.gaussian_filter(histyz,smo)
	histxz_s = ndimage.gaussian_filter(histxz,smo)
	# make a plot
	plt3_hist(x1,y1,histxy_s,y2,z2,histyz_s,x3,z3,histxz_s,fname);


# make a 3x1 plot of different heats from a single file
# bin is pixel size (dx)
# pwidth is distance from midplane for each pixel (dz)
# reading in file and calling all above routines
# fname is output figure name
def plt3heat(file,bin,pwidth,fname):
	#file = 'a1_000100_heat.txt'
	iarr,xarr,yarr,zarr,parr,xrot,yrot = np.loadtxt(file, skiprows=2, unpack=True)
	zrot=zarr
	#pwidth = 0.2; #bin= 0.1;
	# add up energys in width
	x1,y1,histxy=heatxy(xrot,yrot,zrot,parr,bin,pwidth);
	y2,z2,histyz=heatxy(yrot,zrot,xrot,parr,bin,pwidth);
	x3,z3,histxz=heatxy(xrot,zrot,yrot,parr,bin,pwidth);
	# make a plot
	plt3_hist(x1,y1,histxy,y2,z2,histyz,x3,z3,histxz,fname);


# for authomatically making heat file names
def indexfname(froot,index):
	num1 = '{0:d}'.format(index)
	ff = froot + '_'
	if (index < 100000):
		ff = ff + '0'
	if (index < 10000):
		ff = ff + '0'
	if (index < 1000):
		ff = ff + '0'
	if (index < 100):
		ff = ff + '0'
	if (index < 10):
		ff = ff + '0'
	ff = ff + num1 + '_heat.txt'
	return ff

# make plots but using a sum over files
# using xyz positions from a comparison file with index icomp 
#   imin is first file index
#   imax is last file index
#   di is distance between indices for files
#   froot is file root for creating heatfile names
#   bin is xy pixel size, pwidth is distance to midplane pixel size
#   output figure filename is automatically created
def plt3sumheat(froot,imin,imax,di,icomp,bin,pwidth):
	xtot = []
	ytot = []
	ztot = []
	ptot = []
	nfiles = (imax-imin)/di
	fac = 1.0/nfiles;  # for normalization
	print(fac)
	# read in the comparison file from timestep icomp 
	file = indexfname(froot,icomp);
	iarr0,xarr0,yarr0,zarr0,parr0,xrot0,yrot0 = np.loadtxt(file,skiprows=2,unpack=True)
	xcomp=xrot0;  # these are the coordinates used for body frame!
	ycomp=yrot0;
	zcomp=zarr0;
	for ii in range(imin,imax,di):
		file = indexfname(froot,ii);
		iarr,xarr,yarr,zarr,parr,xrot,yrot = np.loadtxt(file,skiprows=2,unpack=True)
		zrot=zarr
		xtot = np.concatenate([xtot,xcomp])   # using comparison file positions
		ytot = np.concatenate([ytot,ycomp])
		ztot = np.concatenate([ztot,zcomp])
		ptot = np.concatenate([ptot,parr])
		print(file)

	# add up energys in width
	x1,y1,histxy=heatxy(xtot,ytot,ztot,ptot,bin,pwidth);
	y2,z2,histyz=heatxy(ytot,ztot,xtot,ptot,bin,pwidth);
	x3,z3,histxz=heatxy(xtot,ztot,ytot,ptot,bin,pwidth);
	#img = ndimage.gaussian_filter(img, sigma=(5, 5));  # if we want to smooth
	fname = froot + '_3heat.png'
	# make a plot
	plt3_hist(x1,y1,histxy*fac,y2,z2,histyz*fac,x3,z3,histxz*fac,fname);
		
		
# add up power per unit volume as a function of radius
# return a radial array and a power per volume array
# dr is pixel size for radial array
# rarr is radial coordinates
#        rarr = np.sqrt(xarr*xarr + yarr*yarr + zarr*zarr)
# parr is power array
def heatdr(rarr,parr,dr):
	rmin = 0.0;
	poff = 0.0;
	n = np.size(rarr);
	rarrh= np.arange(0,1.0+2*dr,dr);
	histr = rarrh*0.0;
	nr = np.size(rarrh);
	for i in range (0,n):
		ri = np.int((rarr[i] - rmin)/dr + poff);
		if ((ri >=0) and (ri<nr)):
			histr[ri] = histr[ri] + parr[i];

	for i in range (0,nr):
		r1 = dr*i;
		r2 = r1 + dr; 
		volume = (4.0*np.pi/3.0)*(r2**3 - r1**3);  # is never zero
		histr[i] = histr[i]/volume; # normalize for volume

	return rarrh,histr;


# make a radial plot summing heats and using a comparison file
# return power average and radius arrays so a comparison plot can be made
def pltsumheatdr(froot,imin,imax,di,icomp,dr):
	rtot = []
	ptot = []
	nfiles = (imax-imin)/di
	fac = 1.0/nfiles;  # for normalization
	print(fac)
	# read in the comparison file from timestep icomp 
	file = indexfname(froot,icomp);
	iarr0,xarr0,yarr0,zarr0,parr0,xrot0,yrot0 = np.loadtxt(file,skiprows=2,unpack=True);
	rarr0 = np.sqrt(xrot0*xrot0 + yrot0*yrot0 + zarr0*zarr0); # radius
	for ii in range(imin,imax,di):
		file = indexfname(froot,ii);
		iarr,xarr,yarr,zarr,parr,xrot,yrot = np.loadtxt(file,skiprows=2,unpack=True)
		rtot = np.concatenate([rtot,rarr0])   # using comparison file positions
		ptot = np.concatenate([ptot,parr])

	rarrh,histr = heatdr(rtot,ptot,dr); # get the radial heat array
	histr = histr*fac;   # normalize it

	fname = froot + '_heatr.png'
	# make a plot
	matplotlib.rcParams.update({'font.size': 20})
	plt.figure();
	plt.plot(rarrh,histr,'ro:');
	plt.xlabel("radius");
	plt.ylabel("power/volume");
	plt.savefig(fname);
	return rarrh,histr;  # useful if you want to plot a bunch together
        
       
# plot 5 heat radii plots
def plt5r(fname,r,h1,h2,h3,h4,h5):
	matplotlib.rcParams.update({'font.size': 20});
	f=plt.figure()
	ax = f.add_subplot(1, 1, 1);
	#f.subplots_adjust(hspace=0);
	#f.subplots_adjust(wspace=0.05);
	f.subplots_adjust(top=0.95)
	f.subplots_adjust(bottom=0.15)
	f.subplots_adjust(left=0.15)
	f.subplots_adjust(right=0.95)
	#ax.yaxis.labelpad = 1; 
	ax.set_xlabel("radius");
	ax.set_ylabel("power/volume", labelpad=22);
	ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2));
	plt.plot(r,h1,"ko:",label="boring");
	plt.plot(r,h2,"bo:",label="ecc");
	plt.plot(r,h3,"go:",label="soft-shell");
	plt.plot(r,h4,"ro:",label="lopsided");
	plt.plot(r,h5,"mo:",label="hard-shell");
	ax.legend(loc='upper right');
	plt.savefig(fname);

