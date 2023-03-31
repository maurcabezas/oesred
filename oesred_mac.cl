procedure oesred (input)
### CONTACT: mauricio.cabezas@asu.cas.cz
string	input		{"e2022******.fit",prompt="Spectrum target to reduce(.fit)"}
string	output		{"target01",prompt="Output filename"}
string	idtarget	{"Target 01",prompt="Target name on header"}
int	napertures		{49,prompt="Number of apertures to be found\n\n# CALIBRATION PARAMETERS\n"}
###################### CALIBRATION
bool 	orgfile 		{yes,prompt="do you want organize files?"}
bool 	zerocomb 		{no,prompt="Combine zero level images?"}
bool 	trimcal 		{no,prompt="Trim flat and comp?"}
bool 	iftrimc 		{no,prompt="Use trim flat & comp?"}
bool 	zerocorcal 		{no,prompt="Apply zero level correction to flat & comp?"}
bool 	compcomb 		{no,prompt="Combine comparison lamp images?"}
bool 	flatcomb 		{no,prompt="Combine flat field images?"}
bool 	flatapall 		{no,prompt="Extract flat apertures?"}
bool 	compapall 		{no,prompt="Extract comparison apertures?"}
bool 	iddatabase 		{no,prompt="Use database folder for identification?"}
file 	idfolder 		{"idcomp",prompt="folder name with identification database"}
bool 	idencomp 		{no,prompt="Identify features in spectrum for dispersion solution?\n\n# OBJECT PARAMETERS\n"}
###################### OBJECT
bool 	trimob 			{no,prompt="Trim object?"}
bool 	iftrimo			{no,prompt="Use trim object?"}
bool		zerocorob		{no,prompt="Apply zero level correction to object?"}
bool 	crays 			{no,prompt="Remove cosmic rays?"}
bool		ifcrays			{no,prompt="Use object with cosmic rays extraction?"}
bool 	objectapall 	{no,prompt="Extract object apertures?"}
bool 	flatcor 		{no,prompt="Apply flat correction to object?"}
bool 	heliocor 		{no,prompt="calculate JD + RV-helio?"}
bool 	idref 			{no,prompt="refer database identification to images?"}
bool 	combine 		{no,prompt="combine NON-normalized spectra?"}
bool 	rvcorr 		{no,prompt="Apply heliocentric correction to NON-normalized spectra?"}
bool 	norm 			{no,prompt="normalize spectra?"}
bool 	ncombine 		{no,prompt="combine normalized spectra?"}
bool 	nrvcorr 		{no,prompt="Apply heliocentric correction to normalized spectra?\n\n# TASK PARAMETERS\n"}
###################### PARAMETERS
string	nfunction		{"legendre",prompt="Continuum fitting function",enum="legendre|spline3|chebyshev"}
int	norder		{5,prompt="Order of continuum fitting function"}
string	t_funct		{"spline3",prompt="Trace apertures fitting fucntion",enum="legendre|spline3|chebyshev"}
int	t_order		{5,prompt="Order of apertures fitting function"}
bool		edit_o	{no,prompt="Edit object apertures?"}
bool		review_o	{no,prompt="Review object apertures?"}
#######
begin
#bool norm
string spec,idname,oname,ecname,oap,refname,wap,fap,nnap,idlinefold
string iddir, inflat, incomp, inobject, hid, iddat, utmid, cfunction
string pathfolder, tfunct, name, email
int nap
real utmidhr, utmidmin, utmidsec, utstart, expt
####### LACOS variables
char inputCR, outputCR, outmaskCR
real gainCR, readnCR, gainh,readnh
int xorderCR,yorderCR
real sigclipCR,sigfracCR,objlimCR
int niterCR
bool verboseCR
char blk,lapla,deriv2,med5,sub5,sub5abs,med3,med7,imstatout,inputmask
char noise,sigmap,firstsel,starreject,gfirstsel,finalsel
char oldoutput,sigima,galaxy,skymod
char kernel,gkernel
real midpt,skylev,sig,sigcliplow
char l1,l2
int i,corder,torder
bool stop
int previous,npix
####### define variables
name="Mauricio Cabezas"
email='mauricio.cabezas@asu.cas.cz'
spec=input
idname=idtarget
oname=output
nap=napertures
idlinefold=idfolder
#print substr(input,0,-30)
cfunction=nfunction
corder=norder
tfunct=t_funct
torder=t_order
!basename "$PWD" > folder.tmp
cat folder.tmp | scan (pathfolder)
!rm folder.tmp
##########package
noao
imred
ccdred
kpnoslit
onedspec
rv
system
################################################## CALIBRATION
#- Create list for each type of calibration and scicence spec
#check object  name in header
# imhead spec1.fit l+ | page
# hselect e*.fit $I,object yes
if (orgfile==yes){
	!ls *.fit > all.spec
	hselect (images="@all.spec",fields="$I,OBJECT", exp=yes)
	hselect (images="@all.spec",fields="$I", exp='OBJECT=="flat"' ,> "flat.dat")
	hselect (images="@all.spec",fields="$I", exp='OBJECT=="zero"' ,> "zero.dat")
	hselect (images="@all.spec",fields="$I", exp='OBJECT=="comp"' ,> "comp.dat")
	hselect (images="@all.spec", fields="$I", exp="OBJECT!=""flat"" && OBJECT!=""comp"" && OBJECT!=""zero""", > "other.dat")
	###### Create folders
	mkdir (newdir="zero")
	mkdir (newdir="flat")
	mkdir (newdir="comp")
	mkdir (newdir="other")
# system inlcude in
# 	###### move files
	move (files="@zero.dat,zero.dat",newdir="zero/")
	move (files="@flat.dat,flat.dat",newdir="flat/")
	move (files="@comp.dat,comp.dat",newdir="comp/")
	move (files="@other.dat,other.dat",newdir="other/")
	cd "other/"
	hselect (images="@other.dat",fields="$I,OBJECT",exp=yes ,> "tmp.list")
	!awk '{print "oesred (input=\""$1"\",output="}' tmp.list > tmp1.list
	!awk '{for (i=2; i<=NF;i++) {printf("%s",$i)} print ""}' tmp.list > tmp2.list
	!awk '{for(i=2;i<=NF;i++) printf $i" ";print ")"}' tmp.list | awk '{print ",idtarget="$0}' > tmp3.list
	!paste -d'\0' tmp1.list tmp2.list tmp3.list > ../red.list
	!rm tmp*.list
	cd "../"
	#copy (input="other/"//spec,output=spec)

}

###### ZEROCOMBINE
if (access("zero/Zero.fit")){
	zerocomb=no
}
if (zerocomb==yes){
	cd "zero/"
	unlearn zerocombine
	zerocombine.reject="minmax"
	zerocombine.rdnoise= "READNOIS"
	zerocombine.gain   = "GAIN"
	zerocombine (input="@zero.dat",output="Zero.fit")
	cd "../"
}


######TRIM FLAT & COMP
# cd "flat/"
# list = "flat.dat"
# while (scan (s1, $<"flat/flat.dat") != EOF)
# 	print (s1)

if (trimcal==yes){
	unlearn ccdproc
	ccdproc.trimsec = "[2:2035,*]"
#	ccdproc.trimsec = "[5:2025,800:1500]"
	ccdproc.trim = yes
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.darkcor= no
	ccdproc.zerocor=no
	ccdproc.flatcor=no
	#
	cd "flat/"
	ccdproc.ccdtype = "flat"
	ccdproc (images="@flat.dat",output="T@flat.dat")
	cd "../"
	#
	cd "comp/"
	ccdproc.ccdtype = "comp"
	ccdproc (images="@comp.dat",output="T@comp.dat")
	cd "../"

}

###### SUBTRACT zero
if (zerocorcal==yes){
	unlearn ccdproc
	ccdproc.ccdtype="zero"
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.darkcor= no
	ccdproc.zerocor=no
	ccdproc.flatcor=no
	#
	cd "flat/"
	ccdproc.ccdtype="flat"
	ccdproc.zerocor=yes
	ccdproc.zero="../zero/Zero.fit"
	!sed 's/e2/Te2/g' flat.dat > tflat.dat
		if (iftrimc==yes){
		ccdproc (images="@tflat.dat",output="Z@tflat.dat")
		} else {
		ccdproc (images="@flat.dat",output="Z@flat.dat")
		}
	cd "../"
	#
	cd "comp/"
	ccdproc.ccdtype="comp"
	ccdproc.zerocor=yes
	ccdproc.zero="../zero/Zero.fit"
	!sed 's/e2/Te2/g' comp.dat > tcomp.dat
			if (iftrimc==yes){
		ccdproc (images="@tcomp.dat",output="Z@tcomp.dat")
		} else {
		ccdproc (images="@comp.dat",output="Z@comp.dat")
		}
	cd "../"
}
############

###### COMBINE - comp/lamp
if (compcomb==yes){
	cd "comp/"
	unlearn imcombine
	imcombine.reject = "none"
	imcombine.lsigma = 3
	imcombine.hsigma = 3
	imcombine.rdnoise= "READNOIS"
	imcombine.gain   = "GAIN"
  	imcombine.scale = "exposure"
  	imcombine.expname="EXPTIME"
#	imcombine (input="@comp.dat",output = "comp.fits")
	!sed 's/Te/ZTe/g' tcomp.dat > ztcomp.dat
		if (iftrimc==yes){
		imcombine (input="@ztcomp.dat",output = "../ZTcomp.fit")
		} else {
		imcombine (input="Z@comp.dat",output = "../Zcomp.fit")
		}
cd "../"
}

###### COMBINE flat
if (flatcomb==yes){
	cd "flat/"
	!sed 's/Te/ZTe/g' tflat.dat > ztflat.dat
		if (iftrimc==yes){
		imcombine (input="@ztflat.dat",output = "../ZTflat.fit")
		} else {
		imcombine (input="Z@flat.dat",output = "../Zflat.fit")
		}
	cd "../"
}

# ###### SUBTRACT zero
# if (zerocorcal==yes){
# 	imarith (operand1="flat/flat.fits", op="/", operand2="zero/Zero.fits", result="Zflat.fits")
# 	imarith (operand1="comp/comp.fits", op="/", operand2="zero/Zero.fits", result="Zcomp.fits")
# }

######IF TRIMMING
if (iftrimc==yes){
	inflat="ZTflat.fit"
	incomp="ZTcomp.fit"
} else {
	inflat="Zflat.fit"
	incomp="Zcomp.fit"
}


#####APERTURES - APALL FLAT
if (flatapall==yes){
	echelle
	unlearn apall
	apall.format = "echelle"

	apall.extras=no
	apall.extract=yes

	apall.nsum=15

	apall.lower=-5
	apall.upper=5
    apall.b_order=3
	apall.b_sample="-10:-6,6:10"

	apall.nfind=nap
	#apall.minsep=10
	apall.minsep=5
	apall.maxsep=1000


	apall.ylevel = 0.04
	apall.bkg=yes
	#apall.bkg=no

	apall.t_nsum = 10
	apall.t_function = tfunct
	apall.t_niter=100
	apall.t_order= torder

	apall.clean=no
	apall.readnoi= 0
	apall.gain   = 1



	#apall.width=9
	#apall.width=5
#	apall.weights = "none"
	apall (input=inflat, output = "A"//inflat)
}

#####APERTURES - COMP
if (compapall==yes){
	apall.referen=inflat
	apall.format = "echelle"
	apall.find=no
	apall.recente=no
	apall.resize=no
	apall.trace=no
	apall.fittrace=no
	apall.extras=no
	apall.ylevel = 0.04
	apall.extract=yes
	apall (input=incomp, output="A"//incomp//".fit")
}

if (iddatabase==yes){
	unlear directory
	directory.sort=yes
	directory "identify/" | scan (iddir)
	if (iddir=="no"){
		mkdir (newdir="identify")
		mkdir (newdir="identify/database")
		copy (input="../"//idlinefold//"/*",output="identify/database/")
	}
	# directory idlinefold//"/" | scan (iddat)
	# print (iddat)
	# if (iddir=="no"){
	#  	mkdir (newdir="identify")
	# }
}

##### IDENTIFY 1D
## :/xflip YES
## :labels
if (idencomp==yes){
	unlearn directory
	unlearn scopy
	directory.sort=yes
	directory "identify/" | scan (iddir)
	if (iddir=="no"){
		mkdir (newdir="identify")
	}
	copy (input="A"//incomp,output="identify/")
	cd "identify/"
	lpar scopy
	scopy.format="onedspec"
	scopy (input="A"//incomp, output="iazcomp")
	print "second"
	unlearn refspectra
	unlearn hedit
	hedit.addonly=yes
	hedit.verify=no
	hedit.show=no
	for (i=1; i <=nap; i+=1) {
		printf ("iazcomp.00%02d.fit\n",(i)) | scan(ecname)
		if (iddatabase==yes){
			printf ("iazcomp.00%02d\n",(i)) | scan(refname)
			hedit (images=ecname, fields="REFSPEC1", value=refname)
		}
		identify.coordli="linelists$thar.dat"
		#lpar identify
		identify (images=ecname)
	}
	cd "../"
}


################################################## OBJECT
###### create obeject folder
###### always
unlear directory
directory.sort=no
directory oname//"/" | scan (iddir)
if (iddir=="no"){
	mkdir (newdir=oname)
	copy (input="other/"//spec,output=oname//"/"//spec)
}
#CALC UTMIDDLE
hselect (images="other/"//spec,fields="TM_START", exp=yes) | scan (utstart)
hselect (images="other/"//spec,fields="EXPTIME", exp=yes) | scan (expt)
utmidhr=int((utstart + expt/2)/3600)
utmidmin=int((((utstart + expt/2)/3600)-utmidhr)*60)
utmidsec=int(((((utstart + expt/2)/3600)-utmidhr)*60 - utmidmin)*60)
utmid = (utmidhr//":"//utmidmin//":"//utmidsec)
printf ("%d:%d:%d\n",utmidhr,utmidmin,utmidsec) | scan (utmid)
############################
if (trimob==yes){
	unlearn ccdproc
	ccdproc.trimsec = "[2:2035,*]"
	ccdproc.trim = yes
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.darkcor= no
	ccdproc.zerocor=no
	ccdproc.flatcor=no
	#
	cd (oname)
	ccdproc.ccdtype = "object"
	ccdproc (images=spec, output="T"//spec)
	cd "../"
}

######IF TRIMMING
if (iftrimo==yes){
	inobject="T"//spec
} else {
	inobject=spec
}

###### ZERO COR

if (zerocorob==yes){
	cd (oname)
	#imarith (operand1=spec, op="/", operand2="../zero/Zero.fits", result="Z"//spec)

	unlearn ccdproc
	ccdproc.ccdtype="zero"
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.darkcor= no
	ccdproc.zerocor=no
	ccdproc.flatcor=no
	#
	ccdproc.ccdtype="object"
	ccdproc.zerocor=yes
	ccdproc.zero="../zero/Zero.fit"
	ccdproc (images=inobject,output="Z"//inobject)
	cd "../"
}

#####COSMIC RAYS - COMP
if (crays==yes){
	stsdas
#######
#read gain
cd (oname)
hselect (images="Z"//inobject,fields="GAIN", exp=yes) | scan (gainh)
hselect (images="Z"//inobject,fields="READNOIS", exp=yes) | scan (readnh)
#print (gainh)
inputCR="Z"//inobject
outputCR="CrZ"//inobject
outmaskCR="MCrZ"//inobject
gainCR = gainh # 2 #3
readnCR = readnh #2
xorderCR = 3
yorderCR = 3
sigclipCR = 4.5
sigfracCR = 0.3
objlimCR = 3
niterCR = 5
verboseCR = yes
##
 if (!deftask("imcalc")) error(2,"Load package stsdas")
 if (gainCR<=0) error(2,"Gain is required")

 if (verboseCR) {
  print("")
  print("_______________________________________________________________________________")
  print("")
  print("                 L.A.Cosmic: Laplacian cosmic ray removal")
  print("")
  print("                   P. van Dokkum, 2001, PASP 113, 1420")
  print("")
  print("                  Spectroscopy version 1.0  (April 2001)")
  print("_______________________________________________________________________________")
  print("")
  }

 # make temporary files

 blk = mktemp("lacos")
 lapla = mktemp("lacos")
 deriv2 = mktemp("lacos")
 kernel = mktemp("lacos")
 gkernel=mktemp("lacos")
 med3 = mktemp("lacos")
 med5 = mktemp("lacos")
 med7 = mktemp("lacos")
 sub5 = mktemp("lacos")
 sub5abs = mktemp("lacos")
 imstatout = mktemp("lacos")
 noise = mktemp("lacos")
 sigmap = mktemp("lacos")
 firstsel = mktemp("lacos")
 starreject=mktemp("lacos")
 gfirstsel = mktemp("lacos")
 finalsel = mktemp("lacos")
 inputmask = mktemp("lacos")
 oldoutput = mktemp("lacos")
 sigima = mktemp("lacos")
 galaxy=mktemp("lacos")
 skymod = mktemp("lacos")

 # set some parameters in standard IRAF tasks
 convolve.xkernel=1
 convolve.ykernel=1
 convolve.bilinear=yes
 convolve.radsym=no

 # create Laplacian kernel

 print("0 -1 0;",>> kernel)
 print("-1 4 -1;",>> kernel)
 print("0 -1 0;",>> kernel)

 # create growth kernel

 print("1 1 1;",>> gkernel)
 print("1 1 1;",>> gkernel)
 print("1 1 1;",>> gkernel)

 # initialize loop

 i=1
 stop=no
 previous=0

 imcopy(inputCR,oldoutput,verb-)
 imcopy(inputCR,outmaskCR,verb-)
 imreplace(outmaskCR,0,upper=INDEF,lower=INDEF)

# subtract object spectra if desired

  if (xorderCR>0) {
   if (verboseCR) {
    print("Subtracting object spectra")
    print("  fit order = "//xorderCR)
    print("")
    }
   fit1d(oldoutput,galaxy,"fit",axis=1,order=xorderCR,func="leg",low=4.,
        high=4.,nav=1,inter-,sample="*",niter=3,grow=0,cursor="")
   imarith(oldoutput,"-",galaxy,oldoutput)
   }

  else {
   imcopy(oldoutput,galaxy,verb-)
   imreplace(galaxy,0,lower=INDEF,upper=INDEF)
   }

  if (yorderCR>0) {
   if (verboseCR) {
    print("Subtracting sky lines")
    print("  fit order = "//yorderCR)
    print("")
    }
   fit1d(oldoutput,skymod,"fit",axis=2,order=yorderCR,func="leg",low=4.,high=4.,
        inter-,sample="*",nav=1,niter=3,grow=0,cursor="")
   imarith(oldoutput,"-",skymod,oldoutput)
   }

  else {
   imcopy(oldoutput,skymod,verb-)
   imreplace(skymod,0,lower=INDEF,upper=INDEF)
   }

  # add object spectra to sky model

  imarith(skymod,"+",galaxy,skymod)

# start iterations

 while(!stop) {

  if (verboseCR) {
   print("")
   if (i<10) print("_______________________________ Iteration "//i//" ___________________________________")
   if (i>9) print("_______________________________ Iteration "//i//"___________________________________")
   print("")
   print("")
   }

  # add median of residuals to sky model

  median(oldoutput,med5,5,5,zlor=INDEF,zhir=INDEF,verb-)
  imarith(skymod,"+",med5,med5)

  # take second-order derivative (Laplacian) of input image
  # kernel is convolved with subsampled image, in order to remove negative
  # pattern around high pixels

  if (verboseCR) {
   print("Convolving with Laplacian kernel")
   print("")
   }
  blkrep(oldoutput,blk,2,2)
  convolve(blk,lapla,kernel,xkernel=1,ykernel=1)
  imreplace(lapla,0,upper=0,lower=INDEF)
  blkavg(lapla,deriv2,2,2,option="average")

  if (verboseCR) {
   print("Creating noise model using:")
   print("  gain = "//gainCR//" electrons/ADU")
   print("  readnoise = "//readnCR//" electrons")
   print("")
   }

  # create noise model

  imcalc(med5,noise,"sqrt(im1*"//gainCR//" + "//readnCR//"**2)/"//gainCR,verb-)
  imreplace(med5,0.00001,upper=0,lower=INDEF)

  # divide Laplacian by noise model

  imarith(deriv2,"/",noise,sigmap)

  # Laplacian of blkreplicated image counts edges twice:

  imarith(sigmap,"/",2.,sigmap)

  # removal of large structure (bright, extended objects)

  imdel(med5)
  median(sigmap,med5,5,5,zlo=INDEF,zhi=INDEF,verb-)
  imarith(sigmap,"-",med5,sigmap)

  # find all candidate cosmic rays
  # this selection includes sharp features such as stars and HII regions

  if (verboseCR) {
   print("Selecting candidate cosmic rays")
   print("  sigma limit = "//sigclipCR)
   print("")
   }
  imcopy(sigmap,firstsel,verb-)
  imreplace(firstsel,0,upper=sigclipCR,lower=INDEF)
  imreplace(firstsel,1,lower=0.1,upper=INDEF)

  # compare candidate CRs to median filtered image
  # this step rejects bright, compact sources from the initial CR list

  if (verboseCR) {
   print("Removing suspected emission lines")
   print("  selecting cosmic rays > "//objlimCR//" times object flux")
   print("")
   }

  # subtract background and smooth component of objects

  median(oldoutput,med3,3,3,zlo=INDEF,zhi=INDEF,verb-)
  median(med3,med7,7,7,zlo=INDEF,zhi=INDEF,verb-)
  imarith(med3,"-",med7,med3)
  imarith(med3,"/",noise,med3)
  imreplace(med3,0.01,upper=0.01,lower=INDEF)

  # compare CR flux to object flux

  imcalc(firstsel//","//sigmap//","//med3,starreject,"(im1*im2)/im3",verb-)

  # discard if CR flux <= objlim * object flux

  imreplace(starreject,0,upper=objlim,lower=INDEF)
  imreplace(starreject,1,lower=0.5,upper=INDEF)
  imarith(firstsel,"*",starreject,firstsel)

  # grow CRs by one pixel and check in original sigma map

  convolve(firstsel,gfirstsel,gkernel)
  imreplace(gfirstsel,1,lower=0.5,upper=INDEF)
  imarith(gfirstsel,"*",sigmap,gfirstsel)
  imreplace(gfirstsel,0,upper=sigclipCR,lower=INDEF)
  imreplace(gfirstsel,1,lower=0.1,upper=INDEF)

  # grow CRs by one pixel and lower detection limit

  sigcliplow = sigfracCR * sigclipCR

  if (verboseCR) {
   print("Finding neighbouring pixels affected by cosmic rays")
   print("  sigma limit = "//sigcliplow)
   print("")
   }

  convolve(gfirstsel,finalsel,gkernel)
  imreplace(finalsel,1,lower=0.5,upper=INDEF)
  imarith(finalsel,"*",sigmap,finalsel)
  imreplace(finalsel,0,upper=sigcliplow,lower=INDEF)
  imreplace(finalsel,1,lower=0.1,upper=INDEF)

  # determine number of CRs found in this iteration

  imdel(gfirstsel)
  imcalc(finalsel//","//outmaskCR,gfirstsel,"(1-im2)*im1",verb-)
  imstat(gfirstsel,fields="npix",lower=0.5,upper=INDEF,for-) | scan(npix)

  # create cleaned output image; use 3x3 median with CRs excluded

  if (verboseCR) {
   print("Creating output:")
   print("  bad pixel mask: "//outmaskCR)
   print("  cleaned image: "//outputCR)
   print("")
   }

  imdel(med5)
  imarith(outmaskCR,"+",finalsel,outmaskCR)
  imreplace(outmaskCR,1,lower=1,upper=INDEF)
  imcalc(outmaskCR//","//oldoutput,inputmask,"(1-10000*im1)*im2",verb-)
  median(inputmask,med5,5,5,zloreject=-9999,zhi=INDEF,verb-)
  imarith(outmaskCR,"*",med5,med5)
  if (i>1) imdel(outputCR)
  imcalc(outmaskCR//","//oldoutput//","//med5,outputCR,"(1-im1)*im2+im3",verb-)

  # add sky and object spectra back in

  imdel(oldoutput)
  imcopy(outputCR,oldoutput,verb-)

  imarith(outputCR,"+",skymod,outputCR)

  # cleanup and get ready for next iteration

  if (verboseCR) {
   print("Cleaning up")
   print("")
   }

   if (verboseCR) {
   print("")
   print(npix//" cosmic rays found in iteration "//i)
   print("")
   }
print(npix//" cosmic rays found in iteration "//i)
  if (npix==0) stop=yes

  i=i+1
  if (i>niterCR) stop=yes

  # delete temp files

  imdel(blk//","//lapla//","//deriv2//","//med5)
  imdel(med3//","//med7//","//noise//","//sigmap)
  imdel(firstsel//","//starreject//","//gfirstsel)
  imdel(finalsel//","//inputmask)

  }

 imdel(oldoutput//","//skymod//","//galaxy)
 delete(kernel//","//gkernel)
###############################################################


cd "../"
}

######IF CRAYS
if (ifcrays==yes){
	inobject="CrZ"//inobject
} else {
	inobject="Z"//inobject
}
#####APERTURES - OBJECT
if (objectapall==yes){
	apall.referen=inflat
	apall.format = "echelle"
	apall.find=no
	apall.recente=no
	apall.resize=no
	apall.trace=no
	apall.fittrace=no
	apall.extras=no
	apall.extract=yes
	apall.edit=edit_o
	apall.review=review_o
	## check database
	unlear directory
	directory.sort=yes
	directory oname//"/database/" | scan (iddir)
	if (iddir=="no"){
		mkdir (newdir=oname//"/database")
		copy (input="database/*",output=oname//"/database/")
	}
	##
	cd (oname)
	apall (input=inobject, output="A"//inobject)
	cd "../"
}

###### FLAT correction
#print ("FA"//inobject)
if (flatcor==yes){
	imarith (operand1=oname//"/A"//inobject, op="/", operand2="A"//inflat, result=oname//"/FA"//inobject)

	#imcopy("FA"//inobject,oldoutput,verb-)
	#move (files="FA"//inobject,newdir=oname//"/")

}

# ##### SETJD + HELIOCOR
if (heliocor==yes){
	cd (oname)
	#
	print (utmid,inobject)
	unlearn hedit
	hedit.addonly=no
	hedit.verify=no
	hedit.show=no
	hedit (images="FA"//inobject, fields="UT", value=utmid)
	##hedit (images="N"//oname//".fit", fields="UT", value=utmid)
	# hedit (images="N"//oname//".fits", fields="UTMIDDLE", value=utmid)
	unlearn setjd
	setjd.observatory="obspars"
	observatory.observatory = "ondrejov"
	observatory.name = "Ondrejov Observatory"
	observatory.longitude = 345:12:59
	observatory.latitude = 49:54:38
	observatory.altitude = 528
	observatory.timezone = -1
	setjd.epoch="epoch"
	setjd.time="UTMIDDLE"
	setjd (images="FA"//inobject)
	rv
	unlearn rvcorrect
	rv.rvcorrect.input=yes
	rv.rvcorrect.imupdate=yes
	rv.rvcorrect.epoch=2000.
	rv.rvcorrect.observatory="ondrejov"
	rv.rvcorrect (images="FA"//inobject)
	##save header
	print ("FA"//inobject)
	#imheader (images="FA"//inobject)
	imhead ("FA"//inobject, l+ ,> oname//".hd")
	cd "../"
}


system
##### REFSPEC 1D
if (idref==yes){
	unlear directory
	directory.sort=no
	directory oname//"/" | scan (iddir)
	if (iddir=="no"){
		mkdir (newdir=oname)
	}
	cd oname//"/"
	unlear onedspec
	scopy.format="onedspec"
	scopy (input="FA"//inobject, output="ap")
	unlearn refspectra
	dispcor.w1=INDEF
	dispcor.w2=INDEF
	dispcor.nw=INDEF
	dispcor.flux=no
	cd "../identify/"
	system
	 for (i=1; i <=nap; i+=1) {
	 	printf ("ap.00%02d.fit\n",(i)) | scan(oap)
	 	printf ("iazcomp.00%02d.fit\n",(i)) | scan(ecname)
		refspectra.sort="epoch"
		refspectra.group="epoch"
		refspectra.answer=yes
		refspectra.confirm=no
		system.move (files="../"//oname//"/"//oap, newdir="../identify/")
		refspectra (input=oap, referen=ecname)
	 	dispcor (input=oap, output="w"//oap)
	 	system.move (files=oap,newdir="../"//oname)
	  	system.move (files="w"//oap,newdir="../"//oname)
	}
	cd "../"
	cd oname//"/"
	wspectext.header=no
	mkdir (newdir="wap_asc")
	for (i=1; i <=nap; i+=1) {
		printf ("ap.00%02d\n",(i)) | scan(oap)
		wspectext (input="w"//oap//".fit",output="w"//oap//".asc")
		system.move (files="w"//oap//".asc",newdir="wap_asc")

	}
	system.move (files=oname//".hd",newdir="wap_asc")
	cd "wap_asc"
	!ls wap*.asc > norm.list
	#### create python script
	cd "../../"
}
####  COMBINE
#ncombine=no
if (combine==yes){
	cd (oname)
	hselect (images=inobject,fields="GAIN", exp=yes) | scan (gainh)
	hselect (images=inobject,fields="READNOIS", exp=yes) | scan (readnh)
	scombine.aperture = "*"
	scombine.group = "apertures"
	scombine.combine = "average"
	scombine.reject = "none"
	scombine.mclip = yes
	scombine.rdnoise = readnh
	scombine.gain = gainh
	scombine.blank=1.
	scombine (input="wap.00*.fit", output="C-"//oname//"_"//pathfolder//".fit")
	#re-set EXPTIME
	##
	cd "../"
}

####  RVcorr helio
if (rvcorr==yes){
	cd (oname)
	dopcor.isvelocity=yes
	dopcor (input="C-"//oname//"_"//pathfolder//".fit", output="DC-"//oname//"_"//pathfolder//".fit", redshift="-vhelio")
	hedit.addonly=yes
	hedit.verify=no
	hedit.show=no
	hedit (images="DC-"//oname//"_"//pathfolder//".fit", fields="COMMENT1", value="Reduced by "//name//"using OESRED script")
	hedit (images="DC-"//oname//"_"//pathfolder//".fit", fields="COMMENT2", value="email: "//email)
	cd "../"
}

#norm=no
##### normalization
if (norm==yes){
	cd (oname)
	unlearn continuum
	unlearn scombine
	continuum.type="fit"
	continuum.function=cfunction
	continuum.order=corder
	continuum.naverage=10
	continuum.markrej=no
	continuum.niterat=2000
	continuum.high_re=2
	continuum.low_re=1.5
	continuum.grow=0
	for (i=1; i <=nap; i+=1) {
	  	printf ("wap.00%02d.fit\n",(i)) | scan(wap)
	  	printf ("fap.00%02d.fit\n",(i)) | scan(fap)
	  	printf ("nap.00%02d.fit\n",(i)) | scan(nnap)
	  	print (wap, fap, nnap)
  	continuum (input=wap, output=fap)
	}
	cd "../"
}

#### N COMBINE
#ncombine=no
if (ncombine==yes){
	cd (oname)
	hselect (images="Z"//inobject,fields="GAIN", exp=yes) | scan (gainh)
	hselect (images="Z"//inobject,fields="READNOIS", exp=yes) | scan (readnh)
	scombine.aperture = "*"
	scombine.group = "all"
	scombine.combine = "sum"
	scombine.reject = "none"
	scombine.mclip = yes
	scombine.rdnoise = readnh
	scombine.gain = gainh
	scombine.blank=1.
	scombine (input="wap.00*.fit", output="wap.fit")
	scombine (input="fap.00*.fit", output="fap.fit")
	imarith (operand1="wap.fit", op="/", operand2="fap.fit", result="CN-"//oname//"_"//pathfolder//".fit")
	#re-set EXPTIME
	unlearn hedit
	hedit.addonly=no
	hedit.verify=no
	hedit.show=no
	hedit (images="N"//oname//pathfolder//".fit", fields="EXPTIME", value=expt)
	##
	cd "../"
}

####  RVcorr helio normalized
if (nrvcorr==yes){
	cd (oname)
	dopcor.isvelocity=yes
	dopcor (input="CN-"//oname//"_"//pathfolder//".fit", output="DCN-"//oname//"_"//pathfolder//".fit", redshift="-vhelio")
#	system.move (files=oname//"-"//pathfolder//".fit",newdir="../")
	##
	cd "../"
}

end
