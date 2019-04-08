pro star_calib_millstone_tim1
 
   ;+
   ; NAME:
   ;       STAR_CALIB_MILLSTONE_TIM1
   ; PURPOSE:
   ;       To locate different stars in all-sky images in order to determine
   ;       the distortion function and of the lens used and the responsivity 
   ;       of the sky
   ; CALLING SEQUENCE:
   ;       STAR_CALIB_MILLSTONE_TIM1
   ; INPUTS:
   ;       Center: The x and y pixel locations of the center of the image
   ;       Rotation Angle: Amount the image needs to be rotated (in degrees)
   ;       Wave: Wavelength of the main image (in Angstroms)
   ;       thedate: Date the images were taken
   ;       deg2pixCoef and pix2degCoef: Initial guesses for the distortion function
   ;                             (determines where the program will look for stars)
   ;       MaxPixRad: Maximum pixel radius from the center inside which the program
   ;                  will attempt to locate stars
   ;       FWHM: Full Width Half Maximum (in Angstroms) of the filter
   ;       TrCo: Peak transmission coefficient of the filter
   ;       TrCoWave: Transmission coefficient at the main wavelength
   ;       SpW: Spectral Width of the filter
   ; OUTPUTS
   ;       Creates two postscript files: One displays the distortion function 
   ;       and its polynomial coefficients, and the other displays the responsivity.
   ;       It also saves the coefficients of the distortion function, a dark image
   ;       for the site, and an outfile, which is a table that contains the star
   ;       number, flux, average dn, zenith distance, zenith distance in pixels,
   ;       Rayleigh-seconds dn and error for each star and image. 
   ; REVISION HISTORY:
   ;       11 Jun 2017: Edited by Timothy MacDonald. Original, unedited code taken 
   ;                    from universal_star_calib2_jeff.pro
   ;-

;============================================
; EDIT THESE VARIABLES FOR EACH DISTINCT RUN
;============================================

  thedate='Nov1118'   ; Set the date for the set of data to be analyzed
  xofst = 254         ; x center of all sky image (assumed to be the zenith)
  yofst = 258         ; y center of all sky image (assumed to be the zenith)
  rot_axis= -6.75      ; Rotation angle (in degrees)
  maxPixRad = 210     ; Pixel radius inside which the program will attempt to locate stars
  
  deg2pixCoef = [0.0186579,3.07649,0.00103662,-7.68886e-005]      ; Initial guess for the unwarping function
  pix2degCoef = [-.0242126,0.338231,-2.97539e-004,2.13789e-006]   ; Initial guess for the unwarping function
  
;===============================================
; MILLSTONE SPECIFIC DATA (Not Usually Changed)
;===============================================

  UberStarList = 'McStandardStarList6300b.txt'                 ; Name of the file with the list of stars to be tested
  thisFolder = 'C:\Program Files\Exelis\IDL82\lib\StarCalib\'  ; Path to the output images 
  wave = '6300'                          ; Main wavelength (in Angstroms)
  FWHM = 15.1                            ; Full Width at Half Maximum of main wavelength Filter (in Angstroms)
  TrCo = .77                             ; Peak Transmission Coefficient of main wavelength Filter
  TrCoWave = TrCo                        ; Transmission coefficient at the main wavelength (if exact number is unknown, assume it is the peak value)
  SpW = FWHM*TrCo                        ; Area of the filter (If an exact number is not known, estimate by multiplying FWHM and TrCo)
  clon =  -71.48444                      ; Longitude
  clat =  42.611667                      ; Latitude
  camh = 0.146                           ; height above sea level
  therot = 2                             ; angle to rotate images by (0=0 deg, 1=90 deg, 2=180 deg, etc.)
  path2img = 'C:\Millstone\'+thedate+'\' ; path to raw images
  curYr = strmid(theDate,5,2)            ; Extracts the current year from the date
  daynum = getdaynumber(thedate)         ; Extracts the current day from the date
  catfile = theDate + '.cat'             ; Creates a string for a name of a .cat file
  ImgXsize = 512                         ; x-dimension of the images (in pixels)
  ImgYsize = 512                         ; y-dimension of the images (in pixels)
    restore,"C:\Program files\Exelis\IDL82\lib\starcalib\Millstonevigarr.sav"  ; Restore and shift the vignette file to the center
    vig=shift(newvig,xofst-256,yofst-256)
  
  ;*************************
  ; Begin Processing Images
  ;*************************
  
  count=0                  ; Sets count equal to 0, to be used later in for loop
  err=fltarr(500)          ; array to hold errors (deltas of predicted vs actual star x,y positions)
  CalibData=fltarr(2000,7) ; make array to hold,starnum,starflux,starAvgdn,zenith dist,zd_pix,raySecDn,err...for each star and image
   
    set_plot,'win'
  
  pix2degprime=[pix2degCoef(1),2.0*pix2degCoef(2),3.0*pix2degCoef(3)]
  
  stfrac = getSidDayFrac(curYr)                ; Gets sidereal time at the beginning of the current year   
           
  ;*** Conversion factor to change deg to radians ***
  q=double(!pi)/180.0
  clat=clat*q
  rot_axis=rot_axis*q
  
  ;*** Read in Star List ***
  NStars =  file_lines(thisFolder + UberStarList)
  junk  = ''
  filex = ''
  NStars =  NStars - 3       ; Subtract by 3, because there are 3 header lines in the file
  ; Create arrays of the correct dimensions and string, integer, or double-precision floating variables
  names   = strarr(NStars)   ; Star Names
  rah   = dblarr(NStars)     ; Right ascension (in hours)
  dec   = dblarr(NStars)     ; Declination (in degrees)
  flx6300 = dblarr(NStars)   ; Flux at 6300 Angstroms (in phot/cm2/A/s)
  width = intarr(NStars)     ; Width of the star (in pixels)
  starnum = intarr(NStars)   ; Number associated with the star
  lnParts = strarr(10)
  StarCount   = 0
  Inter_Flux  = dblarr(NStars)
  
  openr,Runit,thisFolder + UberStarList,/GET_LUN
  readf,Runit,junk      ; ignore first three lines of header info
  readf,Runit,junk      
  readf,Runit,junk
  while not eof(Runit) do begin
    readf,Runit,filex
    if (strmid(filex, 0, 1) eq ';') then continue
    lnParts = STRSPLIT(filex, ' ', /EXTRACT)
    ; Extracts the contents of each column in the UberStarList
    starNum(StarCount)= lnParts(0)            ; Star number
    names(StarCount)  = lnParts(1)            ; Star name
    width(StarCount)  = fix(lnParts(2))       ; Star width
    rah(StarCount)    = double(lnParts(3))    ; Right ascension
    dec(StarCount)    = double(lnParts(4))    ; Declination
    flx6300(StarCount) = double(lnParts(5))   ; Flux @ 6300

    StarCount = StarCount+1                   
  endwhile
  FREE_LUN,Runit
  close, /all
  ;stop
  
  radeg = rah*15.0  ; convert RA from hours to degrees
  dec = dec*q       ; convert DEC from degrees to radians
  
  ;Get a median dark from the night in question to subtract from each image
  exttime = 120.
  AMedianDark = median_image(path2img,thedate,'DARK',exptime=exptime)
  amediandark = rotate(amediandark,therot)
  save,filename='c:\Program Files\Exelis\IDL82\lib\Starcalib\Millstone\Millstonedark.sav',amediandark
 
  ; Set up windows 0 and 1 (they will later display sky images)
  window,0,xs=ImgXsize,ys=ImgYsize
  !p.position=[0,0,1,1]
  window,1,xs=ImgXsize,ys=ImgYsize
  !p.position=[0,0,1,1]
  
    ; Read in images for the night
  line = ''
  openr,mylun,path2img+catfile,/GET_LUN
  imgname=strarr(200)
  element=0
  ; ******** Creates an image array from a .cat file if no list file exists ***********
   while not eof(mylun) do begin
    readf,mylun,line
    temp = (StrSplit(line,/RegEx,/Extract))
    temp = fix_filename(temp)
    thetime = fname2dec(temp(0))
    ;stop
    imagename = temp(0)          ; Reads image name
    ;stop
    if strmid(temp(4),0,4) eq wave then begin
      exptime = fix(temp(5))
    endif
  if file_test(path2img+'list'+wave+'.txt') eq 0 then begin
    if strmid(temp(4),0,4) eq wave then begin                      ; Finds images at main wavelength
      imgname(element)=imagename
      element=element+1
    endif
  endif
 endwhile
 ; ******** Creates an image array from a list file if it exists ************
  if file_test(path2img+'list'+wave+'.txt') eq 1 then begin
    list=path2img+'list'+wave+'.txt'   ; Defines 'list' as the list of pictures taken at the main wavelength

    openr,unit,list,/GET_LUN     ; Open the .txt file
 
    ; Add all of the image names in the file to the array
    while not eof(unit) do begin
      readf,unit,line
      imgname(element)=line
      element=element+1
    endwhile
    close,/all
  endif
  
  ; For loop to load the images and calibrate the stars
  for i=0,element-1 do begin
      ; print,'*********** '+imgname+' ***********************'    ; Print detected images at main wavelength
      load16,path2img+imgname[i],imgXsize,imgYsize,128,therot,img     ; Load detected images at main wavelength
      
      wait, 0.05
      
      img = float(img)           ; Turns img from a Long array to a floating point array
      
      ;***** Fix hot pixels if necessary *****
      ;img = sigma_filter(img,radius=1.0,n_sigma=4,/all,/iterate)
      ;tv,bytscl(img,min=0,max=3000)
      ;stop
      
      ;***** Subtracts dark array from the image *****
      img = img - amediandark
      ;tv,bytscl(img,min=0,max=3000)
      ;stop
      
      img = img/vig           ; Divide by vignette
      img = img/expTime       ; Divide by exposure time
      imgb = img
      imgc = img
      
      scale_img,img,.87,amin,amax,/showimg
      
      testcircle,img,xofst,yofst,maxPixRad,.25,img
      
      lst = getImgSidTime(imgname[i],daynum,stfrac,clon)
      
      ;getHR inputs: radeg,lst,dec,clat,q ; returns: azp,zd
      getHR,radeg,lst,dec,clat,q,azp,zd
      
      raddeg=zd/q                                                                                           ; Radius from center (in degrees)    
      radpix = (deg2pixCoef[0] + deg2pixCoef[1]*raddeg + deg2pixCoef[2]*raddeg^2 + deg2pixCoef[3]*raddeg^3) ; Radius from center (in pixels)
      decl = dec
      
      ; Finds properties of ONLY the stars that are above horizon obstructions
      maxRadIndex = where(radpix lt maxPixRad)
      radpix_p = radpix(maxRadIndex)            ; Pixel distance from center
      azp_p = azp(maxRadIndex)                  ; Azimuth (in radians)
      names_p = names(maxRadIndex)              ; Star names
      starnum_p = starnum(maxRadIndex)          ; Star number in UberStarList
      flx6300_p = flx6300(maxRadIndex)          ; Flux at 6300 Angstroms
      rah_p = rah(maxRadIndex)                  ; Right Ascension (in hours)
      dec_p = decl(maxRadIndex)                 ; Declination (in radians)
      zd_p=zd(maxradindex)                      ; Altitude (in radians)
      
      ;getStarPixLoc inputs: radpix,azp,rot_axis; returns: x,y
      getStarPixLoc,radpix_p,azp_p,rot_axis,xofst,yofst,x,y
      
      plot,[0,imgXsize-1],[0,ImgYsize-1],xstyle=5,ystyle=5,/nodata,/noerase
      
      NStars = n_elements(radpix_p)   ; Determines number of stars found within the maximum pixel radius

      ;********** Star Names ***********
      ;Put circle around standard star & add star names
      oplot,x,y,psym=6,symsize=1.5,color=254
      xcen=intarr(1)
      xcen(0)=xofst
      ycen=intarr(1)
      ycen(0)=yofst
      oplot,xcen,ycen,psym=5
      xyouts,x,y,'  '+names_p,/data,color=254
      xyouts,x,y,'     '+string(starnum_p),/data,color=254
      ;***************************************************
      
      ;********** Find the center of stars ***********
      StarLag = 8
      for k=0, NStars-1 DO BEGIN
        Strimgx = imgc
        xl = x(k)-StarLag
        xh = x(k)+StarLag
        yl = y(k)-StarLag
        yh = y(k)+StarLag
        CntrVal = max(Strimgx(xl:xh,yl:yh))             ; Central value is defined as the brightest point within the test square
        for ix=x(k)-StarLag, x(k)+StarLag do begin
          for iy=y(k)-StarLag, y(k)+StarLag do begin
            if (Strimgx(ix, iy) eq CntrVal) then begin
              x_new = ix
              y_new = iy
              break
            endif
          endfor
        endfor
        ;print, 'Old Coordinates of '+names(k)+' = (' + strcompress(x(k))+ ','+ strcompress(y(k))+')'   ; Prints expected coordinate of star
        ;print, 'New Coordinates of '+names(k)+' = ('+ strcompress(x_new)+ ',' + strcompress(y_new)+')' ; Prints actual coordinate of star in image
        err(k)=((x(k)-x_new)^2.0+(y(k)-y_new)^2.0)^0.5  ; get error distance
        x(k) = x_new      ; Redefines actual x-position of the star
        y(k) = y_new      ; Redefines actual y-position of the star
      endfor
      
      ;******** Plot empirically determined locations of stars ***********
      wset,0
      scale_img,imgb,.97,amin,amax,/showimg
      plot,[0,imgXsize-1],[0,ImgYsize-1],xstyle=5,ystyle=5,/nodata,/noerase
      oplot,x,y,psym=4,symsize=2,color=254
      xyouts,x,y,'  '+names_p,/data,color=254
      xyouts,x-6,y-12,string(starnum_p),/data,color=254
      imgc = imgb
   
      for jj = 0, NStars - 1 do begin
        wset,1
        imgb = imgc
        xc = x[jj]
        yc = y[jj]
        StarFLx = flx6300_p[jj]

        ; alternative method of getting tot dn from star
        
        sub1=Strimgx(xc-4:xc+4,yc-4:yc+4)  ; a 9x9 box around star
        sub2=Strimgx(xc-3:xc+3,yc-3:yc+3)  ; a 7x7 box around star
        sub1big=rebin(sub1,90,90,/sample)  ; take a look to see if star centered
           tv,bytscl(sub1big,min=0,max=300)
        tot1=total(sub1)                   ; sum of all pixels in big box
        tot2=total(sub2)                   ; sum of all pixels in small box
        dif=49.0*(tot1-tot2)/32.0          ; 49pix of background
        
        startot=tot2-dif ; just star [dn],
        StarAvgDN=startot
        
        raysecdn = getRayDn(pix2degCoef,pix2degPrime,xofst,yofst,xc,yc,q,staravgdn,starflx,spw,trcowave,camh)

        zdist_pix=((xofst-xc)^2.0+(yofst-yc)^2.0)^0.5      ; Pixel distance of star from center
        
        CalibData(count,*,*,*,*,*,*,*)=[starnum_p(jj),starflx,StarAvgDn,zd_p(jj)/q,zdist_pix,RaySecDn,err(jj)]
        
       count=count+1
      endfor
    endfor
  save,filename=thisfolder+'Millstone\Millstone_' + wave + '_' + thedate + '_outfile.sav',CalibData  ; Save array of CalibData
  FREE_LUN,mylun
  
  errindex=where(CalibData(*,6) lt 5.0); only use well identified (centered) stars

  ;******** Create responsivity graph in a PostScript file **********
  set_plot,'ps'
  device,/landscape,filename=thisfolder+'Millstone\Millstone_responsivity_'+wave+'_'+thedate+'.ps'
  plot,CalibData(errindex,3),CalibData(errindex,5),psym=2,yra=[0,100],xra=[0,90],ystyle=1,xstyle=1,$
    ytitle='R-Sec/Dn',xtitle='Zenith Distance [Degrees]',title='Millstone Responsivity @ 6300 A'; CalibData(*,5) is the R-Sec/Dn
  xyouts,10,80, thedate
  device,/close
  
  ;********* Calculate Distortion Function ***********
  set_plot,'win'
  zdarrdeg=findgen(91)    ; Creates an array from 0 degrees to 90 degrees
  zdarrpix=findgen(300)   ; Creates an array from 0 pixels to 300 pixels
  coef1=poly_fit(CalibData(errindex,3),CalibData(errindex,4),3)   ; this is the deg2pixCoef
  coef2=poly_fit(CalibData(errindex,4),CalibData(errindex,3),3)   ; this is the pix2degCoef
  newpixarr=coef1(0)+coef1(1)*zdarrdeg+coef1(2)*zdarrdeg^2.0+coef1(3)*zdarrdeg^3.0
  save, filename='C:\Program Files\Exelis\IDL82\lib\StarCalib\Millstone\coef1_'+thedate+'.sav',coef1 ; Save a file of distortion function's coefficients
  save, filename='C:\Program Files\Exelis\IDL82\lib\StarCalib\Millstone\coef2_'+thedate+'.sav',coef2 ; Save a file of distortion function's coefficients
  ; ******** Create distortion function graph in a PostScript file **********
  set_plot,'ps'
  device,/landscape,filename=thisfolder+'Millstone\Millstone_Distortion_function'+thedate+'.ps'
  plot,CalibData(errindex,3),CalibData(errindex,4),psym=1,ystyle=1,xstyle=1,xra=[0,100],yra=[0,300],$
    xtitle='Zenith Distance [Degrees]',ytitle='Zenith Distance [Pixels]',title='Millstone Distortion Function '+wave
  xyouts,10,285,thedate             ; Displays the date
  xyouts,10,265,'deg2pix coef'      ; This and the four lines below display degree to pixel coefficients
    xyouts,10,250,string(coef1(0))
    xyouts,10,235,string(coef1(1))
    xyouts,10,220,string(coef1(2))
    xyouts,10,205,string(coef1(3))
  xyouts,60,145,'pix2deg coef'      ; This and the four lines below display pixel to degree coefficients 
    xyouts,60,130,string(coef2(0))
    xyouts,60,115,string(coef2(1))
    xyouts,60,100,string(coef2(2))
    xyouts,60, 85,string(coef2(3))
  xyouts,60,70,'xcenter='+string(xofst)     ; Displays the x-center pixel
  xyouts,60,55,'ycenter='+string(yofst)     ; Displays the y-center pixel
  xyouts,60,40,'angle='+string(rot_axis/q)  ; Displays rotation angle
  xyouts,60,25,'therot='+string(therot)     ; Displays therot (number of 90 degree turns)
  oplot,zdarrdeg,newpixarr
  device,/close
  ;stop
  print, 'Finished!'
end