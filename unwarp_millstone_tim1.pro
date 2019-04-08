pro unwarp_millstone_tim1

  ;+
  ; NAME:
  ;       UNWARP_MILLSTONE_TIM1
  ; PURPOSE:
  ;       To unwarp all-sky images, subtract images at a background wavelength,
  ;       and map the latitudes and longitudes of each pixel
  ; CALLING SEQUENCE:
  ;       UNWARP_MILLSTONE_TIM1
  ; INPUTS:
  ;       thedate: The date the images were taken
  ;       thisFOV: Field of view to be unwarped (in degrees)
  ;       xofst & yofst: Center of the image
  ;       xNCP and yNCP: Location of the North Celestial Pole
  ;       rot_axis: Rotation angle
  ;       gridh: Height of the atmospheric layer (in kilometers)
  ;       mainwavelength: Wavelength (in Angstroms) of the main images
  ;       lens: Name of distortion function in nnnxy_po.pro
  ;       conjpoint & zenpoint: Strings that determine whether the conjugate and zenith points
  ;                             will be marked on the map
  ;       timegap: Average length (in seconds) between images for the night
  ;       Vignette and VanRijn files: .sav files of Vignette and VanRijn arrays
  ;                                   that must be resotred so they can be divided by later
  ; OUTPUTS
  ;       This program will create a map .sav file, which contains the image of the
  ;       displayed map (as is created in makemap.pro), a latlon .sav file that contains
  ;       the latitude, longitude, x position, and y position of each pixel. For each image
  ;       taken during the night, 4 different files will be saved of the processed image,
  ;       one .gif file, one .map file, one .sav file, and one .png file.
  ; REVISION HISTORY:
  ;       12 Jun 2017: Edited by Timothy MacDonald. Original, unedited code taken
  ;                    from unwarp_millstone.pro
  ;-
  
  ;===================================================
  ; EDIT THE FOLLOWING VARIABLES FOR EACH DISTINCT RUN
  ;===================================================
  
  thedate = 'Sep2717'     ; The date on which the set of data was taken
  thisFOV = 80.0          ; Size of the field of view to be analyzed (in degrees)
  xofst = 258             ; Pixel center of the image (x-coordinate)
  yofst = 260             ; Pixel center of the image (y-coordinate)
  xNCP = 260              ; Pixel location of the North Celestial pole (x-coordinate)
  yNCP = 386              ; Pixel location of the North Celestial Pole (y-coordinate)
  rot_axis = -6.75        ; rotation angle
  gridh = 275.0           ; height of layer (in kilometers) (400 SAR arc)
  mainwavelength = 6300   ; Wavelength of main images to be processed (in Angstroms)
  backwavelength = 6050   ; Wavelength of the background images to be processed (in Angstroms)
  lens = 'millstone_a'    ; distortion function (found in nnnxy_po.pro)
  conjpoint = 'no'       ; If set to 'yes', a square will be plotted on the map marking the field line of the conjugate location
  zenpoint = 'no'        ; If set to 'yes', an x will be plotted on the map at the location of the camera
  timegap = 402           ; Average time (in seconds) between main wavelength images on the night in question  
  restore,"C:\Program Files\Exelis\IDL82\lib\unwarp_stuff\Millstone\Millstonevigarr.sav"        ; Restores vignette file
  restore,"C:\Program Files\Exelis\IDL82\lib\unwarp_stuff\Millstone\MillstonevanRijnarr400.sav" ; Restores van Rijn function
  
  ;===============  GPS and DMSP Data (optional) =================
  ;==== If not desired, comment out the following three lines ====
  ;sat = 'GPS'                     ; Which satellite you are looking at
  ;satfile = 'lpoc354-2015-12-20'  ; Name of .csv file that will be read
  ;satnum = [1]                   ; GPS satellite number(s) to be looked at (GPS Only)
  ;===============================================================
  
  ;=======================================================
  ; Millstone Hill specific variables (not usually edited)
  ;=======================================================
  site = 'Millstone'          ; Name of site
  ix = 512                    ; x-dimension of the image (in pixels)
  iy = 512                    ; y-dimension of the image (in pixels)
  scale = 5.0                 ; Scale to multiply the image by during unwarping
  clat = 42.62                ; center latitude of Millstone
  clon = -71.5                ; center longitude of Millstone
  mlat = 45                   ; Latitude of center of map (usually same as clat, change if need to fit images from multiple sites on same map)
  mlon = -90                  ; Longitude of center of map (usually same as clon, change if need to fit images from multiple sites on same map)
  conlat = -67.5666           ; center latitude of Rothera (Millstone's conjugate observatory)
  conlon = -68.1325           ; center longitude of Rothera (Millstone's conjugate observatory)
  camh = 0.146                ; camera height (in kilometers)
  altc = 90.0                 ; altitude of center - assume pointing straight up
  azc = 0.0                   ; azimuth of center
  glats = 5.0                 ; Latitude interval displayed on map
  glons = 5.0                 ; Longitude interval displayed on map
  xs = 800                    ; x-dimension of the window display (in pixels)
  ys = 800                    ; y-dimension of the window display (in pixels)
  res = 5.0                   ; Responsivity determined from Responsivity graph (in Rsec/dn)
  backfac = 0.5               ; Background scaling factor
  therot = 2                  ; Amount to be rotated by, counterclockwise (0 = 0 deg, 1 = 90 deg, 2 = 180 deg, etc.)
  exptime = 120.              ; Exposure time of main and background images (in seconds)
  type = 'ortho'              ; Type of map projection to be used
  catfile = thedate+'.cat'    ; Name of catalog file
  ;==============================================================
  
  ;*** Read into raw image directory and create output directory ***
  pathToOut = 'C:\Unwarp_Output\Millstone\'+thedate+'\' ; path to output
  check_n_makedir, pathToOut                   ; creates a directory for the output if it does not yet exist
  pathToRaw = 'C:\Millstone\'+thedate+'\'               ; path to raw images
  ;*** Shift vignette and vanRijn arrays to the center of the image
  vig=shift(newvig,xofst-256,yofst-256)                      ; Shifts vignette to empiracally determined center
  vanRijnarr=shift(vanRijnarr400,xofst-256,yofst-256)        ; Shifts van Rijn function to empiracally determined center
  
  theday = strmid(thedate,3,2)                   ; extracts the day number of the date
  themonth = strmid(thedate,0,3)                 ; extracts the month of the date
  theyear = 2000 + fix(strmid(thedate,5,2))      ; extracts the year
  
  if isa(sat) eq 1 then begin
  ; Path to satellite data
  if sat eq 'GPS' then pathToSat = 'C:\Satellite_Data\GPS\Millstone\'+strcompress(theyear,/remove_all)+'\'
  if sat eq 'DMSP' then pathToSat = 'C:\Satellite_Data\DMSP\Millstone\'+strcompress(theyear,/remove_all)+'\'  
  ; Path to satellite outputs            
  if sat eq 'GPS' then pathToOutSat = 'C:\Satellite_Data\GPS\Millstone\'+strcompress(theyear,/remove_all)+'\Outputs\'+thedate+'\' 
  if sat eq 'DMSP' then pathToOutSat = 'C:\Satellite_Data\DMSP\Millstone\'+strcompress(theyear,/remove_all)+'\Outputs\'+thedate+'\' 
  check_n_makedir,pathToOutSat
  endif
  
  ; ******* Read an image list and create arrays of the image names (both for main and background images) *******
  
  ; This procedure will search for a .txt list file and will ONLY make the arrays if the list file is found
  makeimgarray_list,pathToRaw,mainwavelength,imglist,backlist,count1,count2
  
  ; This procedure will ONLY make the arrays if NO list file is found. It will read the image names from the catalog file.
  makeimgarray_cat,pathToRaw,catfile,mainwavelength,backwavelength,imglist,backlist,count1,count2
  
  ;stop
  ;*********************************
 
  ; Remove images that have no correspoding image
  lonefileremove,imglist,backlist,count1,count2,stopcount
  
  thesize = size(imglist)
  print,thesize
  ;stop
  
  ;*********** Create paths for directories of the outputs *************
  ; Embedded in the names are the unwarp height, field of view, and date
  
  unwarpheight = strmid(strcompress(gridh,/remove_all),0,3)
  
  ; Name of output folder
  pathToOut = pathToOut + 'h'+unwarpheight + '_' + strcompress(thisFOV,/remove_all)+ '_' + thedate + '\'
  
  check_n_makedir,pathToOut
  
  ;*** create map and read map information ***
  
   ; Subroutine used to make the map 
  makemap,xs,ys,mlat,mlon,glats,glons,type,latmin,latmax,lonmin,lonmax 
   ; Plots location of conjugate magnetic field point, and satellite paths
  get_map,xs,ys,lat,lon,xx,yy,clat,clon,conlat,conlon,gridh,theday,themonth,theyear,conjpoint,zenpoint
  
  ; save a copy of the created map
  map=tvrd(0,0,xs,ys)
  save,filename=PathToOut+'map.sav',map
  
  ;This is lat/lon/pix info so that I can click on an image and
  ;get lat/lon coordinates to put into L-value code
  save,filename=pathToOut+'latlon.sav',xs,ys,lon,lat,xx,yy
  
  ;*** get unwarp data ***
  ll_to_xy,lat,lon,x,y,scale,xx,yy,fix(thisFOV),xofst,yofst,clat,clon,camh,gridh,altc,azc,lens,rot_axis
  ;stop
  
  ; ******* Create a dark image, to be subtracted later ********
  aMedianDark = imgMedian(pathToRaw,theDate,'DARK',exptime)
  amediandark = rotate(amediandark,therot)
  
  sumcounter = 0
  
  ; Calculate the latitude of each pixel in a vertical cross section down the middle
  vertlat = fltarr(800)
  for h=0,799 do begin
    vertlatitude = convert_coord(400,h,/device,/to_data)
    vertlat[h] = vertlatitude[1]
  endfor

  ; Calculate the longitude of each pixel in a horizontal cross section across the middle
  horlon = fltarr(800)
  for h=0,799 do begin
    horlongitude = convert_coord(h,400,/device,/to_data)
    horlon[h] = horlongitude[0]
  endfor

  ;*******************************
  ; FOR LOOP TO UNWARP THE IMAGES
  ;*******************************
  for ii = 0, stopcount do begin
  ; for ii = 0, 0 do begin
    theimage = fltarr(xs,ys)
    file6300 = pathToRaw + imglist(ii)
    fileback = pathToRaw + backlist(ii)
    ; ******* Find the times of each image to subtract the proper angle *********
    thetime = fname2time(imglist(ii))
    backtime = fname2time(backlist(ii))
    time = float(strmid(thetime,0,2))*3600.+float(strmid(thetime,3,2))*60.+float(strmid(thetime,6,2))
    btime = float(strmid(backtime,0,2))*3600.+float(strmid(backtime,3,2))*60.+float(strmid(backtime,6,2))
    timediff = time - btime
    anglediff = 15./3600.*timediff
    
    ;******** Find necessary width of velogram ********
    if ii eq 0 then begin
      ; ** Determine if sunrise/sunset happen during middle of image set **
      daysplit = float(strmid(imglist,1,2))
      dayindex = where(imglist ne '',daycount)
      hoursplit = fltarr(daycount-1)
      
      for dd = 1, daycount-1 do hoursplit[dd-1] = daysplit[dd] - daysplit[dd-1]
      maxhoursplit = where(hoursplit eq max(hoursplit))
      ; ** If sunrise/sunset occurs during image set, find actual begin/end times
      if max(hoursplit) ge 2 then thebegintime = fname2time(imglist(maxhoursplit[0]+1))
      if max(hoursplit) ge 2 then theendtime = fname2time(imglist(maxhoursplit[0]))
      if max(hoursplit) le 1 then thebegintime = thetime
      if max(hoursplit) le 1 then  theendtime = fname2time(imglist(stopcount))
      
      beginhour = float(strmid(thebegintime,0,2))
      endhour = float(strmid(theendtime,0,2))
      
      ; If endhour is cutting it close, give it an extra hour to be safe
      endseconds = float(strmid(theendtime,3,2))*60 + float(strmid(theendtime,6,2))
      if endseconds ge 3300 then endhour = endhour+1
      
      if endhour ne 23 then endhour = endhour+1
      if endhour eq 23 then endhour = 0
      
      if endhour gt beginhour then nighttimehours = endhour-beginhour
      if endhour lt beginhour then nighttimehours = 24+endhour-beginhour
      ; Create empty array that will become velogram
      columns = round(nighttimehours*3600/timegap)
      ; Create empty arrays for image cross sections
      vertline10 = fltarr(columns,10,800)   ; Ten pixel wide vertical cross section
      horline10 = fltarr(columns,800,10)    ; Ten pixel wide horizontal cross section
      vertline = fltarr(columns,800)        ; One pixel wide vertical cross section (average of the five pixel width)
      horline = fltarr(columns,800)         ; One pixel wide horizontal cross section (average of the five pixel width)
      reltime = fltarr(stopcount+1)     ; Relative time (number of column in which the final image each cross section needs to go)
    endif
    ;**************************************************
    
    load16,file6300,ix,iy,128,therot,img     ; Load the main image
    img=float(img)                           ; Turn image into floating point array
    img = img - amediandark                  ; Subtract the dark image
    img = sigma_filter(img,radius=2.5,n_sigma=3.0,/all,/iterate) ; Apply Sigma Filter
    
    load16,fileback,ix,iy,128,therot,back    ; Load the background image
    back=float(back)                         ; Turn image into floating point array
    back = back - amediandark                ; Subtract the dark image
    back = sigma_filter(back,radius=2.5,n_sigma=3.0,/all,/iterate)  ; Apply Sigma Filter
    back=rot(back,-anglediff,1.0,xNCP,yNCP,/interp,/pivot)          ; Rotate about NCP to align stars
    
    ; ****** Recreate a blank map to continue lat-lon to pixel coordinate conversions ******
    ; ******************** (Only used for satellite observations) **************************
    if ii ne 0 and isa(satfile) eq 1 and isa(sat) eq 1 then begin
      remakemap,clat*180./!pi,clon*180./!pi,latmin,lonmin,latmax,lonmax,xs,ys,type
    endif
    
    ; img(xofst-3:xofst+3,yofst-3:yofst+3)=30000 ; put a white square at zenith
    
    img=img-backfac*back       ; Subtracts the background image
    img=float(img)/vig         ; Divides by the vignette
    img=img/vanRijnarr         ; Divides by the VanRijn
    wset,0
    
    ;tv,bytscl(img,min=0,max=400)
    ;stop
    
    img=res*img/exptime          ; Divide by exposure time and multiply by responsivity
    
    tv,bytscl(img,min=0,max=800)
    imgcal=fix(img)
    ;stop
    
    img=rebin(img,scale*ix,scale*iy)         ; make it large for unwarping
    
    ; ****** Detect if any pixels are saturated ******
    ; sindex = where(img ge 64000,scount)
    ; if scount NE 0 THEN begin
    ;    print,'SATURATION!'
    ;    img(sindex) = -5000
    ; endif
    
    ;testcircle,img,img      ;this is for putting a circle on an image at a particular
    ;zenith distance to see if the code is working right
    
    theimage(xx,yy)=img(x,y) ; remap image
    
    wset,1
    tv,bytscl(theimage,min=0,max=800); unwarped image without overlays
    ;stop
    
    ; ******* Create cross sections of each image ******* 
  
    ; Adjust times so that the earliest time is the beginning hour
    if time ge beginhour*3600 then adjtime=time-beginhour*3600
    if time lt beginhour*3600 then adjtime=time+86400-beginhour*3600
    ; Determine which relative time column to put each image strip in
    reltime(ii) = adjtime/timegap
    reltime(ii) = round(reltime(ii))
    
    ; Reads a vertical cross section five pixels wide down the middle of the image
    vertline10(reltime(ii),*,*) = theimage(395:404,0:799)
    ; Take the average of each five pixel wide element
    for k=0,799 do begin
      if reltime(ii) ne reltime(ii-1) or ii eq 0 then begin           ; Puts slice in an empty column
        vertline(reltime(ii),k) = mean(vertline10(reltime(ii),*,k))
      endif else begin       ; If desired column is already used, take average of both cross sections
        vertline(reltime(ii),k) = (mean(vertline10(reltime(ii),*,k)) + vertline(reltime(ii),k))/2
      endelse
    endfor
    
    ; Reads a horizontal cross section five pixels wide across the middle of the image
    horline10(reltime(ii),*,*) = theimage(0:799,395:404)
    ; Take the average of each five pixel wide element
    for k=0,799 do begin
      if reltime(ii) ne reltime(ii-1) or ii eq 0 then begin           ; Puts slice in an empty column
        horline(reltime(ii),k) = mean(horline10(reltime(ii),k,*))
      endif else begin       ; If desired column is already used, take average of both cross sections
        horline(reltime(ii),k) = (mean(horline10(reltime(ii),k,*)) + horline(reltime(ii),k))/2
      endelse
    endfor
    
    wset,1
    print,ii
    ;******** Create an image with a map overlaid upon it **********
    wh = where(map gt 0)
    bl = where(map le 1,l)
    
    map1 = map
    map1(wh) = 2
    maplinesindex = where(map eq 255)
    mapimage = theimage
    mapimage(maplinesindex) = 10000
    mapimageint=fix(mapimage)
    loadct,0
    
    ;******* Create names for files ********
    giffile=strmid(imglist(ii),0,7)+".gif"
    savfile=strmid(imglist(ii),0,7)+".sav"
    mapfile=strmid(imglist(ii),0,7)+".map"
    pngfile=strmid(imglist(ii),0,7)+".png"
    
    ;***** Create labels at the bottom of the image **********
    xyouts,350,25,strcompress(thisFOV,/remove_all)+' deg FOV',/device
    xyouts,500,25,strcompress(mainwavelength/10.,/remove_all)+' nm',/device
    xyouts,10,10,thedate,/device, charsize=1.0
    xyouts,100,10,thetime+' UT',/device
    xyouts,10,25,'Unwarp Height = '+unwarpheight+'km',/device
    xyouts,200,25,'Site = Millstone',/device
    xyouts,200,10,'Center = ('+strcompress(long(xofst/5),/remove_all)+','+strcompress(long(yofst/5),/remove_all)+')',/device
    xyouts,350,10,'Rotation Angle = '+strcompress(rot_axis*180./!pi,/remove_all)+' degrees',/device
    ; stop
    savim = tvrd()
    if sumcounter eq 0 then sum = theimage
    if sumcounter gt 0 then sum = sum + theimage
    
    ;********* Create and save images ***************
    write_gif,pathToOut + giffile,savim       ;unwarped with map and text
    write_png,pathToOut + pngfile,mapimageint ;calibrated unw format
    ;png is 16 bit in Rayleighs
    ;save,filename=pathToOut+savfile,savim
    save,filename=pathToOut+savfile,theimage  ; just the unwarped image
    save,filename=pathToOut+mapfile,mapimage  ; unwarped image with map
    sumcounter++
    
    ; **************** Calculate GPS or DMSP data *******************
    if isa(sat) eq 1 and isa(satfile) eq 1 then begin 
      masterSat_tim1,sat,satfile,satnum,pathToSat,pathToOutSat,thetime,thedate,imglist[ii]
    endif
    ; ***************************************************************

  endfor
  close,1
  
  ;****************************************************************************
  ; Display the images that contain the vertical and horizontal cross sections
  ;****************************************************************************

  tickname=['22','23','00','01','02','03','04','05','06','07','08','09','10','11','12']   ; Tickmarks that will be put on the graphs
  ; Remove unused hours from tickname array
  beginstring = strtrim(string(fix(beginhour)),1)
  if strlen(beginstring) eq 1 then beginstring = '0'+beginstring
  endstring = strtrim(string(fix(endhour)),1)
  if strlen(endstring) eq 1 then endstring = '0'+endstring
  beginindex = where(tickname eq beginstring)
  endindex = where(tickname eq endstring)
  tickname = tickname[beginindex:endindex]
  maxvalue = 500      ; Maximum value displayed using the bytscl function
  plot_velogram,reltime,vertline,horline,tickname,thedate,vertlat,horlon,pathToOut,columns,maxvalue,site
  ;stop
  ;*************************************************************************
  ; Create a 00README .txt file with information about the input parameters
  ;*************************************************************************
  openw,2,pathToOut+'00README.txt'
  printf,2, 'Date: '+thedate
  printf,2, 'Field of View:'+strcompress(thisFOV)+' degrees'
  printf,2, 'x position of center:'+strcompress(xofst/5)
  printf,2, 'y position of center:'+strcompress(yofst/5)
  printf,2, 'x position of North Celestial Pole:'+strcompress(xNCP)
  printf,2, 'y position of North Celestial Pole:'+strcompress(yNCP)
  printf,2, 'Rotation Angle:'+strcompress(rot_axis*180./!pi)+' degrees'
  printf,2, 'Height of atmospheric layer:'+strcompress(gridh)+' km'
  printf,2, 'Wavelength of images:'+strcompress(mainwavelength)+' Å'
  printf,2, 'Wavelength of background images:'+strcompress(backwavelength)+' Å'
  printf,2, 'Name of distortion function: '+lens
  printf,2, 'x-dimension of image:'+strcompress(ix)
  printf,2, 'y-dimension of image:'+strcompress(iy)
  printf,2, 'Scale by which the image is multiplied:'+strcompress(scale)
  printf,2, 'Latitude of the camera:'+strcompress(clat*180./!pi)
  printf,2, 'Longitude of the camera:'+strcompress(clon*180./!pi)
  printf,2, 'Elevation of the camera:'+strcompress(camh)+' km'
  printf,2, 'Altitude of center pixel:'+strcompress(altc*180./!pi)+' degrees'
  printf,2, 'Azimuth of center pixel:'+strcompress(azc*180./!pi)+' degrees'
  printf,2, 'Latitude interval on map:'+strcompress(glats)+' degrees'
  printf,2, 'Longitude interval on map:'+strcompress(glons)+' degrees'
  printf,2, 'x-dimension of unwarped image:'+strcompress(xs)
  printf,2, 'y-dimension of unwarped image:'+strcompress(ys)
  printf,2, 'Responsivity:'+strcompress(res)
  printf,2, 'Latitude Limits: ['+strcompress(latmin,/remove_all)+','+strcompress(latmax,/remove_all)+']'
  printf,2, 'Longitude Limits: ['+strcompress(lonmin,/remove_all)+','+strcompress(lonmax,/remove_all)+']'
  printf,2, 'Background scaling factor:'+strcompress(backfac)
  printf,2, 'Rotation number:'+strcompress(therot)
  printf,2, 'Exposure time of images:'+strcompress(exptime)+' seconds'
  printf,2, 'Type of map projection: '+type
  printf,2, 'Average time between images:'+strcompress(timegap)+' seconds'
  free_lun,2
  ;stop
  print,"Finished!"
end