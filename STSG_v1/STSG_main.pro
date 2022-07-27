;******************************************************************
;      STSG is used to reconstruct high-quality NDVI time series data (MODIS/SPOT) 
;       
;      - This procedure STSG_v1 is the source code for the first version of STSG.
;      This is a parallel computing code using multiple cpu cores.
;      For STSG, computers with more cpu cores are preferred. 
;              
;      Coded by Ruyin Cao
;      Email: cao.ruyin@uestc.edu.cn or caoruyin119@gmail.com
;      School of Resources and Environment, University of Electronic 
;      Science and Technology of China
;                                          
;      Update history: 2018-06-30: first version released
;                       
;      Reference: Cao,R., Chen, Y., Shen, M.G., Chen, J., Zhou, J., 
;           Wang, C., Yang, W.  A simple method to improve the quality 
;           of NDVI time-series data by integrating spatiotemporal information
;           with the Savitzky-Golay filter, Remote Sensing of Environment, under review
;******************************************************************


pro STSG_main

; Input parameters  
;********************************************************************************
year = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,$
        2012,2013,2014,2015,2016]
        
; the thereshold of correlation coefficient to define similar pixels 
sampcorr = 0.9 

; half of the neighboring window size within which to search similar pixels 
win = 10 

; the path of the NDVI data
NDVI_filepath = 'D:\STSG_v1\test data 2\NDVI images\NDVI_test_'

; the path of the NDVI quality flags (realibility)
reliability_filepath = 'D:\STSG_v1\test data 2\reliability\test_reliability_'

; STSG performs by lines to prevent out of memory. This paramenter is set $
; to determine how many lines are processed by each cup core at a time. 
; For this case, if the computer have 16 cpu cores, 48 lines (i.e. 16*3) are processed in parallel 
cpucore_line = 3

; snow_address indicates whether to deal with snow contamianted NDVI values (1=yes/0=no)
snow_address = 1
;********************************************************************************

dir = file_dirname(routine_filepath('STSG_main')) 
fid_NDVI =  intarr(n_elements(year))
fid_QA =  intarr(n_elements(year))
for yeari = 0, n_elements(year)-1 do begin
  filename = strtrim(NDVI_filepath,2)+strtrim(year[yeari],2)
  envi_open_file, filename, r_fid=fid, /no_realize
  fid_NDVI[yeari] = fid
  if yeari eq 0 then envi_file_query, fid, ns = ns_origin, nl = nl_origin, nb = nb, $
                                      data_type = data_type_NDVI
  
  filename = strtrim(reliability_filepath,2)+strtrim(year[yeari],2)
  envi_open_file, filename, r_fid=fid, /no_realize
  fid_QA[yeari] = fid
  if yeari eq 0 then envi_file_query, fid, data_type = data_type_QA                                    
endfor
map_info = envi_get_map_info(fid=fid_NDVI[0]) 

imgNDVI = make_array(ns_origin+2*win,nl_origin+2*win,nb,type = data_type_NDVI)
imgQA = make_array(ns_origin+2*win,nl_origin+2*win,nb,type = data_type_QA)
dims = [-1,0,ns_origin-1,0,nl_origin-1]
fid_NDVI_buffer = intarr(n_elements(year))
fid_QA_buffer = intarr(n_elements(year))
for yeari = 0, n_elements(year)-1 do begin    
  for k = 0, nb-1 do begin
    imgNDVI[win:(ns_origin+win-1),win:(nl_origin+win-1),k] = envi_get_data(fid = fid_NDVI[yeari],dims = dims,pos=k)
    imgQA[win:(ns_origin+win-1),win:(nl_origin+win-1),k] = envi_get_data(fid = fid_QA[yeari],dims = dims,pos=k)
  endfor
  outfile = strtrim(dir,2)+strtrim('\NDVI_buffer_',1)+strtrim(yeari,1)
  envi_write_envi_file, imgNDVI, out_name=outfile, data_type=data_type_NDVI,$
                        ns=ns_origin+2*win, nl=nl_origin+2*win, nb=nb, /no_open
  envi_open_file, outfile, r_fid= fid
  fid_NDVI_buffer[yeari] = fid
  outfile = strtrim(dir,2)+strtrim('\QA_buffer_',1)+strtrim(yeari,1)
  envi_write_envi_file, imgQA, out_name=outfile, data_type=data_type_QA,$
                        ns=ns_origin+2*win, nl=nl_origin+2*win, nb=nb, /no_open
  envi_open_file, outfile, r_fid= fid
  fid_QA_buffer[yeari] = fid                      
  if yeari eq 0 then begin
    envi_file_query, fid_QA_buffer[yeari], ns=ns, nl=nl, nb = nb
  endif                                       
  envi_file_mng, id=fid_NDVI[yeari], /remove
  envi_file_mng, id=fid_QA[yeari], /remove                     
endfor

imgNDVI = make_array(ns,1,nb,n_elements(year), type = data_type_NDVI)
imgQA = make_array(ns,1,nb,n_elements(year), type = data_type_QA)
NDVI_reference = make_array(ns,nl,nb, type = data_type_NDVI)
print, ''
print, 'Start: Generating NDVI reference time-series data'
t0=systime(1)
for linei = 0, nl-1 do begin
  if linei mod fix(float(nl)/10.) eq 0 then $
  print, 'processing : ', strtrim((linei/fix(float(nl)/10.)+1)*10,2)+strtrim('%',2)
  dims = [-1,0,ns-1,linei,linei]
  for yeari = 0, n_elements(year)-1 do begin    
    for k = 0, nb-1 do begin
    imgNDVI[*,0,k,yeari] = envi_get_data(fid = fid_NDVI_buffer[yeari],dims = dims,pos=k)
    imgQA[*,0,k,yeari] = envi_get_data(fid = fid_QA_buffer[yeari],dims = dims,pos=k)
    endfor
  endfor
  
  for j = 0, ns-1 do begin
    for k = 0, nb-1 do begin
    vector_NDVI = imgNDVI[j,0,k,*]
    vector_QA = imgQA[j,0,k,*]
    res = where(vector_QA eq 1 or vector_QA eq 0, count)
    if count ge 3 then begin
    NDVI_reference[j,linei,k] = mean(vector_NDVI[res])
    endif else begin
    NDVI_reference[j,linei,k] = -1.0
    endelse
    endfor 
    vec_res1 = NDVI_reference[j,linei,*]
    res = where(vec_res1 ne -1.0, count)
    if count lt nb and count ne 0 then begin
    V = vec_res1[res]
    Xout = indgen(nb)
    NDVI_reference[j,linei,*] = interpol(V, res, Xout)
    endif  
  endfor  
endfor

outfile = strtrim(dir,2)+strtrim('\NDVI_reference',2)   
envi_write_envi_file, NDVI_reference, out_name=outfile, data_type=data_type_NDVI,$
                      ns=ns, nl=nl, nb=nb, MAP_INFO=map_info, /no_open
envi_open_file, outfile, r_fid=fid_reference, /no_realize
print, '' 
print, 'Start: STSG'
ncpus = !cpu.hw_ncpu     
stepline = cpucore_line
step = stepline*ncpus
lastline = (nl-2*win) mod step
if lastline ne 0 then begin
mfid = intarr((nl-2*win)/step+1);mosic parameters
mdims = intarr(5,(nl-2*win)/step+1)
mpos = intarr(nb*n_elements(year),(nl-2*win)/step+1)
pos = indgen(nb*n_elements(year))
x0 = intarr((nl-2*win)/step+1)
y0 = intarr((nl-2*win)/step+1)
     
endif else begin
mfid = intarr((nl-2*win)/step);mosic parameters
mdims = intarr(5,(nl-2*win)/step)
mpos = intarr(nb*n_elements(year),(nl-2*win)/step)
pos = indgen(nb*n_elements(year))
x0 = intarr((nl-2*win)/step)
y0 = intarr((nl-2*win)/step)
     
endelse  
  
numfid = 0
numfid_index = 0

for j = win, nl-1-win, step do begin
print,'STSG Finished :', strtrim(fix(float(j-win)/float(nl-1-2*win)*100))+strtrim('%',2)
leftline = (nl-win)-j  
  if leftline ge step then begin ; determine the last multiple lines
  oBridge = objarr(ncpus)
  nls = step
  endif else begin
  ncpus = 1
  oBridge = objarr(ncpus)        
  nls = leftline
  endelse
    
  for i=0, ncpus-1 do begin
    ; read data by multiple lines
    case ncpus of
    1: begin
       dims = [-1,0,ns-1,j-win,nl-1] 
       img_NDVI = fltarr(ns, leftline+2*win, nb, n_elements(year))
       img_QA = fltarr(ns, leftline+2*win, nb, n_elements(year))
       reference_data = fltarr(ns, leftline+2*win, nb)      
       SG_NDVI = fltarr(ns, nls, nb, n_elements(year)) 
       vector_out = fltarr(ns,leftline,nb,n_elements(year))
    end 
    else: begin
       dims = [-1,0,ns-1,j+i*stepline-win,j+(i+1)*stepline-1+win] 
       img_NDVI = fltarr(ns, stepline+2*win, nb, n_elements(year))
       img_QA = fltarr(ns, stepline+2*win, nb, n_elements(year))
       reference_data = fltarr(ns, stepline+2*win, nb)      
       SG_NDVI = fltarr(ns, nls, nb, n_elements(year))
       vector_out = fltarr(ns,stepline,nb,n_elements(year))
    end
    endcase
      
    for yearfill = 0, n_elements(year)-1 do begin
      for k = 0, nb-1 do begin 
        img_NDVI[*,*,k,yearfill] = float(envi_get_data(fid = fid_NDVI_buffer[yearfill],dims = dims, pos = k))/10000.
        img_QA[*,*,k,yearfill] = float(envi_get_data(fid = fid_QA_buffer[yearfill],dims = dims, pos = k)) 
        if yearfill eq 0 then reference_data[*,*,k] = envi_get_data(Fid=fid_reference,dims = dims,pos = k)/10000.      
      endfor
    endfor           
               
    oBridge[i] = obj_new("IDL_IDLBridge", callback = 'demo_bridge_call')
    oBridge[i] ->setvar,"win", win
    oBridge[i] ->setvar,"sampcorr", sampcorr
    oBridge[i] ->setvar,"reference_data", reference_data
    oBridge[i] ->setvar,"img_QA", img_QA
    oBridge[i] ->setvar,"img_NDVI", img_NDVI            
    oBridge[i] ->setvar,"vector_out", vector_out
    oBridge[i] ->setvar,"dims", dims  
    oBridge[i] ->setvar,"snow_address", snow_address

    file_main  = '.compile '+strtrim(dir,1)+strtrim('\STSG_main.pro',2)
    file_STSG = '.compile '+strtrim(dir,2)+strtrim('\STSG_filter.pro',2)
    oBridge[i] ->execute, file_STSG
    oBridge[i] ->execute, file_main    
    if oBridge[i]->status() eq 0 then oBridge[i]->execute, "vector_out = STSG_filter(win,sampcorr,img_NDVI,img_QA,reference_data,dims,snow_address,vector_out)",/nowait 
          
  endfor
         
  notdone = 1
  while notdone eq 1 do begin
    done = 0
    for i=0, n_elements(oBridge)-1 do done = done+oBridge[i]->status()
    if done eq 0 then notdone = done
  endwhile  
  for i=0, n_elements(oBridge)-1 do begin
    vector_out = oBridge[i]->getvar('vector_out')
    case ncpus of 
    1: begin
    SG_NDVI[*,0:(nls-1),*,*] = vector_out
    end
    else: begin
    SG_NDVI[*,(i*stepline):((i+1)*stepline-1),*,*] = vector_out
    end
    endcase   
    obj_destroy, oBridge[i] 
  endfor 
   
  if ncpus ne 1 then outfile = strtrim(dir,2)+strtrim('\tempoutimg_sub',1)+strtrim(j,1)
  if ncpus eq 1 then outfile = strtrim(dir,2)+strtrim('\tempoutimg_sub',1)+strtrim(j+1,1)
  SG_NDVI = reform(SG_NDVI, ns, nls, nb*n_elements(year))

  envi_write_envi_file, SG_NDVI, out_name=outfile, data_type=4,$
                        ns=ns, nl=nls, nb=nb*n_elements(year), /no_open                      
  envi_open_file, outfile, r_fid= sub_fid
  envi_file_query, sub_fid, ns=sub_ns, nl=sub_nl                     
  mfid[numfid] = sub_fid 
  mdims[*,numfid] = [-1,win, sub_ns-1-win,0, sub_nl-1]   
  mpos[*,numfid] = indgen(nb*n_elements(year))
  x0[numfid]= long(0)
  y0[numfid]= long(j)-win
  pixel_size = [1.,1.]
  use_see_through = replicate(0.,ns) 
  numfid = numfid+1 

endfor
    
outresult = strtrim(NDVI_filepath,2)+strtrim('STSG',1)
envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=outresult, georef=0, xsize=ns_origin, $
    ysize=nl_origin, x0=x0, y0=y0,map_info=map_info, $
    out_dt=4, pixel_size=pixel_size, $
    background=0, use_see_through=use_see_through 
    
for i= 0, n_elements(mfid)-1,1 do begin
  envi_file_mng, id=mfid[i], /remove, /delete
endfor
for i = 0, n_elements(year)-1 do begin
  envi_file_mng, id=fid_NDVI_buffer[i], /remove, /delete
  envi_file_mng, id=fid_QA_buffer[i], /remove, /delete
endfor      

envi_file_mng, id=fid_reference, /remove, /delete
print, 'STSG Finished !'
print, 'Time used =', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',floor((systime(1)-t0) mod 60),'s'  
       
end

;==============================================
pro demo_bridge_call, status, error, oBr  
  case status of  
    2: str="Completed"  
    3: str="Warming" 
    4: str=error ; Aborted message  
  endcase  
  print,'CPU core::',str  
end
