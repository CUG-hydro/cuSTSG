function STSG_filter, win, sampcorr, img_NDVI,img_QA, reference_data, dims,snow_address, vector_out            
     
subs = n_elements(img_NDVI[*,0,0,0])              
subl = n_elements(img_NDVI[0,*,0,0])
nb =  n_elements(img_NDVI[0,0,*,0])
yearnum = n_elements(img_NDVI[0,0,0,*])
      
for ii = win, subs-1-win do begin
  for jj = win, subl-1-win do begin
    catch, errorstatus
    if (errorstatus ne 0) then begin
    catch, /cancel
    continue
    endif 
     
    for year = 0, yearnum-1 do begin
      vector_in = reform(img_NDVI[ii,jj,*,year],nb)
      vector_QA = reform(img_QA[ii,jj,*,year],nb)
      res = reverse(sort(vector_in))
      if mean(vector_in[res[0:2]]) gt 0.15 then begin   
          
      num_elements = n_elements(vector_in) 
      indic = 0         
      ;searching similar pixels
      if year eq 0 then begin
      corr_res = fltarr(2*win+1, 2*win+1)
      Slope_res = fltarr(2*win+1, 2*win+1)
      Intercept_res = fltarr(2*win+1, 2*win+1)
      vectorY = reference_data[ii,jj,*]
      new_corr_similar_res = fltarr(2*win+1, 2*win+1)
            
      for si = -1*win, win, 1 do begin
        for sj = -1*win, win, 1 do begin
          vectorX = reference_data[ii+si,jj+sj,*]
          if vectorX[3] ne 0. then begin
          corr_res[si+win,sj+win] = correlate(vectorX, vectorY)
          tempQA1 = img_QA[ii,jj,*,*]
          tempQA2 = img_QA[ii+si,jj+sj,*,*]
          tempNDVI1 = img_NDVI[ii,jj,*,*]
          tempNDVI2 = img_NDVI[ii+si,jj+sj,*,*]
          for kb = 0, nb-1 do begin
            tempNDVI1[0,0,kb,*] = tempNDVI1[0,0,kb,*]/reference_data[ii,jj,kb]
            tempNDVI2[0,0,kb,*] = tempNDVI2[0,0,kb,*]/reference_data[ii+si,jj+sj,kb]
          endfor
          res = where(tempQA1 le 1 and tempQA2 le 1, count) 
          if count ge 30 then begin
             vectorYY = tempNDVI1[res]
             vectorXX = tempNDVI2[res]
             new_corr_similar_res[si+win,sj+win] = correlate(vectorXX, vectorYY)
             slopeintercept = linfit(vectorXX, vectorYY)
             Intercept_res[si+win,sj+win] = slopeintercept[0]
             Slope_res[si+win,sj+win] = slopeintercept[1]
          endif 
          endif  
        endfor
      endfor
 
       res = where(finite(corr_res) eq 0, count)
       if count ne 0 then corr_res[res] = 0.
       res = where(finite(Slope_res) eq 0, count)
       if count ne 0 then Slope_res = 0.
       res = where(corr_res ge sampcorr and Slope_res ne 0., count)
       if count ge 2 then begin
         samp = count-1
         similar_index = indgen(2,samp)
         slope_intercept = fltarr(2,samp)
         corr_similar = fltarr(samp)
              
         new_corr = reverse(sort(corr_res))
         for i = 0, samp-1 do begin
           similar_index[1,i] = fix(new_corr[i+1]/(2*win+1))+jj-win 
           similar_index[0,i] = new_corr[i+1]-fix(new_corr[i+1]/(2*win+1))*(2*win+1)+ii-win
           slope_intercept[1,i] = Slope_res[fix(new_corr[i+1]/(2*win+1)), new_corr[i+1]-fix(new_corr[i+1]/(2*win+1))*(2*win+1)]
           slope_intercept[0,i] = Intercept_res[fix(new_corr[i+1]/(2*win+1)), new_corr[i+1]-fix(new_corr[i+1]/(2*win+1))*(2*win+1)]  
           corr_similar[i] = new_corr_similar_res[fix(new_corr[i+1]/(2*win+1)), new_corr[i+1]-fix(new_corr[i+1]/(2*win+1))*(2*win+1)]   
         endfor
         aap = 1    
       endif else begin
         aap = 0
       endelse
       endif ; if year
            
       ;generate the trend curve
        if aap eq 1 then begin    
          temp_NDVI = fltarr(nb,samp)
          trend_NDVI = fltarr(nb)
          vector_in = reform(img_NDVI[ii,jj,*,year],nb)
          vector_QA = float(reform(img_QA[ii,jj,*,year],nb))   
            
          for i = 0, nb-1 do begin
            for j = 0, samp-1 do begin
              if img_QA[similar_index[0,j],similar_index[1,j],i,year] le 1 then begin
                 new_ratio = img_NDVI[similar_index[0,j],similar_index[1,j],i,year]/reference_data[similar_index[0,j],similar_index[1,j],i]
                 temp_NDVI[i,j] = (slope_intercept[0,j]+new_ratio*slope_intercept[1,j])*reference_data[ii,jj,i]
                 if temp_NDVI[i,j] ge 1 or temp_NDVI[i,j] le -0.2 then temp_NDVI[i,j] = 0.  
              endif       
            endfor
            
            res = where(finite(temp_NDVI[i,*]) eq 0, count)
            if count ne 0 then temp_NDVI[i,res] = 0.
            res = where(temp_NDVI[i,*] ne 0., count)
            if count ne 0 then begin
               new_corr_similar = corr_similar
               new_corr_similar[res] = new_corr_similar[res]/total(new_corr_similar[res])
               trend_NDVI[i] = total(new_corr_similar[res]*temp_NDVI[i,res])
            endif   
          endfor             
            
          ; generating the trend_NDVI
          res = where(finite(trend_NDVI) eq 0, count)
          if count ne 0 then trend_NDVI[res] = 0.
          res = where(trend_NDVI ne 0., count)
          
          if count ge nb/2 then begin
            conres = res[1:count-1]-res[0:count-2]
            continue_index = where(conres ge 3, count)
            if count eq 0 then begin
              V = trend_NDVI[res]
              Xout = indgen(nb)
              trend_NDVI = interpol(V, res, Xout,/spline)
              indic = 1
            endif else begin
            indic = 0
            endelse
          endif else begin
          indic = 0
          endelse
          
        endif else begin ;if aap
        indic = 0 
        endelse                                            
            
        if indic eq 1 then begin ; STSG
          if snow_address eq 1 then begin
          ; processing contaminated NDVI by snow
          snowres = where(vector_QA eq 2, count)
          if count ne 0 then begin
            bv_count = 0.
            bv_total = 0.
            for yeari = 0, yearnum-1 do begin
              res = where(img_QA[ii,jj,0:5,yeari] le 1, count)
              if count ne 0 then begin
                bv_total = bv_total+total(img_NDVI[ii,jj,res,yeari])
                bv_count = bv_count+count
              endif
            endfor
            if bv_count ne 0 then begin
              bv = bv_total/bv_count
              vector_in[snowres] =  bv
              trend_NDVI[snowres] = bv
            endif
          endif
          endif

          rst = trend_NDVI
          num_elements = n_elements(vector_in)
          ; Calculate the weights for each point
          gdis = 0.0
          fl = fltarr(num_elements)
          maxdif = max(abs(vector_in-rst))
          for i = 0,(num_elements-1) do begin
            case vector_QA[i] of 
            0: fl[i] = (0. > vector_in[i]-rst[i])
            1: fl[i] = vector_in[i]-rst[i]
            else: fl[i] = -1.0
            endcase
         endfor
         res = where(fl ne -1.0, count)
         if count ne 0 then meanfl = mean(fl[res])
         res = where(fl eq -1.0, count)
         if count ne 0 then fl[res] = meanfl
         for i =0,(num_elements-1) do begin
           fl[i] = (fl[i]-min(fl))/(max(fl)-min(fl))
           gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])
         endfor
         ra4 = fltarr(num_elements)
         pre = fltarr(num_elements)
         
         ormax = gdis
              
         ress = where(vector_QA eq 0, count)
         if count ne 0 then rst[ress] = vector_in[ress]
         ress = where(vector_QA ne 0 and vector_QA ne 1, count)
         if count ne 0 then vector_in[ress] = rst[ress]
    
         loop_times = 0l
         while (gdis le ormax) && loop_times LT 50 do begin
             loop_times = loop_times +1
             for i =0,(num_elements-1) do begin
               ra4[i] = (vector_in[i] ge rst[i]) ? vector_in[i] : rst[i]
               pre[i] = rst[i]
             endfor
             ; The Savitzky-Golay fitting
             ;savgolFilter = SAVGOL(4, 4, 0, 6) ;set the window width(4,4) and degree (6) for repetition
             savgolFilter = [-0.00543880, 0.0435097, -0.152289, 0.304585, 0.619267, 0.304585, -0.152289, 0.0435097, -0.00543880]
             rst = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
             ormax = gdis
             ; Calculate the fitting-effect index
             gdis = 0.0
             for i =0,(num_elements-1) do begin
                gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])
             endfor
          endwhile
           
          vec_fil = pre
          for smi = 0, nb-5 do begin
            a1 = vec_fil[smi]
            a2 = vec_fil[smi+1]
            a3 = vec_fil[smi+2]
            a4 = vec_fil[smi+3]
            a5 = vec_fil[smi+4]
            if (a1 gt a2) and (a2 lt a3) and (a3 gt a4) and (a4 lt a5) then begin
              pre[smi+1] = (a1+a3)/2.0
              pre[smi+3] = (a3+a5)/2.0
            endif 
          endfor
          vector_out[ii,jj-win,0:(nb-1),year] = pre
        endif
          
        if indic eq 0 then begin ;SG filter
          if snow_address eq 1 then begin
           ; processing contaminated NDVI by snow
          snowres = where(vector_QA eq 2, count)
          if count ne 0 then begin
            bv_count = 0.
            bv_total = 0.
            for yeari = 0, yearnum-1 do begin
              res = where(img_QA[ii,jj,0:5,yeari] le 1, count)
              if count ne 0 then begin
                bv_total = bv_total+total(img_NDVI[ii,jj,res,yeari])
                bv_count = bv_count+count
              endif
            endfor
            if bv_count ne 0 then begin
              bv = bv_total/bv_count
              vector_in[snowres] =  bv
            endif
          endif
          endif
          
          res = where(vector_QA le 2, count)
          if count lt nb then begin
            V = vector_in[res]
            Xout = indgen(nb)
            vector_in = interpol(V, res, Xout)
          endif
          
           num_elements = n_elements(vector_in)
           vector_in=reform(vector_in,num_elements)      
           ;savgolFilter = SAVGOL(4,4,0,2) ;set the window width(4,4) and degree (2) for computing trend curve
           savgolFilter = [-0.0909091, 0.0606061, 0.168831, 0.233766, 0.255411, 0.233766, 0.168831, 0.0606061, -0.0909091]
           rst = CONVOL(vector_in, savgolFilter, /EDGE_TRUNCATE)

           ; Calculate the weights for each point
           gdis = 0.0
           fl = fltarr(num_elements)
           maxdif = max(abs(vector_in-rst))
           for i =0,(num_elements-1) do begin
             fl[i] = (vector_in[i] ge rst[i]) ? 1.0 : (1-abs(vector_in[i]-rst[i])/maxdif)
             gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])
           endfor

           ra4 = fltARR(num_elements)
           pre = fltARR(num_elements)

           ormax = gdis
           num   = 0
    
           loop_times = 0l
           while (gdis le ormax) && loop_times LT 15 do begin
              loop_times = loop_times +1
              for i =0,(num_elements-1) do begin
                ra4[i] = (vector_in[i] ge rst[i]) ? vector_in[i] : rst[i]
                pre[i] = rst[i]
              endfor
              ; The Savitzky-Golay fitting
              ;savgolFilter = SAVGOL(4, 4, 0, 6)  ;set the window width(4,4) and degree (6) for repetition
              savgolFilter = [-0.00543880, 0.0435097, -0.152289, 0.304585, 0.619267, 0.304585, -0.152289, 0.0435097, -0.00543880]
              rst = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
              ormax = gdis
              ; Calculate the fitting-effect index
              gdis = 0.0
              for i =0,(num_elements-1) do begin
                gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])
              endfor
            endwhile
            vector_out[ii,jj-win,0:(nb-1),year] = pre 
          endif
          
        endif;if max(vector_in) gt 0.15 
      endfor;for year
    endfor
 endfor
    
 return, vector_out

    
 


end