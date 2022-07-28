
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
      VI_raw = reform(img_NDVI[ii,jj,*,year],nb)
      vector_QA = reform(img_QA[ii,jj,*,year],nb)
      res = reverse(sort(VI_raw))
      if mean(VI_raw[res[0:2]]) gt 0.15 then begin   
          
      num_elements = n_elements(VI_raw) 
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
           row = fix(new_corr[i+1]/(2*win+1))
           col = new_corr[i+1]-fix(new_corr[i+1]/(2*win+1))*(2*win+1)
           similar_index[1,i] = row + jj - win
           similar_index[0,i] = col + ii - win
           slope_intercept[1,i] = Slope_res[row, col]
           slope_intercept[0,i] = Intercept_res[row, col]  
           corr_similar[i] = new_corr_similar_res[row, col]   
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
          VI_raw = reform(img_NDVI[ii,jj,*,year],nb)
          vector_QA = float(reform(img_QA[ii,jj,*,year],nb))   
            
          for i = 0, nb-1 do begin ; loop for every doy
            for j = 0, samp-1 do begin
              row = similar_index[0,j]
              col = similar_index[1,j]
              if img_QA[row,col,i,year] le 1 then begin
                ; tar: ii,jj
                 new_ratio = img_NDVI[row, col, i, year]/reference_data[row, col, i]
                 temp_NDVI[i,j] = (slope_intercept[0,j]+new_ratio*slope_intercept[1,j])*reference_data[ii,jj,i] ; Eq. 2-3

                 if temp_NDVI[i,j] ge 1 or temp_NDVI[i,j] le -0.2 then temp_NDVI[i,j] = 0.  
              endif
            endfor
            
            res = where(finite(temp_NDVI[i,*]) eq 0, count) ; 查找NA
            if count ne 0 then temp_NDVI[i,res] = 0.
            res = where(temp_NDVI[i,*] ne 0., count)
            if count ne 0 then begin
               new_corr_similar = corr_similar
               ; 此处公式不一致,R_ij并未进行标准化处理
               new_corr_similar[res] = new_corr_similar[res]/total(new_corr_similar[res]) ; Eq. 5
               trend_NDVI[i] = total(new_corr_similar[res]*temp_NDVI[i,res])    ; the initial NDVI time-series synchronized by reference curve  
            endif
          endfor
          
          ; generating the trend_NDVI
          res = where(finite(trend_NDVI) eq 0, count)
          if count ne 0 then trend_NDVI[res] = 0.
          res = where(trend_NDVI ne 0., count)
          
          ; 进行插值处理的操作
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
              VI_raw[snowres] =  bv
              trend_NDVI[snowres] = bv
            endif
          endif
          endif

          VI_init = trend_NDVI
          num_elements = n_elements(VI_raw); 1:nptperyear
          ; Calculate the weights for each point
          fl = fltarr(num_elements)
          maxdif = max(abs(VI_raw-VI_init))
          for i = 0,(num_elements-1) do begin
            case vector_QA[i] of 
            0: fl[i] = (0. > VI_raw[i]-VI_init[i]) ; Eq. 6-1, max(0, VI_raw[i]-VI_init[i])
            1: fl[i] = VI_raw[i]-VI_init[i]        ; Eq. 6-2
            else: fl[i] = -1.0
            endcase
         endfor
         res = where(fl ne -1.0, count)
         if count ne 0 then begin
          ; 由于质量不佳的部分，会被VI_init替代，因此，这部分可以设置一个中等权重
          meanfl = mean(fl[res]) ; Eq. 6-3
          fl[res] = meanfl
         endif

         gdis = 0.0
         for i =0,(num_elements-1) do begin
           fl[i] = (fl[i]-min(fl))/(max(fl)-min(fl)); Eq. 7
           gdis = gdis + fl[i]*abs(VI_raw[i]-VI_init[i])
         endfor
         ra4 = fltarr(num_elements)
         pre = fltarr(num_elements)
         
         ormax = gdis
              
         ress = where(vector_QA eq 0, count)
         if count ne 0 then VI_init[ress] = VI_raw[ress]
         ress = where(vector_QA ne 0 and vector_QA ne 1, count)
         if count ne 0 then VI_raw[ress] = VI_init[ress]
    
         loop_times = 0l
         while (gdis le ormax) && loop_times LT 50 do begin
             loop_times = loop_times +1
             for i =0,(num_elements-1) do begin
               ra4[i] = (VI_raw[i] ge VI_init[i]) ? VI_raw[i] : VI_init[i]
               pre[i] = VI_init[i]
             endfor
             ; The Savitzky-Golay fitting
             ;savgolFilter = SAVGOL(4, 4, 0, 6) ;set the window width(4,4) and degree (6) for repetition
             savgolFilter = [-0.00543880, 0.0435097, -0.152289, 0.304585, 0.619267, 0.304585, -0.152289, 0.0435097, -0.00543880]
            ; 这里并未考虑权重，如何体现加权回归？
             VI_init = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
             ormax = gdis
             ; Calculate the fitting-effect index
             gdis = 0.0
             for i =0,(num_elements-1) do begin
                gdis = gdis + fl[i]*abs(VI_raw[i]-VI_init[i])
             endfor
          endwhile
          
          ; 5点平滑
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
              VI_raw[snowres] =  bv
            endif
          endif
          endif
          
          res = where(vector_QA le 2, count)
          if count lt nb then begin
            V = VI_raw[res]
            Xout = indgen(nb)
            VI_raw = interpol(V, res, Xout)
          endif
          
           num_elements = n_elements(VI_raw)
           VI_raw=reform(VI_raw,num_elements)      
           ;savgolFilter = SAVGOL(4,4,0,2) ;set the window width(4,4) and degree (2) for computing trend curve
           savgolFilter = [-0.0909091, 0.0606061, 0.168831, 0.233766, 0.255411, 0.233766, 0.168831, 0.0606061, -0.0909091]
           VI_init = CONVOL(VI_raw, savgolFilter, /EDGE_TRUNCATE)

           ; Calculate the weights for each point
           fl = fltarr(num_elements)
           maxdif = max(abs(VI_raw-VI_init))
           gdis = 0.0
           for i =0,(num_elements-1) do begin
             fl[i] = (VI_raw[i] ge VI_init[i]) ? 1.0 : (1-abs(VI_raw[i]-VI_init[i])/maxdif)
             gdis = gdis + fl[i]*abs(VI_raw[i]-VI_init[i])
           endfor

           ra4 = fltARR(num_elements)
           pre = fltARR(num_elements)

           ormax = gdis
           num   = 0
    
           loop_times = 0l
           while (gdis le ormax) && loop_times LT 15 do begin
              loop_times = loop_times +1
              for i =0,(num_elements-1) do begin
                ra4[i] = (VI_raw[i] ge VI_init[i]) ? VI_raw[i] : VI_init[i]
                pre[i] = VI_init[i]
              endfor
              ; The Savitzky-Golay fitting
              ;savgolFilter = SAVGOL(4, 4, 0, 6)  ;set the window width(4,4) and degree (6) for repetition
              savgolFilter = [-0.00543880, 0.0435097, -0.152289, 0.304585, 0.619267, 0.304585, -0.152289, 0.0435097, -0.00543880]
              VI_init = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
              ormax = gdis
              ; Calculate the fitting-effect index
              gdis = 0.0
              for i =0,(num_elements-1) do begin
                gdis = gdis + fl[i]*abs(VI_raw[i]-VI_init[i])
              endfor
            endwhile
            vector_out[ii,jj-win,0:(nb-1),year] = pre 
          endif
          
        endif;if max(VI_raw) gt 0.15 
      endfor;for year
    endfor
endfor
  
return, vector_out


end

function cal_gdist, dist_norm, VI_raw, VI_init
; gdis = cal_gdist(fl, VI_raw, VI_init)
  gdist = 0
  num_elements = n_elements(VI_raw) 
  for i =0,(num_elements-1) do begin
    gdist = gdist + dist_norm[i]*abs(VI_raw[i]-VI_init[i])
  endfor
  return, gdist
end
