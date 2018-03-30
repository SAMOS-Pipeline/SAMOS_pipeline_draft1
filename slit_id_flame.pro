pro slitid

  FITS_READ,'LMask1/LMask1master_flat.fits',data,header
  hdu_slits = mrdfits('LDSS3/2017-11-30/ccd8151c1.fits',0,header,/silent)
  print,size(hdu_slits)
	;print, SIZE(data)
	image = data
	approx_edge = 1412
	sz = SIZE(data)
	N_pixel_x = sz[1]
	cutout_size = 30
	binsize = 100
  	starting_pixel = N_pixel_x/2
  	x_edge = []
  	y_edge = []
 
	; starting guess for the ycoord
	previous_ycoord = approx_edge

	; starting from the center (because we assume that's where the approximate values refer to)
	; first do the right half
 	while starting_pixel LT N_pixel_x do begin
                ;print,"starting loop over: old ycoord ",previous_ycoord
    		; extract the bin
   		end_pixel = min([starting_pixel + binsize - 1, N_pixel_x-1])
    		cutout_bin = image[starting_pixel : end_pixel, previous_ycoord - cutout_size/2: previous_ycoord + cutout_size/2]
                ;print,end_pixel,previous_ycoord-cutout_size/2,previous_ycoord+cutout_size/2
    		; spatial profile
                profile = median(cutout_bin, dimension=1)
                ;print,size(cutout_bin)
                ;print,"[",profile,"]"
                ;print,size(profile)

    		; detect the edge
                derivative = shift(profile,1) - profile
                ;print,"shifted: [",shift(profile,1),"]"
                ;print,"not shifted: [",profile,"]"
                ;print,"derivative: [",derivative,"]"
    		derivative[0] = 0
    		derivative[-1] = 0
                peak = max(derivative, peak_location)
               ; print,"derivative is: ", derivative
               ; print,peak,peak_location
               ; print,"my peak ",(previous_ycoord-cutout_size/2)+peak_location

                
    		; add the amount of pixels left out by the cutout
                peak_location += (previous_ycoord - cutout_size/2)
                ;print, "new peak ",peak_location

                ; save x and y coordinate of the detection
                x_edge = [ x_edge, 0.5*(starting_pixel + end_pixel) ]
               
                y_edge = [ y_edge, peak_location ]
                ;print, "( ",0.5*(starting_pixel + end_pixel),peak_location," )"

    		; save the y coordinate for the next bin
    		previous_ycoord = peak_location
                ;print,"new y-coord ",previous_ycoord
                
    		; advance to next bin (to the right)
    		starting_pixel += binsize

             endwhile
        print,x_edge,y_edge

        ; and then, starting from the center, to the left
        starting_pixel = N_pixel_x/2
        previous_ycoord = approx_edge

        while starting_pixel GT 0 do begin

           ; extract the bin
           end_pixel = min([starting_pixel + binsize - 1, N_pixel_x-1])
           cutout_bin = image[starting_pixel : end_pixel, previous_ycoord - cutout_size/2: previous_ycoord + cutout_size/2]

           ; spatial profile
           profile = median(cutout_bin, dimension=1)

           ; detect the edge
           derivative = shift(profile,1) - profile
           derivative[0] = 0
           derivative[-1] = 0
           peak = max(derivative, peak_location)

           ; add the amount of pixels left out by the cutout
           peak_location += previous_ycoord - cutout_size/2

           ; save x and y coordinate of the detection
           x_edge = [ x_edge, 0.5*(starting_pixel + end_pixel) ]
           y_edge = [ y_edge, peak_location ]

           ; save the y coordinate for the next bin
           previous_ycoord = peak_location

           ; advance to next bin (to the right)
           starting_pixel -= binsize

        endwhile

        ; fill in the missing pixels with NaNs
        y_edge_full = fltarr(N_pixel_x) + !values.d_nan
        print,size(y_edge_full)
        print,y_edge_full
        y_edge_full[x_edge] = y_edge

        ; return array with the y-coordinates of the edges
        ;print, y_edge_full

      
     end
