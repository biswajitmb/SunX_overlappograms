pro whereis,image,w,x,y,z,print=print
;+
; NAME: whereis
;
; PURPOSE:
;   Given the 1-d index of a pixel in an array, return the
;   x and y coordinates corresponding to that pixel.
;
;
; NOTES:
;  pro whereis,image,w,x,y
;
; given the index of a pixel in an array return the
; x and y value
;
; jrg/ucb 94/8/16
; 
; if only one pixel, then return scalars rather than vectors
; 4/16/96 MCL
;-
;###########################################################################

if n_params() lt 4 and ~(keyword_set(print)) then begin
	message,'pro whereis, image, w, x, y, [z]'
endif

sz = size(image)


y = floor(w/sz(1))
x =  w - y*sz(1)

if n_elements(w) eq 1 then begin
    x =  x(0)
    y =  y(0)
endif

if keyword_set(print) then print," (X,Y) = "+printcoo(x,y)

end
