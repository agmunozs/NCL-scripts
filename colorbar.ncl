load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
  colors = (/"blue","white","green","white","red"/)

;---Specify number of colors in each range.
  opt                  = True
  opt@NumColorsInRange = (/ 9, 1, 4, 1, 9 /)
  rgb_array = span_named_colors(colors,opt)

  wks = gsn_open_wks("png","colormap")    ; Open a PS file

;----First frame.
  gsn_define_colormap(wks,rgb_array)
  gsn_draw_colormap(wks)

end
