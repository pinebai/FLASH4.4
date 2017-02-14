pro setcolors

  common palette,$
          r_col, g_col, b_col

    r_col = INDGEN(256)
    g_col = INDGEN(256)
    b_col = INDGEN(256)

    a_tmp=256/8*INDGEN(8)
    b_tmp=long(a_tmp)^2/max(a_tmp)
    c_tmp=ROUND(SQRT(a_tmp)*SQRT(max(a_tmp)))
    d_tmp=REVERSE(c_tmp)

    r_col(0)=0    & g_col(0)=0    & b_col(0)=0    ; black
    r_col(1)=100  & g_col(1)=100  & b_col(1)=100  ; dark grey
    r_col(2)=140  & g_col(2)=140  & b_col(2)=140  ;
    r_col(3)=180  & g_col(3)=180  & b_col(3)=180  ; light grey
    r_col(4)=150  & g_col(4)=150  & b_col(4)=150  ; grey
    r_col(5)=0    & g_col(5)=0    & b_col(5)=120  ; dark blue
    r_col(6)=0    & g_col(6)=0    & b_col(6)=160  ;
    r_col(7)=0    & g_col(7)=0    & b_col(7)=200  ;
    r_col(8)=60   & g_col(8)=210  & b_col(8)=200  ; turquoise
    r_col(9)=0    & g_col(9)=120  & b_col(9)=0    ; dark green
    r_col(10)=0   & g_col(10)=160 & b_col(10)=0   ;
    r_col(11)=0   & g_col(11)=200 & b_col(11)=0   ;
    r_col(12)=0   & g_col(12)=240 & b_col(12)=0   ; light green
    r_col(13)=140 & g_col(13)=110 & b_col(13)=0   ; dark orange
    r_col(14)=180 & g_col(14)=100 & b_col(14)=0   ;
    r_col(15)=220 & g_col(15)=90  & b_col(15)=0   ;
    r_col(16)=255 & g_col(16)=215 & b_col(16)=0   ; gold
    r_col(17)=120 & g_col(17)=0   & b_col(17)=0   ; dark red
    r_col(18)=160 & g_col(18)=0   & b_col(18)=0   ;
    r_col(19)=200 & g_col(19)=0   & b_col(19)=0   ;
    r_col(20)=240 & g_col(20)=0   & b_col(20)=0   ; light red
    r_col(21)=127 & g_col(21)=0   & b_col(21)=127 ; dark violet
    r_col(22)=127 & g_col(22)=0   & b_col(22)=167 ;
    r_col(23)=127 & g_col(23)=0   & b_col(23)=207 ;
    r_col(24)=255 & g_col(24)=105 & b_col(24)=180 ; pink
    r_col(25)=200 & g_col(25)=200 & b_col(25)=200 ; very light grey

    r_col(26)=255 & g_col(26)=0   & b_col(26)=0   ; red
    r_col(27)=0   & g_col(27)=255 & b_col(27)=0   ; green
    r_col(28)=0   & g_col(28)=0   & b_col(28)=255 ; blue
    r_col(29)=139 & g_col(29)=69  & b_col(29)=19  ; brown
    r_col(30)=160 & g_col(30)=32  & b_col(30)=240 ; violet
    r_col(31)=255 & g_col(31)=140 & b_col(31)=0   ; orange 
    r_col(32)=0   & g_col(32)=200 & b_col(32)=200 ; turquoise
    r_col(33)=255 & g_col(33)=0   & b_col(33)=255 ; magenta
    r_col(34)=50  & g_col(34)=50  & b_col(34)=50  ; very dark grey
    r_col(35)=225 & g_col(35)=225 & b_col(35)=225  ; very very light grey

    r_col(255)=255 & g_col(255)=255 & b_col(255)=255 ; white

    TVLCT,r_col,g_col,b_col,0

end ;--------------------------------------------------------------

