s11.s22 <- function(TRpar,r, sgg, dist)
  {
   switch(sgg,
   "s11" =  r_moment_gamma_Weib(TRpar,r, dist)*(r_moment_gamma_Weib(TRpar,2-r, dist)-2*r_moment_gamma_Weib(TRpar,1, dist)*r_moment_gamma_Weib(TRpar,1-r, dist)
                    +r_moment_gamma_Weib(TRpar,1, dist)^2*r_moment_gamma_Weib(TRpar,-r, dist)),
   "s22" =  r_moment_gamma_Weib(TRpar,r, dist)*(-4*(2*r_moment_gamma_Weib(TRpar,1, dist)^2-r_moment_gamma_Weib(TRpar,2, dist))*r_moment_gamma_Weib(TRpar,1, dist)*r_moment_gamma_Weib(TRpar,1-r, dist)
                    +(2*r_moment_gamma_Weib(TRpar,1, dist)^2-r_moment_gamma_Weib(TRpar,2, dist))^2*r_moment_gamma_Weib(TRpar,-r, dist)
                    +(8*r_moment_gamma_Weib(TRpar,1, dist)^2-2*r_moment_gamma_Weib(TRpar,2, dist))*r_moment_gamma_Weib(TRpar,2-r, dist)
                    -4*r_moment_gamma_Weib(TRpar,1, dist)*r_moment_gamma_Weib(TRpar,3-r, dist)+r_moment_gamma_Weib(TRpar,4-r, dist)),
  )
}

