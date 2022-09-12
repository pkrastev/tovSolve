# tovSolve
This library provides a set of routines for solving the [Tolman-Oppenheimer-Volkoff](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation) equations to determine the structure of spherically symmetric neutron stars. The provided example main program is set to compute a sequence of cold neutron stars with the [MDI equation of state](https://arxiv.org/abs/nucl-th/0212090) (EOS).

## Compile

To compile the code use the provided <code>Makefile</code>.

## Example Output

```
$ ./tov.x 
         M           R            k2         lambda          I          beta         rhoc
     0.095584   118.842947     0.000174   411.212712     0.088190     0.001188     0.090000
     0.115570    37.401989     0.000485     3.546640     0.063674     0.004563     0.104242
     0.141874    24.443580     0.002786     2.429058     0.076745     0.008572     0.118485
     0.172160    19.555189     0.008374     2.392401     0.094947     0.013002     0.132727
     0.205533    17.109337     0.017010     2.491670     0.117936     0.017741     0.146970
     0.241186    15.700874     0.027527     2.624197     0.144121     0.022686     0.161212
     0.279000    14.805391     0.038857     2.761774     0.173609     0.027830     0.175455
     0.318700    14.200695     0.050045     2.887528     0.206049     0.033144     0.189697
     0.360015    13.774456     0.060483     2.996552     0.241145     0.038599     0.203939
     0.402888    13.462845     0.069870     3.087393     0.279219     0.044196     0.218182
     0.446953    13.230006     0.077968     3.157417     0.319559     0.049892     0.232424
     0.491876    13.052303     0.084789     3.209142     0.362209     0.055654     0.246667
     0.537403    12.913845     0.090169     3.235551     0.406190     0.061458     0.260909
     0.583306    12.803477     0.094525     3.249390     0.452209     0.067282     0.275152
     0.629378    12.713576     0.097675     3.241430     0.499033     0.073110     0.289394
     0.675679    12.638367     0.099919     3.218964     0.547618     0.078955     0.303636
     0.721829    12.574468     0.101268     3.180799     0.596756     0.084777     0.317879
     0.767636    12.518908     0.101800     3.127492     0.645739     0.090557     0.332121
     0.812937    12.469469     0.101867     3.068217     0.695186     0.096281     0.346364
     0.857593    12.424468     0.101332     2.997423     0.745698     0.101938     0.360606
     0.901487    12.382578     0.100261     2.916090     0.793749     0.107518     0.374848
     0.944520    12.342838     0.098932     2.831550     0.842291     0.113013     0.389091
     0.986609    12.304514     0.097316     2.742336     0.889726     0.118416     0.403333
     1.026746    12.267937     0.095532     2.652284     0.935445     0.123601     0.417576
     1.065893    12.231683     0.093543     2.558908     0.980620     0.128694     0.431818
     1.104132    12.195328     0.090900     2.449889     1.024148     0.133708     0.446061
     1.141433    12.158766     0.089184     2.367831     1.066946     0.138641     0.460303
     1.177775    12.121734     0.086873     2.271541     1.108885     0.143492     0.474545
     1.213141    12.084191     0.084566     2.177210     1.149402     0.148260     0.488788
     1.247518    12.046012     0.081757     2.071832     1.188368     0.152945     0.503030
     1.280896    12.007291     0.079813     1.990286     1.226369     0.157543     0.517273
     1.313270    11.967951     0.077421     1.899217     1.263056     0.162056     0.531515
     1.344639    11.927992     0.075084     1.811340     1.297964     0.166483     0.545758
     1.375004    11.887436     0.072718     1.724645     1.331815     0.170823     0.560000
     1.405636    11.844552     0.070368     1.639003     1.365573     0.175261     0.574242
     1.436031    11.800196     0.067980     1.553962     1.398699     0.179724     0.588485
     1.465301    11.755852     0.065622     1.472088     1.430650     0.184079     0.602727
     1.493414    11.711608     0.063325     1.394032     1.460909     0.188319     0.616970
     1.520349    11.667593     0.061190     1.321916     1.489901     0.192439     0.631212
     1.546097    11.623734     0.059119     1.253344     1.517135     0.196437     0.645455
     1.570657    11.580092     0.057108     1.188142     1.542898     0.200309     0.659697
     1.594038    11.536654     0.055187     1.126802     1.566850     0.204056     0.673939
     1.616257    11.493443     0.053344     1.068938     1.589534     0.207679     0.688182
     1.637334    11.450455     0.051583     1.014461     1.610555     0.211177     0.702424
     1.657298    11.407690     0.049766     0.960592     1.629588     0.214553     0.716667
     1.676181    11.365178     0.048181     0.912783     1.647598     0.217809     0.730909
     1.694015    11.322905     0.046791     0.870101     1.664282     0.220948     0.745152
     1.710837    11.280898     0.045341     0.827608     1.679409     0.223973     0.759394
     1.726684    11.239168     0.043955     0.787580     1.693074     0.226887     0.773636
     1.741594    11.197717     0.042620     0.749684     1.705157     0.229694     0.787879
     1.754836    11.158864     0.041419     0.715999     1.714895     0.232246     0.802121
     1.767245    11.120384     0.040299     0.684711     1.724419     0.234697     0.816364
     1.778945    11.081987     0.039205     0.654695     1.732195     0.237070     0.830606
     1.789965    11.043656     0.038159     0.626285     1.739942     0.239366     0.844848
     1.800332    11.005408     0.037141     0.599101     1.745311     0.241589     0.859091
     1.810073    10.967215     0.036167     0.573330     1.750520     0.243742     0.873333
     1.819217    10.929131     0.035238     0.548979     1.754275     0.245828     0.887576
     1.827788    10.891144     0.034360     0.526053     1.758007     0.247847     0.901818
     1.835812    10.853258     0.033425     0.502907     1.761101     0.249804     0.916061
     1.843312    10.815505     0.032649     0.482741     1.762120     0.251700     0.930303
     1.850313    10.777887     0.031871     0.463103     1.763076     0.253538     0.944545
     1.856836    10.740424     0.031080     0.443819     1.763589     0.255319     0.958788
     1.862903    10.703100     0.030360     0.426053     1.762745     0.257047     0.973030
     1.868534    10.665971     0.029647     0.408876     1.761495     0.258721     0.987273
     1.873750    10.629013     0.028956     0.392478     1.759891     0.260346     1.001515
     1.878569    10.592243     0.028300     0.377000     1.757705     0.261921     1.015758
     1.883008    10.555682     0.027673     0.362330     1.755032     0.263449     1.030000
     1.887085    10.519332     0.027042     0.348017     1.751648     0.264932     1.044242
     1.890818    10.483205     0.026449     0.334576     1.747952     0.266371     1.058485
     1.894220    10.447295     0.025883     0.321844     1.743839     0.267768     1.072727
     1.897307    10.411629     0.025332     0.309648     1.739343     0.269123     1.086970
     1.900094    10.376207     0.024761     0.297557     1.734383     0.270438     1.101212
     1.902594    10.341031     0.024245     0.286454     1.729184     0.271715     1.115455
     1.904822    10.306114     0.023797     0.276446     1.723406     0.272955     1.129697
     1.906787    10.271451     0.023314     0.266309     1.717494     0.274159     1.143939
     1.908503    10.237049     0.022859     0.256767     1.711586     0.275327     1.158182
     1.909983    10.202911     0.022405     0.247504     1.704867     0.276463     1.172424
     1.911235    10.169055     0.021975     0.238748     1.698113     0.277565     1.186667
     1.912271    10.135470     0.021565     0.230451     1.691061     0.278636     1.200909
     1.913100    10.102165     0.021123     0.222049     1.683853     0.279676     1.215152
     1.913733    10.069138     0.020763     0.214714     1.676564     0.280686     1.229394
     1.914179    10.036388     0.020360     0.207146     1.669421     0.281667     1.243636
     1.914445    10.003932     0.020018     0.200393     1.661429     0.282621     1.257879
     1.914542     9.971757     0.019677     0.193837     1.653463     0.283547     1.272121
     1.914474     9.939867     0.019331     0.187398     1.645540     0.284446     1.286364
     1.914253     9.908267     0.019000     0.181280     1.637302     0.285321     1.300606
     1.913883     9.876960     0.018681     0.175440     1.629288     0.286170     1.314848
     1.913374     9.845939     0.018363     0.169765     1.620574     0.286995     1.329091
     1.912729     9.815198     0.018065     0.164418     1.612226     0.287797     1.343333
     1.911958     9.784761     0.017778     0.159314     1.603663     0.288576     1.357576
     1.911064     9.754595     0.017489     0.154320     1.595223     0.289333     1.371818
     1.910055     9.724728     0.017190     0.149372     1.586620     0.290068     1.386061
     1.908935     9.695145     0.016952     0.145077     1.577513     0.290782     1.400303
     1.907711     9.665863     0.016700     0.140780     1.568872     0.291476     1.414545
     1.906387     9.636855     0.016449     0.136591     1.560058     0.292151     1.428788
     1.904968     9.608125     0.016185     0.132408     1.551511     0.292806     1.443030
     1.903458     9.579681     0.015969     0.128717     1.542428     0.293443     1.457273
     1.901863     9.551532     0.015742     0.125035     1.533678     0.294061     1.471515
     1.900187     9.523655     0.015516     0.121452     1.524656     0.294662     1.485758
     1.898433     9.496059     0.015306     0.118082     1.515743     0.295246     1.500000
```

## Example main program

```fortran
!=====================================================================
! MAIN PROGRAM
!=====================================================================
program tov2019
  use eos
  use constants
  use terminate
  use profile
  implicit none

  real(8)           :: rhoc
  real(8)           :: rho_start
  real(8)           :: rho_end
  real(8)           :: rho_tmp
  real(8)           :: rho_h
  integer(4)        :: nsteps
  integer(4)        :: i
  character(len=30) :: eos_file

  external tov_derivs, y_deriv, rkqc

  ! +++ Parameters +++
  eos_file  = 'eos_MDI_x0.0.in' ! EOS  
  rho_start = 0.09d0            ! Starting density
  rho_end   = 1.5d0             ! Ending density
  nsteps    = 100               ! Number of density steps
  ! +++ End Parameters +++

  rho_h = (rho_end - rho_start) / dfloat( nsteps - 1 )

  rho_tmp = rho_start

100 format(1x,t10,a,t22,a,t35,a,t46,a,t62a,t73,a,t86,a)
  write(6,100) 'M', 'R', 'k2', 'lambda', 'I', 'beta', 'rhoc'

  ! +++ Loop over central density +++
  do i = 1, nsteps
    rhoc = rho_tmp * fm3cm3
    call solve_tov(eos_file, rhoc)
    rho_tmp = rho_tmp + rho_h
  enddo

end program tov2019
```